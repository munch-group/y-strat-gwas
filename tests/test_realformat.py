#!/usr/bin/env python3
"""Checks the real-data input adapters: a haplogroup file with the lineage in a
differently-named column plus non-I/R males to drop, and a covariate file with
C-prefixed PCs and no age/batch.

Reuses the data built by run_pipeline_test.py (run that first), specifically
tests/work/data/haplogroup_major.txt and raw_covariates.txt.
"""
import os
import subprocess
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "tests", "work", "data")
T = os.path.join(ROOT, "tests", "work", "realfmt")
PIXI = "pixi run --manifest-path %s/pixi.toml" % ROOT
os.makedirs(T, exist_ok=True)
ok = True


def py(args):
    subprocess.run("%s python %s" % (PIXI, args), shell=True, cwd=ROOT,
                   executable="/bin/bash", check=True)


def check(label, passed, detail=""):
    global ok
    ok &= passed
    print("  [%s] %s%s" % ("PASS" if passed else "FAIL", label,
                           ("  -- " + detail) if detail else ""))


def main():
    hapf = "%s/haplogroup_major.txt" % D
    rawc = "%s/raw_covariates.txt" % D
    if not (os.path.exists(hapf) and os.path.exists(rawc)):
        sys.exit("missing real-format data -- run tests/run_pipeline_test.py first")

    maj = pd.read_csv(hapf, sep=r"\s+")
    n_total = len(maj)
    n_ir = int(maj["Major"].str.upper().isin(["I", "R"]).sum())

    # 1. make_strata drops non-I/R and writes keep_IR
    py("scripts/make_strata.py --hap %s --hap-col Major --out-prefix %s" % (hapf, T))
    keep_ir = pd.read_csv("%s/keep_IR.txt" % T, sep=r"\s+", header=None)
    keep_i = pd.read_csv("%s/keep_I.txt" % T, sep=r"\s+", header=None)
    keep_r = pd.read_csv("%s/keep_R.txt" % T, sep=r"\s+", header=None)
    check("make_strata reads 'Major' col and drops non-I/R",
          len(keep_ir) == n_ir and n_ir < n_total,
          "keep_IR=%d of %d (dropped %d)" % (len(keep_ir), n_total, n_total - n_ir))
    check("keep_I + keep_R == keep_IR", len(keep_i) + len(keep_r) == len(keep_ir),
          "%d + %d == %d" % (len(keep_i), len(keep_r), len(keep_ir)))
    # a dropped (non-I/R) individual is in none of the keep lists
    dropped = maj.loc[~maj["Major"].str.upper().isin(["I", "R"]), "IID"].iloc[0]
    in_keep = dropped in set(keep_ir[1].astype(str))
    check("a non-I/R individual is excluded", not in_keep, "IID=%s" % dropped)

    # 2. prep_covar renames C1..C10 -> PC1..PC10 and adds dummy age/batch
    py("scripts/prep_covar.py --raw-covar %s --npc 10 --pc-prefix C "
       "--add-quant age --add-cat batch --out %s/base_covar.txt" % (rawc, T))
    prepped = pd.read_csv("%s/base_covar.txt" % T, sep=r"\s+")
    want = ["FID", "IID"] + ["PC%d" % i for i in range(1, 11)] + ["age", "batch"]
    check("prep_covar produces FID IID PC1..PC10 age batch",
          list(prepped.columns) == want, ",".join(prepped.columns))
    raw = pd.read_csv(rawc, sep=r"\s+")
    check("PC values preserved from C columns",
          bool((abs(prepped["PC1"] - raw["C1"]) < 1e-6).all()), "PC1 == C1")

    # 3. make_interaction_covars uses the prepped covar + 'Major', I/R only
    py("scripts/make_interaction_covars.py --covar %s/base_covar.txt --hap %s "
       "--hap-col Major --npc 10 --out %s/int_covars.txt" % (T, hapf, T))
    intc = pd.read_csv("%s/int_covars.txt" % T, sep=r"\s+")
    check("int_covars restricted to I/R individuals", len(intc) == n_ir,
          "%d rows == %d I/R" % (len(intc), n_ir))
    check("int_covars has Hap (0/1) and PCxHap products",
          set(intc["Hap"].unique()) <= {0, 1} and "PC10_x_Hap" in intc.columns,
          "Hap levels=%s" % sorted(intc["Hap"].unique()))

    # 4. prep_pheno recodes a HEADERLESS PLINK .pheno (FID IID value), 1/2 -> 0/1
    py("scripts/prep_pheno.py --raw-pheno %s/raw_phenotypes.txt --pheno-name autism "
       "--no-header --value-col 3 --case 2 --control 1 --out %s/pheno.txt" % (D, T))
    raw = pd.read_csv("%s/raw_phenotypes.txt" % D, sep=r"\s+", header=None,
                      names=["FID", "IID", "autism_r"])
    check("raw .pheno is headerless (FID IID value)",
          list(raw.columns) == ["FID", "IID", "autism_r"]
          and bool(raw["autism_r"].isin([1, 2]).any()),
          "%d rows" % len(raw))
    m = raw.merge(pd.read_csv("%s/pheno.txt" % T, sep=r"\s+"), on=["FID", "IID"])
    rawv = pd.to_numeric(m["autism_r"], errors="coerce")
    newv = pd.to_numeric(m["autism"], errors="coerce")
    check("prep_pheno: headerless .pheno, case 2->1, control 1->0",
          bool((newv[rawv == 2] == 1).all()) and bool((newv[rawv == 1] == 0).all()),
          "%d cases, %d controls" % (int((newv == 1).sum()), int((newv == 0).sum())))

    # 5. gwf wiring: with RAW_* set, prep tasks exist AND step1_full depends on them
    env = dict(os.environ, YS_RAW_PHENO="%s/raw_phenotypes.txt" % D,
               YS_RAW_COVAR="%s/raw_covariates.txt" % D,
               YS_HAPFILE="%s/haplogroup_major.txt" % D, YS_HAPCOL="Major",
               YS_PC_PREFIX="C")
    code = ("import workflow as w; s=w.gwf.targets['step1_full']; "
            "ins=[i.split('/')[-1] for i in s.inputs]; "
            "print(int('prep_covar' in w.gwf.targets and 'prep_pheno' in w.gwf.targets "
            "and 'pheno_recoded.txt' in ins and 'base_covar.txt' in ins))")
    r = subprocess.run("%s python -c \"%s\"" % (PIXI, code), shell=True, cwd=ROOT,
                       executable="/bin/bash", env=env, capture_output=True, text=True)
    check("gwf orders prep_covar/prep_pheno before step1_full",
          r.stdout.strip().endswith("1"), r.stdout.strip() or r.stderr.strip()[-80:])

    print("\n%s" % ("ALL REAL-FORMAT CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
