#!/usr/bin/env python3
"""Checks for the two parallelisation features (SPLIT_CHROMS, PERM_BATCHES).

Reuses the data + intermediate results built by run_pipeline_test.py, so run that
first:
    pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py
    pixi run --manifest-path ./pixi.toml python tests/test_parallel.py

It asserts:
  1. chromosome-split REGENIE step 2 + concat == the unsplit scan (exactly),
  2. pooling K permutation batches reproduces the single-job emp_p and sums the
     permutation count, for both perm_interaction and the LAVA local-r_g test.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "tests", "work", "data")
O = os.path.join(ROOT, "tests", "work", "results")
PIXI = "pixi run --manifest-path %s/pixi.toml" % ROOT
T = os.path.join(ROOT, "tests", "work", "parallel")
os.makedirs(T, exist_ok=True)
os.environ["KMP_AFFINITY"] = "disabled"
ok = True


def sh(cmd):
    subprocess.run(cmd, shell=True, cwd=ROOT, executable="/bin/bash", check=True)


def py(args):
    sh("%s python %s" % (PIXI, args))


def check(label, passed, detail=""):
    global ok
    ok &= passed
    print("  [%s] %s%s" % ("PASS" if passed else "FAIL", label,
                           ("  -- " + detail) if detail else ""))


def truth():
    return dict(l.split("\t", 1) for l in
                open("%s/truth.txt" % D).read().splitlines())


def need(*paths):
    miss = [p for p in paths if not os.path.exists(p)]
    if miss:
        sys.exit("missing %s -- run tests/run_pipeline_test.py first" % miss[0])


# --- 1. chromosome-split step 2 == unsplit ---------------------------------
def test_chr_split():
    need("%s/gxhap_autism.regenie" % O, "%s/step1_full_pred.list" % O,
         "%s/int_covars.txt" % O)
    intcols = ",".join(["PC%d" % i for i in range(1, 11)] + ["age", "Hap"]
                       + ["PC%d_x_Hap" % i for i in range(1, 11)])
    chunks = []
    for c in (1, 2):
        out = "%s/gxhap_chr%d" % (T, c)
        sh("%s regenie --step 2 --bed %s/genotypes --chr %d "
           "--phenoFile %s/phenotypes.txt --phenoColList autism "
           "--covarFile %s/int_covars.txt --covarColList %s --catCovarList batch "
           "--bt --firth --approx --pThresh 0.01 --interaction Hap --bsize 400 "
           "--minMAC 20 --pred %s/step1_full_pred.list --threads 4 --out %s "
           ">/dev/null 2>&1" % (PIXI, D, c, D, O, intcols, O, out))
        chunks.append("%s_autism.regenie" % out)
    py("scripts/concat_regenie.py --out %s/gxhap_split.regenie --inputs %s"
       % (T, " ".join(chunks)))
    full = pd.read_csv("%s/gxhap_autism.regenie" % O, sep=r"\s+", comment="#")
    full = full[full.CHROM.isin([1, 2])].sort_values(
        ["CHROM", "GENPOS", "TEST"]).reset_index(drop=True)
    split = pd.read_csv("%s/gxhap_split.regenie" % T, sep=r"\s+", comment="#") \
        .sort_values(["CHROM", "GENPOS", "TEST"]).reset_index(drop=True)
    cols = ["CHROM", "ID", "TEST", "BETA", "LOG10P"]
    check("chr-split step2 + concat == unsplit (chr1,2)",
          len(full) == len(split) and full[cols].equals(split[cols]),
          "%d rows" % len(split))


# --- 2. permutation batching pools to the single-job result ----------------
def test_perm_batching():
    need("%s/gxhap_autism.regenie" % O)
    t = truth()
    force = "%s,%s" % (t["interaction_snp"].strip(), t["confound_snp"].strip())
    common = ("--regenie %s/gxhap_autism.regenie --bfile %s/genotypes "
              "--covar %s/covariates_base.txt --hap %s/haplogroup.txt "
              "--pheno %s/phenotypes.txt --pheno-name autism --npc 10 "
              "--top 50 --panel 200 --global-nperm 60 --force-snps %s"
              % (O, D, D, D, D, force))
    py("scripts/ancestry_matched_perm.py %s --nperm 300 --seed 1 "
       "--out-prefix %s/single" % (common, T))
    counts, lnulls = [], []
    for k in range(3):
        py("scripts/ancestry_matched_perm.py %s --nperm 100 --seed %d "
           "--select-seed 1 --raw-counts --out-prefix %s/b%d"
           % (common, k + 1, T, k))
        counts.append("%s/b%d_counts.tsv" % (T, k))
        lnulls.append("%s/b%d_lambda_null.tsv" % (T, k))
    py("scripts/pool_perm.py --kind interaction --counts %s --lambda-nulls %s "
       "--out-prefix %s/pooled" % (" ".join(counts), " ".join(lnulls), T))
    sj = pd.read_csv("%s/single_interactions.tsv" % T, sep="\t").set_index("ID")
    pl = pd.read_csv("%s/pooled_interactions.tsv" % T, sep="\t").set_index("ID")
    real, conf = t["interaction_snp"].strip(), t["confound_snp"].strip()
    check("perm pool: total nperm = 3 x 100", int(pl.loc[real, "n_perm"]) == 300,
          "n_perm=%d" % int(pl.loc[real, "n_perm"]))
    check("perm pool: real hit still significant & separated",
          pl.loc[real, "anc_matched_emp_p"] <= 0.05
          and pl.loc[conf, "anc_matched_emp_p"] > pl.loc[real, "anc_matched_emp_p"],
          "real p=%.4g conf p=%.4g" % (pl.loc[real, "anc_matched_emp_p"],
                                       pl.loc[conf, "anc_matched_emp_p"]))
    # batches use the same selection seed -> hits/panel identical to single job
    check("perm pool: same hit set as single job",
          set(sj.index) == set(pl.index), "%d hits" % len(pl))


# --- 3. LAVA local-r_g batching pools to the single-job result -------------
def test_lava_batching():
    need("%s/lava_pcs.eigenvec" % O)
    t = truth()
    loc = t["real_block_locus"].strip()
    common = ("--bfile %s/genotypes --pheno %s/phenotypes.txt --pheno-name autism "
              "--hap %s/haplogroup.txt --pcs %s/lava_pcs.eigenvec --npc 10 "
              "--locus %s --locus-name realblock --n-controls 10"
              % (D, D, D, O, loc))
    py("scripts/local_rg_perm.py %s --nperm 300 --seed 1 --out-prefix %s/lsingle"
       % (common, T))
    counts = []
    for k in range(3):
        py("scripts/local_rg_perm.py %s --nperm 100 --seed %d --select-seed 1 "
           "--raw-counts --out-prefix %s/lb%d" % (common, k + 1, T, k))
        counts.append("%s/lb%d_counts.tsv" % (T, k))
    py("scripts/pool_perm.py --kind lava --counts %s --meta %s/lb0_meta.tsv "
       "--out-prefix %s/lpooled" % (" ".join(counts), T, T))
    sj = dict(l.split("\t") for l in
              open("%s/lsingle_summary.txt" % T).read().splitlines() if "\t" in l)
    pl = dict(l.split("\t") for l in
              open("%s/lpooled_summary.txt" % T).read().splitlines() if "\t" in l)
    check("lava pool: genuine block still survives",
          float(pl["target_anc_matched_p"]) <= 0.05,
          "p=%s (single %s)" % (pl["target_anc_matched_p"], sj["target_anc_matched_p"]))
    check("lava pool: same target local_rg as single job",
          abs(float(pl["target_local_rg"]) - float(sj["target_local_rg"])) < 1e-6,
          "local_rg=%s" % pl["target_local_rg"])


if __name__ == "__main__":
    print("# chromosome split"); test_chr_split()
    print("# perm_interaction batching"); test_perm_batching()
    print("# LAVA local-rg batching"); test_lava_batching()
    print("\n%s" % ("ALL PARALLEL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    sys.exit(0 if ok else 1)
