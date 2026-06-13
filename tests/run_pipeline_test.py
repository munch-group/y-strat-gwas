#!/usr/bin/env python3
"""End-to-end local smoke test for the GxHaplogroup pipeline.

It does NOT need Slurm. It builds a synthetic dataset + a synthetic LDSC
reference, points the workflow at them via YS_* env vars, then imports
workflow.py and runs *its actual target specs* (the same shell commands gwf
would submit) in dependency order, asserting each output appears.

Run from the repo root:
    pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py

What it proves: the helper scripts, REGENIE step1/step2 (incl. --interaction),
the regenie->munge->LDSC munge/h2/rg chain, and the file formats all line up,
and that workflow.py generates runnable commands. Estimates are noise at this
scale -- correctness here means "every stage runs and emits well-formed output".
"""
import os
import subprocess
import sys

ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK     = os.path.join(ROOT, "tests", "work")
DATA     = os.path.join(WORK, "data")
REF      = os.path.join(WORK, "eur_w_ld_chr")
OUT      = os.path.join(WORK, "results")
TMP      = os.path.join(WORK, "tmp")
PIXI     = "pixi run --manifest-path %s/pixi.toml" % ROOT
NPC      = 10


def _truth():
    return dict(l.split("\t", 1) for l in
                open("%s/truth.txt" % DATA).read().splitlines())


def _force_snps():
    """real autosomal interaction SNP + confound SNP, so both are re-tested."""
    t = _truth()
    return ",".join(t[k].strip() for k in ("interaction_snp", "confound_snp"))


def _lava_loci():
    """real (genuine) + confound local-divergence blocks for the LAVA arm."""
    t = _truth()
    return "realblock=%s,confblock=%s" % (t["real_block_locus"].strip(),
                                          t["conf_block_locus"].strip())


def sh(cmd, **kw):
    """Run a shell command, streaming output; raise on failure."""
    print("\n$ " + cmd)
    subprocess.run(cmd, shell=True, cwd=ROOT, executable="/bin/bash",
                   check=True, **kw)


def pixi_py(script_and_args):
    sh("%s python %s" % (PIXI, script_and_args))


def build_inputs():
    for d in (WORK, DATA, REF, OUT, TMP):
        os.makedirs(d, exist_ok=True)

    # 1. synthetic genotypes + pheno/covar/hap + LDSC reference LD scores
    pixi_py("tests/make_dummy_data.py --out-dir %s --ref-dir %s" % (DATA, REF))

    # 2. VCF -> plink1 bed/bim/fam.  --id-delim _ splits the "FIDk_IIDk" VCF
    #    sample names back into matching FID / IID columns (else regenie sees
    #    every phenotype as missing).
    sh("%s plink2 --vcf %s/geno.vcf --id-delim _ --make-bed --max-alleles 2 "
       "--out %s/genotypes" % (PIXI, DATA, DATA))

    # 2b. chrX fileset (plink codes contig X as chromosome 23), with SEX=1
    #     (male) applied so REGENIE codes the non-PAR region hemizygously.
    sh("%s plink2 --vcf %s/genoX.vcf --id-delim _ --max-alleles 2 "
       "--update-sex %s/sex.txt --make-bed "
       "--out %s/genotypesX" % (PIXI, DATA, DATA, DATA))

    # 3. HapMap3-style merge-alleles list = every SNP in the panel (SNP A1 A2)
    bim = "%s/genotypes.bim" % DATA
    snplist = "%s/w_hm3.snplist" % DATA
    with open(bim) as f, open(snplist, "w") as o:
        o.write("SNP\tA1\tA2\n")
        for line in f:
            c = line.split()
            o.write("%s\t%s\t%s\n" % (c[1], c[4], c[5]))


def run_workflow_targets():
    """Import workflow.py with YS_* overrides and execute each target's spec."""
    env = {
        "YS_ACCOUNT":   "local_test",
        "YS_BFILE":     "%s/genotypes" % DATA,
        "YS_PHENO":     "%s/phenotypes.txt" % DATA,
        "YS_BASECOVAR": "%s/covariates_base.txt" % DATA,
        "YS_HAPFILE":   "%s/haplogroup.txt" % DATA,
        "YS_NPC":       str(NPC),
        "YS_CATCOVAR":  "batch",
        "YS_EUR_LD":    REF,
        "YS_HM3":       "%s/w_hm3.snplist" % DATA,
        "YS_PREV_POP":  "0.05",
        "YS_STRATUM_SPECIFIC_STEP1": "True",
        "YS_XBFILE":    "%s/genotypesX" % DATA,
        "YS_GENOME_BUILD": "hg38",
        # ancestry-matched permutation: small counts for a fast test; force the
        # planted real + confound SNPs so both are always re-tested.
        "YS_PERM_TOP":          "100",
        "YS_PERM_PANEL":        "300",
        "YS_PERM_NPERM":        "300",
        "YS_PERM_GLOBAL_NPERM": "100",
        "YS_PERM_FORCE_SNPS":   _force_snps(),
        # LAVA local-rg arm: enable on the two planted blocks, small counts.
        "YS_LAVA_LOCI":         _lava_loci(),
        "YS_LAVA_N_CONTROLS":   "12",
        "YS_LAVA_NPERM":        "500",
        "YS_OUT":       OUT,
        "YS_TMP":       TMP,
        # MKL/OpenMP thread-affinity binding asserts in some sandboxed/container
        # environments; disabling it is harmless and unrelated to the pipeline.
        "KMP_AFFINITY": "disabled",
        "OMP_PROC_BIND": "FALSE",
        # This in-process smoke test executes targets in a fixed order, so pin the
        # single-job path regardless of the workflow's split defaults (the split /
        # batch fan-out is covered by tests/test_parallel.py + tests/run_via_gwf.sh).
        "YS_SPLIT_CHROMS": "",
        "YS_PERM_BATCHES": "1",
    }
    os.environ.update(env)
    sys.path.insert(0, ROOT)
    import workflow as wf   # noqa: E402  (import after env is set, by design)

    # dependency order (mirrors README target graph)
    order = ["qc", "strata", "int_covars",
             "step1_full", "interaction", "top_int", "perm_interaction",
             "step1_I", "gwas_I", "munge_I", "h2_I",
             "step1_R", "gwas_R", "munge_R", "h2_R",
             "rg",
             # LAVA local-rg arm (only present when LAVA_LOCI is set)
             "lava_pcs", "lava_perm_realblock", "lava_perm_confblock",
             # chrX (only present when XBFILE is set)
             "xqc", "interaction_X", "top_int_X",
             "gwas_X_I", "xforldsc_I", "gwas_X_R", "xforldsc_R"]
    order = [t for t in order if t in wf.gwf.targets]

    for name in order:
        target = wf.gwf.targets[name]
        print("\n" + "=" * 70 + "\n>>> TARGET %s\n" % name + "=" * 70)
        subprocess.run(target.spec, shell=True, cwd=target.working_dir,
                       executable="/bin/bash", check=True)
        for out in target.outputs:
            if not os.path.exists(out):
                raise SystemExit("FAIL: target %s did not produce %s" % (name, out))
            print("  ok output: %s" % out)
    return wf


def validate(wf):
    import numpy as np
    import pandas as pd
    print("\n" + "#" * 70 + "\n# VALIDATION\n" + "#" * 70)
    ok = True

    def check(label, passed, detail=""):
        nonlocal ok
        ok &= passed
        print("  [%s] %s%s" % ("PASS" if passed else "FAIL", label,
                               ("  -- " + detail) if detail else ""))

    # --- Arm A: the planted SNPxHap interaction should rank near the top ----
    truth = dict(l.split("\t", 1) for l in
                 open("%s/truth.txt" % DATA).read().splitlines())
    true_snp = truth["interaction_snp"].strip()
    reg = pd.read_csv("%s/gxhap_%s.regenie" % (OUT, wf.PHENONAME),
                      sep=r"\s+", comment="#")
    inter = reg[reg.TEST == "ADD-INT_SNPxHap"].sort_values(
        "LOG10P", ascending=False).reset_index(drop=True)
    rank = int(inter.index[inter.ID == true_snp][0]) + 1
    pct = 100.0 * rank / len(inter)
    check("interaction scan runs (ADD-INT_SNPxHap rows)", len(inter) > 0,
          "%d variants" % len(inter))
    check("planted interaction SNP %s in top 5%%" % true_snp, pct <= 5.0,
          "rank %d/%d (%.1f%%)" % (rank, len(inter), pct))

    n_int = sum(1 for _ in open("%s/top_interactions.tsv" % OUT)) - 1
    check("top_interactions.tsv populated", n_int > 0, "%d rows" % n_int)

    # --- ancestry-matched permutation: real hit survives, confound collapses --
    perm = pd.read_csv("%s/perm_interactions.tsv" % OUT, sep="\t").set_index("ID")
    conf_snp = truth["confound_snp"].strip()
    p_real = float(perm.loc[true_snp, "anc_matched_emp_p"])
    p_conf = float(perm.loc[conf_snp, "anc_matched_emp_p"])
    check("real interaction SNP survives ancestry-matched perm", p_real <= 0.10,
          "%s emp_p=%.3g" % (true_snp, p_real))
    check("confound SNP separated from real (p_conf > p_real)", p_conf > p_real,
          "%s emp_p=%.3g vs real %.3g" % (conf_snp, p_conf, p_real))
    sel = dict(l.split("\t") for l in
               open("%s/perm_lambda.txt" % OUT).read().splitlines() if "\t" in l)
    check("strata decoupled Hap from PCs (within-stratum AUC ~0.5)",
          float(sel["within_stratum_HapPC_AUC"]) <= 0.60,
          "AUC=%s, eff_N=%s" % (sel["within_stratum_HapPC_AUC"], sel.get("eff_N")))

    # --- LAVA: genuine local divergence survives, confound collapses ---------
    def _summ(name):
        return dict(l.split("\t") for l in
                    open("%s/lava_%s_summary.txt" % (OUT, name)).read().splitlines()
                    if "\t" in l)
    real_s, conf_s = _summ("realblock"), _summ("confblock")
    rp = float(real_s["target_anc_matched_p"])
    cp = float(conf_s["target_anc_matched_p"])
    check("LAVA: genuine local divergence survives permutation", rp <= 0.05,
          "realblock local_rg=%s p=%.3g" % (real_s["target_local_rg"], rp))
    check("LAVA: confound block separated from genuine (p_conf > p_real)", cp > rp,
          "confblock local_rg=%s p=%.3g vs real %.3g"
          % (conf_s["target_local_rg"], cp, rp))
    check("LAVA: genuine block stands out vs negative controls",
          float(real_s["target_vs_controls_tail_frac"]) <= 0.20,
          "tail-fraction=%s" % real_s["target_vs_controls_tail_frac"])

    # --- chrX: hemizygous step-2 pass should recover the planted X interaction
    if "interaction_snp_X" in truth:
        true_x = truth["interaction_snp_X"].strip()
        xreg = pd.read_csv("%s/gxhapX_%s.regenie" % (OUT, wf.PHENONAME),
                           sep=r"\s+", comment="#")
        xinter = xreg[xreg.TEST == "ADD-INT_SNPxHap"].sort_values(
            "LOG10P", ascending=False).reset_index(drop=True)
        xrank = int(xinter.index[xinter.ID == true_x][0]) + 1
        xpct = 100.0 * xrank / len(xinter)
        check("chrX interaction scan runs (hemizygous step-2)", len(xinter) > 0,
              "%d X variants" % len(xinter))
        check("planted chrX interaction SNP %s in top 10%%" % true_x, xpct <= 10.0,
              "rank %d/%d (%.1f%%)" % (xrank, len(xinter), xpct))
        for s in ("I", "R"):
            f = "%s/gwas_X_%s.forldsc.txt" % (OUT, s)
            nrows = sum(1 for _ in open(f)) - 1
            check("chrX per-stratum sumstats ready for LDSC (stratum %s)" % s,
                  nrows > 0, "%d variants" % nrows)

    # --- Arm B: rg parsed, h2 reported on the liability scale ---------------
    rg_txt = open("%s/I_vs_R_rg.log" % OUT).read()
    check("cross-stratum rg summary present",
          "Summary of Genetic Correlation Results" in rg_txt)

    for s in ("I", "R"):
        h2_txt = open("%s/h2_%s.log" % (OUT, s)).read()
        line = next((l for l in h2_txt.splitlines()
                     if "Liability scale h2" in l), "")
        check("Arm B liability-scale h2 reported (stratum %s)" % s,
              bool(line), line.strip())

    print("\n%s" % ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    return ok


if __name__ == "__main__":
    build_inputs()
    wf = run_workflow_targets()
    sys.exit(0 if validate(wf) else 1)
