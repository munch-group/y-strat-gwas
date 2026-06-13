"""GxHaplogroup workflow for GenomeDK (gwf).

Two arms, both starting from plink1 genotypes:

  A. Genome-wide SNP x Y-haplogroup interaction scan (REGENIE --interaction),
     with PC x Hap product terms so ancestry structure (correlated with I vs R)
     cannot leak into the interaction estimate.

  B. Per-stratum (I, R) main-effect GWAS -> LDSC cross-stratum genetic
     correlation. r_g significantly below 1 indicates pervasive G x haplogroup
     interaction (the well-powered "general pattern" readout). Per-stratum
     SNP-h2 is also reported.

Run:
    pixi shell                 # or rely on `pixi run` inside the specs
    gwf -f workflow.py status
    gwf -f workflow.py run
"""
import os
from gwf import Workflow

# ---------------------------------------------------------------------------
# CONFIG  -- edit these (or override any value via the matching YS_* env var,
# e.g. YS_BFILE, YS_NPC; see tests/ for an end-to-end local run that does this)
# ---------------------------------------------------------------------------
def _env(name, default):
    return os.environ.get("YS_" + name, default)

ACCOUNT   = _env("ACCOUNT", None)   # slurm account
ROOT      = os.path.dirname(os.path.abspath(__file__))

# Defaults below point at the bundled SYNTHETIC TEST DATA (tests/work/data) so the
# workflow runs out of the box: build it once with
#   pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py
# then `gwf -b local run`. SWAP EACH for your real data when ready -- NPC and
# ENV_PREFIX below are likewise set for the synthetic data / sandbox.
BFILE     = _env("BFILE", "%s/tests/work/data/genotypes" % ROOT)             # plink1 .bed/.bim/.fam (autosomes)
PHENO     = _env("PHENO", "%s/tests/work/data/phenotypes.txt" % ROOT)        # FID IID autism (1=case,0=control,NA)
BASECOVAR = _env("BASECOVAR", "%s/tests/work/data/covariates_base.txt" % ROOT)  # FID IID PC1..PCk age batch ...
HAPFILE   = _env("HAPFILE", "%s/tests/work/data/haplogroup.txt" % ROOT)      # FID IID Hap (I or R)
PHENONAME = _env("PHENONAME", "autism")
NPC       = int(_env("NPC", "10"))                         # synthetic data has 10 PCs (set yours for real data)
CATCOVAR  = _env("CATCOVAR", "batch")                      # comma-sep categorical covars, or ""

# Parallelise the REGENIE step-2 scans (interaction, gwas_{I,R}) by chromosome:
# a fan-out job per chromosome + a gather. Accepts ranges/lists, e.g. "1-22" or
# "1,2,5-9"; "" = single job. Results are identical to the single job. Defaults
# to the autosomes so a cluster run fans out; set "" for one job per scan.
SPLIT_CHROMS = _env("SPLIT_CHROMS", "1-22")

# Split the Monte-Carlo permutation targets (perm_interaction, lava_perm_*) into
# this many independent batches (different seeds, pooled) for parallel speedup.
# 1 = single job.
PERM_BATCHES = int(_env("PERM_BATCHES", "4"))

# Extra environment prepended to every command. Defaulted to KMP_AFFINITY=disabled
# for the bundled synthetic data (it dodges an MKL/OpenMP thread-affinity assertion
# seen in some sandboxes/containers, and being baked into the spec strings it also
# reaches a gwf worker whose own environment lacks it). CLEAR IT for a real cluster
# run (set ENV_PREFIX = "" or YS_ENV_PREFIX="") unless you actually need it.
ENV_PREFIX = _env("ENV_PREFIX", "KMP_AFFINITY=disabled")

# ancestry-matched permutation re-test of interaction hits (Arm A, post-hoc).
PERM_TOP          = int(_env("PERM_TOP", "500"))      # top interaction hits to re-test
PERM_PANEL        = int(_env("PERM_PANEL", "5000"))   # random SNPs for the lambda diagnostic
PERM_NPERM        = int(_env("PERM_NPERM", "1000"))   # permutations per hit
PERM_GLOBAL_NPERM = int(_env("PERM_GLOBAL_NPERM", "200"))  # permutations for lambda
PERM_FORCE_SNPS   = _env("PERM_FORCE_SNPS", "")       # comma-sep SNP IDs to always re-test

# LAVA local cross-stratum genetic correlation (optional 3rd arm). Set LAVA_LOCI
# to enable; empty disables. Format: "name=chr:start-end,name2=chr:start-end".
LAVA_LOCI       = _env("LAVA_LOCI", "")  # "" disables LAVA; synthetic: realblock=20:9988-2241939,confblock=21:8710-2271328
LAVA_REF        = _env("LAVA_REF", BFILE)             # plink ref for in-sample LD
LAVA_PARTITION  = _env("LAVA_PARTITION", "")          # LAVA blocks file (R headline only)
LAVA_N_CONTROLS = int(_env("LAVA_N_CONTROLS", "300")) # negative-control loci
LAVA_NPERM      = int(_env("LAVA_NPERM", "10000"))    # permutations per locus
RSCRIPT         = _env("RSCRIPT", "Rscript")          # for the LAVA R step

# chrX (optional): plink1 prefix with X variants coded as chromosome 23 and the
# SEX column filled in the .fam (1=male). Empty string disables all X targets.
# REGENIE codes non-PAR males as hemizygous (0/2) automatically from sex + build.
XBFILE       = _env("XBFILE", "")  # "" => skip chrX; synthetic X is tests/work/data/genotypesX
GENOME_BUILD = _env("GENOME_BUILD", "hg38")               # for REGENIE --par-region

# LDSC reference data (download separately; see README). Defaults point at the
# bundled synthetic reference -- swap EUR_LD/HM3 for the real 1000G EUR / HapMap3.
LDSC_DIR  = _env("LDSC_DIR", "%s/ldsc" % ROOT)             # py3 LDSC fork, vendored + pip-installed by pixi
EUR_LD    = _env("EUR_LD", "%s/tests/work/eur_w_ld_chr" % ROOT)   # EUR LD scores dir (synthetic ref)
HM3       = _env("HM3", "%s/tests/work/data/w_hm3.snplist" % ROOT)  # HapMap3 SNP list for --merge-alleles

# liability-scale h2 (edit to your male-specific population prevalence)
PREV_POP  = float(_env("PREV_POP", "0.03"))

# run REGENIE step1 separately within each stratum (cleaner r_g, ~2x compute).
# False = reuse the full-sample step1 predictors (cheaper, slightly conservative
# toward r_g = 1).
STRATUM_SPECIFIC_STEP1 = _env("STRATUM_SPECIFIC_STEP1", "True").lower() in ("1", "true", "yes")

# ---------------------------------------------------------------------------
# derived
# ---------------------------------------------------------------------------
PIXI    = ("%s " % ENV_PREFIX if ENV_PREFIX else "") \
          + "pixi run --manifest-path %s/pixi.toml" % ROOT
LDSCRUN = "%s python" % PIXI                   # py3 LDSC fork, pip-installed by pixi postinstall

OUT = _env("OUT", "%s/results" % ROOT)
TMP = _env("TMP", "%s/tmp" % ROOT)
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

PCS       = ["PC%d" % i for i in range(1, NPC + 1)]
BASE_COLS = ",".join(PCS + ["age"])
INT_COLS  = ",".join(PCS + ["age", "Hap"] + ["%s_x_Hap" % p for p in PCS])
CAT       = ("--catCovarList %s" % CATCOVAR) if CATCOVAR else ""

def _chroms(spec):
    """Parse "1-22" / "1,2,5-9" into a list of ints."""
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out

defaults = {"account": ACCOUNT} if ACCOUNT else {}
gwf = Workflow(defaults=defaults)

# ---------------------------------------------------------------------------
# QC backbone for REGENIE step 1
# ---------------------------------------------------------------------------
gwf.target("qc",
           inputs=["%s.bed" % BFILE],
           outputs=["%s/qc_pass.snplist" % OUT],
           cores=8, memory="16g", walltime="02:00:00") << """
{pixi} plink2 --bfile {bfile} \
  --maf 0.01 --mac 100 --geno 0.05 --hwe 1e-15 \
  --indep-pairwise 1000 100 0.9 \
  --out {out}/qc
mv {out}/qc.prune.in {out}/qc_pass.snplist
""".format(pixi=PIXI, bfile=BFILE, out=OUT)

# ---------------------------------------------------------------------------
# prep: strata keep-lists + interaction covariate file
# ---------------------------------------------------------------------------
gwf.target("strata",
           inputs=[HAPFILE],
           outputs=["%s/keep_I.txt" % OUT, "%s/keep_R.txt" % OUT],
           cores=1, memory="2g", walltime="00:10:00") << """
{pixi} python {root}/scripts/make_strata.py --hap {hap} --out-prefix {out}
""".format(pixi=PIXI, root=ROOT, hap=HAPFILE, out=OUT)

gwf.target("int_covars",
           inputs=[BASECOVAR, HAPFILE],
           outputs=["%s/int_covars.txt" % OUT],
           cores=1, memory="4g", walltime="00:10:00") << """
{pixi} python {root}/scripts/make_interaction_covars.py \
  --covar {covar} --hap {hap} --npc {npc} --out {out}/int_covars.txt
""".format(pixi=PIXI, root=ROOT, covar=BASECOVAR, hap=HAPFILE, npc=NPC, out=OUT)

# ---------------------------------------------------------------------------
# ARM A: full-sample step1 + interaction scan
# ---------------------------------------------------------------------------
gwf.target("step1_full",
           inputs=["%s/qc_pass.snplist" % OUT],
           outputs=["%s/step1_full_pred.list" % OUT],
           cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 1 --bed {bfile} \
  --extract {out}/qc_pass.snplist \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --bsize 1000 --lowmem --lowmem-prefix {tmp}/step1_full \
  --threads 16 --out {out}/step1_full
""".format(pixi=PIXI, bfile=BFILE, out=OUT, pheno=PHENO, ph=PHENONAME,
           covar=BASECOVAR, base=BASE_COLS, cat=CAT, tmp=TMP)

def _interaction_step2(chr_flag, out_prefix):
    return """
{pixi} regenie --step 2 --bed {bfile} {chrflag} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {out}/int_covars.txt --covarColList {intcols} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --interaction Hap \
  --bsize 400 --minMAC 20 \
  --pred {out}/step1_full_pred.list \
  --threads 16 --out {outpref}
""".format(pixi=PIXI, bfile=BFILE, chrflag=chr_flag, out=OUT, pheno=PHENO,
           ph=PHENONAME, intcols=INT_COLS, cat=CAT, outpref=out_prefix)

_GXHAP = "%s/gxhap_%s.regenie" % (OUT, PHENONAME)
if SPLIT_CHROMS:
    _chunks = []
    for _c in _chroms(SPLIT_CHROMS):
        _ch = "%s/gxhap_chr%d_%s.regenie" % (OUT, _c, PHENONAME)
        gwf.target("interaction_chr%d" % _c,
                   inputs=["%s/step1_full_pred.list" % OUT, "%s/int_covars.txt" % OUT],
                   outputs=[_ch], cores=16, memory="32g", walltime="24:00:00") \
            << _interaction_step2("--chr %d" % _c, "%s/gxhap_chr%d" % (OUT, _c))
        _chunks.append(_ch)
    gwf.target("interaction", inputs=_chunks, outputs=[_GXHAP],
               cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/concat_regenie.py --out {gx} --inputs {chunks}
""".format(pixi=PIXI, root=ROOT, gx=_GXHAP, chunks=" ".join(_chunks))
else:
    gwf.target("interaction",
               inputs=["%s/step1_full_pred.list" % OUT, "%s/int_covars.txt" % OUT],
               outputs=[_GXHAP], cores=16, memory="32g", walltime="24:00:00") \
        << _interaction_step2("", "%s/gxhap" % OUT)

gwf.target("top_int",
           inputs=["%s/gxhap_%s.regenie" % (OUT, PHENONAME)],
           outputs=["%s/top_interactions.tsv" % OUT],
           cores=2, memory="8g", walltime="00:30:00") << """
{pixi} python {root}/scripts/top_interactions.py \
  --regenie {out}/gxhap_{ph}.regenie --out {out}/top_interactions.tsv
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME)

# Post-hoc: re-test the interaction hits against an ancestry-matched permutation
# null (permute Hap within PC-space strata). This is the G-side guard REGENIE's
# --interaction cannot provide -- a real I-vs-R hit survives it, an ancestry
# artifact collapses. Controls confounding only to global-PC resolution; residual
# local ancestry at a locus is not removed (flag survivors for local ancestry /
# replication). See scripts/ancestry_matched_perm.py.
_FORCE = ("--force-snps %s" % PERM_FORCE_SNPS) if PERM_FORCE_SNPS else ""

def _perm_cmd(nperm, gnperm, seed, extra, out_prefix):
    return """
{pixi} python {root}/scripts/ancestry_matched_perm.py \
  --regenie {out}/gxhap_{ph}.regenie \
  --bfile {bfile} --covar {covar} --hap {hap} \
  --pheno {pheno} --pheno-name {ph} --npc {npc} \
  --top {top} --panel {panel} --nperm {nperm} --global-nperm {gnperm} {force} \
  --seed {seed} --select-seed 1 {extra} --out-prefix {outpref}
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME, bfile=BFILE,
           covar=BASECOVAR, hap=HAPFILE, pheno=PHENO, npc=NPC, top=PERM_TOP,
           panel=PERM_PANEL, nperm=nperm, gnperm=gnperm, force=_FORCE,
           seed=seed, extra=extra, outpref=out_prefix)

_PERM_INPUTS = ["%s/gxhap_%s.regenie" % (OUT, PHENONAME), "%s.bed" % BFILE]
if PERM_BATCHES > 1:
    _pb_nperm = max(1, PERM_NPERM // PERM_BATCHES)
    _pb_gnperm = max(1, PERM_GLOBAL_NPERM // PERM_BATCHES)
    _pcounts, _plnulls = [], []
    for _k in range(PERM_BATCHES):
        _cp = "%s/perm_b%d_counts.tsv" % (OUT, _k)
        _lp = "%s/perm_b%d_lambda_null.tsv" % (OUT, _k)
        gwf.target("perm_interaction_batch%d" % _k, inputs=_PERM_INPUTS,
                   outputs=[_cp, _lp], cores=8, memory="16g", walltime="08:00:00") \
            << _perm_cmd(_pb_nperm, _pb_gnperm, _k + 1, "--raw-counts",
                         "%s/perm_b%d" % (OUT, _k))
        _pcounts.append(_cp)
        _plnulls.append(_lp)
    gwf.target("perm_interaction", inputs=_pcounts + _plnulls,
               outputs=["%s/perm_interactions.tsv" % OUT, "%s/perm_lambda.txt" % OUT],
               cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/pool_perm.py --kind interaction \
  --counts {counts} --lambda-nulls {lnulls} --out-prefix {out}/perm
""".format(pixi=PIXI, root=ROOT, out=OUT, counts=" ".join(_pcounts),
           lnulls=" ".join(_plnulls))
else:
    gwf.target("perm_interaction", inputs=_PERM_INPUTS,
               outputs=["%s/perm_interactions.tsv" % OUT, "%s/perm_lambda.txt" % OUT,
                        "%s/perm_strata_selection.tsv" % OUT],
               cores=8, memory="16g", walltime="08:00:00") \
        << _perm_cmd(PERM_NPERM, PERM_GLOBAL_NPERM, 1, "", "%s/perm" % OUT)

# ---------------------------------------------------------------------------
# ARM B: per-stratum GWAS -> munge -> rg / h2
# ---------------------------------------------------------------------------
for s in ("I", "R"):
    keep = "%s/keep_%s.txt" % (OUT, s)

    if STRATUM_SPECIFIC_STEP1:
        gwf.target("step1_%s" % s,
                   inputs=[keep, "%s/qc_pass.snplist" % OUT],
                   outputs=["%s/step1_%s_pred.list" % (OUT, s)],
                   cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 1 --bed {bfile} --keep {keep} \
  --extract {out}/qc_pass.snplist \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --bsize 1000 --lowmem --lowmem-prefix {tmp}/step1_{s} \
  --threads 16 --out {out}/step1_{s}
""".format(pixi=PIXI, bfile=BFILE, keep=keep, out=OUT, pheno=PHENO,
           ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT, tmp=TMP, s=s)
        pred = "%s/step1_%s_pred.list" % (OUT, s)
        pred_inputs = [pred]
    else:
        pred = "%s/step1_full_pred.list" % OUT
        pred_inputs = [pred]

    def _gwas_step2(keep, pred, s, chr_flag, out_prefix):
        return """
{pixi} regenie --step 2 --bed {bfile} --keep {keep} {chrflag} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --bsize 400 --minMAC 20 \
  --pred {pred} \
  --threads 16 --out {outpref}
""".format(pixi=PIXI, bfile=BFILE, keep=keep, chrflag=chr_flag, out=OUT,
           pheno=PHENO, ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT,
           pred=pred, outpref=out_prefix)

    _gw = "%s/gwas_%s_%s.regenie" % (OUT, s, PHENONAME)
    if SPLIT_CHROMS:
        _gchunks = []
        for _c in _chroms(SPLIT_CHROMS):
            _gch = "%s/gwas_%s_chr%d_%s.regenie" % (OUT, s, _c, PHENONAME)
            gwf.target("gwas_%s_chr%d" % (s, _c),
                       inputs=[keep] + pred_inputs, outputs=[_gch],
                       cores=16, memory="32g", walltime="24:00:00") \
                << _gwas_step2(keep, pred, s, "--chr %d" % _c,
                               "%s/gwas_%s_chr%d" % (OUT, s, _c))
            _gchunks.append(_gch)
        gwf.target("gwas_%s" % s, inputs=_gchunks, outputs=[_gw],
                   cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/concat_regenie.py --out {gw} --inputs {chunks}
""".format(pixi=PIXI, root=ROOT, gw=_gw, chunks=" ".join(_gchunks))
    else:
        gwf.target("gwas_%s" % s,
                   inputs=[keep] + pred_inputs, outputs=[_gw],
                   cores=16, memory="32g", walltime="24:00:00") \
            << _gwas_step2(keep, pred, s, "", "%s/gwas_%s" % (OUT, s))

    gwf.target("munge_%s" % s,
               inputs=["%s/gwas_%s_%s.regenie" % (OUT, s, PHENONAME)],
               outputs=["%s/munged_%s.sumstats.gz" % (OUT, s)],
               cores=2, memory="8g", walltime="01:00:00") << """
{pixi} python {root}/scripts/regenie_to_munge.py \
  --regenie {out}/gwas_{s}_{ph}.regenie --out {out}/gwas_{s}.forldsc.txt
{ldsc} {ldscdir}/munge_sumstats.py \
  --sumstats {out}/gwas_{s}.forldsc.txt \
  --merge-alleles {hm3} \
  --snp SNP --a1 A1 --a2 A2 --N-col N --p P --signed-sumstats BETA,0 \
  --chunksize 500000 \
  --out {out}/munged_{s}
""".format(pixi=PIXI, root=ROOT, out=OUT, s=s, ph=PHENONAME,
           ldsc=LDSCRUN, ldscdir=LDSC_DIR, hm3=HM3)

    gwf.target("h2_%s" % s,
               inputs=["%s/munged_%s.sumstats.gz" % (OUT, s), keep],
               outputs=["%s/h2_%s.log" % (OUT, s)],
               cores=2, memory="8g", walltime="00:30:00") << """
SAMP=$({pixi} python {root}/scripts/samp_prev.py \
  --pheno {pheno} --pheno-name {ph} --keep {keep})
{ldsc} {ldscdir}/ldsc.py \
  --h2 {out}/munged_{s}.sumstats.gz \
  --ref-ld-chr {eurld}/ --w-ld-chr {eurld}/ \
  --samp-prev $SAMP --pop-prev {pop} \
  --out {out}/h2_{s}
""".format(pixi=PIXI, root=ROOT, pheno=PHENO, ph=PHENONAME, keep=keep,
           ldsc=LDSCRUN, ldscdir=LDSC_DIR, out=OUT, s=s, eurld=EUR_LD, pop=PREV_POP)

gwf.target("rg",
           inputs=["%s/munged_I.sumstats.gz" % OUT, "%s/munged_R.sumstats.gz" % OUT],
           outputs=["%s/I_vs_R_rg.log" % OUT],
           cores=2, memory="8g", walltime="01:00:00") << """
{ldsc} {ldscdir}/ldsc.py \
  --rg {out}/munged_I.sumstats.gz,{out}/munged_R.sumstats.gz \
  --ref-ld-chr {eurld}/ --w-ld-chr {eurld}/ \
  --out {out}/I_vs_R_rg
{pixi} python {root}/scripts/parse_ldsc_rg.py --log {out}/I_vs_R_rg.log
""".format(ldsc=LDSCRUN, ldscdir=LDSC_DIR, out=OUT, eurld=EUR_LD,
           pixi=PIXI, root=ROOT)

# ---------------------------------------------------------------------------
# chrX (optional; only built when XBFILE is set). LDSC is autosomal, so the X
# chromosome -- the mechanistically interesting target (X-linked regulators x
# Hap) -- gets its own REGENIE step-2 pass reusing the autosomal step-1
# predictors. Males in the non-PAR region are coded hemizygously (0/2) by
# REGENIE using the SEX column of XBFILE's .fam plus the --par-region build.
#
#   Arm A:  interaction_X -> top_int_X     (SNP x Hap interaction scan on chrX)
#   Arm B:  gwas_X_{I,R}  -> xforldsc_{I,R} (per-stratum chrX summary stats)
#
# The per-stratum chrX sumstats are written in LDSC-munge-ready form
# (results/gwas_X_{s}.forldsc.txt). Running LDSC munge/h2/rg on them additionally
# needs an X-chromosome LD-score reference and an X-inclusive --merge-alleles
# list (the bundled autosomal eur_w_ld_chr / HapMap3 list do not cover chrX), so
# that final step is left to the user once an X reference is available.
# ---------------------------------------------------------------------------
if XBFILE:
    XHEMI = "--par-region %s --sex-specific male" % GENOME_BUILD

    gwf.target("xqc",
               inputs=["%s.bed" % XBFILE],
               outputs=["%s/x_qc_pass.snplist" % OUT],
               cores=4, memory="8g", walltime="01:00:00") << """
{pixi} plink2 --bfile {xbfile} --chr 23 \
  --maf 0.01 --mac 20 --geno 0.05 \
  --write-snplist \
  --out {out}/x_qc
mv {out}/x_qc.snplist {out}/x_qc_pass.snplist
""".format(pixi=PIXI, xbfile=XBFILE, out=OUT)

    gwf.target("interaction_X",
               inputs=["%s/step1_full_pred.list" % OUT, "%s/int_covars.txt" % OUT,
                       "%s/x_qc_pass.snplist" % OUT],
               outputs=["%s/gxhapX_%s.regenie" % (OUT, PHENONAME)],
               cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 2 --bed {xbfile} \
  --chr 23 --extract {out}/x_qc_pass.snplist {xhemi} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {out}/int_covars.txt --covarColList {intcols} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --interaction Hap \
  --bsize 400 --minMAC 20 \
  --pred {out}/step1_full_pred.list \
  --threads 16 --out {out}/gxhapX
""".format(pixi=PIXI, xbfile=XBFILE, xhemi=XHEMI, out=OUT, pheno=PHENO,
           ph=PHENONAME, intcols=INT_COLS, cat=CAT)

    gwf.target("top_int_X",
               inputs=["%s/gxhapX_%s.regenie" % (OUT, PHENONAME)],
               outputs=["%s/top_interactions_X.tsv" % OUT],
               cores=2, memory="8g", walltime="00:30:00") << """
{pixi} python {root}/scripts/top_interactions.py \
  --regenie {out}/gxhapX_{ph}.regenie --out {out}/top_interactions_X.tsv
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME)

    for s in ("I", "R"):
        keep = "%s/keep_%s.txt" % (OUT, s)
        pred = ("%s/step1_%s_pred.list" % (OUT, s) if STRATUM_SPECIFIC_STEP1
                else "%s/step1_full_pred.list" % OUT)

        gwf.target("gwas_X_%s" % s,
                   inputs=[keep, pred, "%s/x_qc_pass.snplist" % OUT],
                   outputs=["%s/gwas_X_%s_%s.regenie" % (OUT, s, PHENONAME)],
                   cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 2 --bed {xbfile} --keep {keep} \
  --chr 23 --extract {out}/x_qc_pass.snplist {xhemi} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --bsize 400 --minMAC 20 \
  --pred {pred} \
  --threads 16 --out {out}/gwas_X_{s}
""".format(pixi=PIXI, xbfile=XBFILE, xhemi=XHEMI, keep=keep, out=OUT,
           pheno=PHENO, ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS,
           cat=CAT, pred=pred, s=s)

        gwf.target("xforldsc_%s" % s,
                   inputs=["%s/gwas_X_%s_%s.regenie" % (OUT, s, PHENONAME)],
                   outputs=["%s/gwas_X_%s.forldsc.txt" % (OUT, s)],
                   cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/regenie_to_munge.py \
  --regenie {out}/gwas_X_{s}_{ph}.regenie --out {out}/gwas_X_{s}.forldsc.txt
""".format(pixi=PIXI, root=ROOT, out=OUT, s=s, ph=PHENONAME)

# ---------------------------------------------------------------------------
# LAVA arm (optional; only built when LAVA_LOCI is set). Local cross-stratum
# genetic correlation -- the local-resolution version of Arm B's r_g. It is
# corrupted by the SAME ancestry confounding as Arm A (differential LD between
# the I and R backgrounds depresses local r_g even with no epistasis), so the
# headline LAVA estimate is paired with an ancestry-matched permutation test
# plus a negative-control-loci panel. See METHODS.md.
#
#   lava_pcs                 region-excluded PCs (anti-circularity)
#   lava_perm_<name>         per-locus deconfounding: ancestry-matched permutation
#                            + negative controls (pure Python; the testable core)
#   lava_inputs / lava_local optional headline local r_g via the real LAVA R tool
#                            (needs R + the LAVA package + LAVA_PARTITION; see README)
# ---------------------------------------------------------------------------
if LAVA_LOCI:
    _loci = [e.split("=", 1) for e in LAVA_LOCI.split(",") if "=" in e]
    _loci = [(nm.strip(), loc.strip()) for nm, loc in _loci]
    _chroms = ",".join(sorted({loc.split(":")[0] for _, loc in _loci}, key=str))

    # region-excluded PCs: drop the chromosome(s) carrying the test loci so the
    # locus cannot leak into the ancestry features it is matched on.
    gwf.target("lava_pcs",
               inputs=["%s.bed" % BFILE],
               outputs=["%s/lava_pcs.eigenvec" % OUT],
               cores=8, memory="16g", walltime="04:00:00") << """
{pixi} plink2 --bfile {bfile} --not-chr {chroms} \
  --maf 0.01 --indep-pairwise 200 50 0.1 --out {out}/lava_prune
{pixi} plink2 --bfile {bfile} --extract {out}/lava_prune.prune.in \
  --pca {npc} --out {out}/lava_pcs
""".format(pixi=PIXI, bfile=BFILE, chroms=_chroms, out=OUT, npc=NPC)

    def _lava_cmd(nm, loc, nperm, seed, extra, out_prefix):
        return """
{pixi} python {root}/scripts/local_rg_perm.py \
  --bfile {bfile} --pheno {pheno} --pheno-name {ph} --hap {hap} \
  --pcs {out}/lava_pcs.eigenvec --npc {npc} \
  --locus {loc} --locus-name {nm} \
  --n-controls {nc} --nperm {nperm} --seed {seed} --select-seed 1 {extra} \
  --out-prefix {outpref}
""".format(pixi=PIXI, root=ROOT, bfile=BFILE, pheno=PHENO, ph=PHENONAME,
           hap=HAPFILE, out=OUT, npc=NPC, loc=loc, nm=nm, nc=LAVA_N_CONTROLS,
           nperm=nperm, seed=seed, extra=extra, outpref=out_prefix)

    for nm, loc in _loci:
        _li = ["%s/lava_pcs.eigenvec" % OUT, "%s.bed" % BFILE]
        _lout = ["%s/lava_%s_loci.tsv" % (OUT, nm), "%s/lava_%s_summary.txt" % (OUT, nm)]
        if PERM_BATCHES > 1:
            _lb_nperm = max(1, LAVA_NPERM // PERM_BATCHES)
            _lcounts = []
            for _k in range(PERM_BATCHES):
                _cp = "%s/lava_%s_b%d_counts.tsv" % (OUT, nm, _k)
                gwf.target("lava_perm_%s_batch%d" % (nm, _k), inputs=_li,
                           outputs=[_cp, "%s/lava_%s_b%d_meta.tsv" % (OUT, nm, _k)],
                           cores=8, memory="16g", walltime="12:00:00") \
                    << _lava_cmd(nm, loc, _lb_nperm, _k + 1, "--raw-counts",
                                 "%s/lava_%s_b%d" % (OUT, nm, _k))
                _lcounts.append(_cp)
            gwf.target("lava_perm_%s" % nm, inputs=_lcounts, outputs=_lout,
                       cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/pool_perm.py --kind lava \
  --counts {counts} --meta {out}/lava_{nm}_b0_meta.tsv --out-prefix {out}/lava_{nm}
""".format(pixi=PIXI, root=ROOT, out=OUT, nm=nm, counts=" ".join(_lcounts))
        else:
            gwf.target("lava_perm_%s" % nm, inputs=_li, outputs=_lout,
                       cores=8, memory="16g", walltime="12:00:00") \
                << _lava_cmd(nm, loc, LAVA_NPERM, 1, "", "%s/lava_%s" % (OUT, nm))

    # Optional headline: real LAVA local r_g via the R package. Built only when a
    # LAVA partition file is supplied; the R step runs outside the pixi env.
    if LAVA_PARTITION:
        gwf.target("lava_inputs",
                   inputs=["%s/munged_I.sumstats.gz" % OUT,
                           "%s/munged_R.sumstats.gz" % OUT],
                   outputs=["%s/lava_input.info" % OUT,
                            "%s/lava_sample_overlap.txt" % OUT,
                            "%s/lava_loci.tsv" % OUT],
                   cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/lava_inputs.py \
  --pheno {pheno} --pheno-name {ph} \
  --keep-i {out}/keep_I.txt --keep-r {out}/keep_R.txt \
  --sumstats-i {out}/gwas_I.forldsc.txt --sumstats-r {out}/gwas_R.forldsc.txt \
  --loci "{loci}" --out-prefix {out}/lava
""".format(pixi=PIXI, root=ROOT, pheno=PHENO, ph=PHENONAME, out=OUT,
           loci=LAVA_LOCI)

        gwf.target("lava_local",
                   inputs=["%s/lava_input.info" % OUT, "%s/lava_loci.tsv" % OUT],
                   outputs=["%s/lava_local_rg.tsv" % OUT],
                   cores=4, memory="16g", walltime="04:00:00") << """
{rscript} {root}/scripts/lava_local.R \
  {out}/lava_input.info {out}/lava_sample_overlap.txt {ref} \
  {out}/lava_loci.tsv {out}/lava_local_rg.tsv
""".format(rscript=RSCRIPT, root=ROOT, out=OUT, ref=LAVA_REF)
