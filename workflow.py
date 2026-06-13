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

ACCOUNT   = _env("ACCOUNT", "your_genomedk_project")   # slurm account
ROOT      = os.path.dirname(os.path.abspath(__file__))

BFILE     = _env("BFILE", "/path/to/genotypes")            # plink1 prefix .bed/.bim/.fam (autosomes)
PHENO     = _env("PHENO", "/path/to/phenotypes.txt")       # FID IID autism   (1=case, 0=control, NA)
BASECOVAR = _env("BASECOVAR", "/path/to/covariates_base.txt")  # FID IID PC1..PCk age batch ...
HAPFILE   = _env("HAPFILE", "/path/to/haplogroup.txt")     # FID IID Hap      (I or R)
PHENONAME = _env("PHENONAME", "autism")
NPC       = int(_env("NPC", "20"))
CATCOVAR  = _env("CATCOVAR", "batch")                      # comma-sep categorical covars, or ""

# ancestry-matched permutation re-test of interaction hits (Arm A, post-hoc).
PERM_TOP          = int(_env("PERM_TOP", "500"))      # top interaction hits to re-test
PERM_PANEL        = int(_env("PERM_PANEL", "5000"))   # random SNPs for the lambda diagnostic
PERM_NPERM        = int(_env("PERM_NPERM", "1000"))   # permutations per hit
PERM_GLOBAL_NPERM = int(_env("PERM_GLOBAL_NPERM", "200"))  # permutations for lambda
PERM_FORCE_SNPS   = _env("PERM_FORCE_SNPS", "")       # comma-sep SNP IDs to always re-test

# LAVA local cross-stratum genetic correlation (optional 3rd arm). Set LAVA_LOCI
# to enable; empty disables. Format: "name=chr:start-end,name2=chr:start-end".
LAVA_LOCI       = _env("LAVA_LOCI", "")               # "" disables the LAVA arm
LAVA_REF        = _env("LAVA_REF", BFILE)             # plink ref for in-sample LD
LAVA_PARTITION  = _env("LAVA_PARTITION", "")          # LAVA blocks file (R headline only)
LAVA_N_CONTROLS = int(_env("LAVA_N_CONTROLS", "300")) # negative-control loci
LAVA_NPERM      = int(_env("LAVA_NPERM", "10000"))    # permutations per locus
RSCRIPT         = _env("RSCRIPT", "Rscript")          # for the LAVA R step

# chrX (optional): plink1 prefix with X variants coded as chromosome 23 and the
# SEX column filled in the .fam (1=male). Empty string disables all X targets.
# REGENIE codes non-PAR males as hemizygous (0/2) automatically from sex + build.
XBFILE       = _env("XBFILE", "")                          # "" => skip chrX entirely
GENOME_BUILD = _env("GENOME_BUILD", "hg38")               # for REGENIE --par-region

# LDSC reference data (download separately; see README)
LDSC_DIR  = _env("LDSC_DIR", "%s/ldsc" % ROOT)             # py3 LDSC fork, vendored + pip-installed by pixi
EUR_LD    = _env("EUR_LD", "/path/to/eur_w_ld_chr")        # 1000G EUR LD scores (directory)
HM3       = _env("HM3", "/path/to/w_hm3.snplist")          # HapMap3 SNP list for --merge-alleles

# liability-scale h2 (edit to your male-specific population prevalence)
PREV_POP  = float(_env("PREV_POP", "0.03"))

# run REGENIE step1 separately within each stratum (cleaner r_g, ~2x compute).
# False = reuse the full-sample step1 predictors (cheaper, slightly conservative
# toward r_g = 1).
STRATUM_SPECIFIC_STEP1 = _env("STRATUM_SPECIFIC_STEP1", "True").lower() in ("1", "true", "yes")

# ---------------------------------------------------------------------------
# derived
# ---------------------------------------------------------------------------
PIXI    = "pixi run --manifest-path %s/pixi.toml" % ROOT
LDSCRUN = "%s python" % PIXI                   # py3 LDSC fork, pip-installed by pixi postinstall

OUT = _env("OUT", "%s/results" % ROOT)
TMP = _env("TMP", "%s/tmp" % ROOT)
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

PCS       = ["PC%d" % i for i in range(1, NPC + 1)]
BASE_COLS = ",".join(PCS + ["age"])
INT_COLS  = ",".join(PCS + ["age", "Hap"] + ["%s_x_Hap" % p for p in PCS])
CAT       = ("--catCovarList %s" % CATCOVAR) if CATCOVAR else ""

gwf = Workflow(defaults={"account": ACCOUNT})

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

gwf.target("interaction",
           inputs=["%s/step1_full_pred.list" % OUT, "%s/int_covars.txt" % OUT],
           outputs=["%s/gxhap_%s.regenie" % (OUT, PHENONAME)],
           cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 2 --bed {bfile} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {out}/int_covars.txt --covarColList {intcols} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --interaction Hap \
  --bsize 400 --minMAC 20 \
  --pred {out}/step1_full_pred.list \
  --threads 16 --out {out}/gxhap
""".format(pixi=PIXI, bfile=BFILE, out=OUT, pheno=PHENO, ph=PHENONAME,
           intcols=INT_COLS, cat=CAT)

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
gwf.target("perm_interaction",
           inputs=["%s/gxhap_%s.regenie" % (OUT, PHENONAME), "%s.bed" % BFILE],
           outputs=["%s/perm_interactions.tsv" % OUT,
                    "%s/perm_lambda.txt" % OUT,
                    "%s/perm_strata_selection.tsv" % OUT],
           cores=8, memory="16g", walltime="08:00:00") << """
{pixi} python {root}/scripts/ancestry_matched_perm.py \
  --regenie {out}/gxhap_{ph}.regenie \
  --bfile {bfile} --covar {covar} --hap {hap} \
  --pheno {pheno} --pheno-name {ph} --npc {npc} \
  --top {top} --panel {panel} --nperm {nperm} --global-nperm {gnperm} {force} \
  --seed 1 --out-prefix {out}/perm
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME, bfile=BFILE,
           covar=BASECOVAR, hap=HAPFILE, pheno=PHENO, npc=NPC,
           top=PERM_TOP, panel=PERM_PANEL, nperm=PERM_NPERM,
           gnperm=PERM_GLOBAL_NPERM, force=_FORCE)

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

    gwf.target("gwas_%s" % s,
               inputs=[keep] + pred_inputs,
               outputs=["%s/gwas_%s_%s.regenie" % (OUT, s, PHENONAME)],
               cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 2 --bed {bfile} --keep {keep} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --bsize 400 --minMAC 20 \
  --pred {pred} \
  --threads 16 --out {out}/gwas_{s}
""".format(pixi=PIXI, bfile=BFILE, keep=keep, out=OUT, pheno=PHENO,
           ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT,
           pred=pred, s=s)

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

    for nm, loc in _loci:
        gwf.target("lava_perm_%s" % nm,
                   inputs=["%s/lava_pcs.eigenvec" % OUT, "%s.bed" % BFILE],
                   outputs=["%s/lava_%s_loci.tsv" % (OUT, nm),
                            "%s/lava_%s_summary.txt" % (OUT, nm)],
                   cores=8, memory="16g", walltime="12:00:00") << """
{pixi} python {root}/scripts/local_rg_perm.py \
  --bfile {bfile} --pheno {pheno} --pheno-name {ph} --hap {hap} \
  --pcs {out}/lava_pcs.eigenvec --npc {npc} \
  --locus {loc} --locus-name {nm} \
  --n-controls {nc} --nperm {nperm} --seed 1 --out-prefix {out}/lava_{nm}
""".format(pixi=PIXI, root=ROOT, bfile=BFILE, pheno=PHENO, ph=PHENONAME,
           hap=HAPFILE, out=OUT, npc=NPC, loc=loc, nm=nm,
           nc=LAVA_N_CONTROLS, nperm=LAVA_NPERM)

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
