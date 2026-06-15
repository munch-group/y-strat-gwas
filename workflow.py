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

ACCOUNT   = _env("ChrXh2", None)   # slurm account
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
# To run on the real iPSYCH cohort, keep these synthetic defaults and uncomment the
# "REAL DATA" override block at the end of this CONFIG section -- it sets every
# real path/parameter in one place (and routes the raw pheno/covar through the
# prep tasks). The real file layouts: Hap file has the lineage in col 6 named
# "Major" and includes non-I/R males (dropped); the PC file is a plink MDS cov
# "FID IID SOL C1..C10 st1" (no age/batch -- prep_covar renames C->PC and adds
# dummies); the .pheno files are headerless "FID IID value", 2=case/1=control.


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
LAVA_PARTITION  = _env("LAVA_PARTITION", "lava_meta/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")          # LAVA blocks file (R headline only)
LAVA_N_CONTROLS = int(_env("LAVA_N_CONTROLS", "300")) # negative-control loci
LAVA_NPERM      = int(_env("LAVA_NPERM", "10000"))    # permutations per locus
RSCRIPT         = _env("RSCRIPT", "Rscript")          # for the LAVA R step

# chrX (optional): plink1 prefix with X variants coded as chromosome 23 and the
# SEX column filled in the .fam (1=male). Empty string disables all X targets.
# REGENIE codes non-PAR males as hemizygous (0/2) automatically from sex + build.
XBFILE       = _env("XBFILE", "")  # "" => skip chrX; synthetic X is tests/work/data/genotypesX
GENOME_BUILD = _env("GENOME_BUILD", "hg19")               # for REGENIE --par-region

# LDSC reference data (download separately; see README). Defaults point at the
# bundled synthetic reference -- swap EUR_LD/HM3 for the real 1000G EUR / HapMap3.
LDSC_DIR  = _env("LDSC_DIR", "%s/ldsc" % ROOT)             # py3 LDSC fork, vendored + pip-installed by pixi
# EUR_LD    = _env("EUR_LD", "%s/tests/work/eur_w_ld_chr" % ROOT)   # EUR LD scores dir (synthetic ref)
EUR_LD    = _env("EUR_LD", "%s/ldsc/eur_w_ld_chr" % ROOT)   # EUR LD scores dir (synthetic ref) DOWNLOADED FROM https://zenodo.org/records/8182036
# HM3       = _env("HM3", "%s/tests/work/data/w_hm3.snplist" % ROOT)  # HapMap3 SNP list for --merge-alleles
HM3       = _env("HM3", "%s/ldsc/w_hm3.snplist" % ROOT)  # HapMap3 SNP list for --merge-alleles DOWNLOADED FROM https://zenodo.org/records/7773502

# liability-scale h2 (edit to your male-specific population prevalence)
PREV_POP  = float(_env("PREV_POP", "0.03"))

# run REGENIE step1 separately within each stratum (cleaner r_g, ~2x compute).
# False = reuse the full-sample step1 predictors (cheaper, slightly conservative
# toward r_g = 1).
STRATUM_SPECIFIC_STEP1 = _env("STRATUM_SPECIFIC_STEP1", "True").lower() in ("1", "true", "yes")

# Heritability extras.
#  H2_POOLED  : also run a pooled (non-stratified) full-sample GWAS -> h2_full,
#               and tabulate it against the per-stratum h2_{I,R} (h2_by_stratification.tsv).
#  H2_PER_CHR : per-chromosome h2 from the pooled GWAS (h2_by_chromosome.tsv), to
#               weigh each chromosome's contribution; chrX is included only if an
#               X-chromosome LD reference (EUR_LD_X + HM3_X) is supplied.
def _flag(name, default):
    return _env(name, default).lower() in ("1", "true", "yes")
H2_POOLED  = _flag("H2_POOLED", "True")
H2_PER_CHR = _flag("H2_PER_CHR", "True")
EUR_LD_X   = _env("EUR_LD_X", "")     # chrX LD-score reference dir/prefix (for chrX h2)
HM3_X      = _env("HM3_X", "")        # X-inclusive --merge-alleles list (for chrX munge)

# Females-as-negative-control. Females carry no Y, so a genuine Y-haplogroup x
# autosome interaction cannot exist in them; re-testing the male interaction hits
# for SNP x pseudo-Hap in females (pseudo-Hap = the ancestry-matched split,
# assign_pseudo_hap.py) exposes ancestry artefacts. Like the male inputs, these
# default to the bundled SYNTHETIC FEMALE DATA so the arm runs out of the box --
# SWAP for your real female files (set FBFILE="" to skip the arm entirely). The
# female PCs must live in the SAME space as the male PCs (FBASECOVAR), e.g. PCs
# computed jointly across sexes or females projected onto the male PCA.
FBFILE     = _env("FBFILE", "%s/tests/work/data/females" % ROOT)                    # female plink1 prefix
FPHENO     = _env("FPHENO", "%s/tests/work/data/female_phenotypes.txt" % ROOT)      # female FID IID <pheno>
FBASECOVAR = _env("FBASECOVAR", "%s/tests/work/data/female_covariates_base.txt" % ROOT)  # female PCs (male space)
FNEG_FORCE_SNPS = _env("FNEG_FORCE_SNPS", "")  # comma-sep SNP IDs to always re-test in females

# --- Real-data input adapters -------------------------------------------------
# Haplogroup column name (real files may call it e.g. "Major"; rows whose value is
# not I or R are DROPPED from the whole analysis -- a keep_IR list restricts the
# full-sample targets, the per-stratum keep lists restrict the rest).
HAPCOL = _env("HAPCOL", "Hap")

# Covariate prep. If RAW_COVAR is set, a prep_covar task renames the PC columns from
# PC_PREFIX (e.g. "C") to PC1..PCk and adds any model covariates the raw file lacks
# (the EXTRA_COVAR quantitative cols + the CATCOVAR categorical cols) as DUMMY
# placeholders, writing the BASECOVAR the rest of the pipeline consumes. Same for
# RAW_FCOVAR -> FBASECOVAR (females). Leave empty to use BASECOVAR/FBASECOVAR as-is.
RAW_COVAR   = _env("RAW_COVAR", "")      # raw male covar (FID IID ... <PC_PREFIX>1..k ...)
RAW_FCOVAR  = _env("RAW_FCOVAR", "")     # raw female covar
PC_PREFIX   = _env("PC_PREFIX", "PC")    # PC column prefix in the raw file (real: "C")
EXTRA_COVAR = _env("EXTRA_COVAR", "age")  # quantitative covars beyond PCs (comma-sep, or "")

# Phenotype recode. REGENIE --bt wants control=0, case=1; real files are often in
# PLINK 1/2 coding. If RAW_PHENO is set, a prep_pheno task maps PHENO_CASE->1 /
# PHENO_CONTROL->0 (else NA), writing the PHENO the pipeline consumes. Same for
# RAW_FPHENO -> FPHENO (females). Leave empty to use PHENO/FPHENO as-is.
RAW_PHENO    = _env("RAW_PHENO", "")     # raw male pheno (FID IID ... <pheno> ...)
RAW_FPHENO   = _env("RAW_FPHENO", "")    # raw female pheno
PHENO_RAW_COL = _env("PHENO_RAW_COL", PHENONAME)  # phenotype column (header case)
PHENO_CASE    = _env("PHENO_CASE", "2")  # raw value meaning case
PHENO_CONTROL = _env("PHENO_CONTROL", "1")  # raw value meaning control
# Real PLINK .pheno files are usually HEADERLESS ("FID IID value"). Keep
# PHENO_HAS_HEADER False to read positionally (phenotype value = the 1-based
# PHENO_VALUE_COL, default 3). Set True if the raw pheno has a header naming the
# column (then PHENO_RAW_COL selects it).
PHENO_HAS_HEADER = _env("PHENO_HAS_HEADER", "False").lower() in ("1", "true", "yes")
PHENO_VALUE_COL  = int(_env("PHENO_VALUE_COL", "3"))

# =============================================================================
# REAL DATA  -- uncomment this whole block to run on the iPSYCH cohort instead of
# the synthetic test data. It overrides the relevant defaults above in one place;
# paths are from file_paths.txt. The raw pheno/covar go through RAW_PHENO /
# RAW_COVAR (the prep tasks fix their format) rather than PHENO / BASECOVAR.
# Each value can also be set from the shell via its YS_* env var instead.
# =============================================================================
ACCOUNT    = "ChrXh2"          # your slurm account (was None -> no account)
ENV_PREFIX = ""                # clear the synthetic-data sandbox KMP_AFFINITY flag
_DATA = "/faststorage/jail/project/ChrXh2/data"
_HAP  = "/faststorage/jail/project/ChrXh2/shannon/y_haplo/results"

# --- autosomes (males; non-I/R males dropped everywhere via keep_IR) ---------
# NB: plink/REGENIE read <prefix>.bed/.bim/.fam for this prefix. If your files
# are named .bam/.bin, symlink them to .bed/.bim first.
BFILE   = "%s/autosomes/iPSYCH2015_HRC_2020-merge.hg19.ch.fl.bgn" % _DATA
HAPFILE = "%s/haplogroup_combined_haplogroup_assignments.txt" % _HAP
HAPCOL  = "Major"              # haplogroup in col 6, named "Major"; FID/IID in cols 1-2
NPC     = 10

# --- phenotype: headerless PLINK .pheno, value in col 3, 2=case / 1=control --
RAW_PHENO     = "%s/pheno/asdGWAS2015.pheno" % _DATA
PHENO_CASE    = "2"
PHENO_CONTROL = "1"            # PHENO_HAS_HEADER stays False (value = PHENO_VALUE_COL=3)

# --- covariates: plink MDS cov "FID IID SOL C1..C10 st1" (has a header) -------
RAW_COVAR   = "%s/pca/ASD2015_pca2_ell_dk-dim3_8_8_8.menv.mds_cov" % _DATA
PC_PREFIX   = "C"              # rename C1..C10 -> PC1..PC10
EXTRA_COVAR = "age"            # not in the raw file -> added as a dummy (or "" to drop)
CATCOVAR    = "batch"          # not in the raw file -> added as a dummy (or "" to drop)

# # --- chrX (hemizygous males, build hg19) -------------------------------------
XBFILE       = "%s/chrX/M/iPSYCH2015_HRC_ChrX_M_2020-merge.hg19.ch.fl.bgn" % _DATA
GENOME_BUILD = "hg19"

# # --- females negative control (set FBFILE="" to skip the arm) ----------------
# # FBFILE may be the shared autosome fileset: females are picked out by the
# # FEMALE-ONLY pheno/covar below (males in the fileset fall out of the sample
# # intersection). RAW_FCOVAR must be female-only AND in the SAME PC space as the
# # males (restrict the joint .mds_cov to females -- it already is the joint
# # autosomal PCA). RAW_FPHENO is the headerless female .pheno (same 2/1 coding).
# FBFILE     = BFILE
FBFILE     = ""
# RAW_FPHENO = "%s/pheno/asdGWAS2015_females.pheno" % _DATA
# RAW_FCOVAR = "<female-only .mds_cov in the male PC space>"   # restrict the joint mds_cov to females

# ---------------------------------------------------------------------------
# derived
# ---------------------------------------------------------------------------
# LAVA_REF defaults to BFILE, but BFILE may have been reassigned by the REAL DATA
# block above *after* the LAVA_REF default was first bound -- re-bind it to the
# final BFILE unless the user set it explicitly (YS_LAVA_REF).
if not os.environ.get("YS_LAVA_REF"):
    LAVA_REF = BFILE

PIXI    = ("%s " % ENV_PREFIX if ENV_PREFIX else "") \
          + "pixi run --manifest-path %s/pixi.toml" % ROOT
LDSCRUN = "%s python" % PIXI                   # py3 LDSC fork, pip-installed by pixi postinstall

OUT = _env("OUT", "%s/results" % ROOT)
TMP = _env("TMP", "%s/tmp" % ROOT)
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

PCS       = ["PC%d" % i for i in range(1, NPC + 1)]
_EXTRA    = [c.strip() for c in EXTRA_COVAR.split(",") if c.strip()]
BASE_COLS = ",".join(PCS + _EXTRA)
INT_COLS  = ",".join(PCS + _EXTRA + ["Hap"] + ["%s_x_Hap" % p for p in PCS])
CAT       = ("--catCovarList %s" % CATCOVAR) if CATCOVAR else ""
KEEP_IR   = "%s/keep_IR.txt" % OUT

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
# Covariate prep (only when RAW_COVAR / RAW_FCOVAR is set): rename PCs + add the
# dummy covariate columns the model needs but the raw file lacks. The prepped
# files become BASECOVAR / FBASECOVAR for everything downstream.
# ---------------------------------------------------------------------------
def _prep_covar_cmd(raw, out):
    return """
{pixi} python {root}/scripts/prep_covar.py --raw-covar {raw} --npc {npc} \
  --pc-prefix {pcp} --add-quant "{quant}" --add-cat "{cat}" --out {out}
""".format(pixi=PIXI, root=ROOT, raw=raw, npc=NPC, pcp=PC_PREFIX,
           quant=EXTRA_COVAR, cat=CATCOVAR, out=out)

if RAW_COVAR:
    BASECOVAR = "%s/base_covar.txt" % OUT
    gwf.target("prep_covar", inputs=[RAW_COVAR], outputs=[BASECOVAR],
               cores=1, memory="4g", walltime="00:10:00") \
        << _prep_covar_cmd(RAW_COVAR, BASECOVAR)
if RAW_FCOVAR:
    FBASECOVAR = "%s/female_base_covar.txt" % OUT
    gwf.target("prep_fcovar", inputs=[RAW_FCOVAR], outputs=[FBASECOVAR],
               cores=1, memory="4g", walltime="00:10:00") \
        << _prep_covar_cmd(RAW_FCOVAR, FBASECOVAR)

def _prep_pheno_cmd(raw, out):
    sel = ("--raw-col %s" % PHENO_RAW_COL if PHENO_HAS_HEADER
           else "--no-header --value-col %d" % PHENO_VALUE_COL)
    return """
{pixi} python {root}/scripts/prep_pheno.py --raw-pheno {raw} \
  --pheno-name {ph} {sel} --case {case} --control {ctrl} --out {out}
""".format(pixi=PIXI, root=ROOT, raw=raw, ph=PHENONAME, sel=sel,
           case=PHENO_CASE, ctrl=PHENO_CONTROL, out=out)

if RAW_PHENO:
    PHENO = "%s/pheno_recoded.txt" % OUT
    gwf.target("prep_pheno", inputs=[RAW_PHENO], outputs=[PHENO],
               cores=1, memory="4g", walltime="00:10:00") \
        << _prep_pheno_cmd(RAW_PHENO, PHENO)
if RAW_FPHENO:
    FPHENO = "%s/female_pheno_recoded.txt" % OUT
    gwf.target("prep_fpheno", inputs=[RAW_FPHENO], outputs=[FPHENO],
               cores=1, memory="4g", walltime="00:10:00") \
        << _prep_pheno_cmd(RAW_FPHENO, FPHENO)

# ---------------------------------------------------------------------------
# QC backbone for REGENIE step 1
# ---------------------------------------------------------------------------
gwf.target("qc",
           inputs=["%s.bed" % BFILE, KEEP_IR],
           outputs=["%s/qc_pass.snplist" % OUT],
           cores=8, memory="16g", walltime="02:00:00") << """
{pixi} plink2 --bfile {bfile} --keep {keepir} \
  --maf 0.01 --mac 100 --geno 0.05 --hwe 1e-15 \
  --indep-pairwise 1000 100 0.9 \
  --out {out}/qc
mv {out}/qc.prune.in {out}/qc_pass.snplist
""".format(pixi=PIXI, bfile=BFILE, keepir=KEEP_IR, out=OUT)

# ---------------------------------------------------------------------------
# prep: strata keep-lists + interaction covariate file
# ---------------------------------------------------------------------------
gwf.target("strata",
           inputs=[HAPFILE],
           outputs=["%s/keep_I.txt" % OUT, "%s/keep_R.txt" % OUT, KEEP_IR],
           cores=1, memory="2g", walltime="00:10:00") << """
{pixi} python {root}/scripts/make_strata.py --hap {hap} --hap-col {hapcol} --out-prefix {out}
""".format(pixi=PIXI, root=ROOT, hap=HAPFILE, hapcol=HAPCOL, out=OUT)

gwf.target("int_covars",
           inputs=[BASECOVAR, HAPFILE],
           outputs=["%s/int_covars.txt" % OUT],
           cores=1, memory="4g", walltime="00:10:00") << """
{pixi} python {root}/scripts/make_interaction_covars.py \
  --covar {covar} --hap {hap} --hap-col {hapcol} --npc {npc} --out {out}/int_covars.txt
""".format(pixi=PIXI, root=ROOT, covar=BASECOVAR, hap=HAPFILE, hapcol=HAPCOL,
           npc=NPC, out=OUT)

# ---------------------------------------------------------------------------
# ARM A: full-sample step1 + interaction scan
# ---------------------------------------------------------------------------
gwf.target("step1_full",
           inputs=["%s/qc_pass.snplist" % OUT, KEEP_IR, PHENO, BASECOVAR],
           outputs=["%s/step1_full_pred.list" % OUT],
           cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 1 --bed {bfile} --keep {keepir} \
  --extract {out}/qc_pass.snplist \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --force-step1 --bsize 1000 --lowmem --lowmem-prefix {tmp}/step1_full \
  --threads 16 --out {out}/step1_full
""".format(pixi=PIXI, bfile=BFILE, keepir=KEEP_IR, out=OUT, pheno=PHENO,
           ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT, tmp=TMP)

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
           cores=2, memory="64g", walltime="02:00:00") << """
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
  --bfile {bfile} --covar {covar} --hap {hap} --hap-col-name {hap_col_name} \
  --pheno {pheno} --pheno-name {ph} --npc {npc} \
  --top {top} --panel {panel} --nperm {nperm} --global-nperm {gnperm} {force} \
  --seed {seed} --select-seed 1 {extra} --out-prefix {outpref}
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME, bfile=BFILE,
           covar=BASECOVAR, hap=HAPFILE, hap_col_name=HAPCOL, pheno=PHENO, npc=NPC, top=PERM_TOP,
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
                   outputs=[_cp, _lp], cores=8, memory="64g", walltime="08:00:00") \
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
                   inputs=[keep, "%s/qc_pass.snplist" % OUT, PHENO, BASECOVAR],
                   outputs=["%s/step1_%s_pred.list" % (OUT, s)],
                   cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 1 --bed {bfile} --keep {keep} \
  --extract {out}/qc_pass.snplist \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --force-step1 --bsize 1000 --lowmem --lowmem-prefix {tmp}/step1_{s} \
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
               cores=2, memory="24g", walltime="01:00:00") << """
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
# Heritability extras: a pooled (non-stratified) full-sample GWAS for the
# stratified-vs-unstratified h2 comparison, and per-chromosome h2 to weigh each
# chromosome (incl. chrX, given an X LD reference) against the others.
# ---------------------------------------------------------------------------
_X_H2 = bool(XBFILE and EUR_LD_X and HM3_X)
if H2_POOLED or H2_PER_CHR:
    def _gwas_full_step2(chr_flag, out_prefix):
        return """
{pixi} regenie --step 2 --bed {bfile} --keep {keepir} {chrflag} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --bsize 400 --minMAC 20 \
  --pred {out}/step1_full_pred.list \
  --threads 16 --out {outpref}
""".format(pixi=PIXI, bfile=BFILE, keepir=KEEP_IR, chrflag=chr_flag, out=OUT,
           pheno=PHENO, ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT,
           outpref=out_prefix)

    _GWF = "%s/gwas_full_%s.regenie" % (OUT, PHENONAME)
    if SPLIT_CHROMS:
        _fchunks = []
        for _c in _chroms(SPLIT_CHROMS):
            _fc = "%s/gwas_full_chr%d_%s.regenie" % (OUT, _c, PHENONAME)
            gwf.target("gwas_full_chr%d" % _c,
                       inputs=["%s/step1_full_pred.list" % OUT, KEEP_IR], outputs=[_fc],
                       cores=16, memory="32g", walltime="24:00:00") \
                << _gwas_full_step2("--chr %d" % _c, "%s/gwas_full_chr%d" % (OUT, _c))
            _fchunks.append(_fc)
        gwf.target("gwas_full", inputs=_fchunks, outputs=[_GWF],
                   cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/concat_regenie.py --out {gw} --inputs {chunks}
""".format(pixi=PIXI, root=ROOT, gw=_GWF, chunks=" ".join(_fchunks))
    else:
        gwf.target("gwas_full", inputs=["%s/step1_full_pred.list" % OUT, KEEP_IR],
                   outputs=[_GWF], cores=16, memory="32g", walltime="24:00:00") \
            << _gwas_full_step2("", "%s/gwas_full" % OUT)

    gwf.target("munge_full",
               inputs=[_GWF], outputs=["%s/munged_full.sumstats.gz" % OUT],
               cores=2, memory="24g", walltime="05:00:00") << """
{pixi} python {root}/scripts/regenie_to_munge.py \
  --regenie {gw} --out {out}/gwas_full.forldsc.txt
{ldsc} {ldscdir}/munge_sumstats.py \
  --sumstats {out}/gwas_full.forldsc.txt --merge-alleles {hm3} \
  --snp SNP --a1 A1 --a2 A2 --N-col N --p P --signed-sumstats BETA,0 \
  --chunksize 500000 --out {out}/munged_full
""".format(pixi=PIXI, root=ROOT, gw=_GWF, out=OUT, ldsc=LDSCRUN,
           ldscdir=LDSC_DIR, hm3=HM3)

if H2_POOLED:
    gwf.target("h2_full",
               inputs=["%s/munged_full.sumstats.gz" % OUT],
               outputs=["%s/h2_full.log" % OUT],
               cores=2, memory="24g", walltime="02:00:00") << """
SAMP=$({pixi} python {root}/scripts/samp_prev.py --pheno {pheno} --pheno-name {ph})
{ldsc} {ldscdir}/ldsc.py \
  --h2 {out}/munged_full.sumstats.gz \
  --ref-ld-chr {eurld}/ --w-ld-chr {eurld}/ \
  --samp-prev $SAMP --pop-prev {pop} \
  --out {out}/h2_full
""".format(pixi=PIXI, root=ROOT, pheno=PHENO, ph=PHENONAME, ldsc=LDSCRUN,
           ldscdir=LDSC_DIR, out=OUT, eurld=EUR_LD, pop=PREV_POP)

    gwf.target("h2_by_stratification",
               inputs=["%s/h2_full.log" % OUT, "%s/h2_I.log" % OUT, "%s/h2_R.log" % OUT],
               outputs=["%s/h2_by_stratification.tsv" % OUT],
               cores=1, memory="25g", walltime="2:00:00") << """
{pixi} python {root}/scripts/collect_h2.py \
  --logs {out}/h2_full.log {out}/h2_I.log {out}/h2_R.log \
  --labels pooled I R --label-col stratification \
  --out {out}/h2_by_stratification.tsv
""".format(pixi=PIXI, root=ROOT, out=OUT)

if H2_PER_CHR:
    _aut = _chroms(SPLIT_CHROMS) if SPLIT_CHROMS else list(range(1, 23))
    _per_logs = ["%s/h2_chr%d.log" % (OUT, c) for c in _aut]
    _per_labels = [str(c) for c in _aut]
    _x_step = _x_logs = ""
    _x_inputs = []
    if _X_H2:
        _per_logs.append("%s/h2_chrX.log" % OUT)
        _per_labels.append("X")
        _x_inputs = ["%s/munged_X_full.sumstats.gz" % OUT]
        _x_step = ("{ldsc} {ldscdir}/ldsc.py --h2 {out}/munged_X_full.sumstats.gz "
                   "--ref-ld-chr {xld}/ --w-ld-chr {xld}/ --samp-prev $SAMP "
                   "--pop-prev {pop} --out {out}/h2_chrX\n").format(
                       ldsc=LDSCRUN, ldscdir=LDSC_DIR, out=OUT, xld=EUR_LD_X, pop=PREV_POP)
    gwf.target("h2_by_chromosome",
               inputs=["%s/munged_full.sumstats.gz" % OUT] + _x_inputs,
               outputs=["%s/h2_by_chromosome.tsv" % OUT],
               cores=2, memory="8g", walltime="02:00:00") << """
SAMP=$({pixi} python {root}/scripts/samp_prev.py --pheno {pheno} --pheno-name {ph})
for c in {chroms}; do
  {ldsc} {ldscdir}/ldsc.py --h2 {out}/munged_full.sumstats.gz \
    --ref-ld {eurld}/$c --w-ld {eurld}/$c \
    --samp-prev $SAMP --pop-prev {pop} --out {out}/h2_chr$c
done
{xstep}{pixi} python {root}/scripts/collect_h2.py \
  --logs {logs} --labels {labels} --label-col chrom \
  --out {out}/h2_by_chromosome.tsv
""".format(pixi=PIXI, root=ROOT, pheno=PHENO, ph=PHENONAME, ldsc=LDSCRUN,
           ldscdir=LDSC_DIR, out=OUT, eurld=EUR_LD, pop=PREV_POP,
           chroms=" ".join(str(c) for c in _aut), xstep=_x_step,
           logs=" ".join(_per_logs), labels=" ".join(_per_labels))

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

    # pooled (all-male) chrX GWAS -> munged, so chrX can enter h2_by_chromosome.
    # Only when an X LD reference + X merge-alleles list are supplied (_X_H2).
    if _X_H2:
        gwf.target("gwas_X_full",
                   inputs=["%s/step1_full_pred.list" % OUT,
                           "%s/x_qc_pass.snplist" % OUT, KEEP_IR],
                   outputs=["%s/gwas_X_full_%s.regenie" % (OUT, PHENONAME)],
                   cores=16, memory="32g", walltime="24:00:00") << """
{pixi} regenie --step 2 --bed {xbfile} --keep {keepir} \
  --chr 23 --extract {out}/x_qc_pass.snplist {xhemi} \
  --phenoFile {pheno} --phenoColList {ph} \
  --covarFile {covar} --covarColList {base} {cat} \
  --bt --firth --approx --pThresh 0.01 \
  --bsize 400 --minMAC 20 \
  --pred {out}/step1_full_pred.list \
  --threads 16 --out {out}/gwas_X_full
""".format(pixi=PIXI, xbfile=XBFILE, xhemi=XHEMI, keepir=KEEP_IR, out=OUT,
           pheno=PHENO, ph=PHENONAME, covar=BASECOVAR, base=BASE_COLS, cat=CAT)

        gwf.target("munge_X_full",
                   inputs=["%s/gwas_X_full_%s.regenie" % (OUT, PHENONAME)],
                   outputs=["%s/munged_X_full.sumstats.gz" % OUT],
                   cores=2, memory="8g", walltime="01:00:00") << """
{pixi} python {root}/scripts/regenie_to_munge.py \
  --regenie {out}/gwas_X_full_{ph}.regenie --out {out}/gwas_X_full.forldsc.txt
{ldsc} {ldscdir}/munge_sumstats.py \
  --sumstats {out}/gwas_X_full.forldsc.txt --merge-alleles {hm3x} \
  --snp SNP --a1 A1 --a2 A2 --N-col N --p P --signed-sumstats BETA,0 \
  --chunksize 500000 --out {out}/munged_X_full
""".format(pixi=PIXI, root=ROOT, out=OUT, ph=PHENONAME, ldsc=LDSCRUN,
           ldscdir=LDSC_DIR, hm3x=HM3_X)

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
        _li = ["%s/lava_pcs.eigenvec" % OUT, "%s.bed" % BFILE, PHENO]
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

# ---------------------------------------------------------------------------
# Headline LAVA local r_g via the real LAVA R package (runs OUTSIDE the pixi env,
# through RSCRIPT). Gated on LAVA_PARTITION: the genome-wide LD-block file
# (LAVA's blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile -- columns
# "LOC CHR START STOP", positions in the GENOTYPE BUILD, hg19 here) is scanned
# directly, so this is a genome-wide screen INDEPENDENT of LAVA_LOCI. Follow up
# any block of interest with its ancestry-matched deconfounding by listing that
# block in LAVA_LOCI (-> lava_perm_<name>). Consumes Arm B's per-stratum
# sumstats, so it runs after munge_{I,R}.
# ---------------------------------------------------------------------------
if LAVA_PARTITION:
    gwf.target("lava_inputs",
               inputs=["%s/munged_I.sumstats.gz" % OUT,
                       "%s/munged_R.sumstats.gz" % OUT],
               outputs=["%s/lava_input.info" % OUT,
                        "%s/lava_sample_overlap.txt" % OUT],
               cores=1, memory="4g", walltime="00:30:00") << """
{pixi} python {root}/scripts/lava_inputs.py \
  --pheno {pheno} --pheno-name {ph} \
  --keep-i {out}/keep_I.txt --keep-r {out}/keep_R.txt \
  --sumstats-i {out}/gwas_I.forldsc.txt --sumstats-r {out}/gwas_R.forldsc.txt \
  --out-prefix {out}/lava
""".format(pixi=PIXI, root=ROOT, pheno=PHENO, ph=PHENONAME, out=OUT)

    gwf.target("lava_local",
               inputs=["%s/lava_input.info" % OUT, LAVA_PARTITION],
               outputs=["%s/lava_local_rg.tsv" % OUT],
               cores=4, memory="32g", walltime="24:00:00") << """
{rscript} {root}/scripts/lava_local.R \
  {out}/lava_input.info {out}/lava_sample_overlap.txt {ref} \
  {partition} {out}/lava_local_rg.tsv
""".format(rscript=RSCRIPT, root=ROOT, out=OUT, ref=LAVA_REF,
           partition=LAVA_PARTITION)

# ---------------------------------------------------------------------------
# Females-as-negative-control (optional; only built when FBFILE is set). Females
# have no Y, so a real Y-haplogroup x autosome interaction can't exist in them.
# We split females along the male I-vs-R ancestry propensity (pseudo-Hap) and
# re-test the male interaction hits for SNP x pseudo-Hap in females: a hit that
# reproduces in females is an ancestry artefact; a real Y-driven hit stays null.
# See METHODS.md; scripts/assign_pseudo_hap.py + scripts/female_negcontrol.py.
# ---------------------------------------------------------------------------
if FBFILE:
    gwf.target("female_pseudohap",
               inputs=[HAPFILE, BASECOVAR, FBASECOVAR, "%s.bed" % FBFILE],
               outputs=["%s/female_pseudohap.txt" % OUT],
               cores=1, memory="4g", walltime="00:10:00") << """
{pixi} python {root}/scripts/assign_pseudo_hap.py \
  --male-hap {hap} --hap-col {hapcol} --male-covar {covar} --female-covar {fcovar} \
  --npc {npc} --out {out}/female_pseudohap.txt
""".format(pixi=PIXI, root=ROOT, hap=HAPFILE, hapcol=HAPCOL, covar=BASECOVAR,
           fcovar=FBASECOVAR, npc=NPC, out=OUT)

    gwf.target("female_negcontrol",
               inputs=["%s/top_interactions.tsv" % OUT,
                       "%s/female_pseudohap.txt" % OUT, "%s.bed" % FBFILE, FPHENO],
               outputs=["%s/female_negative_control.tsv" % OUT,
                        "%s/female_lambda.txt" % OUT],
               cores=8, memory="16g", walltime="04:00:00") << """
{pixi} python {root}/scripts/female_negcontrol.py \
  --male-interactions {out}/top_interactions.tsv \
  --fbfile {fbfile} --fpheno {fpheno} --pheno-name {ph} \
  --fcovar {fcovar} --female-hap {out}/female_pseudohap.txt --npc {npc} {force} \
  --seed 1 --out-prefix {out}/female
""".format(pixi=PIXI, root=ROOT, out=OUT, fbfile=FBFILE, fpheno=FPHENO,
           ph=PHENONAME, fcovar=FBASECOVAR, npc=NPC,
           force=("--force-snps %s" % FNEG_FORCE_SNPS) if FNEG_FORCE_SNPS else "")
