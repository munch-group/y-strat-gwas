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
# CONFIG  -- edit these
# ---------------------------------------------------------------------------
ACCOUNT   = "your_genomedk_project"          # slurm account
ROOT      = os.path.dirname(os.path.abspath(__file__))

BFILE     = "/path/to/genotypes"             # plink1 prefix .bed/.bim/.fam (autosomes)
PHENO     = "/path/to/phenotypes.txt"        # FID IID autism   (1=case, 0=control, NA)
BASECOVAR = "/path/to/covariates_base.txt"   # FID IID PC1..PCk age batch ...
HAPFILE   = "/path/to/haplogroup.txt"        # FID IID Hap      (I or R)
PHENONAME = "autism"
NPC       = 20
CATCOVAR  = "batch"                          # comma-sep categorical covars, or ""

# LDSC reference data (download separately; see README)
LDSC_DIR  = "/path/to/ldsc"                  # cloned bulik/ldsc repo
EUR_LD    = "/path/to/eur_w_ld_chr"          # 1000G EUR LD scores (directory)
HM3       = "/path/to/w_hm3.snplist"         # HapMap3 SNP list for --merge-alleles

# liability-scale h2 (edit to your male-specific population prevalence)
PREV_POP  = 0.03

# run REGENIE step1 separately within each stratum (cleaner r_g, ~2x compute).
# False = reuse the full-sample step1 predictors (cheaper, slightly conservative
# toward r_g = 1).
STRATUM_SPECIFIC_STEP1 = True

# ---------------------------------------------------------------------------
# derived
# ---------------------------------------------------------------------------
PIXI    = "pixi run --manifest-path %s/pixi.toml" % ROOT
LDSCRUN = "conda run -n ldsc python"          # classic py2.7 LDSC env

OUT = "%s/results" % ROOT
TMP = "%s/tmp" % ROOT
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
               inputs=["%s/munged_%s.sumstats.gz" % (OUT, s)],
               outputs=["%s/h2_%s.log" % (OUT, s)],
               cores=2, memory="8g", walltime="00:30:00") << """
{ldsc} {ldscdir}/ldsc.py \
  --h2 {out}/munged_{s}.sumstats.gz \
  --ref-ld-chr {eurld}/ --w-ld-chr {eurld}/ \
  --pop-prev {pop} \
  --out {out}/h2_{s}
# NOTE: for liability-scale h2 also add  --samp-prev <cases/(cases+controls) in stratum {s}>
""".format(ldsc=LDSCRUN, ldscdir=LDSC_DIR, out=OUT, s=s, eurld=EUR_LD, pop=PREV_POP)

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
