# Source this (don't execute it) to point workflow.py at the synthetic test data
# in tests/work/data, so `gwf` works in an interactive shell:
#
#     source tests/env.sh
#     gwf -b local status        # or: gwf -b local run
#
# gwf re-reads workflow.py every invocation and takes its config from these YS_*
# env vars; without them it falls back to the /path/to/... placeholders and errors
# with "File /path/to/genotypes.bed is required by qc but does not exist".
#
# The data must already be built -- run `python tests/run_pipeline_test.py` once
# first (it converts the VCF to plink, builds chrX, w_hm3.snplist and the LDSC ref).

_YS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
_YS_D="$_YS_ROOT/tests/work/data"

export YS_BFILE="$_YS_D/genotypes"
export YS_PHENO="$_YS_D/phenotypes.txt"
export YS_BASECOVAR="$_YS_D/covariates_base.txt"
export YS_HAPFILE="$_YS_D/haplogroup.txt"
export YS_NPC=10                       # synthetic covariates have only PC1..PC10
export YS_CATCOVAR=batch
export YS_EUR_LD="$_YS_ROOT/tests/work/eur_w_ld_chr"
export YS_HM3="$_YS_D/w_hm3.snplist"
export YS_PREV_POP=0.05
export YS_OUT="$_YS_ROOT/tests/work/gwf_results"
export YS_TMP="$_YS_ROOT/tests/work/gwf_tmp"
# bake the OpenMP workaround into specs so it reaches the gwf workers:
export YS_ENV_PREFIX="KMP_AFFINITY=disabled"
# optional arms:
export YS_XBFILE="$_YS_D/genotypesX"
export YS_GENOME_BUILD=hg38
# small counts for toy data (production defaults are huge):
export YS_LAVA_NPERM=200 YS_LAVA_N_CONTROLS=10
export YS_PERM_TOP=50 YS_PERM_NPERM=200 YS_PERM_PANEL=300 YS_PERM_GLOBAL_NPERM=50
# parallelise across tasks: the test data spans chr 1-22, so the step-2 scans fan
# out to one job per chromosome (+ gather), and the permutation targets split into
# 4 pooled batches. Drop these two to run the simple single-job path.
export YS_SPLIT_CHROMS=1-22
export YS_PERM_BATCHES=4
if [ -f "$_YS_D/truth.txt" ]; then
  export YS_LAVA_LOCI="$(awk -F'\t' '
    /^real_block_locus/{r=$2} /^conf_block_locus/{c=$2}
    END{gsub(/ /,"",r); gsub(/ /,"",c); print "realblock="r",confblock="c}' "$_YS_D/truth.txt")"
fi

echo "YS_* set -> $_YS_D  (gwf will use this data; YS_OUT=$YS_OUT)"
