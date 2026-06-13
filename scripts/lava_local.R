#!/usr/bin/env Rscript
# Headline LAVA bivariate local genetic correlation between the two strata
# (ASD-in-I vs ASD-in-R) at the supplied loci.
#
# This is the canonical-tool counterpart to scripts/local_rg_perm.py (which does
# the ancestry-matched permutation + negative-control deconfounding). LAVA is an
# R package and is NOT in the pixi env; install once with, e.g.:
#     R -e 'remotes::install_github("josefin-werme/LAVA")'
# and point RSCRIPT in workflow.py at the R that has it.
#
# Usage:
#   lava_local.R input.info sample_overlap.txt ref_prefix loci.tsv out.tsv
#
# Inputs are produced by scripts/lava_inputs.py; ref_prefix is a plink1 fileset
# (LAVA_REF, the in-sample genotypes). Output: one row per locus with the local
# r_g, its CI, and p. Treat these as the headline estimate, NOT the deconfounded
# result -- pair every locus with its lava_<name>_summary.txt permutation p.

suppressMessages(library(LAVA))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: lava_local.R input.info sample_overlap ref_prefix loci.tsv out.tsv")
}
info.file <- args[1]; overlap.file <- args[2]; ref.prefix <- args[3]
loci.file <- args[4]; out.file <- args[5]

phenos <- c("I", "R")
input <- process.input(input.info.file = info.file,
                       sample.overlap.file = overlap.file,
                       ref.prefix = ref.prefix,
                       phenos = phenos)
loci <- read.loci(loci.file)

rows <- list()
for (i in seq_len(nrow(loci))) {
  locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
  if (is.null(locus) || is.null(locus$K) || locus$K < 1) next
  ub <- tryCatch(run.univ.bivar(locus, phenos = phenos), error = function(e) NULL)
  if (is.null(ub) || is.null(ub$bivar) || nrow(ub$bivar) == 0) next
  b <- ub$bivar
  b$locus <- as.character(loci$LOC[i])
  b$chr <- locus$chr; b$start <- locus$start; b$stop <- locus$stop
  rows[[length(rows) + 1]] <- b
}

if (length(rows) == 0) {
  warning("no locus yielded a bivariate estimate (univariate signal too weak?)")
  out <- data.frame()
} else {
  out <- do.call(rbind, rows)
}
write.table(out, out.file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d local r_g estimate(s) to %s\n", nrow(out), out.file))
