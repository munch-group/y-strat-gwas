#!/usr/bin/env python3
"""Build the input files LAVA needs for a bivariate local-r_g run between the two
strata (ASD-in-I vs ASD-in-R): the input.info table, a sample-overlap matrix, and
a locus-definition file.

The I and R strata are DISJOINT individuals, so their sample overlap is zero
(off-diagonal 0, diagonal 1) -- no overlap correction needed. Per-stratum case
and control counts are read from the phenotype file restricted to each keep-list.
"""
import argparse
import pandas as pd


def counts(pheno, name, keep):
    ph = pd.read_csv(pheno, sep=r"\s+", engine="python")
    k = pd.read_csv(keep, sep=r"\s+", engine="python", header=None, names=["FID", "IID"])
    ph = ph.merge(k, on=["FID", "IID"], how="inner")
    y = pd.to_numeric(ph[name], errors="coerce")
    return int((y == 1).sum()), int((y == 0).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-name", required=True)
    ap.add_argument("--keep-i", required=True)
    ap.add_argument("--keep-r", required=True)
    ap.add_argument("--sumstats-i", required=True)
    ap.add_argument("--sumstats-r", required=True)
    ap.add_argument("--loci", default="",
                    help="name=chr:start-end,name2=...; omit for a genome-wide "
                         "scan that reads a ready-made LAVA locus/partition file")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    ca_i, co_i = counts(a.pheno, a.pheno_name, a.keep_i)
    ca_r, co_r = counts(a.pheno, a.pheno_name, a.keep_r)

    with open(a.out_prefix + "_input.info", "w") as f:
        f.write("phenotype\tcases\tcontrols\tfilename\n")
        f.write("I\t%d\t%d\t%s\n" % (ca_i, co_i, a.sumstats_i))
        f.write("R\t%d\t%d\t%s\n" % (ca_r, co_r, a.sumstats_r))

    # disjoint strata -> zero sample overlap
    with open(a.out_prefix + "_sample_overlap.txt", "w") as f:
        f.write("I R\n1 0\n0 1\n")

    # Only build a locus file from an inline --loci string. For a genome-wide
    # scan the caller passes LAVA's ready-made partition file straight to the R
    # step instead, so no _loci.tsv is needed here.
    if a.loci:
        with open(a.out_prefix + "_loci.tsv", "w") as f:
            f.write("LOC\tCHR\tSTART\tSTOP\n")
            for ent in a.loci.split(","):
                if "=" not in ent:
                    continue
                nm, loc = ent.split("=", 1)
                chrom, rng = loc.split(":")
                lo, hi = rng.split("-")
                f.write("%s\t%s\t%s\t%s\n" % (nm.strip(), chrom.strip(),
                                              lo.strip(), hi.strip()))
        print("wrote LAVA inputs: %s_{input.info,sample_overlap.txt,loci.tsv}"
              % a.out_prefix)
    else:
        print("wrote LAVA inputs: %s_{input.info,sample_overlap.txt}"
              % a.out_prefix)


if __name__ == "__main__":
    main()
