#!/usr/bin/env python3
"""Print the case fraction (sample prevalence) of a binary phenotype, optionally
within a --keep subset. Used to supply LDSC's --samp-prev for liability-scale h2,
which LDSC requires alongside --pop-prev (it errors if only one is given).
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True, help="FID IID <pheno> (1/0/NA)")
    ap.add_argument("--pheno-name", required=True)
    ap.add_argument("--keep", help="FID IID keep-list (no header), optional")
    a = ap.parse_args()

    ph = pd.read_csv(a.pheno, sep=r"\s+", engine="python")
    if a.keep:
        keep = pd.read_csv(a.keep, sep=r"\s+", engine="python",
                           header=None, names=["FID", "IID"])
        ph = ph.merge(keep, on=["FID", "IID"], how="inner")

    y = pd.to_numeric(ph[a.pheno_name], errors="coerce")
    y = y[y.isin([0, 1])]
    if len(y) == 0:
        raise SystemExit("no 0/1 phenotype values found")
    print("%.6f" % (float((y == 1).sum()) / len(y)))


if __name__ == "__main__":
    main()
