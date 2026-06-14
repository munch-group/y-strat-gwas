#!/usr/bin/env python3
"""Build the REGENIE step-2 covariate file for the SNP x Y-haplogroup scan.

Adds a numeric Hap column (I=1, R=0) and Keller-style PC x Hap product terms.
Including PC x Hap as covariates is what keeps autosomal ancestry structure
(which is correlated with haplogroup I vs R) from leaking into the interaction
estimate. See Keller 2014, Biol Psychiatry, on confounding in GxE designs.
"""
import argparse
import sys
import pandas as pd

HAP_CODE = {"I": 1, "R": 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--covar", required=True,
                    help="base covariates: FID IID PC1..PCk age batch ...")
    ap.add_argument("--hap", required=True, help="haplogroup file (FID IID ... <hap-col>)")
    ap.add_argument("--hap-col", default="Hap", help="name of the haplogroup column")
    ap.add_argument("--npc", type=int, required=True,
                    help="number of PCs to interact with Hap")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cov = pd.read_csv(a.covar, sep=r"\s+", engine="python")
    hap = pd.read_csv(a.hap, sep=r"\s+", engine="python")

    # keep only I/R (others -- different lineages / missing -- are dropped)
    hap["Hap"] = hap[a.hap_col].astype(str).str.upper().str.strip().map(HAP_CODE)
    hap = hap.dropna(subset=["Hap"])
    hap["Hap"] = hap["Hap"].astype(int)

    pcs = ["PC%d" % i for i in range(1, a.npc + 1)]
    missing = [c for c in pcs if c not in cov.columns]
    if missing:
        sys.exit("missing PC columns in covar file: %s" % missing)

    m = cov.merge(hap[["FID", "IID", "Hap"]], on=["FID", "IID"], how="inner")
    for pc in pcs:
        m["%s_x_Hap" % pc] = m[pc] * m["Hap"]

    m.to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print("wrote %s: %d samples, %d columns" % (a.out, m.shape[0], m.shape[1]))


if __name__ == "__main__":
    main()
