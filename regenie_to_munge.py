#!/usr/bin/env python3
"""Convert a REGENIE step-2 main-effect file into an LDSC-munge-ready table.

REGENIE reports LOG10P (not P) and uses ALLELE1 as the effect allele.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regenie", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.regenie, sep=r"\s+", engine="python", comment="#")
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()

    df["P"] = np.power(10.0, -df["LOG10P"])
    out = pd.DataFrame({
        "SNP":  df["ID"],
        "A1":   df["ALLELE1"],   # effect allele in REGENIE
        "A2":   df["ALLELE0"],
        "N":    df["N"],
        "BETA": df["BETA"],
        "SE":   df["SE"],
        "P":    df["P"],
    })
    out.to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print("wrote %s: %d variants" % (a.out, len(out)))


if __name__ == "__main__":
    main()
