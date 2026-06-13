#!/usr/bin/env python3
"""Extract the SNP x Hap interaction test from a REGENIE --interaction output,
report genomic inflation (lambda_GC) on the interaction p-values, and dump the
strongest signals.

REGENIE writes several TEST rows per variant; the interaction test is the row
labelled 'ADD-INT_SNPxHap' (verify the exact label in your output header).
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regenie", required=True)
    ap.add_argument("--test", default="ADD-INT_SNPxHap")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=200)
    a = ap.parse_args()

    df = pd.read_csv(a.regenie, sep=r"\s+", engine="python", comment="#")
    if "TEST" not in df.columns:
        raise SystemExit("no TEST column found; is this an --interaction output?")

    sub = df[df["TEST"] == a.test].copy()
    if sub.empty:
        avail = sorted(df["TEST"].unique())
        raise SystemExit("no rows for TEST=%s. available: %s" % (a.test, avail))

    sub["P"] = np.power(10.0, -sub["LOG10P"])
    chisq = stats.chi2.isf(sub["P"].clip(lower=1e-300), df=1)
    lam = float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))

    print("interaction test : %s" % a.test)
    print("variants tested  : %d" % len(sub))
    print("lambda_GC        : %.3f" % lam)
    if lam > 1.10:
        print("WARNING: inflation > 1.10 -- check PCxHap covars and QC before trusting hits")

    top = sub.sort_values("LOG10P", ascending=False).head(a.top)
    top.to_csv(a.out, sep="\t", index=False)
    print("top %d -> %s" % (a.top, a.out))


if __name__ == "__main__":
    main()
