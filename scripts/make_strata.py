#!/usr/bin/env python3
"""Map haplogroup labels (I/R) to per-stratum REGENIE --keep lists."""
import argparse
import sys
import pandas as pd

HAP_CODE = {"I": 1, "R": 0}  # I=1, R=0 (consistent with interaction coding)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hap", required=True, help="FID IID Hap  (Hap in {I,R})")
    ap.add_argument("--out-prefix", default=".")
    a = ap.parse_args()

    df = pd.read_csv(a.hap, sep=r"\s+", engine="python")
    need = {"FID", "IID", "Hap"}
    if not need.issubset(df.columns):
        sys.exit("haplogroup file must have columns %s, got %s"
                 % (need, list(df.columns)))

    df["Hap"] = df["Hap"].astype(str).str.upper().str.strip()
    bad = ~df["Hap"].isin(HAP_CODE)
    if bad.any():
        sys.exit("%d rows have Hap outside {I,R}" % int(bad.sum()))

    for label in ("I", "R"):
        sub = df.loc[df["Hap"] == label, ["FID", "IID"]]
        path = "%s/keep_%s.txt" % (a.out_prefix, label)
        sub.to_csv(path, sep="\t", index=False, header=False)
        print("keep_%s.txt: %d samples" % (label, len(sub)))


if __name__ == "__main__":
    main()
