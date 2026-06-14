#!/usr/bin/env python3
"""Map Y-haplogroup labels to per-stratum REGENIE --keep lists.

The haplogroup column is given by --hap-col (default "Hap"; real data may call it
e.g. "Major"). Rows whose haplogroup is NOT I or R (other lineages, missing) are
DROPPED -- they are excluded from the whole analysis. Writes keep_I.txt, keep_R.txt
and keep_IR.txt (the union, used to restrict the full-sample targets).
"""
import argparse
import sys
import pandas as pd

HAP_CODE = {"I": 1, "R": 0}  # I=1, R=0 (consistent with interaction coding)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hap", required=True, help="haplogroup file (FID IID ... <hap-col> ...)")
    ap.add_argument("--hap-col", default="Hap", help="name of the haplogroup column")
    ap.add_argument("--out-prefix", default=".")
    a = ap.parse_args()

    df = pd.read_csv(a.hap, sep=r"\s+", engine="python")
    need = {"FID", "IID", a.hap_col}
    if not need.issubset(df.columns):
        sys.exit("haplogroup file must have columns %s, got %s"
                 % (need, list(df.columns)))

    h = df[a.hap_col].astype(str).str.upper().str.strip()
    keep = h.isin(HAP_CODE)
    df = df.loc[keep, ["FID", "IID"]].assign(Hap=h[keep].values)
    print("kept %d of %d samples (dropped %d not in {I,R})"
          % (len(df), len(keep), int((~keep).sum())))

    for label in ("I", "R"):
        sub = df.loc[df["Hap"] == label, ["FID", "IID"]]
        sub.to_csv("%s/keep_%s.txt" % (a.out_prefix, label),
                   sep="\t", index=False, header=False)
        print("keep_%s.txt: %d samples" % (label, len(sub)))
    df[["FID", "IID"]].to_csv("%s/keep_IR.txt" % a.out_prefix,
                              sep="\t", index=False, header=False)
    print("keep_IR.txt: %d samples (I + R, used to restrict full-sample targets)" % len(df))


if __name__ == "__main__":
    main()
