#!/usr/bin/env python3
"""Recode a binary phenotype to what REGENIE --bt expects: control=0, case=1,
everything else = NA (missing).

Real phenotype files are often in PLINK 1/2 coding (control=1, case=2). REGENIE
treats 0/1; with 1/2 it mis-reads the trait. This maps --case -> 1, --control -> 0
and writes FID IID <pheno-name>.

Two input layouts are supported:
  * WITH header (default): columns selected by name -- FID, IID and --raw-col.
  * HEADERLESS (--no-header): a bare PLINK .pheno, "FID IID value" with no header
    row, where the phenotype value is the 1-based --value-col (default 3). The
    output always carries a header (FID IID <pheno-name>) so REGENIE can pick the
    column by name.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-pheno", required=True)
    ap.add_argument("--pheno-name", required=True, help="output phenotype column name")
    ap.add_argument("--raw-col", default=None, help="phenotype column in the raw file (default: --pheno-name)")
    ap.add_argument("--no-header", action="store_true",
                    help="raw file has no header row (bare PLINK .pheno: FID IID value)")
    ap.add_argument("--value-col", type=int, default=3,
                    help="with --no-header: 1-based column holding the phenotype value")
    ap.add_argument("--case", default="2", help="raw value meaning case")
    ap.add_argument("--control", default="1", help="raw value meaning control")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.no_header:
        df = pd.read_csv(a.raw_pheno, sep=r"\s+", engine="python", header=None)
        need = max(2, a.value_col)
        if df.shape[1] < need:
            raise SystemExit("headerless phenotype file has %d columns; need >= %d "
                             "(FID, IID, value at col %d)"
                             % (df.shape[1], need, a.value_col))
        fid, iid = df.iloc[:, 0], df.iloc[:, 1]
        v = pd.to_numeric(df.iloc[:, a.value_col - 1], errors="coerce")
    else:
        raw_col = a.raw_col or a.pheno_name
        df = pd.read_csv(a.raw_pheno, sep=r"\s+", engine="python")
        for c in ("FID", "IID", raw_col):
            if c not in df.columns:
                raise SystemExit("phenotype file missing column %r (got %s)"
                                 % (c, list(df.columns)))
        fid, iid = df["FID"], df["IID"]
        v = pd.to_numeric(df[raw_col], errors="coerce")

    out_v = pd.Series(np.nan, index=df.index)
    out_v[v == float(a.case)] = 1
    out_v[v == float(a.control)] = 0

    out = pd.DataFrame({"FID": fid.values, "IID": iid.values})
    out[a.pheno_name] = out_v.values
    out.to_csv(a.out, sep="\t", index=False, na_rep="NA")
    n_case = int((out_v == 1).sum())
    n_con = int((out_v == 0).sum())
    n_na = int(out_v.isna().sum())
    print("wrote %s: %d cases (was %s), %d controls (was %s), %d NA"
          % (a.out, n_case, a.case, n_con, a.control, n_na))


if __name__ == "__main__":
    main()
