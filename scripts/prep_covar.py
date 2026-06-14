#!/usr/bin/env python3
"""Normalise a raw covariate / PC file into the form the workflow expects:
FID IID PC1..PCk <extra quantitative> <extra categorical>.

Real PC files often (a) name the PCs differently (e.g. C1..C10, with extra columns
like SOL / st1) and (b) lack the age / batch columns the covariate model references.
This renames <pc-prefix>1..<pc-prefix>k -> PC1..PCk, keeps FID/IID, and adds any
requested quantitative (--add-quant) or categorical (--add-cat) columns that are
missing, as **dummy placeholders** (deterministic, seeded; small variance so REGENIE
doesn't drop them for zero variance). Columns already present are kept as-is.

The dummies carry no real information -- they exist only so the model's covariate
list resolves. To avoid adjusting for placeholders entirely, instead leave them out
of the model (set EXTRA_COVAR="" / CATCOVAR="" in workflow.py) and drop --add-*.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-covar", required=True)
    ap.add_argument("--npc", type=int, required=True)
    ap.add_argument("--pc-prefix", default="PC", help="PC column prefix in the raw file (e.g. C)")
    ap.add_argument("--add-quant", default="", help="comma-sep quantitative dummy cols to add if missing")
    ap.add_argument("--add-cat", default="", help="comma-sep categorical dummy cols to add if missing")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    df = pd.read_csv(a.raw_covar, sep=r"\s+", engine="python")
    if not {"FID", "IID"}.issubset(df.columns):
        raise SystemExit("raw covar must have FID and IID columns; got %s" % list(df.columns))

    src = ["%s%d" % (a.pc_prefix, i) for i in range(1, a.npc + 1)]
    miss = [c for c in src if c not in df.columns]
    if miss:
        raise SystemExit("missing PC columns %s in %s (got %s)"
                         % (miss, a.raw_covar, list(df.columns)))
    pcn = ["PC%d" % i for i in range(1, a.npc + 1)]
    out = df[["FID", "IID"] + src].rename(columns=dict(zip(src, pcn)))

    n = len(out)
    added = []
    for c in [x.strip() for x in a.add_quant.split(",") if x.strip()]:
        if c in df.columns:
            out[c] = df[c].values
        else:
            out[c] = rng.integers(18, 80, size=n)          # dummy quantitative
            added.append(c)
    for c in [x.strip() for x in a.add_cat.split(",") if x.strip()]:
        if c in df.columns:
            out[c] = df[c].values
        else:
            out[c] = rng.integers(0, 2, size=n)            # dummy 2-level categorical
            added.append(c)

    out.to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print("wrote %s: %d samples, %s%d..%s%d -> PC1..PC%d%s"
          % (a.out, n, a.pc_prefix, 1, a.pc_prefix, a.npc, a.npc,
             ("; added dummy cols: " + ",".join(added)) if added else ""))


if __name__ == "__main__":
    main()
