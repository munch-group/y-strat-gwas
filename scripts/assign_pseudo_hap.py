#!/usr/bin/env python3
"""Assign females a *pseudo-haplogroup* from the male I-vs-R ancestry propensity.

Females carry no Y, so a genuine Y-haplogroup x autosome interaction cannot exist
in them -- but the ANCESTRY artefact can. To expose it we split females along the
exact autosomal-ancestry axis that distinguishes male haplogroup I from R: fit the
Hap propensity P(Hap=I | PCs) on the males, apply it to the females, and label the
females pseudo-I / pseudo-R to match the male I:R proportion. Feeding this pseudo-Hap
to the same interaction scan in females turns them into a negative control (see
scripts/female_negcontrol.py).

Requires the female PCs to live in the SAME ancestry space as the male PCs (PCs
computed jointly, or females projected onto the male PCA).
"""
import argparse
import numpy as np
import pandas as pd

HAP_CODE = {"I": 1, "R": 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--male-hap", required=True, help="male haplogroup file (FID IID ... <hap-col>)")
    ap.add_argument("--hap-col", default="Hap", help="name of the haplogroup column")
    ap.add_argument("--male-covar", required=True, help="male FID IID PC1..PCk ...")
    ap.add_argument("--female-covar", required=True, help="female FID IID PC1..PCk ...")
    ap.add_argument("--npc", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pcn = ["PC%d" % i for i in range(1, a.npc + 1)]

    hap = pd.read_csv(a.male_hap, sep=r"\s+", engine="python")
    hap["Hap"] = hap[a.hap_col].astype(str).str.upper().str.strip().map(HAP_CODE)  # I/R only
    mcov = pd.read_csv(a.male_covar, sep=r"\s+", engine="python")
    m = mcov.merge(hap[["FID", "IID", "Hap"]], on=["FID", "IID"], how="inner").dropna(
        subset=pcn + ["Hap"])
    fcov = pd.read_csv(a.female_covar, sep=r"\s+", engine="python").dropna(subset=pcn)

    # standardise female PCs with the MALE mean/sd so the propensity transfers
    mu = m[pcn].mean().values
    sd = m[pcn].std().values.copy()
    sd[sd == 0] = 1.0
    Xm = np.column_stack([np.ones(len(m)), (m[pcn].values - mu) / sd])
    Xf = np.column_stack([np.ones(len(fcov)), (fcov[pcn].values - mu) / sd])

    beta, *_ = np.linalg.lstsq(Xm, m["Hap"].values.astype(float), rcond=None)
    fscore = Xf @ beta                                   # higher => more "I-like"

    i_frac = float((m["Hap"] == 1).mean())               # match the male I:R ratio
    thresh = np.quantile(fscore, 1.0 - i_frac)
    pseudo = np.where(fscore >= thresh, "I", "R")

    out = fcov[["FID", "IID"]].copy()
    out["Hap"] = pseudo
    out.to_csv(a.out, sep="\t", index=False)
    print("wrote %s: %d females (%d pseudo-I, %d pseudo-R; male I frac=%.3f)"
          % (a.out, len(out), int((pseudo == "I").sum()),
             int((pseudo == "R").sum()), i_frac))


if __name__ == "__main__":
    main()
