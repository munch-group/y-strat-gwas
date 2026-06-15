#!/usr/bin/env python3
"""Female negative control for the ancestry artefact in the SNP x Hap scan.

Females have no Y chromosome, so a genuine Y-haplogroup x autosome interaction
*cannot* exist in them. So we take the male interaction hits and re-test each one
for SNP x pseudo-Hap interaction in the females (pseudo-Hap = the ancestry-matched
split from assign_pseudo_hap.py, same model: g + pseudoHap + g:pseudoHap + PCs +
pseudoHap:PCs). The logic:

  * a male hit that is REAL (driven by the Y) is NULL in females           -> survives
  * a male hit that is an ANCESTRY ARTEFACT reproduces in females          -> exposed

A significant female interaction at a male hit is positive evidence that the male
signal is ancestry, not Y. A panel of random SNPs gives the genome-wide female
interaction inflation (lambda) -- if the female scan is inflated overall, the male
scan's structure is ancestry too.

The statistic is the t^2 of the g:pseudoHap term in a linear-probability model (its
parametric chi^2_1 p-value); no permutation here -- the female test IS the confound
detector, we are not trying to deconfound it.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ancestry_matched_perm as amp   # read_fam/read_bim_index/read_bed/interaction_t2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--male-interactions", required=True,
                    help="results/top_interactions.tsv (male hits + LOG10P)")
    ap.add_argument("--fbfile", required=True, help="female plink1 prefix")
    ap.add_argument("--fpheno", required=True)
    ap.add_argument("--pheno-name", required=True)
    ap.add_argument("--fcovar", required=True, help="female FID IID PC1..PCk age batch")
    ap.add_argument("--female-hap", required=True, help="female pseudo-Hap (FID IID Hap)")
    ap.add_argument("--npc", type=int, required=True)
    ap.add_argument("--panel", type=int, default=2000, help="random SNPs for the female lambda")
    ap.add_argument("--force-snps", default="")
    ap.add_argument("--alpha", type=float, default=0.05, help="artefact-flag threshold")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    male = pd.read_csv(a.male_interactions, sep=r"\s+", engine="python")
    male_log10p = dict(zip(male["ID"].astype(str), male["LOG10P"]))
    hits = list(dict.fromkeys([s for s in a.force_snps.split(",") if s]
                              + list(male["ID"].astype(str))))

    fid, iid = amp.read_fam(a.fbfile + ".fam")
    n = len(fid)
    bim_idx = amp.read_bim_index(a.fbfile + ".bim")
    panel_pool = [s for s in bim_idx if s not in set(hits)]
    n_panel = min(a.panel, len(panel_pool))
    panel = list(rng.choice(panel_pool, size=n_panel, replace=False)) if n_panel else []

    want = [s for s in hits + panel if s in bim_idx]
    cols = amp.read_bed(a.fbfile + ".bed", n, [bim_idx[s] for s in want])
    geno = {s: cols[:, i] for i, s in enumerate(want)}

    key = pd.DataFrame({"FID": fid, "IID": iid})
    cov = pd.read_csv(a.fcovar, sep=r"\s+", engine="python")
    hap = pd.read_csv(a.female_hap, sep=r"\s+", engine="python")
    hap["Hap"] = hap["Hap"].astype(str).str.upper().str.strip().map({"I": 1, "R": 0})
    phe = pd.read_csv(a.fpheno, sep=r"\s+", engine="python")
    phe[a.pheno_name] = pd.to_numeric(phe[a.pheno_name], errors="coerce")
    for _d in (cov, hap, phe):
        amp.norm_ids(_d)
    df = (key.merge(cov, on=["FID", "IID"], how="left")
             .merge(hap[["FID", "IID", "Hap"]], on=["FID", "IID"], how="left")
             .merge(phe[["FID", "IID", a.pheno_name]], on=["FID", "IID"], how="left"))

    pcn = ["PC%d" % i for i in range(1, a.npc + 1)]
    extra = [c for c in ("age", "batch") if c in df.columns]
    keep = (df[a.pheno_name].isin([0, 1]) & df["Hap"].isin([0, 1])
            & df[pcn + extra].notna().all(axis=1)).values
    print("females used: %d / %d" % (int(keep.sum()), n))

    y = df[a.pheno_name].values[keep].astype(float)
    hapv = df["Hap"].values[keep].astype(float)
    pcs = df[pcn].values[keep].astype(float)
    pcs_std = (pcs - pcs.mean(0)) / pcs.std(0)
    cbase = np.column_stack([np.ones(keep.sum())]
                            + [df[c].values[keep].astype(float) for c in pcn + extra])
    for s in want:
        geno[s] = geno[s][keep]

    def fem_p(snp):
        t2 = amp.interaction_t2(y, geno[snp], hapv, cbase, pcs_std)
        return t2, float(stats.chi2.sf(t2, 1))

    rows = []
    for snp in hits:
        if snp not in geno:
            continue
        t2, p = fem_p(snp)
        rows.append({"ID": snp, "male_LOG10P": round(float(male_log10p.get(snp, np.nan)), 4),
                     "female_int_chi2": round(t2, 4), "female_int_p": p,
                     "looks_like_ancestry_artifact": p < a.alpha})
    out = pd.DataFrame(rows).sort_values("female_int_p")
    out.to_csv(a.out_prefix + "_negative_control.tsv", sep="\t", index=False)

    panel_in = [s for s in panel if s in geno]
    with open(a.out_prefix + "_lambda.txt", "w") as f:
        f.write("females_used\t%d\n" % int(keep.sum()))
        f.write("hits_tested\t%d\n" % len(out))
        f.write("hits_flagged_artifact\t%d\n" % int(out["looks_like_ancestry_artifact"].sum()))
        if panel_in:
            chi = np.array([fem_p(s)[0] for s in panel_in])
            f.write("panel_snps\t%d\n" % len(panel_in))
            f.write("female_interaction_lambda\t%.4f\n" % (np.median(chi) / 0.4549))
    print("wrote %s_negative_control.tsv (%d male hits re-tested in females; %d flagged)"
          % (a.out_prefix, len(out), int(out["looks_like_ancestry_artifact"].sum())))


if __name__ == "__main__":
    main()
