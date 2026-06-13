#!/usr/bin/env python3
"""Ancestry-matched permutation test for LOCAL cross-stratum genetic correlation,
with a negative-control-loci panel.

This is the local-r_g analogue of scripts/ancestry_matched_perm.py. LAVA-style
local genetic correlation between ASD-in-I and ASD-in-R is corrupted by the SAME
ancestry confounding as the single-variant scan -- differential LD between the I
and R ancestral backgrounds depresses the cross-stratum effect concordance even
with no real epistasis, and a local r_g below 1 *genome-wide* is exactly the
signature of residual structure. So a focused local-r_g claim at a locus must be
deconfounded the same way, plus calibrated against control loci.

What it does for a target locus (and a panel of SNP-count-matched control loci):

  * statistic  T = (b_I' R b_R) / sqrt( (b_I' R b_I)(b_R' R b_R) )
    an LD(R)-aware cross-stratum correlation of per-stratum marginal betas
    (b_I, b_R) over the block -- the fast, sign-sensitive proxy for LAVA local
    r_g. PC-residualised so global ancestry is regressed out of betas and LD.
  * null: permute Hap WITHIN Hap-propensity strata (between-stratum ancestry
    structure preserved, within-stratum signal destroyed); recompute T.
  * empirical p per locus, and the position of the target within the control
    distribution -- the check that distinguishes "this locus is special" from
    "local r_g is depressed everywhere by residual structure".

Controls confounding only to the resolution of the supplied PCs (ideally
region-excluded; see the lava_pcs workflow target). Residual LOCAL ancestry at a
locus is not removed -- survivors still need local ancestry or replication.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ancestry_matched_perm as amp   # noqa: E402  (sibling helper, reused)


def read_bim(path):
    b = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    return (b[0].astype(str).values, b[1].astype(str).values, b[3].astype(int).values)


def parse_locus(s):
    chrom, rng = s.split(":")
    lo, hi = rng.replace("_", "").split("-")
    return chrom, int(lo), int(hi)


def block_indices(chrom_arr, pos_arr, chrom, lo, hi):
    return np.where((chrom_arr == str(chrom)) & (pos_arr >= lo) & (pos_arr <= hi))[0]


def sample_control_blocks(chrom_arr, pos_arr, target_idx, m, n_ctrl,
                          exclude_chrom, rng):
    """n_ctrl windows of ~m consecutive SNPs, away from the target chromosome,
    matched on SNP count (a contiguous run = realistic LD)."""
    order = np.argsort(pos_arr, kind="stable")
    blocks, tries = [], 0
    elig = np.where(chrom_arr != str(exclude_chrom))[0]
    by_chrom = {}
    for i in elig:
        by_chrom.setdefault(chrom_arr[i], []).append(i)
    chroms = [c for c, v in by_chrom.items() if len(v) >= m]
    while len(blocks) < n_ctrl and tries < n_ctrl * 50 and chroms:
        tries += 1
        c = chroms[int(rng.integers(len(chroms)))]
        idx = np.array(sorted(by_chrom[c], key=lambda k: pos_arr[k]))
        start = int(rng.integers(0, len(idx) - m + 1))
        blocks.append(idx[start:start + m])
    return blocks


def local_rg(betaI, betaR, R):
    cov = float(betaI @ R @ betaR)
    vI = float(betaI @ R @ betaI)
    vR = float(betaR @ R @ betaR)
    if vI <= 0 or vR <= 0:
        return 0.0
    return cov / np.sqrt(vI * vR)


def write_local_rg_summary(out_prefix, df, ssize, sauc):
    """df: one row per locus with local_rg, anc_matched_emp_p, is_target. Writes the
    _summary.txt (target stats + where it sits among the negative controls). Shared
    by the single-job path and the batch pooler."""
    t = df[df.is_target].iloc[0]
    ctrl = df[~df.is_target]
    tail = (1 + int((ctrl.local_rg.abs() >= abs(t.local_rg)).sum())) / (1 + len(ctrl))
    nan = float("nan")
    with open(out_prefix + "_summary.txt", "w") as f:
        f.write("target_locus\t%s\n" % t.locus)
        f.write("target_local_rg\t%.4f\n" % t.local_rg)
        f.write("target_anc_matched_p\t%.4g\n" % t.anc_matched_emp_p)
        f.write("n_controls\t%d\n" % len(ctrl))
        f.write("target_vs_controls_tail_frac\t%.4g\n" % tail)
        f.write("control_local_rg_mean\t%.4f\n" % (ctrl.local_rg.mean() if len(ctrl) else nan))
        f.write("control_p_median\t%.4f\n" % (ctrl.anc_matched_emp_p.median() if len(ctrl) else nan))
        f.write("strata_mean_size\t%d\n" % ssize)
        f.write("within_stratum_HapPC_AUC\t%.4f\n" % sauc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bfile", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-name", required=True)
    ap.add_argument("--hap", required=True)
    ap.add_argument("--pcs", required=True,
                    help="plink .eigenvec OR a covar file with PC1..PCk columns")
    ap.add_argument("--npc", type=int, required=True)
    ap.add_argument("--locus", required=True, help="target block as chr:start-end")
    ap.add_argument("--locus-name", default="target")
    ap.add_argument("--n-controls", type=int, default=300)
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--strata-sizes", default="20,40,80,160")
    ap.add_argument("--strata-auc", type=float, default=0.55)
    ap.add_argument("--seed", type=int, default=1, help="permutation seed (vary per batch)")
    ap.add_argument("--select-seed", type=int, default=1,
                    help="control-block selection seed; keep FIXED across batches so they pool")
    ap.add_argument("--raw-counts", action="store_true",
                    help="emit poolable per-batch raw counts instead of final p-values")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)            # permutations (batch-specific)
    sel_rng = np.random.default_rng(a.select_seed)  # control selection (fixed across batches)

    # ---- align samples on the .fam order ---------------------------------
    fid, iid = amp.read_fam(a.bfile + ".fam")
    n = len(fid)
    chrom_arr, snp_ids, pos_arr = read_bim(a.bfile + ".bim")
    key = pd.DataFrame({"FID": fid, "IID": iid})

    def read_pc_table(path):
        df = pd.read_csv(path, sep=r"\s+", engine="python")
        df.columns = [c.lstrip("#") for c in df.columns]
        return df

    pcs_df = read_pc_table(a.pcs)
    hap = pd.read_csv(a.hap, sep=r"\s+", engine="python")
    hap["Hap"] = hap["Hap"].astype(str).str.upper().str.strip().map({"I": 1, "R": 0})
    phe = pd.read_csv(a.pheno, sep=r"\s+", engine="python")
    phe[a.pheno_name] = pd.to_numeric(phe[a.pheno_name], errors="coerce")

    pcn = ["PC%d" % i for i in range(1, a.npc + 1)]
    df = (key.merge(pcs_df[["FID", "IID"] + pcn], on=["FID", "IID"], how="left")
             .merge(hap[["FID", "IID", "Hap"]], on=["FID", "IID"], how="left")
             .merge(phe[["FID", "IID", a.pheno_name]], on=["FID", "IID"], how="left"))
    keep = (df[a.pheno_name].isin([0, 1]) & df["Hap"].isin([0, 1])
            & df[pcn].notna().all(axis=1)).values
    print("samples used: %d / %d" % (int(keep.sum()), n))

    y = df[a.pheno_name].values[keep].astype(float)
    hapv = df["Hap"].values[keep].astype(float)
    pcs = df[pcn].values[keep].astype(float)
    pcs_std = (pcs - pcs.mean(0)) / pcs.std(0)
    keep_rows = np.where(keep)[0]

    # PC-residualiser (regress 1+PCs out of Y and genotypes once) ----------
    C = np.column_stack([np.ones(len(y)), pcs_std])
    P = C @ np.linalg.pinv(C)                          # hat matrix
    yr = y - P @ y

    # ---- ancestry-matched strata (global; the tightening) ----------------
    sizes = [int(x) for x in a.strata_sizes.split(",")]
    (labels, ssize, sauc, mixed_lab, eff_n), grid = amp.select_granularity(
        pcs_std, hapv, sizes, a.strata_auc, a.seed)
    grid.to_csv(a.out_prefix + "_strata_selection.tsv", sep="\t", index=False)
    mixed_idx = [np.where(labels == s)[0] for s in mixed_lab]
    print("strata: mean size %d, within-stratum Hap~PC AUC=%.3f, eff_N=%d"
          % (ssize, sauc, eff_n))

    def block_stat_perm(idx):
        """Return (T_obs, n_ge, n_snps): n_ge = #perms with |T| >= |T_obs|."""
        G = amp.read_bed(a.bfile + ".bed", n, list(idx))[keep_rows]
        sd = G.std(0)
        sd[sd == 0] = 1.0
        Z = (G - G.mean(0)) / sd
        Zr = Z - P @ Z                                 # PC-residualised genotypes
        R = (Zr.T @ Zr) / len(y)                       # LD-aware weighting

        def stat(hv):
            mI, mR = hv == 1, hv == 0
            bI = (Zr[mI].T @ yr[mI]) / np.maximum((Zr[mI] ** 2).sum(0), 1e-12)
            bR = (Zr[mR].T @ yr[mR]) / np.maximum((Zr[mR] ** 2).sum(0), 1e-12)
            return local_rg(bI, bR, R)

        t_obs = stat(hapv)
        n_ge = int(sum(abs(stat(amp.permute_within(hapv, mixed_idx, rng))) >= abs(t_obs)
                       for _ in range(a.nperm)))
        return t_obs, n_ge, len(idx)

    # ---- target locus + negative-control panel (controls fixed per batch) -
    tchrom, lo, hi = parse_locus(a.locus)
    tidx = block_indices(chrom_arr, pos_arr, tchrom, lo, hi)
    if len(tidx) < 3:
        raise SystemExit("target locus has <3 SNPs (%d)" % len(tidx))
    blocks = [(a.locus_name, tchrom, tidx, True)]
    for bi, idx in enumerate(sample_control_blocks(chrom_arr, pos_arr, tidx, len(tidx),
                                                    a.n_controls, tchrom, sel_rng)):
        blocks.append(("control_%d" % bi, chrom_arr[idx[0]], idx, False))

    rows = []
    for name, ch, idx, is_t in blocks:
        to, n_ge, m = block_stat_perm(idx)
        rows.append({"locus": name, "chr": ch, "n_snps": m,
                     "local_rg": round(to, 4), "n_ge": n_ge, "n_perm": a.nperm,
                     "is_target": is_t})
    df_rows = pd.DataFrame(rows)

    if a.raw_counts:                                  # one batch -> pool later
        df_rows.to_csv(a.out_prefix + "_counts.tsv", sep="\t", index=False)
        with open(a.out_prefix + "_meta.tsv", "w") as f:
            f.write("strata_mean_size\t%d\n" % ssize)
            f.write("within_stratum_HapPC_AUC\t%.4f\n" % sauc)
        print("wrote %s_counts.tsv (batch seed=%d, nperm=%d)"
              % (a.out_prefix, a.seed, a.nperm))
        return

    df_rows["anc_matched_emp_p"] = (1 + df_rows.n_ge) / (1 + df_rows.n_perm)
    df_rows[["locus", "chr", "n_snps", "local_rg", "anc_matched_emp_p",
             "is_target"]].to_csv(a.out_prefix + "_loci.tsv", sep="\t", index=False)
    write_local_rg_summary(a.out_prefix, df_rows, ssize, sauc)
    print("wrote %s_loci.tsv" % a.out_prefix)


if __name__ == "__main__":
    main()
