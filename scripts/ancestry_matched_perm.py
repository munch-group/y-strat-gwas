#!/usr/bin/env python3
"""Ancestry-matched permutation test for the SNP x Hap interaction.

Why: REGENIE's --interaction guards the E-side (PC x Hap) but cannot add the
per-SNP G x PC term, so a SNP whose effect varies along ancestry can masquerade
as a Hap interaction (Hap tags ancestry). This re-tests interaction hits against
a null that *preserves* the Hap<->ancestry structure but breaks any genuine
within-ancestry G x Hap effect.

How: permute Hap **within ancestry-matched strata** (bins of the Hap propensity
score P(Hap|PCs) -- the PC-space direction along which I and R separate) and
recompute the interaction statistic. Within a stratum ancestry is ~constant, so
genotype is independent of Hap; a real g x Hap effect lives within strata and is
destroyed by the permutation (observed > null => significant), while a pure
ancestry artifact lives between strata and is preserved by it (observed inside
null => not significant).

The strata are only as good as their ancestry resolution, so granularity is
chosen by a principled rule: make strata fine enough that, *within a stratum*,
Hap is unpredictable from the PCs (AUC of Hap ~ within-stratum-centred PCs near
0.5), then take the LOOSEST strata that still achieve it (most power). The
selection curve is written out. Stratifying on the 1-D propensity rather than raw
multi-PC distance is what makes this reachable -- Euclidean matching across many
PCs wastes resolution on ancestry-uninformative (noise) PCs and never decouples.

Statistic: the t^2 of the g:Hap term in a linear-probability model
    Y ~ 1 + PCs + age + batch + g + Hap + Hap:PCs + g:Hap
(a fast, closed-form statistic; its permutation null is valid regardless of the
linear approximation). REGENIE's logistic LOG10P is carried through for ref.

Outputs (with --out-prefix P):
  P_strata_selection.tsv  granularity grid: target size, within-stratum Hap~PC AUC, mixed strata, eff. N
  P_lambda.txt            genome-wide: observed interaction lambda vs the within-strata permutation null
  P_interactions.tsv      per top hit: raw LOG10P, observed stat, ancestry-matched empirical p
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata


# ---------------------------------------------------------------------------
# plink1 .bed reader (specific SNPs only; SNP-major)
# ---------------------------------------------------------------------------
_BED_LUT = np.array([2.0, np.nan, 1.0, 0.0])   # codes 00,01,10,11 -> dosage


def read_fam(path):
    fam = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    return fam[0].astype(str).values, fam[1].astype(str).values


def norm_ids(df, cols=("FID", "IID")):
    """Force FID/IID to str before any merge. Real all-numeric IDs get inferred
    as int64 by read_csv in one file and as str (via read_fam) elsewhere, and
    pandas refuses to merge str against int64. Coercing both sides to str fixes
    it (and is a no-op on the synthetic non-numeric IDs)."""
    for c in cols:
        if c in df.columns:
            df[c] = df[c].astype(str)
    return df


def read_bim_index(path):
    bim = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    return {sid: i for i, sid in enumerate(bim[1].astype(str).values)}


def read_bed(bed_path, n_samples, snp_indices):
    bps = (n_samples + 3) // 4
    out = np.empty((n_samples, len(snp_indices)), dtype=np.float64)
    with open(bed_path, "rb") as f:
        magic = f.read(3)
        if not (magic[0] == 0x6C and magic[1] == 0x1B and magic[2] == 0x01):
            raise SystemExit("not a SNP-major plink1 .bed file")
        for col, j in enumerate(snp_indices):
            f.seek(3 + j * bps)
            raw = np.frombuffer(f.read(bps), dtype=np.uint8)
            codes = np.empty(bps * 4, dtype=np.uint8)
            codes[0::4] = raw & 0b11
            codes[1::4] = (raw >> 2) & 0b11
            codes[2::4] = (raw >> 4) & 0b11
            codes[3::4] = (raw >> 6) & 0b11
            g = _BED_LUT[codes[:n_samples]]
            m = np.isnan(g)
            if m.any():                          # mean-impute missing
                g[m] = np.nanmean(g) if (~m).any() else 0.0
            out[:, col] = g
    return out


# ---------------------------------------------------------------------------
# interaction statistic + within-strata permutation
# ---------------------------------------------------------------------------
def interaction_t2(y, g, hap, cbase, pcs):
    """t^2 of the g:Hap coefficient in OLS  y ~ cbase + g + hap + hap:pcs + g:hap."""
    X = np.column_stack([cbase, g, hap, pcs * hap[:, None], g * hap])
    XtX = X.T @ X
    XtX[np.diag_indices_from(XtX)] += 1e-8
    Xty = X.T @ y
    inv = np.linalg.inv(XtX)
    beta = inv @ Xty
    rss = float(y @ y - beta @ Xty)
    dof = max(1, len(y) - X.shape[1])
    var_last = (rss / dof) * inv[-1, -1]
    return float(beta[-1] ** 2 / var_last) if var_last > 0 else 0.0


def auc(score, label):
    """AUC of `score` discriminating binary `label` (Mann-Whitney)."""
    pos = label == 1
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return 0.5
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def within_stratum_hap_pc_auc(hap, pcs, labels):
    """Predict Hap from within-stratum-centred PCs; AUC ~0.5 => Hap decoupled
    from ancestry inside strata (the homogeneity target)."""
    cp = pcs.copy().astype(float)
    ch = hap.astype(float).copy()
    for s in np.unique(labels):
        m = labels == s
        cp[m] -= cp[m].mean(0)
        ch[m] -= ch[m].mean()
    if np.allclose(ch, 0):                       # all strata pure => decoupled
        return 0.5
    b, *_ = np.linalg.lstsq(cp, ch, rcond=None)
    return auc(cp @ b, hap)


def propensity_score(pcs_std, hap):
    """Linear propensity P(Hap | PCs): the single PC-space direction along which
    I and R separate -- i.e. the confounding axis. Matching on it (not raw
    multi-PC distance) is what actually decouples Hap from ancestry, and it
    ignores ancestry-uninformative (noise) PCs automatically."""
    X = np.column_stack([np.ones(len(hap)), pcs_std])
    beta, *_ = np.linalg.lstsq(X, hap, rcond=None)
    return X @ beta


def strata_by_propensity(score, n_bins):
    edges = np.quantile(score, np.linspace(0, 1, n_bins + 1))
    edges[0] -= 1e-9
    edges[-1] += 1e-9
    return np.digitize(score, edges[1:-1])


def _summ(labels, hap):
    mixed = [s for s in np.unique(labels)
             if hap[labels == s].min() != hap[labels == s].max()]
    eff_n = int(sum((labels == s).sum() for s in mixed))
    return mixed, eff_n


def select_granularity(pcs_std, hap, target_sizes, auc_thresh, seed):
    """Bin individuals by Hap-propensity at varying granularity; pick the LOOSEST
    strata (largest mean size) whose within-stratum Hap~PC AUC <= auc_thresh
    (most power that still decouples). Returns (chosen, grid_rows)."""
    n = len(hap)
    score = propensity_score(pcs_std, hap)
    rows, chosen = [], None
    for size in sorted(target_sizes, reverse=True):     # loose -> tight
        nb = max(2, int(round(n / size)))
        labels = strata_by_propensity(score, nb)
        a = within_stratum_hap_pc_auc(hap, pcs_std, labels)
        mixed, eff_n = _summ(labels, hap)
        rows.append({"target_size": size, "n_bins": nb,
                     "within_stratum_HapPC_AUC": round(a, 4),
                     "n_mixed_strata": len(mixed), "eff_N": eff_n})
        if chosen is None and a <= auc_thresh and eff_n > 0:
            chosen = (labels, size, a, mixed, eff_n)
    if chosen is None:                                  # fall back to tightest
        size = min(target_sizes)
        labels = strata_by_propensity(score, max(2, int(round(n / size))))
        a = within_stratum_hap_pc_auc(hap, pcs_std, labels)
        mixed, eff_n = _summ(labels, hap)
        chosen = (labels, size, a, mixed, eff_n)
    return chosen, pd.DataFrame(rows)


def permute_within(hap, mixed_idx, rng):
    hp = hap.copy()
    for idx in mixed_idx:
        hp[idx] = rng.permutation(hp[idx])
    return hp


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regenie", required=True, help="gxhap interaction .regenie")
    ap.add_argument("--test", default="ADD-INT_SNPxHap")
    ap.add_argument("--bfile", required=True, help="plink1 prefix (.bed/.bim/.fam)")
    ap.add_argument("--covar", required=True)
    ap.add_argument("--hap", required=True)
    ap.add_argument("--hap-col-name", default='Hap')
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-name", required=True)
    ap.add_argument("--npc", type=int, required=True)
    ap.add_argument("--top", type=int, default=500, help="top interaction hits to re-test")
    ap.add_argument("--panel", type=int, default=5000, help="random SNPs for the lambda diagnostic")
    ap.add_argument("--nperm", type=int, default=1000, help="permutations per hit")
    ap.add_argument("--global-nperm", type=int, default=200, help="permutations for lambda")
    ap.add_argument("--force-snps", default="", help="comma-sep SNP IDs to always include")
    ap.add_argument("--strata-sizes", default="20,40,80,160",
                    help="candidate mean stratum sizes (loose->tight search)")
    ap.add_argument("--strata-auc", type=float, default=0.55,
                    help="within-stratum Hap~PC AUC at/below which strata count as homogeneous")
    ap.add_argument("--seed", type=int, default=1, help="permutation seed (vary per batch)")
    ap.add_argument("--select-seed", type=int, default=1,
                    help="panel-selection seed; keep FIXED across batches so they pool")
    ap.add_argument("--raw-counts", action="store_true",
                    help="emit poolable per-batch raw counts instead of final p-values")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)            # permutations (batch-specific)
    sel_rng = np.random.default_rng(a.select_seed)  # selection (fixed across batches)

    # ---- interaction results: pick hits + a random panel -----------------
    reg = pd.read_csv(a.regenie, sep=r"\s+", engine="python", comment="#")
    reg = reg[reg["TEST"] == a.test].copy()
    if reg.empty:
        raise SystemExit("no rows for TEST=%s in %s" % (a.test, a.regenie))
    reg = reg.sort_values("LOG10P", ascending=False)
    hits = list(reg["ID"].head(a.top))
    forced = [s for s in a.force_snps.split(",") if s]
    hits = list(dict.fromkeys(forced + hits))           # forced first, dedup
    panel_pool = reg["ID"].iloc[a.top:].values
    n_panel = min(a.panel, len(panel_pool))
    panel = list(sel_rng.choice(panel_pool, size=n_panel, replace=False)) if n_panel else []
    raw_log10p = dict(zip(reg["ID"], reg["LOG10P"]))

    # ---- align genotypes + covariates + hap + pheno on the .fam order ----
    fid, iid = read_fam(a.bfile + ".fam")
    n = len(fid)
    key = pd.DataFrame({"FID": fid, "IID": iid})
    bim_idx = read_bim_index(a.bfile + ".bim")

    want = [s for s in hits + panel if s in bim_idx]
    missing = [s for s in hits if s not in bim_idx]
    if missing:
        print("WARNING: %d requested SNPs not in .bim (skipped): %s"
              % (len(missing), ",".join(missing[:5])))
    snp_cols = read_bed(a.bfile + ".bed", n, [bim_idx[s] for s in want])
    geno = {s: snp_cols[:, i] for i, s in enumerate(want)}

    cov = pd.read_csv(a.covar, sep=r"\s+", engine="python")
    hap = pd.read_csv(a.hap, sep=r"\s+", engine="python")
    hap[a.hap_col_name] = hap[a.hap_col_name].astype(str).str.upper().str.strip().map({"I": 1, "R": 0})
    phe = pd.read_csv(a.pheno, sep=r"\s+", engine="python")
    phe[a.pheno_name] = pd.to_numeric(phe[a.pheno_name], errors="coerce")

    for _d in (cov, hap, phe):
        norm_ids(_d)
    df = (key.merge(cov, on=["FID", "IID"], how="left")
             .merge(hap[["FID", "IID", a.hap_col_name]], on=["FID", "IID"], how="left")
             .merge(phe[["FID", "IID", a.pheno_name]], on=["FID", "IID"], how="left"))

    pcn = ["PC%d" % i for i in range(1, a.npc + 1)]
    extra = [c for c in ("age", "batch") if c in df.columns]
    keep = df[a.pheno_name].isin([0, 1]) & df[a.hap_col_name].isin([0, 1])
    keep &= df[pcn + extra].notna().all(axis=1)
    keep = keep.values
    print("samples used in permutation test: %d / %d" % (int(keep.sum()), n))

    y = df[a.pheno_name].values[keep].astype(float)
    hapv = df[a.hap_col_name].values[keep].astype(float)
    pcs = df[pcn].values[keep].astype(float)
    cbase = np.column_stack([np.ones(keep.sum())]
                            + [df[c].values[keep].astype(float) for c in pcn + extra])
    pcs_std = (pcs - pcs.mean(0)) / pcs.std(0)
    for s in want:
        geno[s] = geno[s][keep]

    # ---- choose stratum granularity (the tightening) ---------------------
    sizes = [int(x) for x in a.strata_sizes.split(",")]
    (labels, chosen_size, chosen_auc, mixed_lab, eff_n), grid = select_granularity(
        pcs_std, hapv, sizes, a.strata_auc, a.seed)
    grid.to_csv(a.out_prefix + "_strata_selection.tsv", sep="\t", index=False)
    mixed_idx = [np.where(labels == s)[0] for s in mixed_lab]
    print("selected mean stratum size %d : within-stratum Hap~PC AUC=%.3f, "
          "%d mixed strata, eff_N=%d" % (chosen_size, chosen_auc, len(mixed_idx), eff_n))
    if chosen_auc > a.strata_auc:
        print("WARNING: could not decouple Hap from PCs within strata "
              "(AUC=%.3f > %.2f); residual ancestry confounding remains."
              % (chosen_auc, a.strata_auc))

    def stat(snp, hv):
        return interaction_t2(y, geno[snp], hv, cbase, pcs_std)

    # ---- genome-wide lambda: observed + within-strata permutation null ----
    panel_in = [s for s in panel if s in geno]
    obs_lambda, null_lambdas = float("nan"), np.array([])
    if panel_in:
        obs = np.array([stat(s, hapv) for s in panel_in])
        obs_lambda = float(np.median(obs) / 0.4549)
        null_lambdas = np.array(
            [np.median(np.array([stat(s, permute_within(hapv, mixed_idx, rng))
                                 for s in panel_in])) / 0.4549
             for _ in range(a.global_nperm)])

    # ---- per-hit: observed stat + #(null >= obs) over nperm permutations ---
    hit_rows = []
    for snp in hits:
        if snp not in geno:
            continue
        o = stat(snp, hapv)
        n_ge = int(sum(stat(snp, permute_within(hapv, mixed_idx, rng)) >= o
                       for _ in range(a.nperm)))
        hit_rows.append({"ID": snp,
                         "raw_LOG10P": round(float(raw_log10p.get(snp, np.nan)), 4),
                         "obs_stat": round(o, 4), "n_ge": n_ge, "n_perm": a.nperm,
                         "forced": snp in forced})

    meta = [("panel_snps", len(panel_in)), ("obs_lambda", "%.6f" % obs_lambda),
            ("selected_stratum_size", chosen_size),
            ("within_stratum_HapPC_AUC", "%.4f" % chosen_auc),
            ("n_mixed_strata", len(mixed_idx)), ("eff_N", eff_n)]

    if a.raw_counts:                                  # one batch -> pool later
        pd.DataFrame(hit_rows).to_csv(a.out_prefix + "_counts.tsv", sep="\t", index=False)
        with open(a.out_prefix + "_lambda_null.tsv", "w") as f:
            for k, v in meta:
                f.write("%s\t%s\n" % (k, v))
            f.write("#NULL\n")
            for v in null_lambdas:
                f.write("%.6f\n" % v)
        print("wrote %s_{counts,lambda_null}.tsv (batch seed=%d, nperm=%d)"
              % (a.out_prefix, a.seed, a.nperm))
        return

    rows = [{"ID": r["ID"], "raw_LOG10P": r["raw_LOG10P"], "obs_stat": r["obs_stat"],
             "n_perm": r["n_perm"], "anc_matched_emp_p": (1 + r["n_ge"]) / (1 + r["n_perm"]),
             "forced": r["forced"]} for r in hit_rows]
    pd.DataFrame(rows).sort_values("anc_matched_emp_p").to_csv(
        a.out_prefix + "_interactions.tsv", sep="\t", index=False)
    with open(a.out_prefix + "_lambda.txt", "w") as f:
        if panel_in:
            p = (1 + int((null_lambdas >= obs_lambda).sum())) / (1 + a.global_nperm)
            f.write("panel_snps\t%d\n" % len(panel_in))
            f.write("observed_lambda\t%.4f\n" % obs_lambda)
            f.write("perm_null_lambda_mean\t%.4f\n" % null_lambdas.mean())
            f.write("perm_null_lambda_2.5pct\t%.4f\n" % np.percentile(null_lambdas, 2.5))
            f.write("perm_null_lambda_97.5pct\t%.4f\n" % np.percentile(null_lambdas, 97.5))
            f.write("lambda_emp_p\t%.4g\n" % p)
        else:
            f.write("panel_snps\t0\n")
        f.write("selected_stratum_size\t%d\n" % chosen_size)
        f.write("within_stratum_HapPC_AUC\t%.4f\n" % chosen_auc)
        f.write("n_mixed_strata\t%d\n" % len(mixed_idx))
        f.write("eff_N\t%d\n" % eff_n)
    print("wrote %s_interactions.tsv (%d hits re-tested)" % (a.out_prefix, len(rows)))


if __name__ == "__main__":
    main()
