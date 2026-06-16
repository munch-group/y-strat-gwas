#!/usr/bin/env python3
"""Cheapest first look: is chrX more sensitive to Y-haplogroup than the autosome?

Compares the *aggregate* inflation of the SNP x Hap interaction statistic on the
autosome vs chrX, using the two REGENIE --interaction scans the pipeline already
produces (autosomal `interaction` -> gxhap_<pheno>.regenie, and chrX
`interaction_X` -> gxhapX_<pheno>.regenie). No new model fitting, no chrX LDSC
reference -- it just reads the per-SNP interaction p-values.

The question "is X *as a class* more sensitive?" is an aggregate, 1-df contrast
(it pools every SNP on the chromosome) and so is far better powered than the
per-SNP interaction scan or the cross-stratum rg. This script is the GO/NO-GO
screen: if chrX shows no excess interaction inflation over the autosome, there is
nothing to chase and you can skip building the calibrated permutation contrast.

It is deliberately *uncalibrated*: the reported Z/p use a naive standard error
that assumes independent SNPs. Real LD makes the effective N smaller, so |Z| is
OPTIMISTIC -- treat it as a screen, not a result. Two nuisance asymmetries it does
NOT remove (hemizygous 0/2 coding on non-PAR X, and differential X LD/ancestry)
are exactly what the follow-up ancestry-matched permutation contrast is for. So:
read this as "is there an X excess worth a calibrated test?", nothing stronger.

Metrics, per chromosome class:
  * n             : interaction tests retained
  * lambda_GC     : median(chi2) / 0.4549  (robust, median-based; less powerful)
  * mean_chi2     : the aggregate-inflation metric (this is the one to compare)
  * mean_chi2_bg  : mean_chi2 after dropping genome-wide-significant SNPs
                    (P < 5e-8) -- the *diffuse/polygenic* background. If the X
                    excess is in mean_chi2 but vanishes in mean_chi2_bg, it is
                    driven by a few loci, not a pervasive X sensitivity.

chi2 is derived from REGENIE's LOG10P exactly as top_interactions.py does
(P = 10**(-LOG10P); chi2 = isf(P, df=1)), so lambda_GC here matches the value
top_interactions.py already prints for each scan.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

GW_SIG_P = 5e-8
CHI2_MED = stats.chi2.ppf(0.5, df=1)  # 0.4549..., the lambda_GC denominator


def _find(results_dir, pheno, prefix, explicit):
    """Resolve a regenie scan path: explicit > <dir>/<prefix>_<pheno>.regenie > glob."""
    if explicit:
        return explicit
    cand = os.path.join(results_dir, "%s_%s.regenie" % (prefix, pheno))
    if os.path.exists(cand):
        return cand
    hits = sorted(glob.glob(os.path.join(results_dir, "%s_*.regenie" % prefix)))
    # gxhap_* also matches gxhapX_* by prefix 'gxhap'; keep only exact stem.
    hits = [h for h in hits
            if os.path.basename(h).split("_")[0] == prefix.split("_")[0]
            and os.path.basename(h).startswith(prefix + "_")]
    return hits[0] if hits else cand  # return cand so the error message is informative


def load_chi2(path, test_label, min_maf, chunksize=2_000_000):
    """Return (chi2 array, n, maf_dropped, p array) for the interaction TEST rows.

    Reads only TEST/LOG10P(/A1FREQ) with the C engine, in chunks, so peak memory
    is bounded by one chunk regardless of how big the genome-wide scan is. (The
    naive full-frame read with engine='python' uses tens of GB on a real scan --
    millions of rows x ~4 TEST labels per variant x 13 object columns.)
    """
    if not os.path.exists(path):
        raise SystemExit(
            "ERROR: scan not found: %s\n"
            "  (autosomal scan = target `interaction`; chrX scan = target "
            "`interaction_X`, which needs XBFILE set and the chrX arm run.)" % path)
    head = pd.read_csv(path, sep=r"\s+", engine="c", comment="#", nrows=0)
    cols = list(head.columns)
    if "TEST" not in cols or "LOG10P" not in cols:
        raise SystemExit(
            "ERROR: %s lacks TEST/LOG10P columns -- is this a REGENIE --interaction output?" % path)
    use_maf = min_maf > 0 and "A1FREQ" in cols
    usecols = ["TEST", "LOG10P"] + (["A1FREQ"] if use_maf else [])

    p_parts, maf_dropped, seen = [], 0, set()
    reader = pd.read_csv(path, sep=r"\s+", engine="c", comment="#",
                         usecols=usecols, chunksize=chunksize)
    for chunk in reader:
        seen.update(chunk["TEST"].dropna().unique().tolist())
        sub = chunk[chunk["TEST"] == test_label]
        if sub.empty:
            continue
        log10p = pd.to_numeric(sub["LOG10P"], errors="coerce")
        ok = log10p.notna()
        sub, log10p = sub[ok], log10p[ok]
        if use_maf:
            f = pd.to_numeric(sub["A1FREQ"], errors="coerce")
            maf = np.minimum(f, 1.0 - f)
            keep = (maf >= min_maf).to_numpy()
            maf_dropped += int((~keep).sum())
            log10p = log10p[keep]
        p_parts.append(np.power(10.0, -log10p.to_numpy()))

    if not p_parts:
        raise SystemExit(
            "ERROR: no rows for TEST=%s in %s.\n  available labels: %s\n"
            "  (the interaction label varies by REGENIE version; pass --test)"
            % (test_label, path, sorted(seen)))
    p = np.clip(np.concatenate(p_parts), 1e-300, 1.0)
    chi2 = stats.chi2.isf(p, df=1)
    return chi2, len(p), maf_dropped, p


def summarize(chi2, p):
    n = len(chi2)
    lam = float(np.median(chi2) / CHI2_MED)
    mean_chi2 = float(np.mean(chi2))
    se_mean = float(np.std(chi2, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    bg = chi2[p >= GW_SIG_P]
    mean_bg = float(np.mean(bg)) if len(bg) else float("nan")
    n_sig = int((p < GW_SIG_P).sum())
    return dict(n=n, lam=lam, mean=mean_chi2, se=se_mean,
                mean_bg=mean_bg, n_sig=n_sig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results",
                    help="directory holding the regenie scans (default: results)")
    ap.add_argument("--pheno", default=os.environ.get("YS_PHENONAME", "autism"),
                    help="phenotype name in the filenames (default: $YS_PHENONAME or 'autism')")
    ap.add_argument("--auto", default="",
                    help="autosomal interaction scan (default: <dir>/gxhap_<pheno>.regenie)")
    ap.add_argument("--x", default="",
                    help="chrX interaction scan (default: <dir>/gxhapX_<pheno>.regenie)")
    ap.add_argument("--test", default="ADD-INT_SNPxHap",
                    help="interaction TEST label (verify in your REGENIE header)")
    ap.add_argument("--min-maf", type=float, default=0.01,
                    help="drop variants below this MAF on BOTH classes for comparability "
                         "(default 0.01; 0 disables). Rare hemizygous-X chi2 is the most "
                         "coding-artefact-prone, so a floor keeps the contrast fair.")
    ap.add_argument("--neff-frac", type=float, default=1.0,
                    help="multiply N by this to approximate LD-reduced effective N for a "
                         "conservative Z (e.g. 0.1). Default 1.0 = naive (optimistic).")
    ap.add_argument("--out", default="",
                    help="optional path to also write the report (default: stdout only)")
    a = ap.parse_args()

    auto_path = _find(a.results_dir, a.pheno, "gxhap", a.auto)
    x_path = _find(a.results_dir, a.pheno, "gxhapX", a.x)

    c_auto, n_auto, drop_auto, p_auto = load_chi2(auto_path, a.test, a.min_maf)
    c_x, n_x, drop_x, p_x = load_chi2(x_path, a.test, a.min_maf)
    A = summarize(c_auto, p_auto)
    X = summarize(c_x, p_x)

    # difference in mean chi2 with a naive (independent-SNP) standard error.
    d = X["mean"] - A["mean"]
    # apply optional LD deflation to the effective N behind each SE.
    se_auto = A["se"] / np.sqrt(a.neff_frac) if a.neff_frac > 0 else float("nan")
    se_x = X["se"] / np.sqrt(a.neff_frac) if a.neff_frac > 0 else float("nan")
    se_d = float(np.sqrt(se_x ** 2 + se_auto ** 2))
    z = d / se_d if se_d > 0 else float("nan")
    pval = float(2.0 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
    lam_ratio = X["lam"] / A["lam"] if A["lam"] else float("nan")

    L = []
    w = L.append
    w("=" * 70)
    w("chrX vs autosome -- aggregate SNP x Hap interaction inflation")
    w("=" * 70)
    w("interaction TEST : %s" % a.test)
    w("MAF floor        : %s%s" % (
        ("%.3g" % a.min_maf) if a.min_maf > 0 else "none",
        ("  (dropped auto=%d, X=%d)" % (drop_auto, drop_x)) if a.min_maf > 0 else ""))
    w("autosomal scan   : %s" % auto_path)
    w("chrX scan        : %s" % x_path)
    w("")
    hdr = "%-12s %10s %10s %12s %14s %8s"
    w(hdr % ("class", "n", "lambda_GC", "mean_chi2", "mean_chi2_bg", "GWsig"))
    w(hdr % ("-" * 12, "-" * 10, "-" * 10, "-" * 12, "-" * 14, "-" * 8))
    row = "%-12s %10d %10.4f %12.4f %14.4f %8d"
    w(row % ("autosome", A["n"], A["lam"], A["mean"], A["mean_bg"], A["n_sig"]))
    w(row % ("chrX", X["n"], X["lam"], X["mean"], X["mean_bg"], X["n_sig"]))
    w("")
    w("contrast (chrX - autosome):")
    w("  delta mean_chi2 : %+.4f" % d)
    w("  lambda_GC ratio : %.3f  (X / autosome)" % lam_ratio)
    if a.neff_frac != 1.0:
        w("  naive SE(delta) : %.4f   [N deflated x %.3g for LD]" % (se_d, a.neff_frac))
    else:
        w("  naive SE(delta) : %.4f   [independent-SNP SE -- OPTIMISTIC]" % se_d)
    w("  Z               : %+.2f" % z)
    w("  two-sided p     : %.3g" % pval)
    w("")

    # verdict
    excess_mean = d > 0
    excess_lam = X["lam"] > A["lam"]
    excess_diffuse = np.isfinite(X["mean_bg"]) and np.isfinite(A["mean_bg"]) and X["mean_bg"] > A["mean_bg"]
    w("VERDICT")
    if not excess_mean and not excess_lam:
        w("  NO chrX excess. mean_chi2 and lambda_GC are both <= the autosome.")
        w("  --> No aggregate X sensitivity to chase; skip the calibrated contrast.")
    else:
        bits = []
        if excess_mean:
            bits.append("mean_chi2 +%.4f" % d)
        if excess_lam:
            bits.append("lambda_GC x%.3f" % lam_ratio)
        w("  chrX shows EXCESS interaction inflation (%s)." % ", ".join(bits))
        if excess_diffuse:
            w("  The excess persists after dropping GW-sig SNPs (mean_chi2_bg) --")
            w("  consistent with a diffuse/polygenic X sensitivity, not a few loci.")
        else:
            w("  But it weakens once GW-sig SNPs are dropped (mean_chi2_bg) -- the")
            w("  excess may be driven by a handful of loci rather than pervasive.")
        w("  --> Worth building the calibrated mean_chi2 + ancestry-matched")
        w("      permutation contrast. Note: naive Z is OPTIMISTIC (LD), and it")
        w("      does NOT yet remove hemizygous-coding or differential-X-LD/ancestry")
        w("      asymmetries -- that is what the permutation step is for.")
    w("=" * 70)

    report = "\n".join(L)
    print(report)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(report + "\n")
        print("\nwrote %s" % a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
