#!/usr/bin/env python3
"""Read the result files of every arm and write a plain-language conclusion of
what each arm does and does NOT support, applying the interpretation rules from
METHODS.md / OUTPUTS.md. Missing files (disabled or unfinished arms) are reported
as "not available" rather than crashing, so this runs at the end of the pipeline
regardless of which optional arms are on.

It states conclusions; it does not re-do statistics. Every claim is hedged the
way the docs require -- in particular, a surviving interaction / local-rg signal
is always flagged provisional (global-PC ceiling: local ancestry is never removed).
"""
import argparse
import os
import re
import math

GW_SIG_LOG10P = 7.30103   # p = 5e-8
ALPHA = 0.05              # permutation survival / female-artefact threshold
H2_Z_MIN = 4.0            # per-stratum h2 Z below which rg is hopeless


# ----------------------------------------------------------------------------- IO
def path(d, name):
    return os.path.join(d, name)


def have(p):
    return p is not None and os.path.isfile(p) and os.path.getsize(p) > 0


def read_kv(p):
    """key<TAB>value file -> dict (values left as strings)."""
    out = {}
    for line in open(p):
        if "\t" in line and not line.startswith("#"):
            k, v = line.rstrip("\n").split("\t", 1)
            out[k] = v
    return out


def read_tsv(p):
    import pandas as pd
    return pd.read_csv(p, sep="\t", engine="python")


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def h2_from_log(p):
    """Pull h2, se, mean_chi2, n_snps from an LDSC --h2 .log."""
    txt = open(p).read()
    m = re.search(r"Total (?:Liability|Observed) scale h2: ([-\d.eE]+) \(([-\d.eE]+)\)", txt)
    h2, se = (fnum(m.group(1)), fnum(m.group(2))) if m else (float("nan"), float("nan"))
    mc = re.search(r"Mean Chi\^2: ([-\d.eE]+)", txt)
    ns = re.search(r"After merging with reference panel LD, (\d+) SNPs remain", txt)
    return {"h2": h2, "se": se,
            "mean_chi2": fnum(mc.group(1)) if mc else float("nan"),
            "n_snps": int(ns.group(1)) if ns else None}


def rg_from_log(p):
    """Pull rg, se, p from the LDSC --rg summary block (p1 p2 rg se z p ...)."""
    lines = open(p).read().splitlines()
    for k, line in enumerate(lines):
        if "Summary of Genetic Correlation Results" in line:
            rows = [ln for ln in lines[k + 1:k + 6] if ln.strip()]
            if len(rows) >= 2:
                f = rows[1].split()
                # columns: p1 p2 rg se z p h2_obs h2_obs_se ...
                if len(f) >= 6:
                    return {"rg": fnum(f[2]), "se": fnum(f[3]), "p": fnum(f[5])}
            break
    return None


def zscore(est, se):
    if se and se == se and se > 0:
        return est / se
    return float("nan")


# ------------------------------------------------------------------- arm reports
def arm_a(d, lines):
    lines.append("=" * 72)
    lines.append("ARM A -- individual SNP x Y-haplogroup interactions")
    lines.append("=" * 72)
    top = path(d, "top_interactions.tsv")
    perm = path(d, "perm_interactions.tsv")
    plam = path(d, "perm_lambda.txt")

    if not have(top):
        lines.append("  [not available] top_interactions.tsv missing -- arm did not run.")
        return
    df = read_tsv(top)
    sig = df[df["LOG10P"] >= GW_SIG_LOG10P] if "LOG10P" in df.columns else df.iloc[0:0]
    n_sig = len(sig)
    lines.append("  Genome-wide significant interaction hits (p < 5e-8): %d" % n_sig)
    if n_sig:
        s = sig.sort_values("LOG10P", ascending=False).head(5)
        for _, r in s.iterrows():
            lines.append("    - %s  LOG10P=%.2f" % (r.get("ID", "?"), r["LOG10P"]))

    # genome-wide inflation vs the within-strata permutation null
    survivors = []
    if have(plam):
        kv = read_kv(plam)
        obs = fnum(kv.get("observed_lambda"))
        hi = fnum(kv.get("perm_null_lambda_97.5pct"))
        ep = fnum(kv.get("lambda_emp_p"))
        if obs == obs:
            lines.append("  Genome-wide interaction lambda = %.3f (within-strata "
                         "permutation null mean %.3f, 97.5%%=%.3f, emp p=%s)"
                         % (obs, fnum(kv.get("perm_null_lambda_mean")), hi,
                            kv.get("lambda_emp_p", "NA")))
            if ep == ep and ep <= ALPHA:
                lines.append("    => the genome-wide interaction signal EXCEEDS the "
                             "ancestry-matched null: evidence for pervasive interaction "
                             "beyond ancestry.")
            else:
                lines.append("    => inflation is consistent with the ancestry null: no "
                             "pervasive interaction beyond ancestry is supported.")
    else:
        lines.append("  [permutation lambda not available]")

    # which significant hits survive the ancestry-matched permutation
    perm_available = False
    if have(perm) and n_sig:
        pm = read_tsv(perm)
        if "ID" in pm.columns and "anc_matched_emp_p" in pm.columns:
            perm_available = True
            sig_ids = set(sig["ID"]) if "ID" in sig.columns else set()
            surv = pm[(pm["ID"].isin(sig_ids)) & (pm["anc_matched_emp_p"] < ALPHA)]
            survivors = list(surv["ID"])
            lines.append("  Of the %d significant hits, %d survive the ancestry-matched "
                         "permutation (anc_matched_emp_p < %.2f):" % (n_sig, len(surv), ALPHA))
            for _, r in surv.sort_values("anc_matched_emp_p").head(10).iterrows():
                lines.append("    - %s  anc_matched_emp_p=%.4g" % (r["ID"], r["anc_matched_emp_p"]))
    elif n_sig:
        lines.append("  [perm_interactions.tsv not available -- hits not deconfounded]")

    lines.append("  CONCLUSION:")
    if n_sig == 0:
        lines.append("    No individual variant reaches genome-wide significance for")
        lines.append("    SNP x haplogroup interaction. No variant-level claim.")
    elif survivors:
        lines.append("    %d variant(s) are genome-wide significant AND survive the" % len(survivors))
        lines.append("    ancestry-matched permutation: %s." % ", ".join(map(str, survivors[:10])))
        lines.append("    PROVISIONAL -- the permutation only controls ancestry to global-PC")
        lines.append("    resolution; local ancestry is not removed. Confirm with local-")
        lines.append("    ancestry calls at the locus or replication.")
    elif not perm_available:
        lines.append("    %d significant hit(s) exist but the ancestry-matched permutation" % n_sig)
        lines.append("    is not available -- they are NOT yet deconfounded. Run perm_interaction.")
    else:
        lines.append("    Significant hit(s) exist but none survive the ancestry-matched")
        lines.append("    permutation -> consistent with ancestry confounding, not a")
        lines.append("    haplogroup interaction.")

    # chrX
    topx = path(d, "top_interactions_X.tsv")
    if have(topx):
        dx = read_tsv(topx)
        nx = int((dx["LOG10P"] >= GW_SIG_LOG10P).sum()) if "LOG10P" in dx.columns else 0
        lines.append("  chrX: %d genome-wide significant SNP x Hap interaction hit(s)." % nx)


def arm_b(d, lines):
    lines.append("")
    lines.append("=" * 72)
    lines.append("ARM B -- cross-stratum genetic correlation & heritability")
    lines.append("=" * 72)
    hI, hR = path(d, "h2_I.log"), path(d, "h2_R.log")
    rgp = path(d, "I_vs_R_rg.log")

    hi = h2_from_log(hI) if have(hI) else None
    hr = h2_from_log(hR) if have(hR) else None
    if hi:
        lines.append("  h2_I (liability) = %.4f (%.4f), Z=%.2f, mean chi2=%.3f"
                     % (hi["h2"], hi["se"], zscore(hi["h2"], hi["se"]), hi["mean_chi2"]))
    if hr:
        lines.append("  h2_R (liability) = %.4f (%.4f), Z=%.2f, mean chi2=%.3f"
                     % (hr["h2"], hr["se"], zscore(hr["h2"], hr["se"]), hr["mean_chi2"]))

    # --- cross-stratum rg ---
    rg = rg_from_log(rgp) if have(rgp) else None
    lines.append("  CONCLUSION (rg, the pervasive-interaction readout):")
    if not rg or rg["rg"] != rg["rg"]:
        lines.append("    rg not estimable (LDSC returned nan -- h2 out of bounds / no signal).")
    else:
        rgv, se = rg["rg"], rg["se"]
        zI = zscore(hi["h2"], hi["se"]) if hi else float("nan")
        zR = zscore(hr["h2"], hr["se"]) if hr else float("nan")
        weak_h2 = (zI == zI and zI < H2_Z_MIN) or (zR == zR and zR < H2_Z_MIN)
        hi_ci = rgv + 1.96 * se if se == se else float("nan")
        lines.append("    rg = %.3f (se %.3f, p=%s)" % (rgv, se, rg["p"]))
        if se != se or se > 0.25 or weak_h2:
            why = []
            if se != se or se > 0.25:
                why.append("se >> 0 (the cross-stratum covariance is too weak)")
            if weak_h2:
                why.append("a per-stratum h2 Z < %g (weak denominator)" % H2_Z_MIN)
            lines.append("    => UNDERPOWERED / uninformative: %s." % "; ".join(why))
            lines.append("    Cannot distinguish rg < 1 (interaction) from rg = 1 (none).")
            lines.append("    Lead Arm B with the h2 comparison below, not with rg.")
        elif hi_ci == hi_ci and hi_ci < 1.0:
            lines.append("    => rg significantly BELOW 1: evidence for pervasive G x")
            lines.append("    haplogroup interaction (the well-powered general-pattern readout).")
        else:
            lines.append("    => rg consistent with 1 and well estimated: no evidence for")
            lines.append("    pervasive interaction (and powered to have seen it).")

    # --- h2 difference between strata ---
    if hi and hr and hi["se"] == hi["se"] and hr["se"] == hr["se"]:
        diff = hr["h2"] - hi["h2"]
        sed = math.sqrt(hi["se"] ** 2 + hr["se"] ** 2)
        zd = zscore(diff, sed)
        lines.append("  Heritability by stratification: h2_R - h2_I = %.4f (se %.4f, Z=%.2f)"
                     % (diff, sed, zd))
        verdict = ("differ" if abs(zd) >= 1.96 else
                   "are suggestive but not significantly different" if abs(zd) >= 1.64 else
                   "do not significantly differ")
        lines.append("    => per-stratum heritabilities %s." % verdict)
        lines.append("    CAVEAT: an h2 difference is itself ancestry-confounded (I and R have")
        lines.append("    different LD, fit against one EUR reference) and depends on per-stratum")
        lines.append("    case counts -- treat as a flag, not a result.")

    if have(path(d, "h2_by_chromosome.tsv")):
        lines.append("  Per-chromosome h2: see h2_by_chromosome.tsv (chrX vs autosomes).")


def females(d, lines):
    lines.append("")
    lines.append("=" * 72)
    lines.append("FEMALES negative control (no Y -> a real Y interaction cannot exist)")
    lines.append("=" * 72)
    nc, lam = path(d, "female_negative_control.tsv"), path(d, "female_lambda.txt")
    if not have(nc):
        lines.append("  [not run] (FBFILE unset) -- no female negative control.")
        return
    df = read_tsv(nc)
    flagged = df[df["looks_like_ancestry_artifact"] == True] if \
        "looks_like_ancestry_artifact" in df.columns else df.iloc[0:0]
    n = len(df)
    lines.append("  Male interaction hits re-tested in females: %d; flagged as ancestry "
                 "artefact: %d" % (n, len(flagged)))
    lines.append("  CONCLUSION:")
    if len(flagged):
        ids = list(flagged["ID"]) if "ID" in flagged.columns else []
        lines.append("    %d hit(s) REPRODUCE in females (no Y): %s"
                     % (len(flagged), ", ".join(map(str, ids[:10]))))
        lines.append("    => those are ancestry artefacts, NOT Y-driven -- kill them.")
    else:
        lines.append("    No male hit reproduces in females. This is *consistent with* a")
        lines.append("    Y-driven interaction but is NOT proof (autism is sex-differential,")
        lines.append("    so a hit can be null in females for sex reasons).")


def lava(d, lines):
    lines.append("")
    lines.append("=" * 72)
    lines.append("LAVA -- local cross-stratum genetic correlation")
    lines.append("=" * 72)
    screen = path(d, "lava_local_rg.tsv")
    if have(screen):
        df = read_tsv(screen)
        if "rho.upper" in df.columns and "rho" in df.columns:
            cand = df[df["rho.upper"] < 1.0]
            lines.append("  Genome-wide screen: %d block(s) of %d have rho.upper < 1 "
                         "(local divergence candidates)." % (len(cand), len(df)))
            for _, r in cand.sort_values("rho").head(5).iterrows():
                loc = r.get("locus", "%s:%s-%s" % (r.get("chr"), r.get("start"), r.get("stop")))
                lines.append("    - %s  rho=%.3f  (CI %.3f..%.3f)"
                             % (loc, r["rho"], r.get("rho.lower", float("nan")), r["rho.upper"]))
            lines.append("    NOTE: this screen is blind to ancestry -- nominees must be")
            lines.append("    deconfounded (list them in LAVA_LOCI) before being believed.")
        else:
            lines.append("  lava_local_rg.tsv present but no rho columns (LAVA produced no")
            lines.append("  bivariate estimate -- univariate signal likely too weak).")
    else:
        lines.append("  [genome-wide screen not run] (LAVA_PARTITION unset or R/LAVA absent).")

    # per-block deconfounding summaries
    import glob
    summaries = sorted(glob.glob(path(d, "lava_*_summary.txt")))
    if summaries:
        lines.append("  Per-block deconfounding:")
        for s in summaries:
            kv = read_kv(s)
            nm = os.path.basename(s)[len("lava_"):-len("_summary.txt")]
            rg = fnum(kv.get("target_local_rg"))
            ap_ = fnum(kv.get("target_anc_matched_p"))
            tf = fnum(kv.get("target_vs_controls_tail_frac"))
            ok = (ap_ == ap_ and ap_ < ALPHA) and (tf == tf and tf < ALPHA)
            verdict = ("SURVIVES (provisional -- global-PC ceiling)" if ok else
                       "does NOT stand out -- consistent with genome-wide residual structure")
            lines.append("    - %s: local_rg=%.3f, anc_matched_p=%s, tail_frac=%s => %s"
                         % (nm, rg, kv.get("target_anc_matched_p", "NA"),
                            kv.get("target_vs_controls_tail_frac", "NA"), verdict))
    else:
        lines.append("  [no per-block deconfounding] (LAVA_LOCI unset).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = a.results_dir

    L = []
    L.append("#" * 72)
    L.append("# CONCLUSIONS -- what each arm of the analysis supports")
    L.append("# (auto-generated from the result files; see METHODS.md / OUTPUTS.md)")
    L.append("# Central caveat: haplogroup I vs R tags autosomal ancestry, so every")
    L.append("# positive is provisional to the global-PC resolution of the controls.")
    L.append("#" * 72)
    L.append("")
    arm_a(d, L)
    arm_b(d, L)
    females(d, L)
    lava(d, L)
    L.append("")

    text = "\n".join(L) + "\n"
    with open(a.out, "w") as f:
        f.write(text)
    print(text)
    print("wrote %s" % a.out)


if __name__ == "__main__":
    main()
