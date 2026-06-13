#!/usr/bin/env python3
"""Pool the raw counts from independent permutation batches into the final
results, for either the interaction re-test or the LAVA local-r_g test.

Each batch (run with --raw-counts, fixed selection seed, distinct permutation
seed) emits per-unit counts: obs statistic, n_ge = #(perm >= obs), n_perm. Pooling
sums n_ge and n_perm across batches and forms the empirical p; for the genome-wide
lambda it concatenates the per-batch null-lambda samples.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from local_rg_perm import write_local_rg_summary   # noqa: E402


def pool_interaction(counts, lambda_nulls, out_prefix):
    df = pd.concat([pd.read_csv(c, sep="\t") for c in counts], ignore_index=True)
    g = df.groupby("ID", sort=False)
    out = g.agg(raw_LOG10P=("raw_LOG10P", "first"), obs_stat=("obs_stat", "first"),
                n_ge=("n_ge", "sum"), n_perm=("n_perm", "sum"),
                forced=("forced", "first")).reset_index()
    out["anc_matched_emp_p"] = (1 + out.n_ge) / (1 + out.n_perm)
    out = out[["ID", "raw_LOG10P", "obs_stat", "n_perm", "anc_matched_emp_p", "forced"]]
    out.sort_values("anc_matched_emp_p").to_csv(
        out_prefix + "_interactions.tsv", sep="\t", index=False)

    meta, nulls = {}, []
    for k, path in enumerate(lambda_nulls):
        after = False
        for line in open(path):
            line = line.rstrip("\n")
            if line == "#NULL":
                after = True
                continue
            if after:
                if line:
                    nulls.append(float(line))
            elif k == 0 and "\t" in line:
                key, val = line.split("\t", 1)
                meta[key] = val
    nulls = np.array(nulls)
    obs_lambda = float(meta.get("obs_lambda", "nan"))
    with open(out_prefix + "_lambda.txt", "w") as f:
        f.write("panel_snps\t%s\n" % meta.get("panel_snps", "0"))
        if len(nulls) and np.isfinite(obs_lambda):
            p = (1 + int((nulls >= obs_lambda).sum())) / (1 + len(nulls))
            f.write("observed_lambda\t%.4f\n" % obs_lambda)
            f.write("perm_null_lambda_mean\t%.4f\n" % nulls.mean())
            f.write("perm_null_lambda_2.5pct\t%.4f\n" % np.percentile(nulls, 2.5))
            f.write("perm_null_lambda_97.5pct\t%.4f\n" % np.percentile(nulls, 97.5))
            f.write("lambda_emp_p\t%.4g\n" % p)
        f.write("perm_null_lambda_samples\t%d\n" % len(nulls))
        for key in ("selected_stratum_size", "within_stratum_HapPC_AUC",
                    "n_mixed_strata", "eff_N"):
            if key in meta:
                f.write("%s\t%s\n" % (key, meta[key]))
    print("pooled %d batches -> %s_interactions.tsv (%d hits, %d total perms each)"
          % (len(counts), out_prefix, len(out), int(out.n_perm.iloc[0]) if len(out) else 0))


def pool_lava(counts, meta_file, out_prefix):
    df = pd.concat([pd.read_csv(c, sep="\t") for c in counts], ignore_index=True)
    g = df.groupby("locus", sort=False)
    out = g.agg(chr=("chr", "first"), n_snps=("n_snps", "first"),
                local_rg=("local_rg", "first"), n_ge=("n_ge", "sum"),
                n_perm=("n_perm", "sum"), is_target=("is_target", "first")).reset_index()
    out["anc_matched_emp_p"] = (1 + out.n_ge) / (1 + out.n_perm)
    out[["locus", "chr", "n_snps", "local_rg", "anc_matched_emp_p",
         "is_target"]].to_csv(out_prefix + "_loci.tsv", sep="\t", index=False)
    m = dict(l.split("\t", 1) for l in open(meta_file).read().splitlines() if "\t" in l)
    write_local_rg_summary(out_prefix, out, int(m["strata_mean_size"]),
                           float(m["within_stratum_HapPC_AUC"]))
    print("pooled %d batches -> %s_loci.tsv (%d total perms/locus)"
          % (len(counts), out_prefix, int(out.n_perm.iloc[0]) if len(out) else 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["interaction", "lava"])
    ap.add_argument("--counts", nargs="+", required=True, help="per-batch *_counts.tsv")
    ap.add_argument("--lambda-nulls", nargs="*", default=[],
                    help="per-batch *_lambda_null.tsv (interaction only)")
    ap.add_argument("--meta", help="one *_meta.tsv (lava only)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    if a.kind == "interaction":
        pool_interaction(a.counts, a.lambda_nulls, a.out_prefix)
    else:
        pool_lava(a.counts, a.meta, a.out_prefix)


if __name__ == "__main__":
    main()
