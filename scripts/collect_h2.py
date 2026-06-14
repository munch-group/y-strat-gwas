#!/usr/bin/env python3
"""Parse one or more LDSC --h2 .log files into a single tidy table, for comparing
heritability across groups (e.g. pooled vs per-Y-haplogroup) or across
chromosomes (e.g. chrX vs the autosomes).
"""
import argparse
import re

_PATTERNS = {
    "h2": r"Total (?:Liability|Observed) scale h2: ([-\d.eE]+) \(([-\d.eE]+)\)",
    "intercept": r"Intercept: ([-\d.eE]+) \(([-\d.eE]+)\)",
    "lambda_gc": r"Lambda GC: ([-\d.eE]+)",
    "mean_chi2": r"Mean Chi\^2: ([-\d.eE]+)",
    "n_snps": r"After merging with reference panel LD, (\d+) SNPs remain",
}


def parse_log(path):
    txt = open(path).read()
    row = {}
    m = re.search(_PATTERNS["h2"], txt)
    row["h2"], row["h2_se"] = (m.group(1), m.group(2)) if m else ("NA", "NA")
    row["scale"] = "liability" if "Liability scale" in txt else (
        "observed" if "Observed scale" in txt else "NA")
    m = re.search(_PATTERNS["intercept"], txt)
    row["intercept"], row["intercept_se"] = (m.group(1), m.group(2)) if m else ("NA", "NA")
    for key in ("lambda_gc", "mean_chi2", "n_snps"):
        m = re.search(_PATTERNS[key], txt)
        row[key] = m.group(1) if m else "NA"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True,
                    help="one label per log (e.g. full I R, or 1 2 .. 22 X)")
    ap.add_argument("--label-col", default="group", help="name of the label column")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if len(a.logs) != len(a.labels):
        raise SystemExit("need one --labels entry per --logs file")

    cols = [a.label_col, "scale", "h2", "h2_se", "intercept", "intercept_se",
            "lambda_gc", "mean_chi2", "n_snps"]
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for label, log in zip(a.labels, a.logs):
            r = parse_log(log)
            r[a.label_col] = label
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote %s: %d groups" % (a.out, len(a.logs)))


if __name__ == "__main__":
    main()
