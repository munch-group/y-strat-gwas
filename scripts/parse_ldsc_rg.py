#!/usr/bin/env python3
"""Print the genetic-correlation summary block from an ldsc --rg .log file."""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    a = ap.parse_args()

    lines = open(a.log).read().splitlines()
    idx = None
    for k, line in enumerate(lines):
        if "Summary of Genetic Correlation Results" in line:
            idx = k
            break
    if idx is None:
        sys.exit("no rg summary found in %s" % a.log)

    # header line + one data line per trait pair follow the banner
    for line in lines[idx + 1:idx + 5]:
        if line.strip():
            print(line)


if __name__ == "__main__":
    main()
