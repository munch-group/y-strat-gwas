#!/usr/bin/env python3
"""Concatenate per-chromosome REGENIE step-2 outputs into one file, keeping the
comment (#...) block and header from the first chunk only. Used to gather a
chromosome-split scan back into the single file the downstream targets expect.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--inputs", nargs="+", required=True)
    a = ap.parse_args()

    n_data = 0
    with open(a.out, "w") as o:
        for fi, path in enumerate(a.inputs):
            seen_header = False
            with open(path) as f:
                for line in f:
                    if line.startswith("#"):           # comment block
                        if fi == 0:
                            o.write(line)
                        continue
                    if not seen_header:                # header line
                        seen_header = True
                        if fi == 0:
                            o.write(line)
                        continue
                    o.write(line)                      # data row
                    n_data += 1
    print("wrote %s: %d data rows from %d chunks" % (a.out, n_data, len(a.inputs)))


if __name__ == "__main__":
    main()
