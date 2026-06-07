#!/usr/bin/env python3
"""Normalize raw *_frequency_all.txt hotspot files to 0-1 frequencies.

Use this only if the raw frequency files were generated with the final paper cutoff.
For the paper hotspot maps, the recommended denominator is the number of analyzed
residue observations used for each residue class. If all copies/protomers/frames are
pooled uniformly, provide that value with --denominator.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

def parse_line(line: str):
    key, val = line.split(":", 1)
    value = float(val.strip())
    key = key.strip()
    if "|" in key:
        resname, rest = key.split("|", 1)
    else:
        resname, rest = "", key
    if "_" in rest:
        segid, resid = rest.split("_", 1)
    else:
        segid, resid = rest, ""
    return resname, segid, resid, value

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--denominator", type=float, required=True, help="Number of analyzed frame/copy observations used for normalization")
    args = parser.parse_args()

    inp = Path(args.input_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for file in sorted(inp.glob("*_frequency_all.txt")):
        component = file.name.replace("_frequency_all.txt", "")
        out_csv = out / f"{component}_normalized_contact_frequency.csv"
        with file.open() as fh, out_csv.open("w", newline="") as outfh:
            writer = csv.writer(outfh)
            writer.writerow(["Component", "Resname", "Segid", "Resid", "Raw contact count", "Denominator", "Normalized frequency"])
            for line in fh:
                if not line.strip() or ":" not in line:
                    continue
                resname, segid, resid, count = parse_line(line)
                writer.writerow([component, resname, segid, resid, count, args.denominator, count / args.denominator])
        print(f"Wrote {out_csv}")

if __name__ == "__main__":
    main()
