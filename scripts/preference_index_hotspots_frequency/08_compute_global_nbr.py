#!/usr/bin/env python3
"""Compute the pooled global normalized binding ratio (global NBR)."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

from rav_components import COMPONENT_TOTALS, DEFAULT_PI_COMPONENTS, write_component_totals_csv

def read_stats(path: Path) -> pd.DataFrame:
    with path.open() as fh:
        first = fh.readline()
    if "," in first and "Component" in first:
        return pd.read_csv(path)
    return pd.read_csv(path, skiprows=1)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-files", nargs="+", required=True, help="Overall stats CSVs for all reference systems")
    parser.add_argument("--output", default="global_nbr_summary.csv")
    parser.add_argument("--component-totals-output", default="component_totals.csv")
    parser.add_argument("--components", default=",".join(DEFAULT_PI_COMPONENTS), help="Comma-separated components to include. Default excludes D/Q/T mucin subtypes to avoid double counting.")
    args = parser.parse_args()

    include = {c.strip() for c in args.components.split(',') if c.strip()}
    rows = []
    total_mean_sum = 0.0
    total_available_sum = 0.0

    for file in args.input_files:
        path = Path(file)
        df = read_stats(path)
        for _, row in df.iterrows():
            comp = str(row["Component"]).strip()
            if comp not in include:
                continue
            if comp not in COMPONENT_TOTALS:
                print(f"WARNING: no total available for {comp}; skipping")
                continue
            mean = float(row["Mean"])
            total = float(COMPONENT_TOTALS[comp])
            nbr = mean / total if total else 0.0
            rows.append({"input_file": path.name, "Component": comp, "Mean": mean, "Total Available": total, "NBR": nbr})
            total_mean_sum += mean
            total_available_sum += total

    if not rows or total_available_sum == 0:
        raise SystemExit("No valid rows found for global NBR calculation")
    global_nbr = total_mean_sum / total_available_sum

    out = Path(args.output)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["input_file", "Component", "Mean", "Total Available", "NBR", "Global NBR Mean"])
        writer.writeheader()
        for r in rows:
            r["Global NBR Mean"] = global_nbr
            writer.writerow(r)
    write_component_totals_csv(args.component_totals_output)
    print(f"Global NBR Mean: {global_nbr:.16g}")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
