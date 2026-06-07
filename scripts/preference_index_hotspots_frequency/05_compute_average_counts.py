#!/usr/bin/env python3
"""Compute per-reference and overall mean/SD contact-shell counts."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

def read_counts_csv(csv_file: Path) -> pd.DataFrame:
    return pd.read_csv(csv_file)

def write_stats(stats: dict[str, dict[str, float]], outfile: Path, title: str = "") -> None:
    with outfile.open("w", newline="") as fh:
        writer = csv.writer(fh)
        if title:
            writer.writerow([title])
        writer.writerow(["Component", "Mean", "SD"])
        for comp, values in stats.items():
            writer.writerow([comp, values["mean"], values["sd"]])

def compute_stats(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    stats = {}
    for col in df.columns:
        if col == "Frame":
            continue
        values = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(dtype=float)
        stats[col] = {"mean": float(np.mean(values)), "sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0}
    return stats

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--folder-regex", required=True, help="Regex matching reference folders, e.g. '^spike\\d{2}$' or '^R\\w{3}$'")
    parser.add_argument("--output-prefix", required=True, help="Output prefix, e.g. Spikes, Albumins, D_Mucins")
    parser.add_argument("--interactomes", nargs="+", default=["6A", "20A"])
    args = parser.parse_args()

    base = Path(args.base_dir)
    folder_re = re.compile(args.folder_regex)
    folders = sorted([p for p in base.iterdir() if p.is_dir() and folder_re.match(p.name)])
    if not folders:
        raise SystemExit(f"No folders found for regex {args.folder_regex}")

    for interactome in args.interactomes:
        per_folder_means: dict[str, dict[str, float]] = {}
        for folder in folders:
            counts_file = folder / f"{folder.name}_{interactome}_molecule_counts.csv"
            if not counts_file.exists():
                print(f"WARNING: missing {counts_file}")
                continue
            df = read_counts_csv(counts_file)
            stats = compute_stats(df)
            out = folder / f"{folder.name}_{interactome}_stats.csv"
            write_stats(stats, out, f"{interactome} statistics for {folder.name}")
            per_folder_means[folder.name] = {c: s["mean"] for c, s in stats.items()}

        if not per_folder_means:
            continue
        components = list(next(iter(per_folder_means.values())).keys())
        overall = {}
        for comp in components:
            vals = np.array([d[comp] for d in per_folder_means.values() if comp in d], dtype=float)
            overall[comp] = {"mean": float(np.mean(vals)), "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0}
        out = base / f"{args.output_prefix}_{interactome}_overall_stats.csv"
        write_stats(overall, out, f"Overall {interactome} statistics for {args.output_prefix}")
        print(f"Wrote {out}")

if __name__ == "__main__":
    main()
