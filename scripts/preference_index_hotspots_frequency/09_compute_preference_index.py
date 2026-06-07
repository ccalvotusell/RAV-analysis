#!/usr/bin/env python3
"""Compute NBR and preference index for one reference-system overall-stats file."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rav_components import COMPONENT_TOTALS, DEFAULT_PI_COMPONENTS

def read_stats(path: Path) -> pd.DataFrame:
    with path.open() as fh:
        first = fh.readline()
    if "," in first and "Component" in first:
        return pd.read_csv(path)
    return pd.read_csv(path, skiprows=1)

def read_global_nbr(path: Path) -> float:
    df = pd.read_csv(path)
    if "Global NBR Mean" not in df.columns:
        raise ValueError(f"{path} must contain a 'Global NBR Mean' column")
    return float(df["Global NBR Mean"].dropna().iloc[0])

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Overall stats CSV, e.g. Spikes_6A_overall_stats.csv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--global-nbr-summary", default="global_nbr_summary.csv")
    parser.add_argument("--global-nbr-value", type=float, default=None)
    parser.add_argument("--components", default=",".join(DEFAULT_PI_COMPONENTS), help="Comma-separated components to include in final PI table")
    parser.add_argument("--keep-all", action="store_true", help="Keep rows not in --components if totals exist")
    args = parser.parse_args()

    include = {c.strip() for c in args.components.split(',') if c.strip()}
    global_nbr = args.global_nbr_value if args.global_nbr_value is not None else read_global_nbr(Path(args.global_nbr_summary))

    df = read_stats(Path(args.input))
    if not args.keep_all:
        df = df[df["Component"].isin(include)].copy()
    df["Total Available"] = df["Component"].map(COMPONENT_TOTALS)
    missing = df[df["Total Available"].isna()]["Component"].tolist()
    if missing:
        raise SystemExit(f"Missing totals for components: {missing}")

    df["NBR (Component)"] = df["Mean"].astype(float) / df["Total Available"].astype(float)
    df["Global NBR"] = global_nbr
    df["Pref Index (Component)"] = df["NBR (Component)"] / df["NBR (Component)"].mean()
    df["Pref Index (Global)"] = df["NBR (Component)"] / global_nbr

    total_observed = df["Mean"].sum()
    total_available = df["Total Available"].sum()
    df["Expected"] = (df["Total Available"] / total_available) * total_observed
    df["Enrichment Ratio"] = df["Mean"] / df["Expected"]
    df["SNR"] = df.apply(lambda r: r["Mean"] / r["SD"] if r["SD"] > 0 else np.nan, axis=1)
    df["CV"] = df.apply(lambda r: r["SD"] / r["Mean"] if r["Mean"] > 0 else np.nan, axis=1)
    df["Rel Contribution"] = (df["Mean"] / total_observed) / (df["Total Available"] / total_available)

    cols = ["Component", "Mean", "SD", "Total Available", "NBR (Component)", "Global NBR", "Pref Index (Component)", "Pref Index (Global)", "Expected", "Enrichment Ratio", "SNR", "CV", "Rel Contribution"]
    df[cols].to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
