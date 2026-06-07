#!/usr/bin/env python3
"""Plot preference-index bar charts using the paper-consistent global PI column."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from rav_components import COMPONENT_COLORS, DEFAULT_PI_COMPONENTS

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True)
    parser.add_argument("--output-dir", default="preference_index_plots")
    parser.add_argument("--pi-column", default="Pref Index (Global)")
    parser.add_argument("--ylim", default=None, help="Optional y-limits, e.g. 0,110")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    order = DEFAULT_PI_COMPONENTS
    ylim = tuple(float(x) for x in args.ylim.split(',')) if args.ylim else None

    for file in args.input:
        df = pd.read_csv(file)
        pi_col = args.pi_column if args.pi_column in df.columns else "Pref Index"
        df = df[df["Component"].isin(order)].copy()
        df["Component"] = pd.Categorical(df["Component"], categories=order, ordered=True)
        df = df.sort_values("Component")
        colors = [COMPONENT_COLORS.get(c, (0.5, 0.5, 0.5)) for c in df["Component"]]

        plt.figure(figsize=(10, 5))
        plt.bar(df["Component"].astype(str), df[pi_col], color=colors)
        plt.ylabel("Preference Index")
        plt.xticks(rotation=45, ha="right")
        if ylim:
            plt.ylim(*ylim)
        plt.tight_layout()
        out = outdir / f"{Path(file).stem}_preference_index.png"
        plt.savefig(out, dpi=300)
        plt.close()
        print(f"Wrote {out}")

if __name__ == "__main__":
    main()
