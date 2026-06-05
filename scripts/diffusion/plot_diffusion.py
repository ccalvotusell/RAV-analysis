#!/usr/bin/env python3
"""
plot_diffusion.py

Plot diffusion coefficients computed by compute_diffusion.py.

Example:
    python plot_diffusion.py \
        --input results/diffusion/diffusion_coefficients.csv \
        --output results/diffusion/diffusion_kde.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot diffusion coefficient distributions.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--column", default="D_A2_per_ns")
    parser.add_argument("--group-filter", default=None)
    parser.add_argument("--plot-type", choices=["kde", "hist"], default="kde")
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--title", default="Diffusion coefficient distribution")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.input)

    if args.group_filter is not None:
        data = data[data["group_id"].astype(str).str.contains(args.group_filter)]

    values = data[args.column].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()

    if len(values) == 0:
        raise ValueError("No valid diffusion coefficients found for plotting.")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 5))

    if args.plot_type == "kde" and len(values) >= 2 and np.std(values) > 0:
        kde = gaussian_kde(values)
        x_min, x_max = values.min(), values.max()
        padding = 0.05 * (x_max - x_min) if x_max > x_min else 1.0
        x = np.linspace(x_min - padding, x_max + padding, 500)
        y = kde(x)
        plt.plot(x, y, lw=2)
        plt.fill_between(x, y, alpha=0.25)
        plt.ylabel("Density")
    else:
        plt.hist(values, bins=args.bins, alpha=0.75)
        plt.ylabel("Count")

    plt.xlabel("Diffusion coefficient D (Å²/ns)")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()

    print(f"Figure written to: {output}")


if __name__ == "__main__":
    main()
