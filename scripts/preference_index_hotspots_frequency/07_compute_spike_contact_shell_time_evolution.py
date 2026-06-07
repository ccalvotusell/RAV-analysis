#!/usr/bin/env python3
"""Combine spike count chunks and compute time evolution mean/SD across spikes."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from rav_components import COUNT_COLUMNS, COMPONENT_COLORS

COMPONENTS = [c for c in COUNT_COLUMNS if c != "Frame"]

def parse_ids(text: str):
    if not text:
        return []
    ids = []
    for part in text.split(','):
        part = part.strip()
        if '-' in part:
            a, b = part.split('-', 1)
            ids.extend([f"{i:02d}" for i in range(int(a), int(b) + 1)])
        elif part:
            ids.append(f"{int(part):02d}")
    return ids

def read_and_concat(spike_id: str, chunk_dirs: list[Path], stride: float, interactome: str):
    dfs = []
    frame_offset = 0
    for chunk in chunk_dirs:
        f = chunk / f"spike{spike_id}" / f"spike{spike_id}_{interactome}_molecule_counts.csv"
        if not f.exists():
            print(f"WARNING: missing {f}")
            continue
        df = pd.read_csv(f).sort_values("Frame").reset_index(drop=True)
        df["Frame"] = pd.to_numeric(df["Frame"], errors="coerce").astype(int) + frame_offset
        frame_offset = int(df["Frame"].max()) + 1
        dfs.append(df)
    if not dfs:
        return None
    out = pd.concat(dfs, ignore_index=True)
    out["Simulation_Time"] = out["Frame"] * stride
    return out[["Frame", "Simulation_Time"] + [c for c in COMPONENTS if c in out.columns]]

def group_stats(data: dict[str, pd.DataFrame]):
    if not data:
        raise ValueError("No data available for group")
    df = pd.concat(data.values(), keys=data.keys(), names=["SpikeID", "Row"]).reset_index(level="SpikeID")
    grouped = df.groupby("Frame")
    numeric = [c for c in COMPONENTS if c in df.columns]
    mean = grouped[numeric].mean().reset_index()
    std = grouped[numeric].std().reset_index()
    times = grouped["Simulation_Time"].first().reset_index()
    mean = mean.merge(times, on="Frame")
    std = std.merge(times, on="Frame")
    return mean, std

def plot_all(mean, std, output_dir: Path, prefix: str):
    time = mean["Simulation_Time"].to_numpy()
    plt.figure(figsize=(9, 6))
    for comp in COMPONENTS:
        if comp not in mean.columns:
            continue
        color = COMPONENT_COLORS.get(comp, (0.3, 0.3, 0.3))
        y = mean[comp].to_numpy(float)
        s = std[comp].fillna(0).to_numpy(float)
        plt.plot(time, y, label=comp, color=color)
        plt.fill_between(time, y - s, y + s, alpha=0.2, color=color)
    plt.xlabel("Simulation time")
    plt.ylabel("Molecule count in 6 Å shell")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(output_dir / f"{prefix}_time_evolution_mean_sd.png", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-dirs", nargs="+", required=True, help="Directories containing spikeXX count folders, in time order")
    parser.add_argument("--output-dir", default="time_evolution_all")
    parser.add_argument("--spike-ids", default="00-29")
    parser.add_argument("--exclude-ids", default="06")
    parser.add_argument("--open-ids", default="00-09")
    parser.add_argument("--closed-ids", default="10-29")
    parser.add_argument("--stride", type=float, default=25.0)
    parser.add_argument("--interactome", default="6A")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_dirs = [Path(p) for p in args.chunk_dirs]
    spike_ids = [x for x in parse_ids(args.spike_ids) if x not in set(parse_ids(args.exclude_ids))]
    open_ids = set(parse_ids(args.open_ids))
    closed_ids = set(parse_ids(args.closed_ids))

    all_data, open_data, closed_data = {}, {}, {}
    for sid in spike_ids:
        df = read_and_concat(sid, chunk_dirs, args.stride, args.interactome)
        if df is None:
            continue
        all_data[sid] = df
        df.to_csv(output_dir / f"spike{sid}_{args.interactome}_molecule_counts_all_simulation_time.csv", index=False)
        if sid in open_ids:
            open_data[sid] = df
        if sid in closed_ids:
            closed_data[sid] = df

    mean_all, std_all = group_stats(all_data)
    mean_all.to_csv(output_dir / "all_spikes_mean.csv", index=False)
    std_all.to_csv(output_dir / "all_spikes_std.csv", index=False)
    plot_all(mean_all, std_all, output_dir, "all_spikes")

    if open_data and closed_data:
        mean_open, std_open = group_stats(open_data)
        mean_closed, std_closed = group_stats(closed_data)
        mean_open.to_csv(output_dir / "open_spikes_mean.csv", index=False)
        std_open.to_csv(output_dir / "open_spikes_std.csv", index=False)
        mean_closed.to_csv(output_dir / "closed_spikes_mean.csv", index=False)
        std_closed.to_csv(output_dir / "closed_spikes_std.csv", index=False)
        for comp in COMPONENTS:
            if comp not in mean_open.columns or comp not in mean_closed.columns:
                continue
            time = mean_open["Simulation_Time"].to_numpy()
            plt.figure(figsize=(6, 4))
            for label, mean, std, color in [("Open", mean_open, std_open, "blue"), ("Closed", mean_closed, std_closed, "red")]:
                y = mean[comp].to_numpy(float)
                s = std[comp].fillna(0).to_numpy(float)
                plt.plot(time, y, label=label, color=color)
                plt.fill_between(time, y - s, y + s, alpha=0.2, color=color)
            plt.xlabel("Simulation time")
            plt.ylabel("Molecule count in 6 Å shell")
            plt.title(comp)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_dir / f"{comp}_open_vs_closed.png", dpi=300)
            plt.close()

if __name__ == "__main__":
    main()
