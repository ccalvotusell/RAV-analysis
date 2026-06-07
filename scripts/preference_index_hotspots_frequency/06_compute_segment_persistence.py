#!/usr/bin/env python3
"""Compute segment persistence at 25/50/75% from segment_names files."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

def parse_segment_names_file(filename: Path):
    frames = []
    current = {}
    current_category = None
    with filename.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("Frame:"):
                if current:
                    frames.append(current)
                current = {}
                current_category = None
                continue
            if set(line.strip()) in [set("-"), set("=")]:
                continue
            if line.endswith(":"):
                current_category = line[:-1].strip()
                current.setdefault(current_category, [])
            elif current_category:
                current[current_category].append(line.strip())
    if current:
        frames.append(current)
    return frames

def compute_persistence(frames):
    persistence = {}
    for frame in frames:
        for category, segs in frame.items():
            persistence.setdefault(category, {})
            for seg in set(segs):
                persistence[category][seg] = persistence[category].get(seg, 0) + 1
    return persistence, len(frames)

def write_outputs(folder: Path, interactome: str, persistence, total_frames: int):
    thresholds = [0.25, 0.50, 0.75]
    summary = {}
    selected_by_threshold = {int(thr * 100): {} for thr in thresholds}

    for category, counts in persistence.items():
        total_unique = len(counts)
        summary[category] = {"total": total_unique}
        for thr in thresholds:
            label = int(thr * 100)
            selected = sorted([seg for seg, count in counts.items() if total_frames and count / total_frames >= thr])
            selected_by_threshold[label][category] = selected
            summary[category][f"count_{label}"] = len(selected)
            summary[category][f"pct_{label}"] = (len(selected) / total_unique * 100) if total_unique else 0.0

    for label, category_map in selected_by_threshold.items():
        out = folder / f"{folder.name}_{interactome}_persistence_{label}.txt"
        with out.open("w") as fh:
            fh.write(f"Persistence >= {label}% (Total frames: {total_frames})\n\n")
            for category, selected in category_map.items():
                fh.write(f"{category}:\n")
                for seg in selected:
                    fh.write(f"  {seg}\n")
                fh.write("\n")

    out_summary = folder / f"{folder.name}_{interactome}_persistence_summary.csv"
    with out_summary.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Category", "Total Unique", "Persistent_25 (count)", "Persistent_25 (%)", "Persistent_50 (count)", "Persistent_50 (%)", "Persistent_75 (count)", "Persistent_75 (%)"])
        for category, s in summary.items():
            writer.writerow([category, s["total"], s["count_25"], f"{s['pct_25']:.1f}", s["count_50"], f"{s['pct_50']:.1f}", s["count_75"], f"{s['pct_75']:.1f}"])
    return summary

def write_overall(overall, out_file: Path):
    cats = sorted({cat for folder in overall.values() for cat in folder})
    with out_file.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Category", "Mean_25", "SD_25", "Mean_50", "SD_50", "Mean_75", "SD_75"])
        for cat in cats:
            row = [cat]
            for thr in [25, 50, 75]:
                vals = [d[cat][f"pct_{thr}"] for d in overall.values() if cat in d]
                row += [f"{np.mean(vals):.1f}", f"{np.std(vals, ddof=1) if len(vals) > 1 else 0.0:.1f}"]
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--folder-regex", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--interactomes", nargs="+", default=["6A", "20A"])
    args = parser.parse_args()

    base = Path(args.base_dir)
    folder_re = re.compile(args.folder_regex)
    folders = sorted([p for p in base.iterdir() if p.is_dir() and folder_re.match(p.name)])
    for interactome in args.interactomes:
        overall = {}
        for folder in folders:
            seg_file = folder / f"{folder.name}_{interactome}_segment_names.txt"
            if not seg_file.exists():
                print(f"WARNING: missing {seg_file}")
                continue
            frames = parse_segment_names_file(seg_file)
            persistence, total_frames = compute_persistence(frames)
            overall[folder.name] = write_outputs(folder, interactome, persistence, total_frames)
        if overall:
            out = base / f"{args.output_prefix}_{interactome}_overall_persistence.csv"
            write_overall(overall, out)
            print(f"Wrote {out}")

if __name__ == "__main__":
    main()
