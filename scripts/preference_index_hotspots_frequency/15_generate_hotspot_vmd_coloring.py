#!/usr/bin/env python3
"""Generate a VMD/Tcl coloring script from a hotspot frequency file."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from rav_components import COMPONENT_COLORS

def read_frequency_file(path: Path, value_column: str | None = None):
    data = []
    if path.suffix == ".csv":
        import pandas as pd
        df = pd.read_csv(path)
        if value_column is None:
            value_column = "Frequency" if "Frequency" in df.columns else "Normalized frequency"
        for _, row in df.iterrows():
            segid = str(row.get("Segid", ""))
            resid = str(row.get("Resid", ""))
            value = float(row[value_column])
            if segid and resid:
                data.append((segid, resid, value))
    else:
        for line in path.read_text().splitlines():
            if not line.strip() or ":" not in line:
                continue
            key, val = line.split(":", 1)
            value = float(val.strip())
            selection = key.split("|", 1)[-1]
            segid, resid = selection.split("_", 1)
            data.append((segid, resid, value))
    return data

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Frequency TXT or CSV")
    parser.add_argument("--output", required=True)
    parser.add_argument("--pdb-reference", required=True)
    parser.add_argument("--component", required=True)
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--value-column", default=None)
    args = parser.parse_args()

    data = read_frequency_file(Path(args.input), args.value_column)
    if not data:
        raise SystemExit("No data found")
    values = np.array([x[2] for x in data], dtype=float)
    min_v, max_v = values.min(), values.max()
    edges = np.linspace(min_v, max_v, args.bins + 1) if max_v != min_v else None
    bin_data = {i: [] for i in range(args.bins)}
    for segid, resid, value in data:
        if edges is None:
            idx = args.bins - 1
        elif value == max_v:
            idx = args.bins - 1
        else:
            idx = int(np.searchsorted(edges, value, side="right") - 1)
            idx = min(max(idx, 0), args.bins - 1)
        bin_data[idx].append((segid, resid))

    base = np.array(COMPONENT_COLORS.get(args.component, (0.5, 0.5, 0.5)))
    white = np.array([1.0, 1.0, 1.0])
    lines = [f"# Hotspot coloring generated for {args.component}", f"mol new {args.pdb_reference}", "display projection Orthographic", ""]
    for i in range(args.bins):
        frac = i / (args.bins - 1) if args.bins > 1 else 1.0
        color = white * (1 - frac) + base * frac
        lines.append(f"color change rgb {i} {color[0]:.3f} {color[1]:.3f} {color[2]:.3f}")
    lines.append("")
    for i, residues in bin_data.items():
        if not residues:
            continue
        selection = " or ".join([f"(segid {s} and resid {r})" for s, r in residues])
        lines += [f"mol selection {{{selection}}}", "mol material Opaque", f"mol color ColorID {i}", "mol representation VDW", "mol addrep 0", ""]
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
