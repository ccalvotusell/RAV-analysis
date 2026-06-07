#!/usr/bin/env python3
"""Combine per-spike hotspot unique/frequency files across all spike folders."""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

def read_unique(path: Path):
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}

def read_frequency(path: Path):
    data = {}
    for line in path.read_text().splitlines():
        if not line.strip() or ":" not in line:
            continue
        key, val = line.split(":", 1)
        data[key.strip()] = int(float(val.strip()))
    return data

def rename_segment(res_id: str) -> str:
    try:
        resname, rest = res_id.split("|")
        segment, resnum = rest.split("_")
    except ValueError:
        return res_id
    m = re.match(r"^([ABC])([12])(\d{2})$", segment)
    if m:
        letter, subunit, _ = m.groups()
        segment = f"{letter}{'100' if subunit == '1' else '200'}"
    else:
        g = re.match(r"^(\d{2})(\w{2})$", segment)
        if g:
            segment = f"00{g.group(2)}"
    return f"{resname}|{segment}_{resnum}"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--folder-regex", default=r"^spike\d{2}$")
    parser.add_argument("--analysis-dir", default="analysis_results_75")
    parser.add_argument("--output-dir", default="hotspots/summary_75")
    args = parser.parse_args()

    base = Path(args.base_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    folder_re = re.compile(args.folder_regex)
    combined_unique = defaultdict(set)
    combined_freq = defaultdict(lambda: defaultdict(int))

    for folder in sorted([p for p in base.iterdir() if p.is_dir() and folder_re.match(p.name)]):
        adir = folder / args.analysis_dir
        if not adir.exists():
            print(f"WARNING: missing {adir}")
            continue
        for f in adir.glob("*_unique.txt"):
            comp = f.name.replace("_unique.txt", "")
            combined_unique[comp].update(read_unique(f))
        for f in adir.glob("*_frequency.txt"):
            comp = f.name.replace("_frequency.txt", "")
            for key, val in read_frequency(f).items():
                combined_freq[comp][rename_segment(key)] += val

    for comp, values in combined_unique.items():
        (outdir / f"{comp}_unique_all.txt").write_text("\n".join(sorted(values)) + "\n")
    for comp, freq in combined_freq.items():
        with (outdir / f"{comp}_frequency_all.txt").open("w") as fh:
            for key, val in sorted(freq.items()):
                fh.write(f"{key}: {val}\n")
    print(f"Wrote hotspot summary files in {outdir}")

if __name__ == "__main__":
    main()
