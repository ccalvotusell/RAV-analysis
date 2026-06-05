#!/usr/bin/env python3
"""
compute_diffusion.py

Compute mean-squared displacement (MSD) curves and diffusion coefficients
for one or more molecular groups in an MD trajectory using MDAnalysis.

The script expects a CSV file with two columns:

    group_id,selection

where `selection` is a valid MDAnalysis atom selection string.

Example:
    python compute_diffusion.py \
        --topology example_data/system.psf \
        --trajectory example_data/traj.dcd \
        --groups scripts/diffusion/example_groups.csv \
        --frame-interval-ns 0.5 \
        --fit-start-ns 20 \
        --output results/diffusion

Notes:
    - Diffusion coefficients are estimated from MSD = 2*d*D*t.
    - By default d = 3, so MSD = 6Dt.
    - Input trajectories should be imaged/aligned consistently before analysis.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Dict, Tuple

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.stats import linregress


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute MSD and diffusion coefficients for MDAnalysis selections."
    )
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True, nargs="+")
    parser.add_argument("--groups", required=True)
    parser.add_argument(
        "--frame-interval-ns",
        type=float,
        default=None,
        help="Time between saved frames in ns. If omitted, trajectory time is used.",
    )
    parser.add_argument("--fit-start-ns", type=float, default=20.0)
    parser.add_argument("--fit-end-ns", type=float, default=None)
    parser.add_argument("--dimensions", type=int, default=3, choices=[1, 2, 3])
    parser.add_argument("--output", required=True)
    parser.add_argument("--save-per-group-msd", action="store_true")
    return parser.parse_args()


def read_groups(groups_csv: str | Path) -> List[Dict[str, str]]:
    groups = []
    with open(groups_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"group_id", "selection"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                f"Group CSV must contain columns {required}. Found: {reader.fieldnames}"
            )
        for row in reader:
            group_id = row["group_id"].strip()
            selection = row["selection"].strip()
            if group_id and selection:
                groups.append({"group_id": group_id, "selection": selection})
    if not groups:
        raise ValueError(f"No valid groups found in {groups_csv}")
    return groups


def get_times_ns(universe: mda.Universe, frame_interval_ns: float | None) -> np.ndarray:
    n_frames = len(universe.trajectory)

    if frame_interval_ns is not None:
        return np.arange(n_frames, dtype=float) * frame_interval_ns

    times = []
    for ts in universe.trajectory:
        times.append(float(ts.time))
    universe.trajectory[0]

    times = np.asarray(times, dtype=float)

    # MDAnalysis often reports times in ps. Convert to ns if values look like ps.
    if np.nanmax(times) > 1000:
        times = times / 1000.0

    return times


def compute_msd_for_selection(
    universe: mda.Universe,
    selection: str,
) -> Tuple[np.ndarray, int]:
    atoms = universe.select_atoms(selection)

    if len(atoms) == 0:
        raise ValueError(f"Selection returned zero atoms: {selection}")

    universe.trajectory[0]
    reference = atoms.positions.copy()
    msd = np.zeros(len(universe.trajectory), dtype=float)

    for i, _ts in enumerate(universe.trajectory):
        displacement = atoms.positions - reference
        squared_displacement = np.sum(displacement * displacement, axis=1)
        msd[i] = np.mean(squared_displacement)

    universe.trajectory[0]
    return msd, len(atoms)


def fit_diffusion(
    times_ns: np.ndarray,
    msd_a2: np.ndarray,
    fit_start_ns: float,
    fit_end_ns: float | None,
    dimensions: int,
) -> Tuple[float, float, float, float, int]:
    if fit_end_ns is None:
        fit_end_ns = float(times_ns[-1])

    mask = (times_ns >= fit_start_ns) & (times_ns <= fit_end_ns)

    if np.count_nonzero(mask) < 3:
        raise ValueError(
            "Not enough points for linear regression. "
            f"fit_start_ns={fit_start_ns}, fit_end_ns={fit_end_ns}, "
            f"selected_points={np.count_nonzero(mask)}"
        )

    model = linregress(times_ns[mask], msd_a2[mask])
    slope = model.slope
    slope_stderr = model.stderr if model.stderr is not None else np.nan

    diffusion = slope / (2.0 * dimensions)
    diffusion_error = slope_stderr / (2.0 * dimensions)

    return diffusion, diffusion_error, float(fit_end_ns), float(model.rvalue), int(np.count_nonzero(mask))


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    groups = read_groups(args.groups)

    print("Loading trajectory...")
    universe = mda.Universe(args.topology, *args.trajectory)

    print("Reading time axis...")
    times_ns = get_times_ns(universe, args.frame_interval_ns)

    all_msd_records = []
    summary_records = []

    for group in groups:
        group_id = group["group_id"]
        selection = group["selection"]

        print(f"Computing MSD for {group_id}: {selection}")

        try:
            msd, n_atoms = compute_msd_for_selection(universe, selection)
            D, D_err, fit_end_ns, r_value, n_fit_points = fit_diffusion(
                times_ns=times_ns,
                msd_a2=msd,
                fit_start_ns=args.fit_start_ns,
                fit_end_ns=args.fit_end_ns,
                dimensions=args.dimensions,
            )

            summary_records.append(
                {
                    "group_id": group_id,
                    "selection": selection,
                    "n_atoms": n_atoms,
                    "D_A2_per_ns": D,
                    "D_error_A2_per_ns": D_err,
                    "fit_start_ns": args.fit_start_ns,
                    "fit_end_ns": fit_end_ns,
                    "fit_r_value": r_value,
                    "n_fit_points": n_fit_points,
                    "dimensions": args.dimensions,
                    "error": "",
                }
            )

            for t, value in zip(times_ns, msd):
                all_msd_records.append(
                    {"group_id": group_id, "time_ns": t, "MSD_A2": value}
                )

            if args.save_per_group_msd:
                pd.DataFrame({"time_ns": times_ns, "MSD_A2": msd}).to_csv(
                    output_dir / f"{group_id}_msd.csv", index=False
                )

        except Exception as exc:
            print(f"WARNING: failed for group {group_id}: {exc}")
            summary_records.append(
                {
                    "group_id": group_id,
                    "selection": selection,
                    "n_atoms": 0,
                    "D_A2_per_ns": np.nan,
                    "D_error_A2_per_ns": np.nan,
                    "fit_start_ns": args.fit_start_ns,
                    "fit_end_ns": args.fit_end_ns if args.fit_end_ns is not None else np.nan,
                    "fit_r_value": np.nan,
                    "n_fit_points": 0,
                    "dimensions": args.dimensions,
                    "error": str(exc),
                }
            )

    pd.DataFrame(summary_records).to_csv(
        output_dir / "diffusion_coefficients.csv", index=False
    )
    pd.DataFrame(all_msd_records).to_csv(
        output_dir / "msd_timeseries.csv", index=False
    )

    print(f"Done. Results written to: {output_dir}")


if __name__ == "__main__":
    main()
