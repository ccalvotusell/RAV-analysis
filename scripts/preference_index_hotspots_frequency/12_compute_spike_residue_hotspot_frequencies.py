#!/usr/bin/env python3
"""Compute normalized spike residue-level contact frequencies from spike contact-shell PDBs.

This script implements the Methods-style hotspot calculation: any heavy atom of a
component within --cutoff Å of a target spike residue marks that residue as contacted
by that component for that frame. Frequencies are normalized by the number of
residue observations, so values are in [0, 1].
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from rav_components import classify_segment, is_heavy_atom, load_mucin_pairs, spike_group, SPIKE_PROTEIN_RE, SPIKE_GLYCAN_RE

DEFAULT_COMPONENTS = ["DPPG", "DPPC", "CHL", "CAL", "MG", "CLA", "SOD", "POT", "Mucin", "Albumin"]

def target_residue_ok(segid: str, selected_group: str, target: str) -> bool:
    seg = segid.upper()
    if spike_group(seg) != selected_group:
        return False
    if target == "protein":
        return bool(SPIKE_PROTEIN_RE.match(seg))
    if target == "glycans":
        return bool(SPIKE_GLYCAN_RE.match(seg))
    return bool(SPIKE_PROTEIN_RE.match(seg) or SPIKE_GLYCAN_RE.match(seg))

def residue_key(res):
    segid = res.segid.strip().upper()
    protomer = segid[0] if SPIKE_PROTEIN_RE.match(segid) else "glycan"
    return (res.resname.strip(), segid, str(res.resid), protomer)

def process_pdb(pdb_file: Path, selected_group: str, mucin_pairs, components, cutoff: float, target: str):
    u = mda.Universe(str(pdb_file))
    target_atoms = []
    target_reskeys = []
    observed_keys = set()

    for res in u.residues:
        if target_residue_ok(res.segid.strip(), selected_group, target):
            key = residue_key(res)
            heavy = [a for a in res.atoms if is_heavy_atom(a.name)]
            if heavy:
                observed_keys.add(key)
                for atom in heavy:
                    target_atoms.append(atom)
                    target_reskeys.append(key)

    if not target_atoms:
        return observed_keys, {c: set() for c in components}

    target_coords = np.array([a.position for a in target_atoms], dtype=float)
    component_atoms = {c: [] for c in components}

    for atom in u.atoms:
        if not is_heavy_atom(atom.name):
            continue
        segid = atom.segid.strip().upper()
        # Skip target spike atoms.
        if spike_group(segid) == selected_group:
            continue
        comp, _ = classify_segment(segid, mucin_pairs)
        if comp in component_atoms:
            component_atoms[comp].append(atom)

    contacts = {c: set() for c in components}
    for comp, atoms in component_atoms.items():
        if not atoms:
            continue
        coords = np.array([a.position for a in atoms], dtype=float)
        d = distance_array(target_coords, coords)
        contacted_atom_rows = np.where((d <= cutoff).any(axis=1))[0]
        for row in contacted_atom_rows:
            contacts[comp].add(target_reskeys[row])
    return observed_keys, contacts

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--output-dir", default="hotspots/normalized_frequency")
    parser.add_argument("--folder-regex", default=r"^spike\d{2}$")
    parser.add_argument("--mucin-pairs", nargs="+", default=["inputs/mucin_pairs_T.txt", "inputs/mucin_pairs_D.txt", "inputs/mucin_pairs_Q.txt"])
    parser.add_argument("--cutoff", type=float, default=6.0)
    parser.add_argument("--target", choices=["protein", "glycans", "protein_and_glycans"], default="protein")
    parser.add_argument("--components", default=",".join(DEFAULT_COMPONENTS))
    args = parser.parse_args()

    base = Path(args.base_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    mucin_pairs = load_mucin_pairs(args.mucin_pairs)
    folder_re = re.compile(args.folder_regex)
    components = [c.strip() for c in args.components.split(',') if c.strip()]

    detail_counts = {c: Counter() for c in components}
    aggregate_counts = {c: Counter() for c in components}
    detail_observed = Counter()
    aggregate_observed = Counter()
    n_pdbs = 0

    for folder in sorted([p for p in base.iterdir() if p.is_dir() and folder_re.match(p.name)]):
        selected_group = folder.name.replace("spike", "")
        pdb_files = sorted([p for p in folder.glob("*.pdb") if "20A" not in p.name])
        for pdb in pdb_files:
            observed, contacts = process_pdb(pdb, selected_group, mucin_pairs, components, args.cutoff, args.target)
            n_pdbs += 1
            for key in observed:
                detail_observed[key] += 1
                resname, segid, resid, protomer = key
                agg_key = (resname, resid, protomer)
                aggregate_observed[agg_key] += 1
            for comp, keys in contacts.items():
                for key in keys:
                    detail_counts[comp][key] += 1
                    resname, segid, resid, protomer = key
                    aggregate_counts[comp][(resname, resid, protomer)] += 1

    detail_out = outdir / "spike_residue_contact_frequency_by_copy.csv"
    with detail_out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Component", "Resname", "Segid", "Resid", "Protomer", "Contact frame count", "Analyzed residue observations", "Frequency"])
        for comp in components:
            all_keys = set(detail_observed) | set(detail_counts[comp])
            for key in sorted(all_keys):
                denom = detail_observed.get(key, 0)
                if denom == 0:
                    continue
                count = detail_counts[comp].get(key, 0)
                writer.writerow([comp, *key, count, denom, count / denom])

    aggregate_out = outdir / "spike_residue_contact_frequency_aggregated.csv"
    with aggregate_out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Component", "Resname", "Resid", "Protomer", "Contact frame count", "Analyzed residue observations", "Frequency"])
        for comp in components:
            all_keys = set(aggregate_observed) | set(aggregate_counts[comp])
            for key in sorted(all_keys, key=lambda x: (x[2], int(x[1]) if str(x[1]).isdigit() else str(x[1]))):
                denom = aggregate_observed.get(key, 0)
                if denom == 0:
                    continue
                count = aggregate_counts[comp].get(key, 0)
                writer.writerow([comp, *key, count, denom, count / denom])

    print(f"Processed {n_pdbs} 6A PDB files")
    print(f"Wrote {detail_out}")
    print(f"Wrote {aggregate_out}")

if __name__ == "__main__":
    main()
