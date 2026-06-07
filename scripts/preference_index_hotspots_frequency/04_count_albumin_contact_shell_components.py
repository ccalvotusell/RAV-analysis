#!/usr/bin/env python3
"""Count component composition in albumin-centered contact-shell PDB files."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from rav_components import COUNT_COLUMNS, classify_segment, load_mucin_pairs, mucin_subtype

RESIDUE_COMPONENTS = {"DPPG", "DPPC", "CHL", "CAL", "CLA", "SOD", "MG", "POT"}

def process_pdb(pdb_file: Path, selected_albumin: str, mucin_pairs: dict[str, str]):
    molecule_segments = {"Spike": set(), "Albumin": set(), "M Prot": set(), "E Prot": set(), "Mucin": set()}
    molecule_residues = {comp: set() for comp in RESIDUE_COMPONENTS}

    with pdb_file.open() as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            segid = line[72:76].strip().upper()
            if segid == selected_albumin:
                continue
            chain_id = line[21].strip()
            resid = line[22:26].strip()

            comp, unit_id = classify_segment(segid, mucin_pairs)
            if comp is None:
                continue
            if comp == "Mucin" and unit_id:
                molecule_segments["Mucin"].add(unit_id)
            elif comp in {"Spike", "Albumin", "M Prot", "E Prot"} and unit_id:
                molecule_segments[comp].add(unit_id)
            elif comp in RESIDUE_COMPONENTS:
                molecule_residues[comp].add((segid, chain_id, resid))

    d_mucin = {m for m in molecule_segments["Mucin"] if mucin_subtype(m) == "D Mucin"}
    q_mucin = {m for m in molecule_segments["Mucin"] if mucin_subtype(m) == "Q Mucin"}
    t_mucin = {m for m in molecule_segments["Mucin"] if mucin_subtype(m) == "T Mucin"}

    counts = {
        "Spike": len(molecule_segments["Spike"]),
        "Mucin": len(d_mucin) + len(q_mucin) + len(t_mucin),
        "D Mucin": len(d_mucin),
        "Q Mucin": len(q_mucin),
        "T Mucin": len(t_mucin),
        "Albumin": len(molecule_segments["Albumin"]),
        "M Prot": len(molecule_segments["M Prot"]),
        "E Prot": len(molecule_segments["E Prot"]),
    }
    for comp in RESIDUE_COMPONENTS:
        counts[comp] = len(molecule_residues[comp])

    segment_blocks = {
        "Spike (interactors)": sorted(molecule_segments["Spike"]),
        "Mucin (protein|glycan)": [f"{m} ({mucin_pairs.get(m, 'NA')})" for m in sorted(molecule_segments["Mucin"])],
        "D Mucin": [f"{m} ({mucin_pairs.get(m, 'NA')})" for m in sorted(d_mucin)],
        "Q Mucin": [f"{m} ({mucin_pairs.get(m, 'NA')})" for m in sorted(q_mucin)],
        "T Mucin": [f"{m} ({mucin_pairs.get(m, 'NA')})" for m in sorted(t_mucin)],
        "Albumin (interactors)": sorted(molecule_segments["Albumin"]),
        "M Prot": sorted(molecule_segments["M Prot"]),
        "E Prot": sorted(molecule_segments["E Prot"]),
    }
    for comp in sorted(RESIDUE_COMPONENTS):
        segment_blocks[f"{comp} (unique residues)"] = sorted({r[0] for r in molecule_residues[comp]})
    return counts, segment_blocks

def write_interactome(folder: Path, interactome: str, pdb_files: list[Path], mucin_pairs: dict[str, str]):
    selected_albumin = folder.name.upper()
    counts_txt = folder / f"{folder.name}_{interactome}_molecule_counts.txt"
    counts_csv = folder / f"{folder.name}_{interactome}_molecule_counts.csv"
    segments_file = folder / f"{folder.name}_{interactome}_segment_names.txt"

    with counts_txt.open("w") as out_txt, counts_csv.open("w", newline="") as out_csv, segments_file.open("w") as out_segments:
        out_txt.write("\t".join(COUNT_COLUMNS) + "\n")
        writer = csv.DictWriter(out_csv, fieldnames=COUNT_COLUMNS)
        writer.writeheader()
        for pdb_file in pdb_files:
            frame = pdb_file.stem.split("_")[-1]
            counts, segment_blocks = process_pdb(pdb_file, selected_albumin, mucin_pairs)
            row = {"Frame": frame, **{c: counts.get(c, 0) for c in COUNT_COLUMNS if c != "Frame"}}
            out_txt.write("\t".join(str(row[col]) for col in COUNT_COLUMNS) + "\n")
            writer.writerow(row)
            out_segments.write(f"Frame: {frame}\n")
            for title, values in segment_blocks.items():
                out_segments.write(f"{title}:\n")
                for value in values:
                    out_segments.write(f"  {value}\n")
                out_segments.write("\n")
            out_segments.write("-" * 40 + "\n")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--mucin-pairs", nargs="+", default=["inputs/mucin_pairs_T.txt", "inputs/mucin_pairs_D.txt", "inputs/mucin_pairs_Q.txt"])
    parser.add_argument("--folder-regex", default=r"^R\w{3}$")
    args = parser.parse_args()

    mucin_pairs = load_mucin_pairs(args.mucin_pairs)
    folder_re = re.compile(args.folder_regex)
    folders = sorted([p for p in Path(args.base_dir).iterdir() if p.is_dir() and folder_re.match(p.name)])
    if not folders:
        raise SystemExit("No albumin folders found")

    for folder in folders:
        pdbs = sorted(folder.glob("*.pdb"))
        pdbs_6a = [p for p in pdbs if "20A" not in p.name]
        pdbs_20a = [p for p in pdbs if "20A" in p.name]
        print(f"Processing {folder.name}: {len(pdbs_6a)} 6A PDBs, {len(pdbs_20a)} 20A PDBs")
        write_interactome(folder, "6A", pdbs_6a, mucin_pairs)
        write_interactome(folder, "20A", pdbs_20a, mucin_pairs)

if __name__ == "__main__":
    main()
