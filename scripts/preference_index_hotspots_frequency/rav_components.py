#!/usr/bin/env python3
"""Shared component definitions for RAV interactome/preference-index analyses."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple

COMPONENT_TOTALS = {
    "Albumin": 550,
    "Spike": 29,
    "Mucin": 147 + 101 + 97,
    "D Mucin": 147,
    "Q Mucin": 101,
    "T Mucin": 97,
    "M Prot": 720,
    "E Prot": 20,
    "DPPG": 24807,
    "DPPC": 236055,
    "CHL": 21213,
    "CAL": 44827,
    "CLA": 4885636,
    "SOD": 4627417,
    "MG": 27281,
    "POT": 179254,
}

# Default components used for PI tables. This excludes D/Q/T mucin subtypes to avoid
# double-counting them together with the combined "Mucin" category.
DEFAULT_PI_COMPONENTS = [
    "Spike", "Mucin", "Albumin", "M Prot", "E Prot",
    "DPPG", "DPPC", "CHL", "CAL", "CLA", "SOD", "MG", "POT",
]

COUNT_COLUMNS = [
    "Frame", "Spike", "Mucin", "D Mucin", "Q Mucin", "T Mucin",
    "Albumin", "M Prot", "E Prot", "DPPG", "DPPC", "CHL", "CAL", "CLA", "SOD", "MG", "POT",
]

COMPONENT_COLORS = {
    "Mucin": (0.992, 0.306, 0.357),
    "D Mucin": (0.992, 0.306, 0.357),
    "Q Mucin": (0.992, 0.306, 0.357),
    "T Mucin": (0.992, 0.306, 0.357),
    "DPPC": (0.925, 0.353, 0.004),
    "M Prot": (0.55, 0.55, 0.55),
    "Albumin": (0.133, 0.545, 0.133),
    "POT": (0.325, 0.173, 0.51),
    "CHL": (0.9, 0.4, 0.7),
    "MG": (0.5, 0.3, 0.0),
    "SOD": (0.5, 0.5, 0.75),
    "E Prot": (1.0, 0.85, 0.0),
    "CAL": (1.0, 0.85, 0.0),
    "CLA": (0.0, 0.9, 0.5),
    "Spike": (0.0, 0.88, 1.0),
    "DPPG": (0.96, 0.72, 0.0),
}

SPIKE_PROTEIN_RE = re.compile(r"^[ABC][12]\d{2}$")
SPIKE_GLYCAN_RE = re.compile(r"^\d{2}\w{2}$")
ALBUMIN_RE = re.compile(r"^R\w{3}$")
D_MUCIN_RE = re.compile(r"^[D4]\w{3}$")
Q_MUCIN_RE = re.compile(r"^[Q5]\w{3}$")
T_MUCIN_RE = re.compile(r"^[T6]\w{3}$")
MPROT_RE = re.compile(r"^[WV]\w{3}$")
EPROT_RE = re.compile(r"^E[1-4][1-4]\w$")
DPPG_RE = re.compile(r"^G\w{3}$")
DPPC_RE = re.compile(r"^U\w{3}$")
CHL_RE = re.compile(r"^I\w{3}$")
CAL_RE = re.compile(r"^L\w{3}$")
CLA_RE = re.compile(r"^F\w{3}$")
SOD_RE = re.compile(r"^N\w{3}$")
MG_RE = re.compile(r"^M\w{3}$")
POT_RE = re.compile(r"^K\w{3}$")

ION_LIPID_COMPONENTS = {
    "DPPG": DPPG_RE,
    "DPPC": DPPC_RE,
    "CHL": CHL_RE,
    "CAL": CAL_RE,
    "CLA": CLA_RE,
    "SOD": SOD_RE,
    "MG": MG_RE,
    "POT": POT_RE,
}

def load_mucin_pairs(paths: Iterable[str | Path]) -> Dict[str, str]:
    pairs: Dict[str, str] = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Missing mucin-pair file: {p}")
        with p.open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) == 2:
                    pairs[parts[0].upper()] = parts[1].upper()
    return pairs

def mucin_subtype(protein_seg: str) -> Optional[str]:
    seg = protein_seg.upper()
    if D_MUCIN_RE.match(seg):
        return "D Mucin"
    if Q_MUCIN_RE.match(seg):
        return "Q Mucin"
    if T_MUCIN_RE.match(seg):
        return "T Mucin"
    return None

def spike_group(segid: str) -> Optional[str]:
    seg = segid.upper()
    if SPIKE_PROTEIN_RE.match(seg):
        return seg[2:4]
    if SPIKE_GLYCAN_RE.match(seg):
        return seg[:2]
    return None

def is_spike_segment(segid: str) -> bool:
    return spike_group(segid) is not None

def is_heavy_atom(atom_name: str) -> bool:
    # PDB atom names can start with a digit; strip leading digits before testing H.
    name = atom_name.strip().lstrip("0123456789").upper()
    return not name.startswith("H")

def classify_segment(segid: str, mucin_pairs: Optional[Mapping[str, str]] = None) -> Tuple[Optional[str], Optional[str]]:
    """Return (component, unit_id) for a segment ID.

    unit_id is the unique molecule identifier used for molecule-level counting.
    For lipids/ions this function returns the component only; use a residue key
    such as (segid, chain, resid) as unit_id when counting unique residues.
    """
    seg = segid.strip().upper()
    if is_spike_segment(seg):
        return "Spike", spike_group(seg)
    if mucin_pairs:
        if seg in mucin_pairs:
            return "Mucin", seg
        # If the hit is a mucin glycan segment, map it back to the protein segment.
        for prot, gly in mucin_pairs.items():
            if seg == gly:
                return "Mucin", prot
    if ALBUMIN_RE.match(seg):
        return "Albumin", seg
    if MPROT_RE.match(seg):
        return "M Prot", seg
    if EPROT_RE.match(seg):
        return "E Prot", seg
    for comp, regex in ION_LIPID_COMPONENTS.items():
        if regex.match(seg):
            return comp, None
    return None, None

def write_component_totals_csv(path: str | Path) -> None:
    p = Path(path)
    with p.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Component", "Total Available", "Include in default PI"])
        for comp, total in COMPONENT_TOTALS.items():
            writer.writerow([comp, total, "yes" if comp in DEFAULT_PI_COMPONENTS else "no"])
