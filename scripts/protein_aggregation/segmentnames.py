import MDAnalysis as mda
import numpy as np

# Input files
psf_file = "/net/gpfs-amarolab/nwauer/GORDONBELL/RAV_v2.7/4_PRODUCTION/dry/dry.psf"
dcd_file = "/net/gpfs-amarolab/nwauer/GORDONBELL/RAV_v2.7/4_PRODUCTION/dry/dry.50.dcd"

# Load the trajectory
u = mda.Universe(psf_file, dcd_file)

# Define segment selections based on original script
sel_spike = "segid A2* and name CA and resid 1016"
sel_albumin = "segid R* and name CA and resid 194"
sel_mucin_D = "segid D* and name CA and resid 376"
sel_mucin_Q = "segid Q* and name CA and resid 122"
sel_mucin_T = "segid T* and name CA and resid 118"

# Store segment metadata
segment_order = []

def add_segment_metadata(selection, comp_type, subtype=None):
    group = u.select_atoms(selection)
    if len(group) == 0:
        print(f"Warning: No atoms found for selection: {selection}")
        return
    
    sorted_atoms = sorted(group, key=lambda a: a.segid)
    for atom in sorted_atoms:
        segment_order.append(f"{atom.segid},{comp_type},{subtype if subtype else 'None'}")

# Collect segments in the same order as processed before
add_segment_metadata(sel_spike, "spike")
add_segment_metadata(sel_albumin, "albumin")
add_segment_metadata(sel_mucin_D, "mucin", "D")
add_segment_metadata(sel_mucin_Q, "mucin", "Q")
add_segment_metadata(sel_mucin_T, "mucin", "T")

# Save segment order to file
output_file = "segment_order_mapping.txt"
with open(output_file, "w") as f:
    f.write("SegmentID,ComponentType,Subtype\n")
    for entry in segment_order:
        f.write(entry + "\n")

print(f"Segment order saved to {output_file}")

