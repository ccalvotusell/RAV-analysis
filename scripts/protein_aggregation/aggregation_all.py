#!/usr/bin/env python
"""
Mega-script to calculate and visualize distance matrices between multiple components
in a large MD simulation. The components are:

  1. Spikes:  segid starting with "A2" (e.g. A2XX) with representative atom: CA of resid 1016.
  2. Albumins: segid starting with "R" (e.g. Rxxx) with representative atom: CA of resid 194.
  3. Mucins:   There are three subtypes:
      - Type D: segid starting with "D" (Dxxx) with representative atom: CA of resid 376.
      - Type Q: segid starting with "Q" (Qxxx) with representative atom: CA of resid 122.
      - Type T: segid starting with "T" (Txxx) with representative atom: CA of resid 118.

The simulation topology is provided by a PSF file and the trajectory is split among multiple DCD files.
The script computes the average pairwise distance matrix (averaged over frames) among all the
selected representative atoms, plots an initial heatmap (with and without quadrant contours), then
reorders the matrix via hierarchical clustering and produces:
   - a reordered viridis heatmap (with quadrant contours), and
   - a final heatmap in which each interaction’s color is determined by its component pair type.

Author: Your Name
Date: YYYY-MM-DD
"""

import os
import glob
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib import rcParams, colors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

#############################
# User parameters & paths
#############################

# Topology and Trajectory paths:
TOPOLOGY_PATH = "/net/gpfs-amarolab/nwauer/GORDONBELL/RAV_v2.7/4_PRODUCTION/dry/dry.psf"
TRAJECTORY_DIR = "/net/gpfs-amarolab/nwauer/GORDONBELL/RAV_v2.7/4_PRODUCTION/dry/"

# All DCD files in the trajectory directory (sorted for consistency)
dcd_files = sorted(glob.glob(os.path.join(TRAJECTORY_DIR, "*.dcd")))
if len(dcd_files) == 0:
    raise IOError("No DCD files found in the specified trajectory directory.")

print("Found {} DCD files.".format(len(dcd_files)))

# Plot parameters
figsize = (12, 10)
dpi = 300
fontsize_labels = 14
fontsize_ticks = 10
title_fontsize = 16

# Font: try to use Arial, falling back to DejaVu Sans.
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

#############################
# Define representative selections
#############################
# For each component type we define a selection string. We assume the following:
#   - Spikes: segid starts with "A2", CA of resid 1016.
#   - Albumins: segid starts with "R", CA of resid 194.
#   - Mucins:
#         Type D: segid starts with "D", CA of resid 376.
#         Type Q: segid starts with "Q", CA of resid 122.
#         Type T: segid starts with "T", CA of resid 118.
#
# (Adjust the selection strings if your MD system uses a different naming convention.)

sel_spike = "segid A2* and name CA and resid 1016"
sel_albumin = "segid R* and name CA and resid 194"
sel_mucin_D = "segid D* and name CA and resid 376"
sel_mucin_Q = "segid Q* and name CA and resid 122"
sel_mucin_T = "segid T* and name CA and resid 118"

#############################
# Load Universe
#############################
print("Loading Universe with topology and trajectories...")
u = mda.Universe(TOPOLOGY_PATH, *dcd_files)

#############################
# Select representative atoms and store metadata
#############################
# For each selected atom we will store a dictionary with:
#   'atom_index'  : index of the atom in the universe
#   'label'       : a string label (e.g. the segid)
#   'comp_type'   : "spike", "albumin", or "mucin"
#   'subtype'     : for mucins only (D, Q, or T), otherwise None

selected_atoms = []  # list of atom indices (order matters)
metadata = []        # list of dicts corresponding to each atom

def add_atoms(selection, comp_type, subtype=None):
    group = u.select_atoms(selection)
    if len(group) == 0:
        print("Warning: No atoms found for selection:", selection)
        return
    # sort by segid (or by atom index) to have a consistent ordering
    # Here we sort by segid string.
    sorted_atoms = sorted(group, key=lambda a: a.segid)
    for atom in sorted_atoms:
        selected_atoms.append(atom.index)
        metadata.append({
            'label': atom.segid,  # you may include additional info if desired
            'comp_type': comp_type,
            'subtype': subtype
        })

# Note: the ordering here will define the quadrant boundaries for the initial heatmap.
# We choose to order by component type: spikes, then albumins, then mucins (D, Q, T).
add_atoms(sel_spike, comp_type="spike")
add_atoms(sel_albumin, comp_type="albumin")
add_atoms(sel_mucin_D, comp_type="mucin", subtype="D")
add_atoms(sel_mucin_Q, comp_type="mucin", subtype="Q")
add_atoms(sel_mucin_T, comp_type="mucin", subtype="T")

N = len(selected_atoms)
if N == 0:
    raise ValueError("No atoms were selected. Check your selection strings.")
print("Total selected representative atoms: {}".format(N))

# For convenience, extract lists of labels and component types:
labels = [md['label'] for md in metadata]
comp_types = [md['comp_type'] for md in metadata]

# Also record the boundaries for each component type for drawing quadrant lines.
# (Assuming the ordering we built above)
n_spike = sum(1 for md in metadata if md['comp_type'] == "spike")
n_albumin = sum(1 for md in metadata if md['comp_type'] == "albumin")
n_mucin = N - (n_spike + n_albumin)
boundaries = [n_spike, n_spike+n_albumin]  # vertical/horizontal lines at these indices

#############################
# Compute average distance matrix over all frames
#############################
print("Computing average distance matrix over {} frames...".format(u.trajectory.n_frames))
dist_sum = np.zeros((N, N), dtype=np.float64)
frame_count = 0

# Loop over frames (vectorized within each frame)
for ts in u.trajectory:
    pos = u.atoms.positions[selected_atoms]  # shape (N,3)
    diff = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]  # (N,N,3)
    dists = np.sqrt(np.sum(diff**2, axis=2))  # (N,N)
    dist_sum += dists
    frame_count += 1
    if frame_count % 50 == 0:
        print("Processed {} frames...".format(frame_count))

avg_dist_matrix = dist_sum / frame_count
print("Finished processing {} frames.".format(frame_count))

#############################
# Save the distance matrix to a file (optional)
#############################
np.savetxt("combined_distance_matrix.txt", avg_dist_matrix, fmt="%.3f")
print("Distance matrix saved to combined_distance_matrix.txt")

#############################
# Plot 1: Initial Heatmap using Viridis (with and without quadrant contour lines)
#############################

def plot_heatmap(matrix, tick_labels, output_file, title="", draw_contours=False, boundaries=None):
    plt.figure(figsize=figsize, dpi=dpi)
    im = plt.imshow(matrix, cmap='viridis', origin='lower')
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    ticks = np.arange(len(tick_labels))
    plt.xticks(ticks=ticks, labels=tick_labels, fontsize=fontsize_ticks, rotation=45)
    plt.yticks(ticks=ticks, labels=tick_labels, fontsize=fontsize_ticks)
    plt.xlabel("Component", fontsize=fontsize_labels)
    plt.ylabel("Component", fontsize=fontsize_labels)
    plt.title(title, fontsize=title_fontsize)
    if draw_contours and boundaries is not None:
        # Draw horizontal and vertical lines at the boundaries.
        for b in boundaries:
            plt.axhline(b - 0.5, color='black', linewidth=1.5)
            plt.axvline(b - 0.5, color='black', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi)
    print("Saved heatmap to", output_file)
    plt.close()

# Plot without contours:
plot_heatmap(avg_dist_matrix, labels, "initial_heatmap_no_contours.png",
             title="Average Distance Matrix (Viridis) - No Contours", 
             draw_contours=False)

# Plot with contour lines (quadrant boundaries for spikes, albumins, mucins):
plot_heatmap(avg_dist_matrix, labels, "initial_heatmap_with_contours.png",
             title="Average Distance Matrix (Viridis) - With Contours", 
             draw_contours=True, boundaries=boundaries)

#############################
# Reorder (Cluster) the Distance Matrix
#############################
print("Reordering distance matrix by hierarchical clustering...")
# Convert to condensed form
condensed = squareform(avg_dist_matrix)
Z = linkage(condensed, method='average')
order = leaves_list(Z)

# Reorder the matrix and metadata
reordered_matrix = avg_dist_matrix[order][:, order]
reordered_labels = [labels[i] for i in order]
reordered_comp_types = [comp_types[i] for i in order]

# Save reordered heatmap (using viridis) with contour lines for the new grouping (boundaries not drawn here)
plot_heatmap(reordered_matrix, reordered_labels, "reordered_heatmap.png",
             title="Reordered Distance Matrix (Viridis)", draw_contours=False)

#############################
# Plot 3: Final Custom-Colored Heatmap Based on Interaction Type
#############################
print("Generating final custom-colored heatmap...")

# Define interaction categories based on the pair of component types.
def interaction_category(type1, type2):
    # types: "spike", "albumin", "mucin"
    if type1 == "spike" and type2 == "spike":
        return "spike-spike"
    elif (type1 == "spike" and type2 == "albumin") or (type1 == "albumin" and type2 == "spike"):
        return "spike-albumin"
    elif (type1 == "spike" and type2 == "mucin") or (type1 == "mucin" and type2 == "spike"):
        return "spike-mucin"
    elif type1 == "albumin" and type2 == "albumin":
        return "albumin-albumin"
    elif type1 == "mucin" and type2 == "mucin":
        return "mucin-mucin"
    elif (type1 == "albumin" and type2 == "mucin") or (type1 == "mucin" and type2 == "albumin"):
        return "albumin-mucin"
    else:
        return "other"

# Define a colormap for each interaction category.
# Each colormap goes from white (for high distances) to a dark color (for low distances).
custom_cmaps = {
    'spike-spike': colors.LinearSegmentedColormap.from_list('spike_spike', ['white', 'purple']),
    'spike-albumin': colors.LinearSegmentedColormap.from_list('spike_albumin', ['white', 'gray']),
    'spike-mucin': colors.LinearSegmentedColormap.from_list('spike_mucin', ['white', 'blue']),
    'albumin-albumin': colors.LinearSegmentedColormap.from_list('albumin_albumin', ['white', 'red']),
    'mucin-mucin': colors.LinearSegmentedColormap.from_list('mucin_mucin', ['white', 'green']),
    'albumin-mucin': colors.LinearSegmentedColormap.from_list('albumin_mucin', ['white', 'brown'])
}

# For each interaction category, we will determine the local min and max distances over the reordered matrix.
# We do this by building a mask (a boolean array) for each category.
N_reord = reordered_matrix.shape[0]
# Create an array to store the category for each (i,j) pair
interaction_cat = np.empty((N_reord, N_reord), dtype=object)
for i in range(N_reord):
    for j in range(N_reord):
        interaction_cat[i, j] = interaction_category(reordered_comp_types[i], reordered_comp_types[j])

# For each category, find the min and max of the distances in that category.
cat_minmax = {}
for cat in custom_cmaps.keys():
    mask = (interaction_cat == cat)
    if np.any(mask):
        cat_min = np.min(reordered_matrix[mask])
        cat_max = np.max(reordered_matrix[mask])
        # To avoid division by zero later:
        if cat_max - cat_min < 1e-6:
            cat_max = cat_min + 1e-6
        cat_minmax[cat] = (cat_min, cat_max)
    else:
        cat_minmax[cat] = (0, 1)  # dummy values

# Build the final RGB image (RGBA) from the reordered matrix using the custom colormaps.
final_image = np.zeros((N_reord, N_reord, 4))  # RGBA

# Loop over each pixel (vectorization per category is possible if needed)
for i in range(N_reord):
    for j in range(N_reord):
        cat = interaction_cat[i, j]
        cmap = custom_cmaps.get(cat, plt.cm.viridis)
        dmin, dmax = cat_minmax[cat]
        # Normalize so that low distance gives 1 (dark color) and high distance gives 0 (white).
        norm_val = (reordered_matrix[i, j] - dmin) / (dmax - dmin)
        norm_val = np.clip(norm_val, 0, 1)
        # Invert: lower distances -> higher value for the colormap
        color = cmap(1 - norm_val)
        final_image[i, j, :] = color

# Plot and save the final custom-colored heatmap.
plt.figure(figsize=figsize, dpi=dpi)
plt.imshow(final_image, origin='lower')
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.viridis), fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=fontsize_ticks)
plt.xticks(ticks=np.arange(N_reord), labels=reordered_labels, fontsize=fontsize_ticks, rotation=45)
plt.yticks(ticks=np.arange(N_reord), labels=reordered_labels, fontsize=fontsize_ticks)
plt.xlabel("Component", fontsize=fontsize_labels)
plt.ylabel("Component", fontsize=fontsize_labels)
plt.title("Custom-Colored Distance Matrix\n(Reordered by Clustering)", fontsize=title_fontsize)
plt.tight_layout()
plt.savefig("final_custom_colored_heatmap.png", dpi=dpi)
print("Saved final custom-colored heatmap to final_custom_colored_heatmap.png")
plt.close()

print("All processing and plotting done!")

