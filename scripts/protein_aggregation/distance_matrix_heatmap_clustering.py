#!/usr/bin/env python
"""
Script to visualize and reorder the distance matrix from an MD simulation.
It reads a precomputed distance matrix from "combined_distance_matrix.txt" 
and generates heatmaps with a custom "Plasma-White" colormap.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

#############################
# Load Precomputed Distance Matrix
#############################
distance_matrix_file = "combined_distance_matrix.txt"
print(f"Loading distance matrix from {distance_matrix_file}...")
distance_matrix = np.loadtxt(distance_matrix_file)

# Ensure it's a square matrix
N = distance_matrix.shape[0]
assert distance_matrix.shape[1] == N, "Error: Distance matrix is not square."

# Apply the color cap: everything above 500 will be set to 500
color_cap = 600
distance_matrix_clipped = np.clip(distance_matrix, 0, color_cap)

print(f"Loaded distance matrix of size {N} x {N} with values capped at {color_cap}")

#############################
# Define Custom Plasma-White Colormap
#############################
plasma_colors = [
    (0.0, "#0D0887"),  # Dark purple (low distances)
    (0.25, "#9C179E"),  # Red-purple
    (0.5, "#ED7953"),  # Orange-red
    (0.75, "#F7D13D"),  # Yellow
    (1.0, "#FFFFFF")   # White (high distances)
]
custom_plasma_white = colors.LinearSegmentedColormap.from_list("plasma_white", plasma_colors)

#############################
# Plotting Function with Custom Colormap
#############################
figsize = (12, 10)
dpi = 300
fontsize_labels = 14
fontsize_ticks = 10
title_fontsize = 16

def plot_heatmap(matrix, output_file, title="", cmap=custom_plasma_white, vmin=0, vmax=600):
    plt.figure(figsize=figsize, dpi=dpi)
    im = plt.imshow(matrix, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel("Component", fontsize=fontsize_labels)
    plt.ylabel("Component", fontsize=fontsize_labels)
    plt.title(title, fontsize=title_fontsize)
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi)
    print("Saved heatmap to", output_file)
    plt.close()

# Initial heatmap (custom plasma-white)
plot_heatmap(distance_matrix_clipped, "initial_heatmap.png", title="Initial Distance Matrix (Plasma-White)")

#############################
# Reorder (Cluster) the Distance Matrix
#############################
print("Reordering distance matrix by hierarchical clustering...")

# Convert to condensed form
condensed = squareform(distance_matrix_clipped)
Z = linkage(condensed, method='average')
order = leaves_list(Z)

# Reorder the matrix
reordered_matrix = distance_matrix_clipped[order][:, order]

# Save reordered heatmap (Plasma-White)
plot_heatmap(reordered_matrix, "reordered_heatmap.png", title="Reordered Distance Matrix (Plasma-White)")

print("All processing and plotting done!")

