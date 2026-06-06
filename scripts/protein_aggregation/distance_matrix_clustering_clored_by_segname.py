#!/usr/bin/env python
"""
Final custom heatmap with recoloring based on interaction type.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import pandas as pd

#############################
# Load Data
#############################

# Load the distance matrix
distance_matrix_file = "combined_distance_matrix.txt"
print(f"Loading distance matrix from {distance_matrix_file}...")
distance_matrix = np.loadtxt(distance_matrix_file)

# Cap distance at 500
color_cap = 500
distance_matrix_clipped = np.clip(distance_matrix, 0, color_cap)

# Load segment order mapping
segment_mapping_file = "segment_order_mapping.txt"
print(f"Loading segment order mapping from {segment_mapping_file}...")
segment_data = pd.read_csv(segment_mapping_file)

# Ensure matrix is square and mapping is correct
N = distance_matrix.shape[0]
assert distance_matrix.shape[1] == N, "Error: Distance matrix is not square."
assert len(segment_data) == N, "Error: Segment mapping size does not match matrix dimensions."

# Extract segment IDs and component types
segment_ids = segment_data["SegmentID"].values
component_types = segment_data["ComponentType"].values

# Create a dictionary for quick lookup
segment_dict = dict(zip(segment_ids, component_types))

#############################
# Reorder (Cluster) the Distance Matrix
#############################
print("Reordering distance matrix by hierarchical clustering...")

# Convert to condensed form
condensed = squareform(distance_matrix_clipped)
Z = linkage(condensed, method='average')
order = leaves_list(Z)

# Reorder the matrix and segment metadata
reordered_matrix = distance_matrix_clipped[order][:, order]
reordered_segments = segment_ids[order]
reordered_types = component_types[order]

#############################
# Define Interaction-Based Colormaps
#############################

# Colormaps for each interaction type
interaction_colormaps = {
    "spike-spike": colors.LinearSegmentedColormap.from_list("spike_spike", ["white", "blue"]),
    "spike-albumin": colors.LinearSegmentedColormap.from_list("spike_albumin", ["white", "cyan"]),
    "spike-mucin": colors.LinearSegmentedColormap.from_list("spike_mucin", ["white", "purple"]),
    "albumin-albumin": colors.LinearSegmentedColormap.from_list("albumin_albumin", ["white", "green"]),
    "mucin-mucin": colors.LinearSegmentedColormap.from_list("mucin_mucin", ["white", "red"]),
    "mucin-albumin": colors.LinearSegmentedColormap.from_list("mucin_albumin", ["white", "brown"])
}

# Function to determine interaction type
def get_interaction_type(type1, type2):
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
        return "mucin-albumin"
    else:
        return None

#############################
# Generate Custom-Colored Heatmap
#############################

# Create an empty RGBA image
N_reord = reordered_matrix.shape[0]
final_image = np.zeros((N_reord, N_reord, 4))  # RGBA format

# Loop over each element in the reordered matrix
for i in range(N_reord):
    for j in range(N_reord):
        type1 = reordered_types[i]
        type2 = reordered_types[j]
        interaction_type = get_interaction_type(type1, type2)

        # Select the corresponding colormap
        cmap = interaction_colormaps.get(interaction_type, plt.cm.viridis)

        # Normalize distance to 0-1 scale within capped range
        norm_val = reordered_matrix[i, j] / color_cap  # Scale between 0 and 1
        norm_val = np.clip(norm_val, 0, 1)  # Ensure values stay within 0-1
        color = cmap(1 - norm_val)  # Invert color mapping (shorter distances darker)

        final_image[i, j, :] = color  # Assign color to image

#############################
# Plot Final Heatmap
#############################
figsize = (12, 10)
dpi = 1200
fontsize_labels = 14
fontsize_ticks = 10
title_fontsize = 16

plt.figure(figsize=figsize, dpi=dpi)
plt.imshow(final_image, origin='lower')
plt.xticks(ticks=np.arange(0, N_reord, 50), labels=[""] * len(np.arange(0, N_reord, 50)))
plt.yticks(ticks=np.arange(0, N_reord, 50), labels=[""] * len(np.arange(0, N_reord, 50)))
plt.xlabel("Component", fontsize=fontsize_labels)
plt.ylabel("Component", fontsize=fontsize_labels)
plt.title("Final Recolored Distance Matrix", fontsize=title_fontsize)

# Save the plot
output_file = "final_interaction_heatmap.png"
plt.tight_layout()
plt.savefig(output_file, dpi=dpi)
print(f"Saved final heatmap to {output_file}")
plt.close()

