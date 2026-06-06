# Protein aggregation and network analysis

This folder documents the protein aggregation and proximity-network analysis used for Figure 3 of the RAV manuscript.

The workflow quantifies mesoscale spatial organization between the major protein classes in the respiratory aerosol virus (RAV) system:

- SARS-CoV-2 spike proteins
- human serum albumin
- mucins

The scripts are provided in their original production-analysis form where possible. Some scripts contain project-specific absolute paths and were designed to run on the original analysis filesystem. To run them elsewhere, update the input/output paths and provide the corresponding topology and trajectory files.

## Repository placement

Recommended organization:

```text
RAV-analysis/
├── scripts/
│   └── protein_aggregation/
│       ├── aggregation_all.py
│       ├── segmentnames.py
│       ├── distance_matrix_heatmap_clustering.py
│       ├── distance_matrix_clustering_clored_by_segname.py
│       ├── graph_bycolor_aggregation_100.py
│       ├── environment_full_gt_python_3-13.yml
│       └── README.md
│
├── notebooks/
│   └── interactive-graph-aggregation.ipynb
│
└── data/
    └── figure3_protein_aggregation/
        ├── combined_distance_matrix.txt
        ├── segment_order_mapping.txt
        ├── updated_segment_order_mapping.txt
        └── README.md
```

## Input data

### `combined_distance_matrix.txt`

Processed pairwise distance matrix between protein representatives in the RAV system.

Rows and columns follow the order defined in:

```text
segment_order_mapping.txt
```

This matrix is used as the processed input for the heatmap, hierarchical clustering, and graph analyses.

### `segment_order_mapping.txt`

Metadata file mapping each row/column in the distance matrix to a protein segment and component class:

```text
SegmentID,ComponentType,Subtype
```

Example component classes:

```text
spike
albumin
mucin
```

Mucin entries can include subtype labels:

```text
D
Q
T
```

### `updated_segment_order_mapping.txt`

Extended mapping file used by graph-visualization workflows. This file contains the segment identifiers and component classes, and may also include additional per-node metadata such as diffusion coefficients.

## Scripts

### `aggregation_all.py`

Main production script for generating the protein-protein distance matrix and heatmap visualizations from trajectory data.

The script:

1. loads the RAV topology and trajectory files;
2. selects representative Cα atoms for spikes, albumins, and mucin subtypes;
3. stores segment metadata for each selected representative atom;
4. loops over trajectory frames;
5. computes the pairwise Euclidean distance matrix between all selected representatives for each frame;
6. averages the distance matrix over all frames;
7. writes the processed matrix to:

```text
combined_distance_matrix.txt
```

The representative selections are:

```python
sel_spike   = "segid A2* and name CA and resid 1016"
sel_albumin = "segid R* and name CA and resid 194"
sel_mucin_D = "segid D* and name CA and resid 376"
sel_mucin_Q = "segid Q* and name CA and resid 122"
sel_mucin_T = "segid T* and name CA and resid 118"
```

The script also generates initial and clustered heatmaps:

```text
initial_heatmap_no_contours.png
initial_heatmap_with_contours.png
reordered_heatmap.png
final_custom_colored_heatmap.png
```

### `segmentnames.py`

Utility script for generating the segment-to-component mapping file:

```text
segment_order_mapping.txt
```

The script uses the same representative selections as `aggregation_all.py` and writes:

```text
SegmentID,ComponentType,Subtype
```

This file defines the ordering of segments used to interpret rows and columns of the distance matrix.

### `distance_matrix_heatmap_clustering.py`

Visualization script for plotting the raw and hierarchically reordered distance matrix.

Typical outputs:

```text
initial_heatmap.png
reordered_heatmap.png
```

This script starts from the already processed matrix:

```text
combined_distance_matrix.txt
```

### `distance_matrix_clustering_clored_by_segname.py`

Visualization script that combines the distance matrix with segment metadata to generate an interaction-colored clustered heatmap.

Typical output:

```text
final_interaction_heatmap.png
```

Interaction classes include combinations such as:

```text
spike-spike
spike-albumin
spike-mucin
albumin-albumin
albumin-mucin
mucin-mucin
```

### `graph_bycolor_aggregation_100.py`

Graph-based aggregation analysis using a 100 Å distance cutoff.

The script:

1. loads `combined_distance_matrix.txt`;
2. loads `updated_segment_order_mapping.txt`;
3. creates one node per protein representative;
4. creates an edge between pairs separated by less than the distance threshold;
5. colors/shapes nodes by component class;
6. analyzes connected components and their composition.

Typical outputs:

```text
diffusion_colormap_graph_100.png
diffusion_colormap_graph_100_expanded.png
cluster_analysis.txt
```

The manuscript describes graph analyses at:

```text
100 Å
150 Å
200 Å
```

The included script is the original 100 Å version. Equivalent 150 Å and 200 Å analyses can be generated by changing the distance threshold and output filenames, or by adding the corresponding original scripts if available.

### `interactive-graph-aggregation.ipynb`

Interactive notebook for exploring protein aggregation networks.

This notebook should be placed under:

```text
notebooks/
```

The notebook is intended for interactive visualization and exploration of the processed aggregation graph. It should be documented as a visualization aid rather than the primary production script that computes the distance matrix.

## Workflow summary

### Step 1 — Generate segment metadata

```bash
python scripts/protein_aggregation/segmentnames.py
```

Output:

```text
segment_order_mapping.txt
```

### Step 2 — Generate the average distance matrix

```bash
python scripts/protein_aggregation/aggregation_all.py
```

Output:

```text
combined_distance_matrix.txt
```

and heatmap images.

### Step 3 — Generate matrix visualizations

```bash
python scripts/protein_aggregation/distance_matrix_heatmap_clustering.py
python scripts/protein_aggregation/distance_matrix_clustering_clored_by_segname.py
```

### Step 4 — Generate graph-based aggregation networks

```bash
python scripts/protein_aggregation/graph_bycolor_aggregation_100.py
```

For 150 Å and 200 Å graph analyses, update the distance threshold in the graph script or include separate graph scripts for those cutoffs.

### Step 5 — Optional interactive graph exploration

Open:

```text
notebooks/interactive-graph-aggregation.ipynb
```

## Dependencies

Core dependencies:

```text
Python
NumPy
pandas
matplotlib
MDAnalysis
SciPy
```

Graph workflow dependencies:

```text
graph-tool
pyvis
```

The graph-tool dependency can be difficult to install through pip, so the original environment file is included:

```text
environment_full_gt_python_3-13.yml
```

## Important notes

### Large-system limitation

The full RAV trajectory is not included in this repository because of file-size limitations. The processed matrix file is provided so that downstream heatmap and graph visualizations can be regenerated without distributing the full billion-atom trajectory.

### Path assumptions

Several scripts contain hard-coded paths to the original RAV analysis filesystem. These must be updated before running the scripts on a different machine.

### Representative-atom versus COM definition

The current `aggregation_all.py` workflow uses representative Cα atoms for each protein segment, not a full center-of-mass calculation. This should be checked against the final manuscript Methods text. If the manuscript states that center-of-mass distances were used, either:

1. replace this script with the final COM-based script, or
2. update the Methods text to describe representative Cα atom distances.

### Hierarchical clustering method

The uploaded scripts use average-linkage hierarchical clustering. If the final manuscript states Ward’s method, this should be checked for consistency before final publication.

## Data availability note

This folder contains processed data sufficient to regenerate the matrix-based and graph-based visualizations. Full trajectory-level reproduction requires access to the RAV topology and trajectory files.
