# RAV-analysis

Analysis scripts and supporting workflows for the manuscript:

**Billion-atom simulations reveal condensate-like organization in virus-laden respiratory aerosols**

This repository contains analysis scripts used to characterize molecular transport, lipid self-organization, protein aggregation, spike contact-shell composition, preference indices, residue-level spike hotspots, and related source-data workflows for the Respiratory Aerosol Virus (RAV) molecular dynamics simulation.

The repository is designed to make the manuscript analysis transparent and reproducible while keeping very large simulation files outside GitHub. Full trajectories, topologies, and original coordinate files are provided through external data archives.

## Data availability

Large RAV simulation files are not stored in this GitHub repository because the full all-atom system and trajectory-derived intermediate files are too large for version control.

Original PDB/topology/trajectory data are available from:

- Amaro Lab data portal: <https://amarolab.ucsd.edu/data.php>
- Zenodo archive: **[ADD ZENODO DOI / LINK HERE]**

This GitHub repository provides the analysis code, lightweight configuration files, component-mapping files, and manuscript source-data tables where possible. Full reproduction of the trajectory-dependent analyses requires downloading the large input files from the data archives above.

## Repository contents

```text
RAV-analysis/
├── README.md
├── LICENSE
├── environment.yml
├── scripts/
│   ├── diffusion/
│   ├── lipid_clustering/
│   ├── preference_index_hotspots_frequency/
│   ├── protein_aggregation/
│   └── spike_dynamics/
```

The main analysis modules are:

| Module | Purpose | Large external data needed? |
|---|---|---:|
| `scripts/diffusion/` | Water, ion, lipid, and protein diffusion analyses | Yes |
| `scripts/lipid_clustering/` | Respiratory lipid DBSCAN clustering and surface-coverage analysis | Yes |
| `scripts/protein_aggregation/` | Protein COM distance matrices, hierarchical clustering, and aggregation networks | Yes for full rerun; no for protein aggregation analysis from provided source files |
| `scripts/preference_index_hotspots_frequency/` | 6 Å/20 Å contact-shell counts, time evolution, NBR/preference index, and spike residue hotspots | Yes for full rerun; no for PI from provided source tables |
| `scripts/spike_dynamics/` | Local helpers/links for spike dynamics; full workflow is external | External spike-dynamics repository |

## What is included and what is not included

### Included

- Source code and analysis scripts.
- Conda environment file with Python dependencies.
- License file.
- Component mapping files and lightweight source-data tables.
- Documentation describing the workflow, inputs, outputs, and limitations.

### Not included directly in GitHub

- Full RAV trajectory files.
- Full topology/coordinate files for the billion-atom simulation.
- Frame-wise extracted PDB files for every spike, albumin, mucin, or time point.
- Large temporary/intermediate files generated from trajectory processing.

These files are excluded from GitHub because they are large simulation data products. They are made available through the Amaro Lab data portal and Zenodo archive listed above and upon request.

## System requirements

### Operating systems tested

The scripts were developed for Linux-based HPC environments and Unix-like local systems. The lightweight post-processing scripts should also run on local computers with a working conda installation.

### Python requirements

Install the conda environment provided in this repository:

```bash
conda env create -f environment.yml
conda activate rav-analysis
```

The core Python dependencies include:

- Python 3.10
- NumPy
- pandas
- SciPy
- Matplotlib
- scikit-learn
- MDAnalysis
- NetworkX
- Jupyter/IPython tools for notebooks

### External software

Some analyses require software that is not installed through `environment.yml`:

- **VMD 1.9.x or later** for Tcl/VMD trajectory and structure workflows.
- **Blender / Molecular Nodes** for final rendering workflows, where applicable.
- **Mayavi** for selected 3D streamline visualizations, where applicable.

These external visualization tools are only required for the relevant visualization workflows, not for all Python source-data calculations.

### Non-standard hardware

No non-standard hardware is required to run the small demo or post-process lightweight CSV/source-data files.

Full reproduction from the RAV trajectory requires large-memory HPC resources and access to the external trajectory/topology files.

## Installation guide

Clone the repository:

```bash
git clone https://github.com/ccalvotusell/RAV-analysis.git
cd RAV-analysis
```

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate rav-analysis
```

Typical installation time on a normal desktop or laptop is approximately 5–20 minutes, depending on conda solver speed and internet connection.

Check that the environment works:

```bash
python -c "import MDAnalysis, numpy, pandas, scipy, sklearn, networkx; print('RAV-analysis environment OK')"
```

## Analysis workflows

### 1. Diffusion analysis

Location:

```text
scripts/diffusion/
```

Purpose:

- Calculate spatially resolved water diffusion.
- Estimate diffusion coefficients for ions, lipids, and proteins.
- Generate source data for diffusion distributions and spatial diffusion maps.

Typical inputs:

```text
RAV topology file
RAV trajectory file or trajectory subset
component selection definitions
voxel/grid definitions, where applicable
```

Typical outputs:

```text
water_diffusion_voxels.csv
ion_diffusion_coefficients.csv
lipid_diffusion_coefficients.csv
protein_diffusion_coefficients.csv
spatial_diffusion_maps/
```

Full diffusion calculations require external RAV trajectory/topology files and may require HPC resources.

### 2. Respiratory lipid clustering

Location:

```text
scripts/lipid_clustering/
```

Purpose:

- Identify respiratory lipid clusters using DBSCAN.
- Track cluster number, mean cluster size, and lipid surface coverage over time.
- Generate source data for lipid clustering and surface organization figures.

Typical inputs:

```text
RAV topology file
RAV trajectory file or trajectory subset
respiratory lipid selections
surface/interior classification parameters
```

Typical outputs:

```text
lipid_cluster_counts_over_time.csv
lipid_cluster_size_distributions.csv
surface_lipid_clusters.csv
surface_coverage_over_time.csv
```

### 3. Protein aggregation analysis

Location:

```text
scripts/protein_aggregation/
```

Purpose:

- Compute centers of mass for mucins, albumins, and spike proteins.
- Generate pairwise protein COM distance matrices.
- Apply hierarchical clustering to reorder distance matrices.
- Construct proximity graphs using distance cutoffs such as 100 Å, 150 Å, and 200 Å.
- Classify aggregate composition and network connectivity.

Typical inputs:

```text
RAV topology file
RAV trajectory file or trajectory subset
protein segment definitions
protein class labels: spike, mucin, albumin
```

Typical outputs:

```text
protein_com_distance_matrices/
clustered_distance_matrices/
protein_network_edges_100A.csv
protein_network_edges_150A.csv
protein_network_edges_200A.csv
connected_components_summary.csv
aggregate_composition_summary.csv
```

### 4. Preference index, contact-shell composition, and hotspot analysis

Location:

```text
scripts/preference_index_hotspots_frequency/
```

Purpose:

- Extract and count molecules within 6 Å and 20 Å contact shells around reference components.
- Generate spike-, albumin-, and mucin-centered contact-shell composition tables.
- Compute normalized binding ratios (NBRs).
- Compute preference indices using a pooled global NBR.
- Generate residue-level spike contact-frequency hotspot source data.

Important definitions:

- The primary interaction shell used for preference indices is **6 Å**.
- The 20 Å shell is used as an auxiliary broader-neighborhood comparison.
- Preference index values are calculated from NBR values normalized by a single global NBR computed across reference systems.
- For final manuscript-level preference-index tables, use the global preference index column.

Required lightweight mapping inputs:

```text
mucin_pairs_D.txt
mucin_pairs_Q.txt
mucin_pairs_T.txt
albumin_segs.txt
component_totals.csv
```

Typical source-data inputs:

```text
Spikes_6A_overall_stats.csv
Albumins_6A_overall_stats.csv
D_Mucins_6A_overall_stats.csv
Q_Mucins_6A_overall_stats.csv
T_Mucins_6A_overall_stats.csv
```

Typical outputs:

```text
global_nbr_summary.csv
Spikes_6A_preference_analysis_summary.csv
Albumins_6A_preference_analysis_summary.csv
Mucins_6A_preference_analysis_summary.csv
Spikes_6A_preference_indices.csv
Albumins_6A_preference_indices.csv
Mucins_6A_preference_indices.csv
normalized_spike_residue_contact_frequencies.csv
```

Example preference-index workflow:

```bash
conda activate rav-analysis
cd scripts/preference_index_hotspots_frequency

python 08_compute_global_nbr.py \
  --component-totals ../../source_data/preference_index_hotspots_frequency/component_totals.csv \
  --inputs \
    ../../source_data/preference_index_hotspots_frequency/Spikes_6A_overall_stats.csv \
    ../../source_data/preference_index_hotspots_frequency/Albumins_6A_overall_stats.csv \
    ../../source_data/preference_index_hotspots_frequency/D_Mucins_6A_overall_stats.csv \
    ../../source_data/preference_index_hotspots_frequency/Q_Mucins_6A_overall_stats.csv \
    ../../source_data/preference_index_hotspots_frequency/T_Mucins_6A_overall_stats.csv \
  --output ../../source_data/preference_index_hotspots_frequency/global_nbr_summary.csv

python 09_compute_preference_index.py \
  --input ../../source_data/preference_index_hotspots_frequency/Spikes_6A_overall_stats.csv \
  --component-totals ../../source_data/preference_index_hotspots_frequency/component_totals.csv \
  --global-nbr ../../source_data/preference_index_hotspots_frequency/global_nbr_summary.csv \
  --output ../../source_data/preference_index_hotspots_frequency/Spikes_6A_preference_analysis_summary.csv
```

Expected runtime from source-data CSV files is less than a few minutes on a normal desktop or laptop.

Full contact-shell extraction from trajectory-derived structures may take substantially longer and may require HPC resources.

### 5. Spike conformational dynamics

The complete workflow for spike hip, knee, and ankle tilting; NTD triangle area; and RBD–central-helix distance is maintained at:

<https://github.com/lcasalino/spike-conformational-dynamics>

Use that repository to reproduce the conformational dynamics analyses associated with the RAV and single-spike reference simulations.

## Reproducing manuscript source data

To reproduce all quantitative source data, users need:

1. This repository.
2. Large RAV topology, structure, and trajectory files from the Amaro Lab data portal and/or Zenodo.
3. The dedicated spike conformational dynamics repository.
4. The component-mapping files and source-data inputs listed in each module.

Because the full RAV system is extremely large, some analyses are best reproduced from intermediate/source-data CSV files rather than from raw trajectories. Where possible, those lightweight source-data files should be included under `source_data/`.

## Manuscript-to-repository mapping

| Manuscript analysis | Repository location | Notes |
|---|---|---|
| Water diffusion and spatial heterogeneity | `scripts/diffusion/` | Requires full trajectory/topology for full rerun |
| Ion/lipid/protein diffusion | `scripts/diffusion/` | Requires full trajectory/topology for full rerun |
| Respiratory lipid clustering | `scripts/lipid_clustering/` | Uses DBSCAN-style clustering and surface analysis |
| Protein aggregation networks | `scripts/protein_aggregation/` | COM distance matrices and network cutoffs |
| Spike contact-shell time evolution | `scripts/preference_index_hotspots_frequency/` | Uses 6 Å shell counts around spikes |
| Preference-index analysis | `scripts/preference_index_hotspots_frequency/` | Uses global NBR across reference systems |
| Spike residue-level hotspots | `scripts/preference_index_hotspots_frequency/` | Requires normalized residue contact-frequency outputs |
| Spike hip/knee/ankle, NTD area, RBD distance | External repository | See `spike-conformational-dynamics` repository |

## Known limitations and justifications

- Full trajectories and topology files are not included in GitHub due to size. They are available through external data repositories.
- Large extracted PDB files are not included because they are trajectory-derived intermediates and can be regenerated from the archived data.
- Spike conformational dynamics scripts are maintained in a separate public repository.
- Some visualization workflows require external software such as VMD, Blender, or Mayavi and may not be fully reproducible from the Python environment alone.

## License

This repository is released under the **Apache License 2.0**. See the `LICENSE` file for details.

## Citation

If you use this repository, please cite the manuscript:

**Billion-atom simulations reveal condensate-like organization in virus-laden respiratory aerosols**
Wauer, N. A.‡; Calvó-Tusell, C.‡; Dommer, A. C.; Casalino, L.; Kearns, F. L.; Caparotta, M.; Rosenfeld, M. A.; Morris, C. K.; Amaro, R. E. (https://www.biorxiv.org/content/10.64898/2026.04.30.721971v1)


## Contact

For questions about the RAV analysis scripts, please contact the manuscript corresponding author.

For questions about the full simulation data, please refer to the manuscript data availability statement, the Amaro Lab data portal, and the Zenodo archive.
