# RAV spike/albumin interactome, hotspot, and preference-index analysis

This folder contains scripts for the RAV contact-shell analyses used to quantify spike, mucin and albumin interactomes, contact-shell time evolution, residue-level spike hotspots, and preference indices.
## 1. Environment

```bash
conda create -n rav_interactome python=3.10
conda activate rav_interactome
pip install -r requirements.txt
```

VMD is required for the Tcl extraction scripts.

## 2. Spike-centered contact-shell workflow

### 2.1 Extract 6 Å and 20 Å spike contact-shell PDBs

Edit paths in:

```text
scripts/01_extract_spike_contact_shells.tcl
```

Then run:

```bash
vmd -dispdev text -e scripts/01_extract_spike_contact_shells.tcl
```

Expected output:

```text
spike00/RAV_v27_s00_frame_0000.pdb
spike00/RAV_v27_s00_20A_frame_0000.pdb
spike01/...
```

### 2.2 Count molecules in each spike contact shell

```bash
python scripts/02_count_spike_contact_shell_components.py \
  --base-dir . \
  --mucin-pairs inputs/mucin_pairs_T.txt inputs/mucin_pairs_D.txt inputs/mucin_pairs_Q.txt
```

Expected output per spike:

```text
spikeXX/spikeXX_6A_molecule_counts.csv
spikeXX/spikeXX_6A_molecule_counts.txt
spikeXX/spikeXX_6A_segment_names.txt
spikeXX/spikeXX_20A_molecule_counts.csv
spikeXX/spikeXX_20A_molecule_counts.txt
spikeXX/spikeXX_20A_segment_names.txt
```

### 2.3 Compute average spike contact-shell statistics

```bash
python scripts/05_compute_average_counts.py \
  --base-dir . \
  --folder-regex '^spike\d{2}$' \
  --output-prefix Spikes
```

Expected output:

```text
Spikes_6A_overall_stats.csv
Spikes_20A_overall_stats.csv
spikeXX/spikeXX_6A_stats.csv
spikeXX/spikeXX_20A_stats.csv
```

### 2.4 Compute segment persistence

```bash
python scripts/06_compute_segment_persistence.py \
  --base-dir . \
  --folder-regex '^spike\d{2}$' \
  --output-prefix Spikes
```

Expected output:

```text
spikeXX/spikeXX_6A_persistence_25.txt
spikeXX/spikeXX_6A_persistence_50.txt
spikeXX/spikeXX_6A_persistence_75.txt
spikeXX/spikeXX_6A_persistence_summary.csv
Spikes_6A_overall_persistence.csv
```

## 3. Albumin-centered contact-shell workflow

### 3.1 Extract albumin contact-shell PDBs

Edit paths in:

```text
scripts/03_extract_albumin_contact_shells.tcl
```

Then run:

```bash
vmd -dispdev text -e scripts/03_extract_albumin_contact_shells.tcl
```

Expected output:

```text
R000/Albumin_R000_frame_0000.pdb
R000/Albumin_R000_20A_frame_0000.pdb
R001/...
```

### 3.2 Count molecules in each albumin shell

```bash
python scripts/04_count_albumin_contact_shell_components.py \
  --base-dir . \
  --mucin-pairs inputs/mucin_pairs_T.txt inputs/mucin_pairs_D.txt inputs/mucin_pairs_Q.txt
```

### 3.3 Compute average albumin contact-shell statistics

```bash
python scripts/05_compute_average_counts.py \
  --base-dir . \
  --folder-regex '^R\w{3}$' \
  --output-prefix Albumins
```

Expected output:

```text
Albumins_6A_overall_stats.csv
Albumins_20A_overall_stats.csv
```

## 4. Contact-shell time evolution

```bash
python scripts/07_compute_spike_contact_shell_time_evolution.py \
  --chunk-dirs /path/to/chunk1 /path/to/chunk2 \
  --output-dir time_evolution_all \
  --stride 25 \
  --spike-ids 00-29 \
  --exclude-ids 06 \
  --open-ids 00-09 \
  --closed-ids 10-29
```

Expected output:

```text
time_evolution_all/all_spikes_mean.csv
time_evolution_all/all_spikes_std.csv
time_evolution_all/open_spikes_mean.csv
time_evolution_all/open_spikes_std.csv
time_evolution_all/closed_spikes_mean.csv
time_evolution_all/closed_spikes_std.csv
time_evolution_all/*_open_vs_closed.png
```

## 5. Preference-index workflow

The paper reports the global preference index, where NBR values are normalized by a pooled global NBR calculated across reference interactomes.

### 5.1 Compute global NBR

```bash
python scripts/08_compute_global_nbr.py \
  --input-files \
    Spikes_6A_overall_stats.csv \
    Albumins_6A_overall_stats.csv \
    D_Mucins_6A_overall_stats.csv \
    Q_Mucins_6A_overall_stats.csv \
    T_Mucins_6A_overall_stats.csv \
  --output global_nbr_summary.csv \
  --component-totals-output component_totals.csv
```

Expected output:

```text
global_nbr_summary.csv
component_totals.csv
```

### 5.2 Compute PI for each reference system

```bash
python scripts/09_compute_preference_index.py \
  --input Spikes_6A_overall_stats.csv \
  --global-nbr-summary global_nbr_summary.csv \
  --output Spikes_6A_preference_analysis_summary.csv

python scripts/09_compute_preference_index.py \
  --input Albumins_6A_overall_stats.csv \
  --global-nbr-summary global_nbr_summary.csv \
  --output Albumins_6A_preference_analysis_summary.csv

python scripts/09_compute_preference_index.py \
  --input Mucins_6A_overall_stats.csv \
  --global-nbr-summary global_nbr_summary.csv \
  --output Mucins_6A_preference_analysis_summary.csv
```


### 5.3 Extract compact PI tables and plot

```bash
python scripts/10_extract_preference_columns.py \
  --input Spikes_6A_preference_analysis_summary.csv \
  --output Spikes_6A_preference_indices.csv

python scripts/11_plot_preference_index.py \
  --input Spikes_6A_preference_analysis_summary.csv Albumins_6A_preference_analysis_summary.csv Mucins_6A_preference_analysis_summary.csv \
  --output-dir preference_index_plots
```

## 6. Spike residue-level hotspot frequencies

```bash
python scripts/12_compute_spike_residue_hotspot_frequencies.py \
  --base-dir . \
  --output-dir hotspots/normalized_frequency \
  --cutoff 6.0 \
  --target protein \
  --mucin-pairs inputs/mucin_pairs_T.txt inputs/mucin_pairs_D.txt inputs/mucin_pairs_Q.txt
```

Expected output:

```text
hotspots/normalized_frequency/spike_residue_contact_frequency_by_copy.csv
hotspots/normalized_frequency/spike_residue_contact_frequency_aggregated.csv
```

## 7. Hotspot visualization in VMD

```bash
python scripts/15_generate_hotspot_vmd_coloring.py \
  --input hotspots/normalized_frequency/spike_residue_contact_frequency_aggregated.csv \
  --component CAL \
  --pdb-reference reference_spike.pdb \
  --output hotspots/visualization/color_CAL_hotspots.tcl
```

