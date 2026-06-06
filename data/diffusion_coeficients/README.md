# Figure 2 diffusion processed data

This folder contains processed diffusion coefficients used for the diffusion-related analyses and plots associated with Figure 2 and related supplementary figures.

Diffusion coefficients are reported in Å²/ns and were obtained from linear fits to mean-squared displacement (MSD) curves.
## Files

- 'spike_diffusion_coefficients.csv': diffusion coefficients for SARS-CoV-2 spike proteins.
- 'protein_diffusion_coefficients.csv': diffusion coefficients for albumin and M-protein components.
- 'mucin_diffusion_coefficients.csv': diffusion coefficients for all mucin chains.
- 'ion_diffusion_coefficients.csv': diffusion coefficients for CAL, CLA, MG, POT and SOD ions.
- 'lipid_diffusion_coefficients.csv': diffusion coefficients for DPPC, DPPG, CHL1 and filtered viral membrane lipids.

## Notes

These are processed data tables, not raw molecular dynamics trajectories. The full billion-atom trajectory is not included because of file-size limitations.


## Large processed tables

The full per-ion diffusion table contains 9,764,415 rows. To keep the repository below standard GitHub file-size limits, it is provided as: ion_diffusion_coefficients.csv.gz

This file can be decompressed as follows: gunzip -k ion_diffusion_coefficients.csv.gz
