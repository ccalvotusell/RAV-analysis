# Diffusion analysis

This module computes mean-squared displacement (MSD) curves and diffusion coefficients for selected molecular groups in molecular dynamics trajectories. The reported units are Å²/ns.
More information can be found here: https://docs.mdanalysis.org/2.7.0/documentation_pages/analysis/msd.html

## Files

compute_diffusion.py      Compute MSD curves and diffusion coefficients.
plot_diffusion.py         Plot diffusion coefficient distributions.
example_groups.csv        Example group definitions using MDAnalysis selections.

## Requirements

- Python >= 3.10
- MDAnalysis
- NumPy
- pandas
- SciPy
- matplotlib

## Input

The following script requires:

1. A topology file, for example PSF, PDB, GRO, or PRMTOP.
2. One or more trajectory files, for example DCD, XTC, or TRR.
3. A CSV file defining molecular groups and MDAnalysis atom selections.

See example 'groups.csv' file.

## Example to run command

python scripts/diffusion/compute_diffusion.py \
  --topology example_data/system.psf \
  --trajectory example_data/traj.dcd \
  --groups scripts/diffusion/example_groups.csv \
  --frame-interval-ns 0.5 \
  --fit-start-ns 20 \
  --output results/diffusion

## Outputs

The script writes:

results/diffusion/diffusion_coefficients.csv
results/diffusion/msd_timeseries.csv


## Plotting

python scripts/diffusion/plot_diffusion.py \
  --input results/diffusion/diffusion_coefficients.csv \
  --output results/diffusion/diffusion_kde.png

To plot only groups whose ID contains 'spike':

python scripts/diffusion/plot_diffusion.py \
  --input results/diffusion/diffusion_coefficients.csv \
  --output results/diffusion/spike_diffusion_kde.png \
  --group-filter spike


## Notes

Input trajectories should be preprocessed to remove periodic imaging artifacts. If diffusion is interpreted relative to the virion or aerosol reference frame, trajectories should be aligned consistently before running this analysis.

The example dataset used for testing is intended only to demonstrate code functionality. The full billion-atom production trajectory is not included because of file-size limitations.
