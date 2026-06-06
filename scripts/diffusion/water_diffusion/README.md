# Water diffusion voxel analysis

This folder documents the original water-diffusion voxel workflow used to compute spatially resolved water diffusion coefficients in the respiratory aerosol virus (RAV) simulations.

The scripts are intentionally provided in their original form to document the production workflow. They contain project-specific absolute paths and were designed to run on the original analysis filesystem. To run them elsewhere, users must update their directory and provide the corresponding voxel-level PDB/PSF/DCD files.

## Files

scripts/water_diffusion/
├── getVoxel.tcl
├── Voxel_msd_script.py
├── Voxel_mass-charge_script.py
└── Voxel_msd_plotting_script.py
└── pdb_name_list.dat

## Workflow overview

### 1. Generate voxel trajectory files

'getVoxel.tcl' is a VMD/Tcl script that loops over a three-dimensional grid and extracts local cubic regions from the full RAV trajectory.

For each grid center, the script selects atoms within approximately ±100 Å in x, y and z. If the selected region contains more than 10,000 atoms, it writes:

pdbs/<voxel_name>.pdb
pdbs/<voxel_name>.psf
pdbs/<voxel_name>.dcd

where voxel names follow the convention:

x1000_y1200_z1400

The list of voxel names used for downstream processing is stored in the file 'pdb_name_list.dat'

### 2. Compute water diffusion within each voxel

'Voxel_msd_script.py' computes water MSDs and diffusion coefficients for one voxel at a time.

Example usage on the original workflow:

python Voxel_msd_script.py x1000_y1200_z1400

The script expects:

${maindir}/pdbs/x1000_y1200_z1400.psf
${maindir}/pdbs/x1000_y1200_z1400.dcd

For each voxel, the script divides the region into eight sub-voxels using greater-than/less-than selections relative to the x/y/z coordinates encoded in the filename.

Water selection: resname TIP3

The MSD is calculated with MDAnalysis 'EinsteinMSD' using 'msd_type='xyz''.

The diffusion coefficient is estimated from the linear MSD regime.

Output:

${maindir}/Diff_raw2/<voxel_name>.dat

Each output '.dat' file contains four rows across the eight sub-voxels:

row 1: diffusion coefficients
row 2: MSD slopes
row 3: regression r values
row 4: number of particles

### 3. Compute local mass and charge composition

'Voxel_mass-charge_script.py' computes mass and charge composition for the same eight sub-voxels.

Example usage:

python Voxel_mass-charge_script.py x1000_y1200_z1400

The script expects:

${maindir}/pdbs/x1000_y1200_z1400.psf
${maindir}/pdbs/x1000_y1200_z1400.pdb

It reports mass and charge for the following molecular categories:

total
water
respiratory lipids
viral membrane
albumin
spike
mucin
divalent ions
ions

Outputs:

${maindir}/Mass_raw/<voxel_name>.dat
${maindir}/Charge_raw/<voxel_name>.dat

### 4. Assemble and plot voxel diffusion maps

'Voxel_msd_plotting_script.py' reads the voxel list from:

${maindir}/pdb_name_list.dat

and combines the eight sub-voxel values from each 'Diff_raw2/<voxel_name>.dat' file into a dataframe with columns:

Diff Coeff
Slope
Fit Error
Particles
x
y
z

The script then plots two-dimensional contour slices along the z direction.

These scripts document the original production analysis. The full voxel-level inputs are derived from the billion-atom trajectory and are too large be included.

