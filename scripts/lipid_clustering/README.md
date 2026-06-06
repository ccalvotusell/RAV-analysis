# Respiratory lipid clustering and surface coverage analysis

This folder documents the original respiratory-lipid clustering workflow used to quantify lipid self-organization and surface-associated lipid coverage in the respiratory aerosol virus (RAV) simulations.

The scripts are intentionally provided in their original form to document the production workflow. They contain project-specific paths and were designed to run on the original RAV analysis filesystem. To run them elsewhere, users must update the input/output directories and provide the corresponding lipid-only PSF/DCD files.

## Files

scripts/lipid_clustering/
├── Resp_lipid_clustering_script.py
└── Resp_lipid_clustering_surface_coverage_script.py

## Workflow overview

### 1. Respiratory lipid clustering

'Resp_lipid_clustering_script.py' analyzes respiratory lipid self-organization using the headgroup atoms of respiratory lipids.

Main atom selection:

name P or (resname CHL1 and name O3)

This selection captures lipid phosphate atoms and the cholesterol O3 atom as headgroup/reference atoms.

For each trajectory frame, the script applies DBSCAN clustering to the lipid headgroup coordinates:

DBSCAN(eps=20, min_samples=5)

The script reports:

clusters.npy        number of clusters per frame
cluster_sizes.npy   cluster-size distribution per frame
unclustered.npy     number of unclustered/noise lipids per frame


### 2. Surface-associated lipid clustering and coverage

'Resp_lipid_clustering_surface_coverage_script.py' analyzes the subset of lipids located at the aerosol surface and estimates their surface coverage.

Surface-lipid selection:

same residue as not point 1500 1500 1500 1300

This selects lipids outside a sphere of radius 1300 Å centered at `(1500, 1500, 1500)`, corresponding to the approximate aerosol surface region in the production coordinate frame.

The script then:

1. calculates lipid orientation vectors;
2. assigns lipids to inner/outer leaflets based on the lipid vector orientation relative to the aerosol center;
3. keeps outer/surface-facing lipids;
4. clusters surface lipid headgroups with DBSCAN using: DBSCAN(eps=20, min_samples=5)
5. estimates fractional surface coverage using a triangulation-based estimate on a unit sphere.

The script reports:

surface_clusters.npy          number of surface clusters per frame
surface_coverage.npy          estimated surface coverage fraction per frame
surface_clusters_sizes.npy    surface cluster-size distribution per frame


