#!/usr/bin/env python
# coding: utf-8


import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import freud
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull, SphericalVoronoi, geometric_slerp
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.neighbors import kneighbors_graph

filedir1 = 'RAV_v2/4_PRODUCTION/dry/Slices/'
filedir2 = 'RAV_v2.7/4_PRODUCTION/dry/Slices/'
output1 = 'RAV_v2/Analysis/Resp_lipid_clustering/RAV-v2/'
output2 = 'RAV_v2/Analysis/Resp_lipid_clustering/RAV-v2.7/'

def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = "" 
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,3), status)
    sys.stdout.write(text)
    sys.stdout.flush()
    
resp1 = mda.Universe(filedir1+'resp_lipids_all.psf',filedir1+'resp_lipids_all.dcd')

hgs = resp1.select_atoms('name P or (resname CHL1 and name O3)')
print(hgs.n_residues)
print(len(hgs.universe.trajectory))


clusters = []
cluster_sizes = []
unclustered = []
for i in np.arange(0,len(resp1.trajectory)):
    update_progress(i/len(resp1.trajectory))
    hgs.universe.trajectory[i]
    hgs_pos = hgs.positions
    X = hgs_pos
    db = DBSCAN(eps=20,min_samples=5).fit(X)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    size = []
    labels_arr = np.array(labels)
    
    for j in unique_labels:
        size.append(np.sum(labels_arr == j))
    sizes_count = []
    for j in np.arange(1,201):
        count = 0
        for k in size:
            if k == j:
                count += 1
        sizes_count.append(count)
    cluster_sizes.append(sizes_count)
    clusters.append(n_clusters_)
    unclustered.append(n_noise_)
    
clusters = np.array(clusters)
cluster_sizes = np.array(cluster_sizes)
unclustered = np.array(unclustered)

np.save(output1+'clusters.npy',clusters)
np.save(output1+'cluster_sizes.npy',cluster_sizes)
np.save(output1+'unclustered.npy',unclustered)

