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

# Define lipid atom dictionaries

lipid_res = ['CHL1','POPC','POPE','POPS','POPI','DPPC','DPPG']
lipid_hg_dict = {'CHL1':'O3','POPC':'P','POPE':'P','POPS':'P','POPI':'P','DPPC':'P','DPPG':'P'}
lipid_tail_dict = {
    'CHL1':[
        'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17'
    ],
    'POPC':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216','C217','C218',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ],
    'POPE':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216','C217','C218',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ],
    'POPS':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216','C217','C218',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ],
    'POPI':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216','C217','C218',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ],
    'DPPC':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ],
    'DPPG':[
        'C21','C22','C23','C24','C25','C26','C27','C28','C29','C210','C211','C212','C213','C214','C215','C216',
        'C31','C32','C33','C34','C35','C36','C37','C38','C39','C310','C311','C312','C313','C314','C315','C316'
        ]
    }
lipid_hg_index_dict = {'CHL1': 2,'POPC': 19,'POPE': 10,'POPS': 12,'POPI': 22,'DPPC':19,'DPPG':12}
lipid_tail_index_dict = {
    'CHL1':[0, 4, 7, 8, 10, 13, 15, 17, 20, 23, 25, 30, 33, 36, 38, 43, 46],
    'POPC':[30, 32, 39, 41, 44, 47, 50, 53, 56, 59, 62, 64, 66, 69, 72, 75, 78, 81, 84, 87, 91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121, 124, 127, 130],
    'POPE':[21, 23, 30, 32, 35, 38, 41, 44, 47, 50, 53, 55, 57, 60, 63, 66, 69, 72, 75, 78, 82, 85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121],
    'POPS':[23, 25, 32, 34, 37, 40, 43, 46, 49, 52, 55, 57, 59, 62, 65, 68, 71, 74, 77, 80, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123],
    'POPI':[33, 35, 42, 44, 47, 50, 53, 56, 59, 62, 65, 67, 69, 72, 75, 78, 81, 84, 87, 90, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121, 124, 127, 130, 133],
    'DPPC':[30, 32, 39, 41, 44, 47, 50, 53, 56, 59, 62, 65, 67, 10, 73, 76, 79, 82, 85, 89, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128],
    'DPPG':[23, 25, 32, 34, 37, 40, 43, 46, 49, 52, 55, 57, 60, 63, 66, 69, 72, 75, 78, 82, 85, 88, 91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121]
}
chl1_atom_num = 74
popc_atom_num = 134
pope_atom_num = 125
pops_atom_num = 127
popi_atom_num = 137
dppc_atom_num = 130
dppg_atom_num = 123


def sphere_coverage(points):
    """
    Estimates the surface coverage of points on a unit sphere using Delaunay triangulation.

    Args:
        points (numpy.ndarray): An array of 3D points (shape: (n_points, 3)).

    Returns:
        float: The estimated surface area covered by the points.
    """
    polar_coords = appendSpherical_np(points)
    tri = Delaunay(polar_coords[:,4:])
    #tri = Delaunay(points)
    surface_area = 0.0
    for simplex in tri.simplices:
        triangle = points[simplex]
        v1 = triangle[1] - triangle[0]
        v2 = triangle[2] - triangle[0]
        cross_product = np.cross(v1, v2)
        area = 0.5 * np.linalg.norm(cross_product)
        surface_area += area
    return surface_area

def appendSpherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew

### Define Functions ###
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

def CalcVector(xyz1,xyz2):
    #print(xyz1)
    #print(xyz2)
    vect = [xyz1[0]-xyz2[0],xyz1[1]-xyz2[1],xyz1[2]-xyz2[2]]
    magn = np.sqrt((vect[0]**2)+(vect[1]**2)+(vect[2]**2))
    unit_vect = [vect[0]/magn,vect[1]/magn,vect[2]/magn]
    return(unit_vect)

def CalcVectAngle(vec1,vec2):
    #dot = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])
    dot = np.dot(vec1,vec2)
    magn1 = np.sqrt((vec1[0]**2)+(vec1[1]**2)+(vec1[2]**2))
    magn2 = np.sqrt((vec2[0]**2)+(vec2[1]**2)+(vec2[2]**2))
    angle = np.arccos(dot/(magn1*magn2))
    return(angle)

def CalcLipidVectorsPerFrame(atomgroup,frame_index):#For MDAanlysis
    #update_progress(frame_index/(atomgroup.universe.trajectory.n_frames))
    atomgroup.universe.trajectory[frame_index]
    COM = atomgroup.center_of_mass()
    #print('Calculating lipid vectors')
    for num,i in enumerate(atomgroup.residues):
        resname = i.resname
        hg_index = lipid_hg_index_dict[resname]
        hg_xyz = i.atoms.positions[hg_index]
        res_COM = i.atoms.center_of_mass()
        cent_vect = CalcVector(hg_xyz,COM) #HG to COM of virion
        lipid_vector = CalcVector(hg_xyz,res_COM)
        cent_vect_angle = CalcVectAngle(lipid_vector,cent_vect)
        if cent_vect_angle >= np.pi/2 and cent_vect_angle <= 3*np.pi/2:
            leaflet = 1
        else:
            leaflet = 0
        vect_xyz = [lipid_vector[0],lipid_vector[1],lipid_vector[2],hg_xyz[0],hg_xyz[1],hg_xyz[2],cent_vect[0],cent_vect[1],cent_vect[2],cent_vect_angle,resname,leaflet]
        if num == 0:
            vectors_xyz = [vect_xyz]
        else:
            vectors_xyz.append(vect_xyz)
        
    return(vectors_xyz)

resp2 = mda.Universe(filedir2+'resp_lipids_all.psf',filedir2+'resp_lipids_all.dcd')
surface_clusters = []
surface_coverage = []
surface_cluster_sizes = []
for frame in np.arange(0,len(resp2.trajectory)):
    update_progress(frame/len(resp2.trajectory))
    resp2.trajectory[frame]
    surface_lipids = resp2.select_atoms('(same residue as not point 1500 1500 1500 1300) and not name H*')
    surface_hg = resp2.select_atoms('(same residue as not point 1500 1500 1500 1300) and (name P or (resname CHL1 and name O3))')

    vects = CalcLipidVectorsPerFrame(surface_lipids,frame)
    hgs_pos = surface_hg.positions

    labels = ['lipid vector X','lipid vector Y','lipid vector Z','headgroup X','headgroup Y','headgroup Z','cent vect X','cent vect Y','cent vect Z','Vectors Angle','Resname','Leaflet ID']
    df = pd.DataFrame(vects,columns=labels)
    #print(np.mean(df['Leaflet ID'].to_numpy()))
    surface_df = df[df['Leaflet ID'] == 1]


    #X = hgs_pos
    x_coor = surface_df['headgroup X'].to_numpy()
    y_coor = surface_df['headgroup Y'].to_numpy()
    z_coor = surface_df['headgroup Z'].to_numpy()

    hg_coords = np.column_stack((x_coor, y_coor, z_coor))
    #print(np.shape(hg_coords))
    #print(hg_coords[:5])

    X = hg_coords
    db = DBSCAN(eps=20,min_samples=5).fit(X)
    labels = db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    data = []
    for i in np.arange(0,len(labels)-1):
        app = [hg_coords[i][0],hg_coords[i][1],hg_coords[i][2],labels[i]]
        #print(app)
        data.append(app)
    data_re = np.reshape(data,(-1,4))
    cluster_df = pd.DataFrame(data_re,columns=['x_coor','y_coor','z_coor','clusterid'])
    surface_clusters.append(n_clusters_)
    #print("Number of total lipids: %d" % len(X))
    #print("Estimated number of clusters: %d" % n_clusters_)
    #print("Estimated number of noise points: %d" % n_noise_)
    #print(cluster_df.head())
    labels_arr = np.array(labels)
    unique_labels = set(labels)
    size=[]
    for j in unique_labels:
        size.append(np.sum(labels_arr == j))
    sizes_count = []
    for j in np.arange(1,201):
        count = 0
        for k in size:
            if k == j:
                count += 1
        sizes_count.append(count)
    surface_cluster_sizes.append(sizes_count)
    cluster_coverage = []
    if len(np.unique(labels)) > 1:
        for cluster in np.unique(labels)[1:]:
            #update_progress(cluster/len(np.unique(labels[1:])))
            cluster_df_tmp = cluster_df[cluster_df['clusterid']==cluster]
            #print(len(cluster_df_tmp))
            x_pos = cluster_df_tmp['x_coor'].to_numpy()
            y_pos = cluster_df_tmp['y_coor'].to_numpy()
            z_pos = cluster_df_tmp['z_coor'].to_numpy()
            center = [1500, 1500, 1500]
            dist = []
            sph_vect = []
            #print(x_pos)

            if size[cluster] >= 4:
                for i in np.arange(0,len(x_pos)):
                    #print(x_pos)
                    vect = np.array([x_pos[i],y_pos[i],z_pos[i]])-center
                    magn = np.sqrt((vect[0]**2)+(vect[1]**2)+(vect[2]**2))
                    unit_vect = vect/magn
                    #new_vect = (unit_vect*radius)+center
                    sph_vect.append(unit_vect)
                sph_vect = np.array(sph_vect)
                #print(sph_vect)
                coverage = sphere_coverage(sph_vect)
                sphere_surface_area = 4 * np.pi  # Surface area of a unit sphere
                coverage_ratio = (coverage / sphere_surface_area)
            else:
                coverage_ratio = 0
            cluster_coverage.append(coverage_ratio)

        surface_coverage.append(sum(cluster_coverage))
    else:
        surface_coverage.append(0)
        
np.save(output2+'surface_clusters.npy',surface_clusters)
np.save(output2+'surface_coverage.npy',surface_coverage)
np.save(output2+'surface_clusters_sizes.npy',surface_cluster_sizes)

