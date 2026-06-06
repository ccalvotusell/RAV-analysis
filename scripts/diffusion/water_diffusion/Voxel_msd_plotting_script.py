#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis as mda
import re
import ipympl

maindir = '/.../Analysis/Mass_voxel/'
filenames_path = f'{maindir}pdb_name_list.dat'
with open(filenames_path) as f:
    filenames = f.read().splitlines()


x = []
y = []
z = []
diff_coeff = []
slope = []
fit_reg = []
particles = []
for i in filenames:
    data = np.genfromtxt(f'{maindir}Diff_raw2/{i}.dat')
    for m in range(0,8):
        diff_coeff.append(data[0,m])
        slope.append(data[1,m])
        fit_reg.append(data[2,m])
        particles.append(data[3,m])

    xyz = re.sub(r"[xyz]","",i).split('_',-1)
    for l in [1,-1]:
        for j in [1,-1]:
            for k in [1,-1]:
                x.append(int(xyz[0]) + l*50)
                y.append(int(xyz[1]) + j*50)
                z.append(int(xyz[2]) + k*50)
            
#print(x,y,z)
d = {'Diff Coeff': diff_coeff,'Slope': slope,'Fit Error': fit_reg,'Particles': particles,'x': x, 'y': y,'z':z}
df = pd.DataFrame(data=d)
df.head()

#Plot slices in the Z direction
plt.rcParams['font.size'] = 16
plt.rcParams['text.color'] = 'w'
plt.rcParams['axes.edgecolor'] = 'w'
plt.rcParams['xtick.color'] = 'w'
plt.rcParams['ytick.color'] = 'w'
plt.rcParams['axes.grid'] = False
for i in pd.unique(df_good['z']):
    df_slice = df_good[df_good['z'] == i]
    X = df_slice['x'].to_numpy()
    Y = df_slice['y'].to_numpy()
    #Z = df_slice['z'].to_numpy()
    Diff = df_slice['Diff Coeff'].to_numpy()
    grid_x, grid_y = np.mgrid[0:3000:1000j,0:3000:1000j]
    grid_diff = griddata((X,Y),Diff,(grid_x,grid_y),method='cubic')
    fig, ax = plt.subplots(layout='constrained',figsize=(10,8))
    CS = ax.contourf(grid_x,grid_y,grid_diff,cmap='viridis',vmin=150,vmax=220,levels=500)
    #cbar = fig.colorbar(CS,extend='both',ticks=np.arange(140,230,10))
    #cbar.ax.set_ylabel('Diffusion Coefficient')
    ax.set(xlim=(0,3000),ylim=(0,3000))
    ax.set_aspect('equal','box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.grid(None)
    plt.grid(None)
    #plt.show()
    plt.savefig(f'{maindir}/z{int(i)}_diffusion.png',dpi=300,transparent=True)
    plt.close()

