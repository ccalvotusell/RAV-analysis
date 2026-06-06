import time, sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
from scipy.stats import linregress
import re

def getLinearRegression(MSD,msd_timeseries,lagtimes,start,end):
    linear_model = linregress(lagtimes[start:end],msd_timeseries[start:end])
    slope = linear_model.slope
    error = linear_model.rvalue
    D = slope * 1/(2*MSD.dim_fac)
    return(D, slope, error)

ns_to_slice = 10
num_slices = round(ns_to_slice/10)
ns_per_frame = 0.2
timestep = ns_per_frame
#slice every 10 ns or 167 frames
last_frame = -1
slice_length = round(ns_to_slice/(num_slices*ns_per_frame))
start_frames = []
end_frames = []
for i in np.arange(0,num_slices):
    frame  = last_frame - (i * slice_length)
    end_frames.append(frame)

for i in end_frames:
    start_frames.append(i-slice_length)

start_frame = start_frames[0]
end_frame = end_frames[0]

maindir = '/.../Analysis/Mass_voxel/'

file_name = sys.argv[1]
#file_name = "x800_y800_z800"
xyz = re.sub(r"[xyz]","",file_name).split('_',-1)

u = mda.Universe(f'{maindir}pdbs/{file_name}.psf', f'{maindir}pdbs/{file_name}.dcd')

#print('%d frames loaded into memory' %(len(u.trajectory)))
if len(u.trajectory) < (ns_to_slice/ns_per_frame):
    raise Exception('Too few frames loaded. %d frames loaded, %d frames required' %(len(u.trajectory),ns_to_slice/ns_per_frame))


u.trajectory[end_frame]
#print('Slice %d' % (i))

great_less = ['>','<']
D_wat_arr = []
slope_wat_arr = []
error_wat_arr = []
particles_arr = []
for i in great_less:
    for j in great_less:
        for k in great_less:
            u.trajectory[end_frame]
            MSD_wat = msd.EinsteinMSD(u, select=f'resname TIP3 and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}', msd_type='xyz')
            MSD_wat.run(start=start_frame,stop=end_frame)
            #print(0)
            msd_wat = MSD_wat.results.timeseries
            #np.savetxt('MSD_all_water_slice'+str(i)+'.csv',msd_wat,delimiter=' ')
            #print(MSD_wat.n_particles)
            
            start_time = 1
            start_index = int(start_time/timestep)
            end_time = 9
            end_index = int(end_time/timestep)
            particles = MSD_wat.n_particles
            nframes = MSD_wat.n_frames
            #timestep = 0.06 # 6 ps between frames
            lagtimes = np.arange(nframes)*timestep # make the lag-time axis
            D_wat, slope_wat, error_wat = getLinearRegression(MSD_wat,msd_wat,lagtimes,start_index,end_index)
            D_wat_arr.append(D_wat)
            slope_wat_arr.append(slope_wat)
            error_wat_arr.append(error_wat)
            particles_arr.append(particles)
            #print(D_wat, slope_wat, error_wat)
np.savetxt(f'{maindir}Diff_raw2/{file_name}.dat',[D_wat_arr,slope_wat_arr,error_wat_arr,particles_arr])

