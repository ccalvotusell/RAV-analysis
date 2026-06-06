import time, sys
import numpy as np
import re
import MDAnalysis as mda

maindir = '/.../Analysis/Mass_voxel/'

file_name = sys.argv[1]
#file_name = "x1200_y1200_z1600"
xyz = re.sub(r"[xyz]","",file_name).split('_',-1)

u = mda.Universe(f'{maindir}pdbs/{file_name}.psf', f'{maindir}pdbs/{file_name}.pdb')

great_less = ['>','<']

total_mass = []
total_charge = []
wat_mass = []
wat_charge = []
resp_lipid_mass = []
resp_lipid_charge = []
memb_mass = []
memb_charge = []
alb_mass = []
alb_charge = []
spike_mass = []
spike_charge = []
mucin_mass = []
mucin_charge = []
divalent_mass = []
divalent_charge = []
ions_mass = []
ions_charge = []


for i in great_less:
    for j in great_less:
        for k in great_less:
            total = u.select_atoms(f'prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            total_mass.append(total.accumulate('masses'))
            total_charge.append(total.accumulate('charges'))
            wat = u.select_atoms(f'resname TIP3 and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            wat_mass.append(wat.accumulate('masses'))
            wat_charge.append(wat.accumulate('charges'))
            resp_lipid = u.select_atoms(f'segid [G,U,I]* and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            resp_lipid_mass.append(resp_lipid.accumulate('masses'))
            resp_lipid_charge.append(resp_lipid.accumulate('charges'))
            memb = u.select_atoms(f'(resname POP* CHL1 and not segid I*) and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            memb_mass.append(memb.accumulate('masses'))
            memb_charge.append(memb.accumulate('charges'))
            alb = u.select_atoms(f'segid R* and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            alb_mass.append(alb.accumulate('masses'))
            alb_charge.append(alb.accumulate('charges'))
            spike = u.select_atoms(f'segid [A,B,C]* [0-9][0-9][0-9][0-9] and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            spike_mass.append(spike.accumulate('masses'))
            spike_charge.append(spike.accumulate('charges'))
            mucin = u.select_atoms(f'segid [D,Q,T,4,5,6]* and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            mucin_mass.append(mucin.accumulate('masses'))
            mucin_charge.append(mucin.accumulate('charges'))
            divalent = u.select_atoms(f'resname CAL MG and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            divalent_mass.append(divalent.accumulate('masses'))
            divalent_charge.append(divalent.accumulate('charges'))
            ions = u.select_atoms(f'resname CAL MG POT SOD and prop x {i} {xyz[0]} and prop y {j} {xyz[1]} and prop z {k} {xyz[2]}')
            ions_mass.append(ions.accumulate('masses'))
            ions_charge.append(ions.accumulate('charges'))

np.savetxt(f'{maindir}Mass_raw/{file_name}.dat',[total_mass,wat_mass,resp_lipid_mass,memb_mass,alb_mass,spike_mass,mucin_mass,divalent_mass,ions_mass])
np.savetxt(f'{maindir}Charge_raw/{file_name}.dat',[total_charge,wat_charge,resp_lipid_charge,memb_charge,alb_charge,spike_charge,mucin_charge,divalent_charge,ions_charge])
