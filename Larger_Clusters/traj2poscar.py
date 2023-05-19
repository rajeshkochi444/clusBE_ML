import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory


adsorb_list = ['NH2NH2', 'NH2NH', 'NHNH2', 'NH2N', 'NNH2', 'NNH', 'NHN', 'NHNH', 'NH2', 'NH', 'N', 'N2']

clus_size=20
angle = 90
folder_path1 = 'Adsorbates'
folder_path2 = 'Adsorbates_angle_'+str(angle)

def generate_poscar(folder_path, angle=None):
    print(folder_path)
    for adsorb in adsorb_list:
        for i in range(clus_size):
            try:
                if angle ==None:
                    fname = folder_path + '/adsorb_' + adsorb + '_' + str(i)
                else:
                     fname = folder_path + '/adsorb_' +adsorb + '_' + str(i) + '_' + str(angle)
                traj_file = fname +'.traj'
                #print(i,traj_file)
                traj = Trajectory(traj_file)

                img = traj[0]
                #print(adsorb,traj_file, img)
                #print(img.get_positions())
                #img.set_cell([14.0, 14.0, 14.0])
                img.pbc = (True, True, True)
                img.center(vacuum=10)
                #print(img)
                #print(img.get_positions())
                #print('\n')

                write(fname+'_POSCAR', img, format='vasp')
            except:
                print("File not found", i, traj_file)

generate_poscar(folder_path1)
print('\n')
generate_poscar(folder_path2, angle=angle)
