import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import seaborn as sns
import matplotlib.pyplot as plt

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']

for adsorb in adsorb_list:
    traj_file = 'filtered_traj_'+adsorb+'.traj'
    traj = Trajectory(traj_file)
    print(adsorb)
    print(len(traj))
    ene_list = []
    image_list = []
    for i, image in enumerate(traj):
        if i % 5000 == 0:
            print(i)
        ene = image.get_potential_energy()
        ene_list.append(ene)
        image_list.append(image)

    print(len(ene_list))
    print(len(image_list))
    lowest_geom_ix = np.argmin(ene_list)
    lowest_energy_image = image_list[lowest_geom_ix]
    write('lowest_energy_'+adsorb+ '.cif', lowest_energy_image, format='cif')
    write('POSCAR_'+adsorb, lowest_energy_image, format='vasp')

    highest_idx = np.argmax(ene_list)
    highest_energy_image = image_list[highest_idx]
    write('highest_energy_'+adsorb+ '.cif', highest_energy_image, format='cif')
