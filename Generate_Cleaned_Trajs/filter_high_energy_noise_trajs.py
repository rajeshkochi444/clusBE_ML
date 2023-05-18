import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import seaborn as sns
import matplotlib.pyplot as plt
import random
import itertools
import copy
from collections import Counter

threshold = 2.0 # with respect to minimum energy configuration

def find_min_energy(traj_file):

    traj = Trajectory(traj_file)

    ene_list = []
    for i, image in enumerate(traj):
        ene = image.get_potential_energy()
        ene_list.append(ene)


    min_ene = min(ene_list)
    max_ene = max(ene_list)
    mean_ene = np.mean(ene_list)

    return min_ene,  max_ene, mean_ene

def filter_high_energy_configs(traj_file):


    traj = Trajectory(traj_file)
    print("total images:", len(traj))

    min_ene,  max_ene, mean_ene = find_min_energy(traj_file)
    print('min, max and mean energy before filtering', min_ene,  max_ene, mean_ene)

    filtered_ene_list = []
    filtered_traj = []
    for i, image in enumerate(traj):
        ene = image.get_potential_energy()
        #if ene < min_ene + threshold:
        if (ene >= min_ene) and (ene <=  mean_ene + threshold):
            filtered_traj.append(image)
            filtered_ene_list.append(ene)

    filter_fname = 'filtered_traj_' +  adsorb +  '.traj'
    write(filter_fname,filtered_traj )



    filter_energy_y = np.array(filtered_ene_list)
    save_yene_fname = 'filtered_y_' + adsorb + '.npy'
    np.save(save_yene_fname, filter_energy_y)

    print(len(Trajectory(filter_fname)))
    print(filter_energy_y.shape)

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']

for adsorb in adsorb_list:
    traj_file = 'cleaned_overlap_'+adsorb+'.traj'
    print(traj_file)
    filter_high_energy_configs(traj_file)

    min_ene,  max_ene, mean_ene = find_min_energy('filtered_traj_'  + adsorb +  '.traj')
    print('min, max and mean energy after filtering',min_ene,  max_ene, mean_ene)
