import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import matplotlib.pyplot as plt

select_frac = 0.7
max_num = 3000

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']
#adsorb_list = ['NHNH2']
for adsorb in adsorb_list:
    img_list_adsorb = []
    img_list_bareclus = []

    traj_file_adsorb = 'random_traj_adsorb_final_'+adsorb+'.traj'
    traj_file_bareclus = 'random_traj_bareclus_final_'+adsorb+'.traj'
    traj_adsorb = Trajectory(traj_file_adsorb)
    traj_bareclus = Trajectory(traj_file_bareclus)
    tot_img_adsorb = len(traj_adsorb)
    tot_img_bareclus = len(traj_bareclus)
    print(adsorb, tot_img_adsorb, tot_img_bareclus)

    for i in range(2000):
        img_list_adsorb.append(traj_adsorb[i])
        img_list_bareclus.append(traj_bareclus[i])

    fname_adsorb = 'random_traj_adsorb_2k_'+adsorb+'.traj'
    fname_bareclus = 'random_traj_bareclus_2k_'+adsorb+'.traj'
    write(fname_adsorb, img_list_adsorb)
    write(fname_bareclus, img_list_bareclus)

    new_traj_adsorb = Trajectory(fname_adsorb)
    print(len(new_traj_adsorb))
    new_traj_bareclus = Trajectory(fname_bareclus)
    print(len(new_traj_bareclus))
    print(new_traj_adsorb[0].get_potential_energy())
    print(new_traj_bareclus[0].get_potential_energy())

