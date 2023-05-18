import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import seaborn as sns
import matplotlib.pyplot as plt

select_frac = 0.7
max_num = 3000

def random_traj_create(adsorb):
    img_list = []

    traj_file = '../Trajs_Filtered_1eV/filtered_traj_1eV_' + adsorb + '.traj'
    traj = Trajectory(traj_file)
    tot_img = len(traj)
    print(traj_file, tot_img)

    np.random.seed(42)
    permute_idx = np.random.permutation(len(traj))
    print(permute_idx)


    tot_img_selected = round(tot_img * 0.7)
    print(tot_img_selected)



    if tot_img_selected <= max_num:
        random_idx = permute_idx[:tot_img_selected]
    else:
        random_idx = permute_idx[:max_num]

    for i in range(len(random_idx)):
        k = random_idx[i]
        img = traj[k]
        img_list.append(img)

    fname = 'random_traj_'+adsorb+'.traj'
    write(fname, img_list)

    new_traj = Trajectory(fname)
    print(len(new_traj))

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']
#adsorb_list = ['NHNH2']
for adsorb in adsorb_list:
    random_traj_create(adsorb)
