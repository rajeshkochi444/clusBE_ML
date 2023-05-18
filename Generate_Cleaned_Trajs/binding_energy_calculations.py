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

metal = 'Cu'

GM_ene_dict = { 'Co' :  -51.38823794, 
                'Cr' :  -73.82860343,
                'Cu' :  -24.98181639,
                'Fe' :  -63.70054494,
                'Mn' :	-70.90382205,
                'Ni' :	-37.60661532,
                'Sc' :	-45.99787217,
                'Ti' :	-59.46573482,
                'V'  :	-70.57994066,
             }



def binding_energy(adsorb, n_NH3, n_H, metal):
    traj_file = 'filtered_traj_' + adsorb + '.traj' 
    traj = Trajectory(traj_file)
    print("total images:", len(traj))
    
    ene_GM      = GM_ene_dict[metal]
    print("metal, GM_energy: ", metal, ene_GM)
    ene_mol_H2  = -6.77255467
    ene_mol_N2  = -16.6419041
    ene_mol_NH3 = -19.54417378
    ene_mol_H = 0.5 * ene_mol_H2

    ene_list = []
    binding_ene_list = []
    for i, image in enumerate(traj):
        ene = image.get_potential_energy()
        ene_list.append(ene)
        binding_ene = ene + (n_NH3 * ene_mol_NH3) - ene_GM - ene_mol_N2 - (n_H * ene_mol_H)
        binding_ene_list.append(binding_ene)
        
    y = np.array(binding_ene_list)
    print(y.shape)
        
    save_y_fname = 'y_' + adsorb + '.npy' 
    np.save(save_y_fname, y)


# * + N2 --> N2*
# * + N2 + 1H --> NNH* 1H = 1/2H2
# * + N2 + 2H --> NNH2* 
# * + N2 + 2H --> NHNH* 
# * + N2 + 3H --> N* + NH3(g) 3H = 3/2H2
# * + N2 + 3H --> NHNH2*
# * + N2 + 4H --> NH* + NH3(g)
# * + N2 + 4H --> NH2NH2*
# * + N2 + 5H --> NH2* + NH3(g) 5H = 5/2 H2
# * + N2 + 6H --> 2NH3(g)


binding_energy_dict = {'N2'   : (0, 0),
                      'NNH'   : (0, 1),
                      'NNH2'  : (0, 2),
                      'NHNH'  : (0, 2),
                      'N'     : (1, 3),
                      'NHNH2' : (0, 3),
                      'NH'    : (1, 4),
                      'NH2NH2': (0, 4),
                      'NH2'   : (1, 5)
                      }


#traj_dir = './Au/All_Full_Trajs/'
#metal_list = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']

for adsorb in adsorb_list:
    
    n_NH3 = binding_energy_dict[adsorb][0]
    n_H = binding_energy_dict[adsorb][1]
    print(adsorb, n_NH3, n_H)
    binding_energy(adsorb, n_NH3, n_H, metal)



for adsorb in adsorb_list:
    fname = 'y_' + adsorb + '.npy'
    if adsorb == 'N2':
        a = np.load(fname)
        print(a.shape)
    else:
        b = np.load(fname)
        print(b.shape)
        a = np.concatenate((a,b))
        print("stacked a", a.shape)



np.save('y_all_adsorb.npy', a)


np.load('y_all_adsorb.npy').shape


