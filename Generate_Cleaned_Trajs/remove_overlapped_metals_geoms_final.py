#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import seaborn as sns
import matplotlib.pyplot as plt
import copy
from ase.data import *


# In[2]:


adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']
#adsorb_list = ['NNH2']


# In[3]:


for adsorb in adsorb_list:
    
    traj_file = 'cleaned_'+adsorb+'.traj'
    traj = Trajectory(traj_file)
    print(adsorb)
    print(len(traj))
    
    cleaned_traj_list = []
    removed_list = []
    
    for i, image in enumerate(traj):
        list_range = list(range(10))
        min_dist_list = []
        for i in list_range:
            list1 = copy.deepcopy(list_range)
            list1.remove(i)
            #print(i, l1)
            dist1 = image.get_distances(i, list1)
            #print(dist1)
            min_dist1 = min(dist1)
            #print(min_dist1)
            min_dist_list.append(min_dist1)
    
        #print(min_dist_list)
        shortest_dist = min(min_dist_list)
        #print(shortest_dist)

        elements = list(set(image.get_chemical_symbols()))

        for item in elements:
            if 'N' in elements:
                elements.remove('N')
            if 'H' in elements:
                elements.remove('H')
        metal_covalent_radii = covalent_radii[atomic_numbers[elements[0]]]
        #print(metal_covalent_radii)
        
        if shortest_dist < 0.8 * metal_covalent_radii:
            print("metal overlap", shortest_dist, image.get_potential_energy())
            removed_list.append(shortest_dist)
        else:
            cleaned_traj_list.append(image)
            
 
    print(len(cleaned_traj_list))
    write('cleaned_overlap_'+adsorb+'.traj', cleaned_traj_list)
    print(len(removed_list))


# In[5]:


for adsorb in adsorb_list:
    new_traj = Trajectory('cleaned_overlap_'+adsorb+'.traj')
    print(adsorb, len(new_traj))


# In[ ]:




