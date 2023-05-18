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


# In[2]:


def remove_dissociated_adsorbates(traj_file, adsorb):
   
    traj = Trajectory(traj_file)
    len_adsorb = adsorb_tot_atoms[adsorb]
    tot_N = tot_N_list[adsorb]
    tot_H = tot_H_list[adsorb]
    print(len(traj))
    print(len_adsorb,tot_N, tot_H )
    
    if tot_N ==2:
        cleaned_traj_list = []
        
        for i, image in enumerate(traj):

            dist1_cluster = image.get_distances(10, list(range(10)))
            dist2_cluster = image.get_distances(11, list(range(10)))
            min_dist1 = min(dist1_cluster)
            min_dist2 = min(dist2_cluster)
            min_distNM = min(min_dist1, min_dist2)
            
            if min_distNM < clus_dist_threshold: 

                list_range = list(range(10, 10 + len_adsorb))
                list1 = [i for i in list_range if i != 10]
                list2 = [i for i in list_range if i != 11]
               
    
                dist1 = image.get_distances(10, list1)
                dist2 = image.get_distances(11, list2)
                ene1 = image.get_potential_energy()

    
                if dist1[0] < bond_distance_NN:
                    
                    H_bonded_list1 = [ dist1[i] for i in range(1,len(dist1)) if dist1[i] < bond_distance_NH]
                    H_bonded_list2 = [ dist2[i] for i in range(1,len(dist2)) if dist2[i] < bond_distance_NH]               
            
                    if adsorb == 'NH2NH2':
                        if len(H_bonded_list1) == 2 and len(H_bonded_list2) == 2:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)                            
                    
                    elif adsorb == 'NHNH':    
                        if len(H_bonded_list1) == 1 and len(H_bonded_list2) == 1:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)                            
                    
                    elif adsorb == 'N2':           
                        if len(H_bonded_list1) == 0 and len(H_bonded_list2) == 0:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)
                    
                    elif adsorb == 'NHNH2':            
                        if len(H_bonded_list1) == 1 and len(H_bonded_list2) == 2:
                            cleaned_traj_list.append(image)                            
                        elif len(H_bonded_list1) == 2 and len(H_bonded_list2) == 1:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)                            
                        
                    elif adsorb == 'NNH':    
                        if len(H_bonded_list1) == 0 and len(H_bonded_list2) == 1:
                            cleaned_traj_list.append(image)                            
                        elif len(H_bonded_list1) == 1 and len(H_bonded_list2) == 0:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)                            
                        
                    elif adsorb == 'NNH2':             
                        if len(H_bonded_list1) == 0 and len(H_bonded_list2) == 2:
                            cleaned_traj_list.append(image)                            
                        elif len(H_bonded_list1) == 2 and len(H_bonded_list2) == 0:
                            cleaned_traj_list.append(image)                            
                        else:
                            print("Dissocated", H_bonded_list1, H_bonded_list2)
                    else:
                        print("Some issue. Check trajectory")
                        
                else:
                    print("N-N Bond Broken", dist1[0])
                    
            else:
                print("NonBonded to Cluster", min_dist1, min_dist2)
                
    elif tot_N ==1 and tot_H != 0 :
       
        cleaned_traj_list = []        
        for i, image in enumerate(traj):

            dist1_cluster = image.get_distances(10, list(range(10)))
            min_dist1 = min(dist1_cluster)            
           
            if min_dist1 < clus_dist_threshold: 

                list_range = list(range(10, 10 + len_adsorb))
                list1 = [i for i in list_range if i != 10]                
    
                dist1 = image.get_distances(10, list1)
                ene1 = image.get_potential_energy()
                
                H_bonded_list1 = [ dist1[i] for i in range(len(dist1)) if dist1[i] < bond_distance_NH]
                    
                if adsorb == 'NH':
                    if len(H_bonded_list1) == 1:
                        cleaned_traj_list.append(image)                        
                    else:
                        print("Dissocated", H_bonded_list1, dist1)
                    
                elif adsorb == 'NH2':    
                    if len(H_bonded_list1) == 2:
                        cleaned_traj_list.append(image)                        
                    else:
                        print("Dissocated", H_bonded_list1, dist1)
                
                else:
                    print("Some issue. Check trajectory")
            else:
                print('Non Bonded to Cluster',min_dist1 )
                    
    elif tot_N ==1 and tot_H ==0:
       
        cleaned_traj_list = []        
        for i, image in enumerate(traj):

            dist1_cluster = image.get_distances(10, list(range(10)))
            min_dist1 = min(dist1_cluster)            
            
            if min_dist1 < clus_dist_threshold: 
                cleaned_traj_list.append(image)                 
            else:
                print("Some issue. Check trajectory") 
    else:
        print("Something wrong in the code")

    print("Cleaned Traj images")
    print(len(cleaned_traj_list))
    write('cleaned_'+adsorb+'.traj', cleaned_traj_list)
    print('\n')


# In[3]:


#Parameters
bond_distance_NN = 1.62
bond_distance_NH = 1.20
clus_dist_threshold = 4.0

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']
adsorb_tot_atoms = {'N2': 2, 'NH2NH2': 6, 'NHNH': 4, 'NHNH2': 5, 'NNH': 3, 'NNH2':4, 'N': 1, 'NH':2, 'NH2': 3 }
tot_N_list = {'N2': 2, 'NH2NH2': 2, 'NHNH': 2, 'NHNH2': 2, 'NNH': 2, 'NNH2':2, 'N': 1, 'NH':1, 'NH2': 1 }
tot_H_list = {'N2': 0, 'NH2NH2': 4, 'NHNH': 2, 'NHNH2': 3, 'NNH': 1, 'NNH2':2, 'N': 0, 'NH':1, 'NH2': 2 }


# In[4]:


#traj_dir = './All_Full_Trajs/'
for adsorb in adsorb_list:
    traj_file =  'single_traj_' + adsorb + '.traj'
    print(traj_file)
    print(adsorb)
    remove_dissociated_adsorbates(traj_file, adsorb)


# In[ ]:





# In[ ]:




