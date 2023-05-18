#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
#import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


ene_mol_H2  = -6.77255467
ene_mol_N2  = -16.6419041
ene_mol_NH3 = -19.54417378
ene_mol_H = 0.5 * ene_mol_H2

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

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']
all_traj_adsorb_img_list = [] 
for adsorb in adsorb_list:
    traj_file_adsorb = 'random_traj_adsorb_2k_' +  adsorb + '.traj'
    traj_file_bareclus = 'random_traj_bareclus_2k_' +  adsorb + '.traj'
    traj_adsorb = Trajectory(traj_file_adsorb)
    traj_bareclus = Trajectory(traj_file_bareclus)
    print(traj_file_adsorb, len(traj_adsorb))
    print(traj_file_bareclus, len(traj_bareclus))
    
    
    ene_list_adsorb = []
    ene_list_bareclus = []   
    binding_ene_list = []
    
    n_NH3 = binding_energy_dict[adsorb][0]
    n_H = binding_energy_dict[adsorb][1]
    print(adsorb, n_NH3, n_H)
    
    for i in range(len(traj_bareclus)): 
        #if i == 2000:
            #break
        img_adsorb = traj_adsorb[i]
        all_traj_adsorb_img_list.append(img_adsorb)
        ene_adsorb = traj_adsorb[i].get_potential_energy()
        ene_list_adsorb.append(ene_adsorb)
        ene_bareclus = traj_bareclus[i].get_potential_energy()
        ene_list_bareclus.append(ene_bareclus)
        
        binding_ene = ene_adsorb + (n_NH3 * ene_mol_NH3) - ene_bareclus - ene_mol_N2 - (n_H * ene_mol_H)
        binding_ene_list.append(binding_ene)
        
    arr_ene_adsorb = np.array(ene_list_adsorb)
    arr_ene_bareclus = np.array(ene_list_bareclus)
    arr_binding_ene = np.array(binding_ene_list)         
    
    print(np.max(arr_ene_adsorb), np.min(arr_ene_adsorb), np.mean(arr_ene_adsorb))
    print(np.max(arr_ene_bareclus), np.min(arr_ene_bareclus), np.mean(arr_ene_bareclus))
    print(np.max(arr_binding_ene), np.min(arr_binding_ene), np.mean(arr_binding_ene))
    
   
    y = np.array(binding_ene_list)
    print(y.shape)

    save_y_fname = 'y_' + adsorb + '.npy'
    np.save(save_y_fname, y)

write('all_adsorb_single_traj_2k.traj', all_traj_adsorb_img_list)
all_img_traj = Trajectory('all_adsorb_single_traj_2k.traj')
print('tot images in all_adsorb_single_traj_2k.traj:', len(all_img_traj))
# In[3]:



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

print("total IE values in y_all_adsorb.npy")
print(np.load('y_all_adsorb.npy').shape)


# In[5]:


#for adsorb in adsorb_list:
    #fname = 'y_' + adsorb + '.npy'
    #be_list = np.load(fname)
    #print(adsorb, len(be_list), np.max(be_list), np.min(be_list), np.mean(be_list))
    #sns.displot(be_list)
    #plt.xlim(min_ene-0.5,min_ene+2.5)
    #plt.show()


# In[6]:

'''
def be_list_data(adsorb):
    fname = 'y_' + adsorb + '.npy'
    be_list = np.load(fname)
    return be_list


# In[7]:


nrows = 3
ncols = 3
fig, axes = plt.subplots(nrows=3,ncols=3, figsize=(20,20))

row_col = [(i,j) for i in range(nrows) for j in range(ncols)]
print(row_col)

for i, adsorb in enumerate(adsorb_list):
    row_num, col_num = row_col[i][0], row_col[i][1]
    be_list = be_list_data(adsorb)
    sns.histplot(ax=axes[row_num, col_num],x = be_list)
    axes[row_num,col_num].set_title(adsorb, fontsize=16)
    axes[row_num,col_num].set_xlabel('B.E(eV)', fontsize=16)
plt.savefig('2k_data_dist.png')

# In[8]:


be_list = np.load('y_all_adsorb.npy')
sns.histplot(be_list, kde=True)
plt.savefig('2k_data_dist_all.png')

# In[ ]:

'''


