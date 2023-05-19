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
from ase import Atoms
from scipy.spatial.transform import Rotation as R
from ase.data import *
from ase.calculators.emt import EMT
from ase.optimize import BFGS
import copy


# In[2]:


def CoM(clus):
    """
    Support function to set the origin of the cluster at the centre of the mass
    """
    (cx, cy, cz) = clus.get_center_of_mass()
    new_xyz = []
    for i, a in enumerate(clus):
        x, y, z = a.position
        x -= cx
        y -= cy
        z -= cz
        new_xyz.append((x, y, z))
    clus.set_positions(new_xyz)
    return clus


# In[3]:


def find_unit_vec(clus,atom_idx1, atom_idx2):
    clus_pos = clus.get_positions()
    p1 = clus_pos[atom_idx1,:]
    p2 = clus_pos[atom_idx2,:]
    v = p2 - p1
    mod_v = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
    unit_v = v/mod_v
    return unit_v


# In[4]:


def append_new_atoms(clus, new_atom, new_atom_pos):
    '''
    Function: Append new atom to a given cluster
    '''
    clus_pos = clus.get_positions()
    new_pos = np.vstack((clus_pos, new_atom_pos))
    orginal_atom_list = clus.get_chemical_symbols()
    new_atom_list = orginal_atom_list + [new_atom]
    new_clus = Atoms(new_atom_list, new_pos)
    return new_clus


# In[5]:


def clus_after_removing_atoms(clus, index_list):
    clus_pos =  clus.get_positions()
    new_pos = np.delete(clus_pos, index_list, axis=0)
    orginal_atom_list = clus.get_chemical_symbols()
    
    for idx in sorted(index_list, reverse=True):
        del orginal_atom_list[idx]
    
    new_atom_list = orginal_atom_list #after removing the atoms at the indices
    new_clus = Atoms(new_atom_list, new_pos)
    return new_clus


# In[6]:


def remove_atoms(atom_pos, index_list):
    new_pos = np.delete(atom_pos, index_list, axis=0)
    return new_pos


# In[7]:


def rotate_coords(theta, p, rot_axis):
    
    Rx = np.matrix([[1,0,0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    Ry = np.matrix([[np.cos(theta), 0, np.sin(theta)], [0,1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    Rz = np.matrix([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    
    p_matrix = np.matrix(p) 
    p_mat_T = p_matrix.transpose()
    
    if rot_axis == 'x':
        R = Rx
    elif rot_axis == 'y':
        R = Ry
    else:
        R = Rz
                   
    rotated_coords = R*p_mat_T               
    rotated_coords_arr = np.array(rotated_coords).flatten()
    new_x, new_y, new_z = rotated_coords_arr
    
    new_coords  = (new_x, new_y, new_z)
    return new_coords


# In[8]:


def generate_clus_at_CoM(traj_fname,i, clus_size):
    traj = Trajectory(traj_fname)
    img = traj[i]
    
    write('initial_image.cif', img)
    
    img = CoM(img)
    pos = img.get_positions()
    atom_list = img.get_chemical_symbols()[:clus_size]
    
    #if not bare_cluster, adsorbate atoms will be removed 
    if len(img) != clus_size:
        index_list = np.arange(len(img))
        #print('index_list', index_list)
        remove_atom_index_list = index_list[clus_size:len(img)]
        #print(remove_atom_index_list)
        new_pos  = np.delete(pos, remove_atom_index_list, axis=0)
    else:
        new_pos = pos
    clus = Atoms(atom_list, new_pos)
    clus_pos = clus.get_positions()
    
    #adding a hypthetical atom at the center of mass; will remove later step
    com = img.get_center_of_mass()
    clus_COM_atom = append_new_atoms(clus,'N', com)
    
    #write('a_clus_COM.com', clus_COM_atom, format='gaussian-in')
    
    return clus_COM_atom


# In[9]:


def dist_from_CoM(clus_COM_atom):
    #Finding the atom positions from center of mass in descending order
    tot_atoms = len(clus_COM_atom) - 1 # -1 for added extra N atom at the center
    
    dist = clus_COM_atom.get_distances(tot_atoms, list(range(tot_atoms)))
    #print('dist', dist)
    dist_order = np.argsort(dist)[::-1]
    return dist_order


# In[10]:


def add_adsorbate_atom(clus, metal_idx, N_idx, d, new_atom):
    '''
   Function: Add at an given two atom positions at a distance 'd' from the second atom
    '''
    unit_vec = find_unit_vec(clus, metal_idx, N_idx)
    clus_pos = clus.get_positions()
 
    q1 =  clus_pos[N_idx,:] + unit_vec * d
    new_pos = np.vstack((clus_pos, q1))
    
    orginal_atom_list = clus.get_chemical_symbols()
    new_atom_list = orginal_atom_list + [new_atom]
    new_clus = Atoms(new_atom_list, new_pos)
    
    return new_clus


# In[11]:


def add_NH_fragment(clus, clus_size, idx_dist_order, d, d_NH, theta, rot_axis, N_idx, H_idx, add_Natom):
    
    com_idx = clus_size #CoM dummy atom position
    
    metal_idx = idx_dist_order # add N atom at the metal atom farthest from CoM. If already created N atom, it will pass
    if add_Natom == True:
        clus_N = add_adsorbate_atom(clus, com_idx, metal_idx, d, 'N')
        #write('a_clus_adsorb_N.com', clus_N, format='gaussian-in')
    else:
        clus_N = clus
    
    
    #add the H atom to the newly added/existing  N atom 
    new_atom = 'H'
    clus_NH = add_adsorbate_atom(clus_N, com_idx, N_idx, d_NH, 'H')
    #write('a_clus_adsorb_NH.com', clus_NH, format='gaussian-in')
 
    
    #Rotate the newly created NH bond
    clus_NH_pos = clus_NH.get_positions()
    p1 = clus_NH_pos[N_idx,:]
    p2 = clus_NH_pos[H_idx,:]
    p = p2 -p1
    
    theta = theta *np.pi/180
    NH_coords_after_rotation  = rotate_coords(theta, p, rot_axis)
    #print(NH_coords_after_rotation)
    new_coord_after_rotation = p2 + NH_coords_after_rotation
    
    #removing existing H atom coords at position and added the rotated H atom 
    clus_NH_removed_H = clus_after_removing_atoms(clus_NH, [H_idx])
    clus_NH_rotated = append_new_atoms(clus_NH_removed_H,'H', new_coord_after_rotation)
    
    #write('clus_rotate_coords.com', clus_NH_rotated, format='gaussian-in')
    
    return clus_NH_rotated


# In[12]:


def find_best_position_NH(clus, clus_size, d, idx_dist_order, rot_axis_list,theta_list,d_NH_list,  N_idx, H_idx, add_Natom ):
    
    order_dist = []
    order_theta = []
    order_axis = []
    order_dist_NH = []
    
    
    for axis in rot_axis_list:
        for theta in theta_list:
            for dist_NH in d_NH_list:
                clus_NH_rotated =  add_NH_fragment(clus, clus_size, idx_dist_order, d, dist_NH, theta, axis, N_idx, H_idx, add_Natom)
                                            
                
                dist = clus_NH_rotated.get_distance( N_idx, H_idx)
                order_dist.append(dist)
                order_theta.append(theta)
                order_axis.append(axis)
                order_dist_NH.append(dist_NH)
    
    diff_order = [abs(value- 1.0) for value in order_dist]
    arg_sort_diff_order = np.argsort(diff_order)
    
    optimal_index  = arg_sort_diff_order[0]
    
    optimal_dist = order_dist[optimal_index]
    optimal_theta = order_theta[optimal_index]
    optimal_axis = order_axis[optimal_index]
    optimal_d_NH = order_dist_NH[optimal_index]  
    
    #print(f' optimal_dist: {optimal_dist}, optimal_d_NH: {optimal_d_NH}, optimal_theta: {optimal_theta}, optimal_axis: {optimal_axis}') 
    
    #generate
    clus_NH_rotated_optimal =  add_NH_fragment(clus, clus_size, idx_dist_order, d, optimal_d_NH, optimal_theta, optimal_axis, N_idx, H_idx, add_Natom)
    
    cluster_composition = str(clus_NH_rotated_optimal.symbols)
    fname = 'clus_' + cluster_composition +'.com'
    #print('fname', fname)
    #write(fname, clus_NH_rotated_optimal, format='gaussian-in')
    
    return clus_NH_rotated_optimal, optimal_axis
   


# In[13]:


def relax(clus):
    clus.calc = EMT()
    dyn = BFGS(clus)
    dyn.run(fmax=0.05, steps=200)
    return clus


# In[14]:


def generate_NH2NH2_adsorbate(traj_filename, img_idx,  clus_size, dist_atom_idx):

    theta_list = [30, 60, 90, 150]
    d_NH_list = [0.6, 0.7, 0.8]
    rot_axis_list = ['x', 'y', 'z']
    d=2.0
    total_images = len(Trajectory(traj_filename))
    
    #print("total images in trajectory:", total_images)

    clus_COM_atom = generate_clus_at_CoM(traj_filename, img_idx, clus_size)

    atom_list = clus_COM_atom.get_chemical_symbols()
    covalent_radii_metal = covalent_radii[atomic_numbers[atom_list[0]]]
    covalent_radii_N = covalent_radii[atomic_numbers['N']]
    covalent_radii_H = covalent_radii[atomic_numbers['H']]
                                      
    d_MN1 = (covalent_radii_metal + covalent_radii_N) * 0.9 
    d_NN = 2* covalent_radii_N * 0.9
    d_MN2 = d_MN1 + d_NN


    dist_order = dist_from_CoM(clus_COM_atom)
    idx_dist_order = dist_order[dist_atom_idx] 
    
    com = clus_size
    clus_atom = idx_dist_order
    N1 = clus_size+1
    H1 = clus_size+2
    H2 = clus_size+3
    N2 = clus_size+4
    H3 = clus_size+5
    H4 = clus_size+6  

    #Generating NH 
    clus_NH_rotated_optimal1, optimal_axis1 = find_best_position_NH(clus_COM_atom, clus_size, d_MN1, idx_dist_order, rot_axis_list,                                                     theta_list, d_NH_list, N_idx=N1, H_idx=H1, add_Natom=True ) 
    
    #Generating NH2
    rot_axis_list.remove(optimal_axis1)
    clus_NH_rotated_optimal2, optimal_axis2 = find_best_position_NH(clus_NH_rotated_optimal1, clus_size, d_MN1, idx_dist_order,                       rot_axis_list, theta_list, d_NH_list, N_idx=N1, H_idx=H2, add_Natom=False ) 

    #Generate NH2NH
    theta_list = [30, 60, 90, 150]
    d_NH_list = [0.6, 0.7, 0.8]
    rot_axis_list = ['x', 'y', 'z']
    clus_NH_rotated_optimal3, optimal_axis3 = find_best_position_NH(clus_NH_rotated_optimal2, clus_size, d_MN2, idx_dist_order,                                                      rot_axis_list, theta_list, d_NH_list, N_idx=N2, H_idx=H3, add_Natom=True ) 

    #Generate NH2NH2
    rot_axis_list.remove(optimal_axis3)
    clus_NH_rotated_optimal4, optimal_axis4 = find_best_position_NH(clus_NH_rotated_optimal3, clus_size, d_MN1, idx_dist_order,                                                      rot_axis_list, theta_list, d_NH_list, N_idx=N2, H_idx=H4, add_Natom=False )
    return clus_NH_rotated_optimal4, dist_order, idx_dist_order


# In[15]:


def NH_bond_angle_rotations(clus_NH_rotated_optimal4, clus_size, idx_dist_order, angle):
    '''
    clus_NH_rotated_optimal4 will have additional N atom at the COM position  which will be removed later
    position of COM will be clus_size based on 0-indexing in python
    For clus_size 10: 
        position of COM N is at 10 (clus_size)
        First N atom position: 11 (clus_size+1)
        H attached to first N atoms positions: 12,13 (Clus_size+2, clus_size+3)
        Second N atom position: 14 (clus_size+4)
        H atoms attached to Second N atom position:15,16 (clus_size+5, clus_size+6)
    '''
    
    clus = copy.deepcopy(clus_NH_rotated_optimal4)
    
    com = clus_size
    clus_atom = idx_dist_order
    N1 = clus_size+1
    H1 = clus_size+2
    H2 = clus_size+3
    N2 = clus_size+4
    H3 = clus_size+5
    H4 = clus_size+6    

    clus.set_distance(H1, H2, 2.0) # moving two H-atoms by 2 Ang
    clus.set_distance(N1, H1, 1.0, fix=0, indices=[H1]) # fixing the NH bond distance to 1.0
    clus.set_distance(N1, H2, 1.0, fix=0, indices=[H2]) # fixing the NH bond distance to 1.0

    
    clus.set_distance(N1, N2, 1.4, fix=0, indices=[N2,H3,H4]) # fixing the N-N bond distance to 1.4

    clus.set_distance(H3,H4, 2.0)
    clus.set_distance(N2,H3, 1.0, fix=0, indices=[H3])
    clus.set_distance(N2,H4, 1.0, fix=0, indices=[H4])
    
    #fname_dist_move = "final_dist_NH2NH2_CoM_" + str(clus_atom) +   ".com"
    #write(fname_dist_move, clus, format='gaussian-in')
    
    
     
    clus_angle_copy = copy.deepcopy(clus)
    clus_angle_copy.set_angle(a1=com,a2=clus_atom, a3=N1, angle = angle, indices=[N1, H1,H2, N2, H3, H4])
    
    #fname_angle_move = "final_angle_NH2NH2_CoM_" + str(clus_atom) +  "_"+ str(angle) + ".com"
    #write(fname_angle_move, clus_angle_copy, format='gaussian-in')
    
    
    
    #delete CoM Extra N atom
    clus_dist_com = copy.deepcopy(clus)
    clus_dist = clus_after_removing_atoms(clus_dist_com, [clus_size])
    #fname_dist_CoM_removed = "final_distance_NH2NH2_" + str(clus_atom) 
    #write(fname_dist_CoM_removed + '.traj', clus_dist)
    #write(fname_dist_CoM_removed + '.com', clus_dist, format='gaussian-in')
    
    generate_adsorbates_from_NH2NH2(clus_dist,clus_size, clus_atom)
    
    clus_angle_com = copy.deepcopy(clus_angle_copy)
    clus_angle = clus_after_removing_atoms(clus_angle_com, [clus_size])
    #fname_angle_CoM_removed = "final_angle_NH2NH2_" + str(clus_atom) +  "_"+ str(angle)
    #write(fname_angle_CoM_removed + '.traj', clus_angle)
    #write(fname_angle_CoM_removed + '.com', clus_angle, format='gaussian-in')
    
    generate_adsorbates_from_NH2NH2_angle(clus_angle,clus_size, clus_atom, angle)

        


# In[16]:


def generate_adsorbates_from_NH2NH2(clus, clus_size, clus_atom_num):
    '''
    Arguments:
        clus: clus without CoM dummay tom
        file_name: base file name

    
    clus_NH_rotated_optimal4 will have additional N atom at the COM position  which will be removed later
    position of COM will be clus_size based on 0-indexing in python
    For clus_size 10: 
        index number will change as CoM dummy atom removed
        First N atom position: 10 (clus_size)
        H attached to first N atoms positions: 11,12 (Clus_size+1, clus_size+2)
        Second N atom position: 13 (clus_size+3)
        H atoms attached to Second N atom position:14,15 (clus_size+4, clus_size+5)
    '''    

    N1 = clus_size+0
    H1 = clus_size+1
    H2 = clus_size+2
    N2 = clus_size+3
    H3 = clus_size+4
    H4 = clus_size+5   
    
      
    folder_path = './Adsorbates/'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    adsorb_list = ['NH2NH2', 'NH2NH', 'NHNH2', 'NH2N', 'NNH2', 'NNH', 'NHN', 'NHNH', 'NH2', 'NH', 'N', 'N2']
    
    clus_NH2NH2 = clus
    #traj_path = os.path.join(folder_path, 'adsorb_NH2NH2.traj')
    #write(traj_path,clus_NH2NH2 )
    
    clus_NH2NH = clus_after_removing_atoms(clus, [H4])
    clus_NHNH2 = clus_after_removing_atoms(clus, [H2])
    
    clus_NH2N = clus_after_removing_atoms(clus, [H3,H4])    
    clus_NNH2 = clus_after_removing_atoms(clus, [H1,H2])
    
    clus_NNH = clus_after_removing_atoms(clus, [H1,H2, H3])
    clus_NHN = clus_after_removing_atoms(clus, [H2,H3, H4])
    
    clus_NHNH = clus_after_removing_atoms(clus, [H1,H3])
 
    clus_NH2 = clus_after_removing_atoms(clus, [N2, H3,H4])
    clus_NH = clus_after_removing_atoms(clus, [H2,N2, H3,H4])
    clus_N = clus_after_removing_atoms(clus, [H1, H2,N2, H3,H4])
    clus_N2 = clus_after_removing_atoms(clus, [H1, H2, H3,H4])
    
    adsorb_dict = {'NH2NH2': clus_NH2NH2,
                  'NH2NH': clus_NH2NH,
                  'NHNH2': clus_NHNH2,
                  'NH2N': clus_NH2N,
                  'NNH2': clus_NNH2,
                  'NNH': clus_NNH,
                  'NHN': clus_NHN,
                  'NHNH': clus_NHNH,
                  'NH2': clus_NH2,
                  'NH': clus_NH,
                  'N': clus_N,
                  'N2': clus_N2}
    
    for adsorb in adsorb_list:
        #print(adsorb)
        fname = 'adsorb_' + adsorb + '_' + str(clus_atom_num)
        ftraj = os.path.join(folder_path, fname + '.traj' )
        fcom = os.path.join(folder_path, fname + '.com' )
        
        clus_obj = adsorb_dict[adsorb]
                             
        write(ftraj, clus_obj)
        write(fcom, clus_obj, format='gaussian-in')
        
    


# In[17]:


def generate_adsorbates_from_NH2NH2_angle(clus, clus_size, clus_atom_num, angle):
    '''
    Arguments:
        clus: clus without CoM dummay tom
        file_name: base file name

    
    clus_NH_rotated_optimal4 will have additional N atom at the COM position  which will be removed later
    position of COM will be clus_size based on 0-indexing in python
    For clus_size 10: 
        index number will change as CoM dummy atom removed
        First N atom position: 10 (clus_size)
        H attached to first N atoms positions: 11,12 (Clus_size+1, clus_size+2)
        Second N atom position: 13 (clus_size+3)
        H atoms attached to Second N atom position:14,15 (clus_size+4, clus_size+5)
    '''    

    N1 = clus_size+0
    H1 = clus_size+1
    H2 = clus_size+2
    N2 = clus_size+3
    H3 = clus_size+4
    H4 = clus_size+5   
    
      
    folder_path = './Adsorbates_angle_'+ str(angle) +'/'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    adsorb_list = ['NH2NH2', 'NH2NH', 'NHNH2', 'NH2N', 'NNH2', 'NNH', 'NHN', 'NHNH', 'NH2', 'NH', 'N', 'N2']
    
    clus_NH2NH2 = copy.deepcopy(clus)

    clus_NH2NH = clus_after_removing_atoms(clus, [H4])
    clus_NHNH2 = clus_after_removing_atoms(clus, [H2])
    
    clus_NH2N = clus_after_removing_atoms(clus, [H3,H4])    
    clus_NNH2 = clus_after_removing_atoms(clus, [H1,H2])
    
    clus_NNH = clus_after_removing_atoms(clus, [H1,H2, H3])
    clus_NHN = clus_after_removing_atoms(clus, [H2,H3, H4])
    
    clus_NHNH = clus_after_removing_atoms(clus, [H1,H3])
 
    clus_NH2 = clus_after_removing_atoms(clus, [N2, H3,H4])
    clus_NH = clus_after_removing_atoms(clus, [H2,N2, H3,H4])
    clus_N = clus_after_removing_atoms(clus, [H1, H2,N2, H3,H4])
    clus_N2 = clus_after_removing_atoms(clus, [H1, H2, H3,H4])
    
    adsorb_dict = {'NH2NH2': clus_NH2NH2,
                  'NH2NH': clus_NH2NH,
                  'NHNH2': clus_NHNH2,
                  'NH2N': clus_NH2N,
                  'NNH2': clus_NNH2,
                  'NNH': clus_NNH,
                  'NHN': clus_NHN,
                  'NHNH': clus_NHNH,
                  'NH2': clus_NH2,
                  'NH': clus_NH,
                  'N': clus_N,
                  'N2': clus_N2}
    
    for adsorb in adsorb_list:
        #print(adsorb)
        fname = 'adsorb_' + adsorb + '_' + str(clus_atom_num) + "_"+ str(angle)
        ftraj = os.path.join(folder_path, fname + '.traj' )
        fcom = os.path.join(folder_path, fname + '.com' )
       
        clus_obj = adsorb_dict[adsorb]
                             
        write(ftraj, clus_obj)
        write(fcom, clus_obj, format='gaussian-in')
        
    


# In[18]:


traj_file_name = 'relaxed_final_geom.traj'
traj = Trajectory(traj_file_name)
print('total number of images in traj:',len(traj))


# In[19]:


range(10)


# In[20]:


clus_size=20
dist_atom_idx_list = [1,10,15,18]
#angle_list = [60,90, 150]
angle_list = [90]
img_idx = 0

for dist_atom_idx in range(clus_size):
#for dist_atom_idx in dist_atom_idx_list :
    for angle in angle_list:
        try:
            clus_NH_rotated_optimal4, dist_order, idx_dist_order = generate_NH2NH2_adsorbate(traj_file_name, img_idx, clus_size, dist_atom_idx)
            print(idx_dist_order)
            NH_bond_angle_rotations(clus_NH_rotated_optimal4, clus_size, idx_dist_order, angle)
        except:
            print("error: ", idx_dist_order)
            pass


print(dist_order)
