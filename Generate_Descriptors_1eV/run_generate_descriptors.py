import numpy as np
from generate_descriptors import locate_descriptor_atoms, generate_acsf_descriptor, generate_soap_descriptor

descriptor_type = 'acsf' #or soap
metal_list = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
n_descriptors = 1
metal_atom = 'Cu'

adsorb_list = ['N2', 'N', 'NH', 'NH2', 'NH2NH2', 'NHNH', 'NHNH2', 'NNH', 'NNH2']

for  adsorb in adsorb_list:

    traj_file = '../Trajs_Filtered_1eV/filtered_traj_1eV_' + adsorb + '.traj'
    print(traj_file)
    atom_pos_descriptors = locate_descriptor_atoms(traj_file, n_descriptors)
    if descriptor_type == 'acsf':
        descriptor_name = descriptor_type+ '_descriptor'
        descriptor_name = generate_acsf_descriptor(traj_file, atom_pos_descriptors, metal_list)
    elif  descriptor_type == 'soap':
        descriptor_name = descriptor_type+ '_descriptor'
        descriptor_name = generate_soap_descriptor(traj_file, atom_pos_descriptors, metal_list)
    
    print(descriptor_name.shape)
    
    save_descriptor_fname = descriptor_type + '_ndes'+ str(n_descriptors) + '_1eV_' + adsorb + '.npy'    
    np.save(save_descriptor_fname, descriptor_name)
print("Completed descriptor generation")

print("Stacking descripor arrays")
for adsorb in adsorb_list:
    fname = descriptor_type + '_ndes'+ str(n_descriptors) + '_1eV_' + adsorb + '.npy'
    if adsorb == 'N2':
        a = np.load(fname)
        print(a.shape)
    else:
        b = np.load(fname)
        print(b.shape)
        a = np.vstack((a,b))
        print("stacked a", a.shape)

all_descriptors_single_npy = descriptor_type + '_ndes'+ str(n_descriptors) +  '_all_adsorb_1eV_' + metal_atom + '.npy'
np.save(all_descriptors_single_npy, a)

print(np.load(all_descriptors_single_npy).shape)
