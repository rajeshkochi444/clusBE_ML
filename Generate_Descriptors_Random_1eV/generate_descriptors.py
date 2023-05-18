import numpy as np
import itertools
from collections import Counter
from dscribe.descriptors import ACSF, SOAP
from ase.io import Trajectory

def locate_descriptor_atoms(traj_file, n_descriptors):
    atomic_descriptors = []
    traj = Trajectory(traj_file)
    #print(len(traj))

    for i, img in enumerate(traj):
        #print(i)
        element_list = img.get_chemical_symbols()
        unique_elements = set(element_list)
        count = dict(Counter(element_list))
    
        if "N" in unique_elements:
            tot_N = count['N']
        else:
            tot_N = 0
            
        if "H" in unique_elements:
            tot_H = count['H']
        else:
            tot_H = 0
            
        if tot_N == 2:
            dist1 = img.get_distances(10, list(range(10)))
            dist2 = img.get_distances(11, list(range(10)))
            
            argsort1 = np.argsort(dist1)
            argsort2 = np.argsort(dist2)
            

            if dist1[argsort1[0]] < dist2[argsort2[0]]:
                descriptor_list = [10, argsort1[0]]
            else:
                descriptor_list = [11, argsort2[0]] # two-atom descriptors
            
            if n_descriptors == 1: # Onlythe closest N atom to metal cluster
                descriptor_list = [descriptor_list[0]]
                
            if n_descriptors > 2: #clostest N atom and the metal atom
                if 10 in descriptor_list:
                    descriptor_list.append(11)
                else:
                    descriptor_list.append(10)
            
            # for n_descriptors =3, we will add the second N atom
            # for n_descriptors > 3, we will add the metal atoms based on argsort list for the closte N atom to the cluster
            if n_descriptors > 3:
                if descriptor_list[0] == 10:
                    reference_argsort = argsort1
                else:
                    reference_argsort = argsort2
                        
                for i in range(1, len(reference_argsort)):
                    if len(descriptor_list) != n_descriptors:
                        descriptor_list.append(reference_argsort[i])

                    
        elif tot_N ==1 or tot_N == 0 :
            
            dist = img.get_distances(10, list(range(10)))
            arg_sort = np.argsort(dist).tolist()
            
            descriptor_list = [10, arg_sort[0]]
            if n_descriptors == 1:
                descriptor_list = [descriptor_list[0]]
            if n_descriptors > 2:  
                for i in range(1, len(arg_sort)):
                    if len(descriptor_list) != n_descriptors:
                        descriptor_list.append(arg_sort[i])

            
        else:
            print("Something Wrong. Check trajectory file")
            
        atomic_descriptors.append(descriptor_list)  
        #print("atomic_descriptors", atomic_descriptors)
    
    return atomic_descriptors



def generate_acsf_descriptor(traj_file,pos_list, metal_list):
    epsilon = [1,2,3,4,5,6,7]
    kappa = [0.5,1.0,1.5,2.0,2.5]
    eta = [0.01,0.03,0.06,0.1,0.2,0.4,1.0,2.5,5.0]
    R_s = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    lamb = [-1, 1]

    g2_params = [list(item) for item in itertools.product(eta, R_s)]
    g3_params = kappa
    g4_params = [list(item) for item in itertools.product(eta, epsilon, lamb)]
    #g5_params = g4_params

    # Set up: Instantiating an ACSF descripto
    species = [ 'H', "N"] + metal_list
    print('species:', species)
    acsf = ACSF(
        species=species,
        rcut=6.5,
        g2_params=g2_params,
        g3_params=g3_params,
        g4_params=g4_params,
        #g5_params=g5_params
    )

    # Create ACSF output for the system
    traj = Trajectory(traj_file)
    pos_list = pos_list
    acsf_traj  = acsf.create(traj, positions=pos_list, n_jobs=-1)
    
    return acsf_traj


def generate_soap_descriptor(traj_file,pos_list, metal_list):
    species = [ 'H', "N"] + metal_list
    r_cut = 6.5
    n_max = 9
    l_max = 10


    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=False,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
    )

    traj = Trajectory(traj_file)

    # Create ACSF output for the system
    pos_list = pos_list
    soap_traj  = soap.create(traj, positions=pos_list, n_jobs=-1)
    #print(soap_traj.shape)

    return soap_traj


