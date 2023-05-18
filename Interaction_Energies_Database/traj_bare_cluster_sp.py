import ase
import os
import glob
import numpy as np
from ase.io import read, write, Trajectory
import subprocess

adsorb = 'NH2NH2'
metal = 'Cu'
calc_folder  = 'VASP_'+adsorb
img_idx = 0

traj_file = 'random_traj_' + adsorb + '.traj'
traj = Trajectory(traj_file)
tot_img = len(traj)
print(traj_file, tot_img)

img =  traj[img_idx]
del img[[atom.index for atom in img if atom.symbol != metal]]
print(img)
write(calc_folder + '/POSCAR', img, format='vasp')
#if not os.path.exists(calc_folder):
#cd  os.makedirs(calc_folder)
#write(calc_folder + '/POSCAR_'+str(img_idx), img, format='vasp')
#potcar_copy  = 'cp /rds/projects/2018/johnston-copper-clusters-rr/Rajesh-2/potpaw_PBE_54/' + metal + '/POTCAR' + ' ' +  calc_folder+'/'
#subprocess.run(potcar_copy, shell=True)

