#!/bin/bash

for adsorb in  N2  NH NH2 NH2NH2 NHNH2 NHNH  NNH2 NNH; do
	cd $adsorb/VASP_$adsorb
	pwd
	cp ~/bin/NRR/random_out2traj_IE.py .
	rm -r Traj_OUTCAR
	python random_out2traj_IE.py > random_out2traj_IE_$adsorb.out
	cd ../../
done



