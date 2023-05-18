#!/bin/bash

metal='Cu'
for adsorb in  N2 N NH NH2 NH2NH2 NHNH2 NHNH  NNH2 NNH; do
	mkdir $adsorb
	mv random_traj_${adsorb}.traj $adsorb
	cp traj_bare_cluster_sp.* $adsorb
	cd $adsorb
		sed -i "s/adsorb='NH2NH2'/adsorb='${adsorb}'/g"	traj_bare_cluster_sp.sh
		sed -i "s/Cu/$metal/g"	traj_bare_cluster_sp.sh
		sbatch traj_bare_cluster_sp.sh
	cd ..
done



