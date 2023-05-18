#!/bin/bash
#SBATCH --job-name=n2ads
#SBATCH --ntasks=24
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu 8G


module purge
module load bluebear


cd "$PBS_O_WORKDIR"
module load VASP/5.4.4-iomkl-2020b


adsorb='NH2NH2'
metal='Cu'

mkdir VASP_$adsorb
mkdir VASP_$adsorb/OUTCARS_List

sed -i "s/Cu/$metal/g" traj_bare_cluster_sp.py 
sed -i "s/NH2NH2/$adsorb/g" traj_bare_cluster_sp.py 

for i in {0..3000}; do
	echo $i 
	sed -i "s/img_idx = 0/img_idx = $i/g" traj_bare_cluster_sp.py
	cp /rds/projects/2018/johnston-copper-clusters-rr/Rajesh-2/potpaw_PBE_54/$metal/POTCAR VASP_$adsorb 
	cp ~/bin/NRR/INCAR_SP VASP_$adsorb/INCAR 
	cp ~/bin/NRR/KPOINTS VASP_$adsorb
        python traj_bare_cluster_sp.py
	cd VASP_$adsorb
		echo "Starting calculation in $adsorb ", `date`
		mpirun vasp_gam > vasp.out
		echo "Finished calculation in $adsorb ", `date`
		mv OUTCAR OUTCARS_List/OUTCAR_$i
	cd ..	
	sed -i "s/img_idx = $i/img_idx = 0/g" traj_bare_cluster_sp.py
done



