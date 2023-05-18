#!/bin/bash
#SBATCH --ntasks=24
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --time=48:00:00
#SBATCH --mem 32G


module purge
module load bluebear


cd "$PBS_O_WORKDIR"

source /rds/projects/2018/johnston-copper-clusters-rr/Rajesh-2/anaconda3/etc/profile.d/conda.sh
conda activate dscribe 


	echo "Starting at "`date`
	python -u pca_nn_hp_tuning.py >> result_pca_nn_hp_tuning.out 
	echo "Ending at "`date`

