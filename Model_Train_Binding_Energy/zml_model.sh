#!/bin/bash
#SBATCH --job-name=Wang
#SBATCH --ntasks=24
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --time=6:00:00
#SBATCH --mem 124G



module purge
module load bluebear


cd "$PBS_O_WORKDIR"

source /rds/projects/2018/johnston-copper-clusters-rr/Rajesh-2/anaconda3/etc/profile.d/conda.sh
conda activate dscribe 


	echo "Starting at "`date`
	python -u neural_network_reg.py >> results_model.out
	echo "Ending at "`date`
