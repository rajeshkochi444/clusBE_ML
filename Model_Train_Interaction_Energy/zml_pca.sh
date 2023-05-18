#!/bin/bash
#SBATCH --job-name=Wang
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --nodes 1
#SBATCH --mem-per-cpu 6900M
#SBATCH --ntasks 24


module purge
module load bluebear


cd "$PBS_O_WORKDIR"

source /rds/projects/2018/johnston-copper-clusters-rr/Rajesh-2/anaconda3/etc/profile.d/conda.sh
conda activate dscribe 


	echo "Starting at "`date`
	#python -u pca.py >> results_pca.out
	python -u pca_999.py >> results_pca_999.out
	echo "Ending at "`date`
