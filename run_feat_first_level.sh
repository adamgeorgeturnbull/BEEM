#!/bin/bash

#SBATCH -J first_level_array
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t 2-
#SBATCH -o first_level_logs/array.%A_%a.log
#SBATCH -e first_level_logs/array.%A_%a.err
#SBATCH --array=0-39
#SBATCH --mail-type=END
#SBATCH --mail-user=adam_turnbull@urmc.rochester.edu

module load fsl/6.0.5.1
fsl-setup

subsid=($(find . -wholename "./BEEM*/bold/*first_level*.fsf" | sort -n))

feat ${subsid[$SLURM_ARRAY_TASK_ID]}
