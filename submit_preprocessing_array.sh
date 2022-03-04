#!/bin/bash

#SBATCH -J preprocessing_array
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t 2-
#SBATCH -o array.%A_%a.log
#SBATCH -e array.%A_%a.err
#SBATCH --array=1-39
#SBATCH --mail-type=END
#SBATCH --mail-user=adam_turnbull@urmc.rochester.edu

subsid=($(find . -wholename "./BEEM*/*preprocessing*.sh" | sort -n))

source ${subsid[$SLURM_ARRAY_TASK_ID]}
