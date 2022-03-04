#!/bin/bash

#SBATCH -J first_level_subjectID
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t 2-
#SBATCH -o first_level/first_level.log
#SBATCH -e first_level/first_level.err
#SBATCH --mail-type=END
#SBATCH --mail-user=adam_turnbull@urmc.rochester.edu

feat subjectID_first_level.fsf
