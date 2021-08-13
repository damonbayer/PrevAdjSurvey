#!/bin/bash

#SBATCH --partition=norm
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mem=4g
#SBATCH --array=120-250
#SBATCH --gres=lscratch:10
#SBATCH --job-name="PrevAdjSurvey"

module load R

Rscript --vanilla fixed_weights_many_arrayjob.R ${SLURM_ARRAY_TASK_ID}


