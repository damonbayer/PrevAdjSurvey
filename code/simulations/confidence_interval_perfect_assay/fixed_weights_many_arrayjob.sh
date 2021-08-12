#!/bin/bash

#SBATCH --partition=norm
#SBATCH --constraint=x2650
#SBATCH --ntasks=16
#SBATCH --ntasks-per-core=1
#SBATCH --time=2:00:00
#SBATCH --mem=32g
#SBATCH --array=1-2
#SBATCH --gres=lscratch:10
#SBATCH --job-name="PrevAdjSurvey"

module load R

Rscript --vanilla fixed_weights_many_arrayjob.R ${SLURM_ARRAY_TASK_ID}


