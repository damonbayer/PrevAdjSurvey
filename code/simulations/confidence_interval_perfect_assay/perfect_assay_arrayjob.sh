#!/bin/bash

#SBATCH --partition=norm,quick
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mem=4g
#SBATCH --array=205
#SBATCH --gres=lscratch:10
#SBATCH --job-name="perfect"

module load R

Rscript --vanilla perfect_assay_arrayjob.R ${SLURM_ARRAY_TASK_ID}


