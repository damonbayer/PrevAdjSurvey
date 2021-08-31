#!/bin/bash

#SBATCH --partition=norm
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mem=4g
#SBATCH --array=1-5
#SBATCH --gres=lscratch:10
#SBATCH --job-name="imperfect"

module load R

Rscript --vanilla imperfect_assay_arrayjob.R ${SLURM_ARRAY_TASK_ID}


