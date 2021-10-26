#!/bin/bash

#SBATCH --partition=norm
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mem=4g
#SBATCH --array=1-1000
#SBATCH --gres=lscratch:10
#SBATCH --job-name="imperfect"

module load R

Rscript --vanilla imperfect_assay_arrayjob_perfect_test.R ${SLURM_ARRAY_TASK_ID}


