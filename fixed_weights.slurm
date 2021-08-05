#!/bin/bash

#SBATCH --partition=quick
#SBATCH --constraint=x2680
#SBATCH --ntasks=28
#SBATCH --ntasks-per-core=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=247g
#SBATCH --job-name="prevAdjSurvey"

module load R

Rscript --vanilla fixed_weights.R
