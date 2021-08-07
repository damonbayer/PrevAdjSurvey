#!/bin/bash

#SBATCH --partition=multinode
#SBATCH --constraint=x2680
#SBATCH --nodes=10
#SBATCH --ntasks=280
#SBATCH --ntasks-per-core=1
#SBATCH --time=1-00:00:00
#SBATCH --exclusive
#SBATCH --mem=247g
#SBATCH --gres=lscratch:10
#SBATCH --job-name="prev100"
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=damonmbayer@gmail.com

module load R

mpiexec Rscript --vanilla fixed_weights_doMPI.R
