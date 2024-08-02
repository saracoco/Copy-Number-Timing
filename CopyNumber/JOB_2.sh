#!/bin/bash
#SBATCH --job-name=BM
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=05:00:00
#SBATCH --partition=THIN
#SBATCH --mem=200gb

module load R

R CMD BATCH scripts/simulate_2.R


module purge
