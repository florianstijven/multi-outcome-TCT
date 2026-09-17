#!/bin/bash
#SBATCH --output=par-%J.out
#SBATCH --ntasks=1 --cpus-per-task=36 --nodes=1
#SBATCH --time=02:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH -A lp_doctoralresearch

ml R

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

Rscript R/simulations/data-generating-mechanism.R
