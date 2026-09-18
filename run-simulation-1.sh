#!/bin/bash
#SBATCH --output=par-%J.out
#SBATCH --ntasks=1 --cpus-per-task=36 --nodes=1
#SBATCH --time=02:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH -A lp_doctoralresearch

module load cluster/wice/batch
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=1

Rscript -e "if(!require('ellmer')) install.packages('ellmer')"
Rscript -e "if(!require('A4LEARN')) install.packages('A4LEARN.tar.gz')"
Rscript -e "renv::restore(exclude = c('A4LEARN'))"



Rscript R/simulations/data-generating-mechanism.R
