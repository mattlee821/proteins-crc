#!/bin/bash

#SBATCH --job-name=004_missing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-5:0:00
#SBATCH --mem=15000M
#SBATCH --partition=low_p

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/

/opt/R/4.4.1/bin/Rscript 004_identify-missing.R



