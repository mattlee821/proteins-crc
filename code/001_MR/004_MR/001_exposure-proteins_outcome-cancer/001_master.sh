#!/bin/bash

#SBATCH --job-name=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-02:0:00
#SBATCH --mem=5000M
#SBATCH --partition=low_p

VAR1=filelist-01

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc/code/001_MR/004_MR/001_exposure-proteins_outcome-cancer/filelist/

/opt/R/4.4.1/bin/Rscript ${VAR1}.R

