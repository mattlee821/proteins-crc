#!/bin/bash

#SBATCH --job-name=master
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-02:0:00
#SBATCH --mem=10000M
#SBATCH --partition=low_p

VAR1=001_master

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc/code/001_MR/003_harmonise/001_exposure-proteins_outcome-cancer/filelist/

/opt/R/4.4.1/bin/Rscript ${VAR1}.R

