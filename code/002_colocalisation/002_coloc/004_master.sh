#!/bin/bash

#SBATCH --job-name=test-coloc-parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00-05:00:00
#SBATCH --mem=10000M
#SBATCH --partition=low_p

VAR1=004_master

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/filelist/

time /opt/R/4.4.1/bin/Rscript ${VAR1}.R

