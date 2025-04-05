#!/bin/bash

#SBATCH --job-name=instruments-cancer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=20000M
#SBATCH --partition=high_p

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc

time /opt/R/4.4.1/bin/Rscript code/001_MR/001_instruments/002_instruments-cancer/001_instruments-cancer.R
