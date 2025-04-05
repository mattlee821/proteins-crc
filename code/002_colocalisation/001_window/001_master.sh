#!/bin/bash

#SBATCH --job-name=window-parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=10-01:00:00
#SBATCH --mem=380000M
#SBATCH --partition=high_p

VAR1=001_master

export TMPDIR=/scratch/leem/temp/ 

cd /data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/001_window/

time /opt/R/4.4.1/bin/Rscript ${VAR1}.R
