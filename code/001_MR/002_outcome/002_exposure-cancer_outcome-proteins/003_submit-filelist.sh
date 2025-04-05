#!/bin/bash

#SBATCH --job-name=submit-jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=5000M
#SBATCH --partition=low_p

export DIRECTORY="/data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/002_exposure-cancer_outcome-proteins/filelist/"
cd ${DIRECTORY}

# Find all .sh files in the directory
sh_files=($DIRECTORY/*.sh)

# Set batch size
batch_size=100

# Get the total number of .sh files
total_files=${#sh_files[@]}

# Loop through the files in batches of 100
for ((i=0; i<total_files; i+=batch_size)); do
    # Submit the batch of 100 files
    for ((j=i; j<i+batch_size && j<total_files; j++)); do
        # Submit the .sh file
        qsub "${sh_files[$j]}"
    done

    # Optionally add a delay to avoid overloading the scheduler (e.g., 1-second pause)
    sleep 1
done
