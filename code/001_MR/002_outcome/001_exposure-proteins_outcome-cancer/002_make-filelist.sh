#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=5000M
#SBATCH --partition=low_p

cd /data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/001_exposure-proteins_outcome-cancer/
mkdir -p filelist
mkdir -p /data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/
mkdir -p /data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/002_outcome/001_exposure-proteins_outcome-cancer/

export PATH_OUT="/data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/001_exposure-proteins_outcome-cancer/filelist/"
export FILE_OUT_PREFIX="filelist-"
find "${PATH_OUT}" -name "filelist-*" -type f -exec rm {} +
find "${PATH_OUT}" -name "slurm-*" -type f -exec rm {} +
# Define the directories
directories=(
"/data/GWAS_data/work/huyghe_2018_PMID30510241/GWAS/"
"/data/GWAS_data/work/laskar_2024_PMID38408508/GWAS/"
"/data/GWAS_data/work/fernandez-rozadilla_2022_PMID36539618/GWAS/"
)
# Initialize an empty array to store file paths
files=()
# Loop through each directory
for dir in "${directories[@]}"; do
    # Use find to get all files with the pattern "gz" and add them to the files array
    files+=($(find "$dir" -type f -name "*.gz"))
done
# Save the file paths to a text file
printf "%s\n" "${files[@]}" > filenames
export FILE_IN="/data/MET_share/work/001_projects/proteins-crc/code/001_MR/002_outcome/001_exposure-proteins_outcome-cancer/filenames"

# Create one file per line in FILE_IN
i=1
while IFS= read -r line; do
    # Format the number with leading zeros
    printf -v num "%02d" "$i"
    
    # Create the file
    echo "$line" > "${PATH_OUT}${FILE_OUT_PREFIX}${num}"
    
    # Print the filename to the screen
    echo "${FILE_OUT_PREFIX}${num}"

    # Increment the counter
    i=$((i + 1))
done < "$FILE_IN"

# Check the number of files is consistent  
wc -l ${FILE_IN}
find ${PATH_OUT} -name "filelist-*" -print0 | xargs -0 wc -l | tail -n 1
find ${PATH_OUT} -name "filelist-*" > filenames
wc -l filenames

# Create multiple .sh scripts based on a master script
while IFS= read -r sh_file; do 
    # Extract the base filename
    base_filename=$(basename "$sh_file" .sh)
    
    # Use awk to process the master SLURM script
    awk -v fname="$base_filename" '{
        if (NR == 3) {
            print "#SBATCH --job-name=protein-crc_" fname;
        } else if (NR == 10) {
            print "VAR1=\"" fname "\"";
        } else {
            print $0;
        }
    }' 001_master.sh > "${sh_file}.sh"

    echo "${base_filename}.sh"

done < filenames

# Generate scripts
while IFS= read -r filepath; do
    # Extract the base filename
    base_filename=$(basename "$filepath")
    
    # Use awk to process the master script
    awk -v fname="$base_filename" '{
        if (NR == 18) {
            print "VAR_filelist=\"" fname "\"";
        } else {
            print $0;
        }
    }' 001_master.R > "${PATH_OUT}${base_filename}.R"

    echo "${base_filename}.R"

done < filenames
