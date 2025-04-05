#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=5000M
#SBATCH --partition=low_p

cd /data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/
mkdir -p filelist

export PATH_OUT="/data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/filelist/"
export FILE_OUT_PREFIX="filelist-"
find "${PATH_OUT}" -name "filelist-*" -type f -exec rm {} +
find "${PATH_OUT}" -name "filelist-*sh.e*" -type f -exec rm {} +
find "${PATH_OUT}" -name "filelist-*sh.o*" -type f -exec rm {} +
find "${PATH_OUT}" -name "slurm-*" -type f -exec rm {} +
find /data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/001_window/ -type f > filenames
export FILE_IN="/data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/filenames"

# focus on missing
DIR_FILES="/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/001_window/"
DIR_DONE="/data/MET_share/work/001_projects/proteins-crc/analysis/002_colocalisation/002_coloc/results/"
# Get list of .txt files (without extension)
mapfile -t LIST_FILES < <(find "$DIR_FILES" -maxdepth 1 -type f -name "*.txt" | sed 's/\.txt$//' | xargs -n 1 basename)
# Get list of directories (folder names only)
mapfile -t LIST_FILES_DONE < <(find "$DIR_DONE" -mindepth 1 -maxdepth 1 -type d | xargs -n 1 basename)
# Create a hash table (associative array) for fast lookups
declare -A DONE_SET
for dir in "${LIST_FILES_DONE[@]}"; do
    DONE_SET["$dir"]=1
done
# Filter files (skip those in DONE_SET)
LIST_FILES_FILTERED=()
for file in "${LIST_FILES[@]}"; do
    if [[ -z "${DONE_SET[$file]}" ]]; then
        LIST_FILES_FILTERED+=("${DIR_FILES}${file}.txt")
    fi
done
printf "%s\n" "${LIST_FILES_FILTERED[@]}" > filenames
export FILE_IN="/data/MET_share/work/001_projects/proteins-crc/code/002_colocalisation/002_coloc/filenames"

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
            print "#SBATCH --job-name=coloc_" fname;
        } else if (NR == 11) {
            print "VAR1=\"" fname "\"";
        } else {
            print $0;
        }
    }' 004_master.sh > "${sh_file}.sh"

    echo "${base_filename}.sh"

done < filenames

# Generate scripts
while IFS= read -r filepath; do
    # Extract the base filename
    base_filename=$(basename "$filepath")
    
    # Use awk to process the master script
    awk -v fname="$base_filename" '{
        if (NR == 19) {
            print "VAR_filelist=\"" fname "\"";
        } else {
            print $0;
        }
    }' 004_master.R > "${PATH_OUT}${base_filename}.R"

    echo "${base_filename}.R"

done < filenames
