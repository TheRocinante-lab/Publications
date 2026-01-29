#!/bin/bash
#SBATCH --job-name=Filtered_Trimmed_FastQC
#SBATCH --output=Filtered_Trimmed_FastQC_%j.log
#SBATCH --error=Filtered_Trimmed_FastQC_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load necessary modules
module load FastQC/0.11.8-Java-1.8
module load Python

# Check if MultiQC is installed, if not, install it
if ! [ -x "$(command -v multiqc)" ]; then
  pip install --user multiqc
fi
export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Navigate to the FastQC directory
cd 2_fastqc

# Create a new directory for the filtered trimmed FastQC reports
mkdir -p FastQC_reports_Trimmed_Filtered

# Define maximum parallel jobs
MAX_JOBS=8

# FastQC analysis on trimmed FASTQ files, excluding specific samples
running_jobs=0
for i in trimmedsamples/*.fq.gz; do
    # Skip files that start with the specified prefixes
    if [[ $i != *"S101_EKDN220047288-1A"* ]]; then
        fastqc "$i" -o FastQC_reports_Trimmed_Filtered &
        ((running_jobs++))
        if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
            while [ $running_jobs -ge $MAX_JOBS ]; do
                sleep 1  # wait for 1 second before checking again
                running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
            done
        fi
    fi
done
wait

# Run MultiQC to compile FastQC reports into a single summary
multiqc FastQC_reports_Trimmed_Filtered -o FastQC_reports_Trimmed_Filtered_Summary

echo "MultiQC Job finished"

cd ..

# Directory containing the original trimmed samples
SOURCE_DIR="2_fastqc/trimmedsamples"

# Directory to store the filtered trimmed samples
DEST_DIR="2_fastqc/filtered_trimmedsamples"
mkdir -p "$DEST_DIR"

# Array of sample prefixes to exclude
EXCLUDE_SAMPLES=("S101_EKDN220047288-1A")

# Function to check if a file name starts with any of the excluded prefixes
function is_excluded {
    for prefix in "${EXCLUDE_SAMPLES[@]}"; do
        if [[ $1 == "$prefix"* ]]; then
            return 0 # 0 indicates true in bash, meaning the file should be excluded
        fi
    done
    return 1 # False, the file should not be excluded
}

# Iterate over all .fq.gz files in the source directory
for file in "$SOURCE_DIR"/*.fq.gz; do
    filename=$(basename "$file")
    if ! is_excluded "$filename"; then
        # Copy the file to the destination directory if it's not excluded
        cp "$file" "$DEST_DIR/"
    fi
done

echo "Filtering and copying of trimmed samples completed."

# Print the date and time and finish the script
date
echo "Job finished"


