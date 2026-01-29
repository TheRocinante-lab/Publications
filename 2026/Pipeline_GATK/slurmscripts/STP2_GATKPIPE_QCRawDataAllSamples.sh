#!/bin/bash
#SBATCH --job-name=QC_report
#SBATCH --output=QC_report_%j.log
#SBATCH --error=QC_report_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# to run the script: 
# loadmod
# sbatch slurmscripts/STP2_GATKPIPE_QCRawDataAllSamples.sh

export _JAVA_OPTIONS="-Xms256m -Xmx60g"


module load FastQC/0.11.8-Java-1.8
module load Python
# Install MultiQC in the user's local Python environment
pip install --user multiqc

# Navigate to the directory where FastQC reports will be generated
cd 2_fastqc
mkdir -p FastQC_reports_Raw

MAX_JOBS=8  # Control the maximum number of parallel jobs
running_jobs=0

# FastQC analysis on extracted FASTQ files
# The fastq files are located inside each sample subdirectory under '01.RawData'
for set_dir in ../1_data/fastq_set1 ../1_data/fastq_set2
do
    for sample_dir in "$set_dir"/*/01.RawData/*
    do
        for fq_file in "$sample_dir"/*.fq.gz
        do
            fastqc "$fq_file" -o FastQC_reports_Raw &
            ((running_jobs++))
            if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
                while [ $running_jobs -ge $MAX_JOBS ]; do
                    sleep 1  # wait for 1 second before checking again
                    running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
                done
            fi
        done
    done
done
wait

# Run MultiQC to compile FastQC reports into a single summary
~/.local/bin/multiqc --export FastQC_reports_Raw -o FastQC_reports_Raw

# print the date and time and finish the script
date
echo "MultiQC Job finished"



# Output CSV file name
output_csv="sample_fastq_files.csv"

# Add CSV header
echo "Set Name,Sample Name,Fastq File 1,Fastq File 2" > "$output_csv"

# Iterate over set directories
for set_dir in ../1_data/fastq_set1 ../1_data/fastq_set2; do
    set_name=$(basename "$set_dir")  # Extract set name

    # Iterate over sample directories
    for sample_dir in "$set_dir"/*/01.RawData/*; do
        sample_name=$(basename "$sample_dir")  # Extract sample name

        # Find fastq files
        fastq_files=($(find "$sample_dir" -type f -name "*.fq.gz"))
        if [ ${#fastq_files[@]} -eq 2 ]; then
            # If there are exactly 2 fastq files, add their names to the CSV
            echo "$set_name,$sample_name,${fastq_files[0]},${fastq_files[1]}" >> "$output_csv"
        fi
    done
done

echo "CSV file created: $output_csv"

