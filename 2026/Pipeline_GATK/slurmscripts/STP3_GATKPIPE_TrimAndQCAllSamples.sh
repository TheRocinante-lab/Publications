#!/bin/bash
#SBATCH --job-name=QC_report
#SBATCH --output=Trimmed_QC_report_%j.log
#SBATCH --error=Trimmed_QC_report_%j.err
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
# sbatch slurmscripts/STP3_GATKPIPE_TrimAndQCAllSamples.sh

# Load necessary modules
module load FastQC/0.11.8-Java-1.8
module load Python
module load Trimmomatic/0.38-Java-1.8

export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Check if MultiQC is installed, if not, install it
if ! [ -x "$(command -v multiqc)" ]; then
  pip install --user multiqc
fi

# Navigate to the directory for FastQC and Trimmomatic processing
cd 2_fastqc
mkdir -p trimmedsamples
mkdir -p FastQC_reports_Trimmed



export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Trimming with Trimmomatic
MAX_JOBS=8
running_jobs=0
for forward in ../1_data/fastq_set1/*/01.RawData/*/*_1.fq.gz ../1_data/fastq_set2/*/01.RawData/*/*_1.fq.gz; do
    reverse="${forward/_1.fq.gz/_2.fq.gz}"
    basename=$(basename "$forward" "_1.fq.gz")
    forward_out="trimmedsamples/${basename}_1_trimmed.fq.gz"
    reverse_out="trimmedsamples/${basename}_2_trimmed.fq.gz"
    forward_unpaired="${forward_out}_unpaired"
    reverse_unpaired="${reverse_out}_unpaired"

    # Ensure the path to Trimmomatic and its adapters is correct
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 4 -phred33 \
    "$forward" "$reverse" \
    "$forward_out" "$forward_unpaired" \
    "$reverse_out" "$reverse_unpaired" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1  # wait for 1 second before checking again
            running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
        done
    fi
done
wait

# FastQC analysis on trimmed FASTQ files
MAX_JOBS=8
running_jobs=0
for i in trimmedsamples/*.fq.gz; do
    fastqc "$i" -o FastQC_reports_Trimmed &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1  # wait for 1 second before checking again
            running_jobs=$(jobs -p | wc -l)  # update the count of running jobs
        done
    fi
done
wait

# Run MultiQC to compile FastQC reports into a single summary
multiqc FastQC_reports_Trimmed -o FastQC_reports_Trimmed

# Print the date and time and finish the script
date
echo "Job finished"
