#!/bin/bash
#SBATCH --job-name=QC_report
#SBATCH --output=MQC_aft_MDUP_report_%j.log
#SBATCH --error=MQC_aft_MDUP_report_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# to run the script: 
# loadmod
# sbatch slurmscripts/STP2_GATKPIPE_QCRawDataAllSamples.sh

module load FastQC/0.11.8-Java-1.8
module load Python
# Install MultiQC in the user's local Python environment
pip install --user multiqc



# Run MultiQC to compile reports
multiqc 4_processing/Metrics

echo "Metrics collection and MultiQC analysis completed."

# Print the date and time and finish the script
date
echo "Job finished"
