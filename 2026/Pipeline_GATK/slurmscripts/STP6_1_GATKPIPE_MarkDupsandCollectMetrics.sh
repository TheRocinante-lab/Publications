#!/bin/bash
#SBATCH --job-name=MarkDuplicates_Metrics
#SBATCH --output=MarkDuplicates_Metrics_%j.out
#SBATCH --error=MarkDuplicates_Metrics_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load necessary modules
module load Java/1.8.0_212
module load picard/2.26.10-Java-15
module load SAMtools/1.16.1-GCC-11.3.0
module load R/4.2.0-foss-2021b
module load Python

# Install Python packages
pip install --user numpy --upgrade

export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Reference genome
ref="/home/deniz/GenomicData/github/GATKPIPE2/0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"

# Directory for processed BAM files
mkdir -p 4_processing/Metrics

# Change directory to where the BAM files are located
cd 3_mapping/aligned_samples

# MarkDuplicates
MAX_JOBS=16
running_jobs=0
for bam in *_sorted.bam; do  # Change here to iterate over sorted BAM files
    name=$(basename "$bam" _sorted.bam)
    sorted_md_bam="${name}_sorted_md.bam"
    md_metrics="../../4_processing/Metrics/${name}_md_metrics.txt"

    # Mark Duplicates
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$bam O=$sorted_md_bam M=$md_metrics

    # Index the final BAM file
    samtools index $sorted_md_bam &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
    fi
done
wait


echo "MarkDuplicates processing completed."

# Collect Metrics
for bam_file in *_sorted_md.bam; do  # Change here to iterate over sorted-marked-duplicate BAM files
    base=$(basename "${bam_file}" "_sorted_md.bam")

    # Collect Alignment Summary Metrics
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
        -R $ref \
        -I "$bam_file" \
        -O "../../4_processing/Metrics/${base}_alignment_metrics.txt"
  
    # Collect Insert Size Metrics
    java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
        -I "$bam_file" \
        -O "../../4_processing/Metrics/${base}_insert_metrics.txt" \
        -H "../../4_processing/Metrics/${base}_insert_size_histogram.pdf"
  
    # Depth of Coverage
    samtools depth -a "$bam_file" > "../../4_processing/Metrics/${base}_depth_out.txt"
done


# Print the date and time and finish the script
date
echo "Job finished"
