#!/bin/bash
#SBATCH --job-name=Align_to_reference
#SBATCH --output=Align_to_reference_%j.out
#SBATCH --error=Align_to_reference_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load necessary modules
module load BWA/0.7.17-GCC-8.2.0-2.31.1
module load SAMtools/1.16.1-GCC-11.3.0
module load Java/1.8.0_212
module load picard/2.26.10-Java-15
export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Directory containing the filtered trimmed samples
SAMPLE_DIR="2_fastqc/filtered_trimmedsamples"

# Reference genome file
REF_GENOME="/home/deniz/GenomicData/github/GATKPIPE2/0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"

# Index the reference genome for BWA
bwa index "$REF_GENOME"
echo "Reference genome indexed with BWA."

# Create sequence dictionary for the reference genome
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R="$REF_GENOME" O="${REF_GENOME%.fa}.dict"
echo "Sequence dictionary created for reference genome."

# Output directory for SAM files
mkdir -p 3_mapping/aligned_samples

MAX_JOBS=16
running_jobs=0

# Align the filtered trimmed reads to the reference genome
cd "$SAMPLE_DIR"
for forward in *_1_trimmed.fq.gz; do
    reverse="${forward/_1_trimmed.fq.gz/_2_trimmed.fq.gz}"
    base=$(basename "$forward" "_1_trimmed.fq.gz")
    sam_output="../../3_mapping/aligned_samples/${base}.sam"
    echo "Aligning: $forward and $reverse to $sam_output"

    bwa mem -K 100000000 -v 3 -t 1 -Y -R "@RG\tID:${base}\tLB:${base}\tPL:Illumina\tPM:machine\tSM:${base}" "$REF_GENOME" "$forward" "$reverse" > "$sam_output" &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
    fi
    echo "Alignment of $forward and $reverse finished"
done
wait

# Convert SAM to BAM, Sort, and Index
cd ../../3_mapping/aligned_samples

MAX_JOBS=16
running_jobs=0

for sam_file in *.sam; do
    bam_file="${sam_file%.sam}.bam"

    # Convert SAM to BAM
    samtools view -bS "$sam_file" > "$bam_file" &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
    fi
done
wait  # Wait for all SAM to BAM conversions to finish

echo "SAM to BAM conversion completed."

MAX_JOBS=16
running_jobs=0

# Sort BAM files
for bam_file in *.bam; do
    sorted_bam_file="${bam_file%.bam}_sorted.bam"

    # Sort the BAM file
    samtools sort "$bam_file" -o "$sorted_bam_file" &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
    fi
done
wait  # Wait for all BAM sorting to finish

echo "BAM sorting completed."

MAX_JOBS=16
running_jobs=0

# Index BAM files
for sorted_bam_file in *_sorted.bam; do
    samtools index "$sorted_bam_file" &
    ((running_jobs++))
    if [ "$running_jobs" -ge "$MAX_JOBS" ]; then
        while [ $running_jobs -ge $MAX_JOBS ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
    fi
done
wait  # Wait for all BAM indexing to finish

echo "BAM indexing completed."

# Print the date and time and finish the script
date
echo "Job finished"
