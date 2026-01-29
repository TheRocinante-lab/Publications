#!/bin/bash
#SBATCH --job-name=PrepareGVCFs
#SBATCH --output=PrepareGVCFs_%j.log
#SBATCH --error=PrepareGVCFs_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Assuming your base directory is the parent directory of '4_processing'
BASE_DIR="/home/deniz/GenomicData/github/GATKPIPE2"  # Update this path if needed
REF_DIR="/home/deniz/GenomicData/github/GATKPIPE2/0_index/referenceIPO323"  # Adjusted to the detailed path
REF_GENOME="$REF_DIR/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
GVCF_DIR="/home/deniz/GenomicData/github/GATKPIPE2/4_processing/GVCF"

# Load modules
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 
module load Python

# Check if the GVCF directory exists
if [ ! -d "$GVCF_DIR" ]; then
    echo "GVCF directory does not exist: $GVCF_DIR"
    exit 1
fi

# Navigate to the GVCF directory
cd "$GVCF_DIR"

# Remove the existing sample map file if it exists
rm -f zymoseptoria.sample_map

# Check if there are any GVCF files
if ls *.g.vcf.gz 1> /dev/null 2>&1; then
    # Create a new sample map file with absolute paths
    for gvcf in *.g.vcf.gz; do
        sample_name=$(basename "$gvcf" .g.vcf.gz)
        echo -e "$sample_name\t$GVCF_DIR/$gvcf" >> zymoseptoria.sample_map
    done
else
    echo "No GVCF files found in $GVCF_DIR"
    exit 1
fi

# Check for the existence of the reference genome index file
if [ ! -f "$REF_GENOME.fai" ]; then
    echo "Reference genome index file not found: $REF_GENOME.fai"
    exit 1
fi

# Create interval list
awk '{print $1 FS "1" FS $2}' "$REF_GENOME.fai" > "$GVCF_DIR/interval.list"

echo "Preparation and GVCF File Generation Completed"
