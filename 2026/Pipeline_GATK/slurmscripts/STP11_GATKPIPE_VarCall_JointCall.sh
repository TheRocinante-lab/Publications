#!/bin/bash
#SBATCH --job-name=ExtractSequences
#SBATCH --output=ExtractSequences_%j.log
#SBATCH --error=ExtractSequences_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load required modules
module load BCFtools/1.9-GCC-8.2.0-2.31.1

# Define variables
VCF_FILE="4_processing/GVCF/final_joint_called.vcf"
OUTPUT_DIR="extracted_sequences"
REGION="3:247434-251178"  # Define the region you are interested in
OUTPUT_FASTA="${OUTPUT_DIR}/extracted_sequences_${REGION}.fasta"

# Create the output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

# Extract the variant sequences from the VCF file for the specified region
echo "Extracting sequences for region ${REGION} from VCF file: ${VCF_FILE}"

# Process the VCF file
bcftools view -r $REGION $VCF_FILE | while read -r line; do
    # Skip headers
    if [[ $line == "#"* ]]; then
        continue
    fi
    
    # Extract fields from VCF line (CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE)
    CHROM=$(echo $line | cut -f 1)
    POS=$(echo $line | cut -f 2)
    REF=$(echo $line | cut -f 4)
    ALT=$(echo $line | cut -f 5)
    GENOTYPE=$(echo $line | cut -f 10 | cut -d':' -f1)
    
    # Print variant information for debugging
    echo "Processing variant at ${CHROM}:${POS} with REF=${REF}, ALT=${ALT}, GENOTYPE=${GENOTYPE}"
    
    # Handle genotype 0 (reference allele)
    if [[ "$GENOTYPE" == "0" ]]; then
        echo "Genotype 0 (reference) at ${POS}: ${REF}"
    
    # Handle genotype 1 (first ALT allele)
    elif [[ "$GENOTYPE" == "1" ]]; then
        ALT1=$(echo $ALT | cut -d',' -f1)
        echo "Genotype 1 (ALT) at ${POS}: ${ALT1}"
    
    # Handle <NON_REF> or missing genotypes
    elif [[ "$GENOTYPE" == "<NON_REF>" || "$GENOTYPE" == "." ]]; then
        echo "Warning: Unrecognized genotype at ${POS}, treating as missing"
    
    else
        echo "Warning: Unrecognized genotype ${GENOTYPE} at ${POS}"
    fi

done > $OUTPUT_FASTA

# Check if the sequence extraction was successful
if [[ -s ${OUTPUT_FASTA} ]]; then
    echo "Sequence extraction completed successfully. Output saved to: ${OUTPUT_FASTA}"
else
    echo "Error: No sequences were extracted. Please check if the region contains variants or if the VCF file is correct."
fi
