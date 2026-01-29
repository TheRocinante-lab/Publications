#!/bin/bash
#SBATCH --job-name=ExtractVCFSequences
#SBATCH --output=ExtractVCFSequences_%j.log
#SBATCH --error=ExtractVCFSequences_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load required modules
module load BCFtools

# Define variables
VCF_DIR="4_processing/GVCF"
OUTPUT_DIR="extracted_sequences"
OUTPUT_DIR_CONSENSUS="${OUTPUT_DIR}/consensus_sequences"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR_CONSENSUS}

# Define the region to extract
REGIONS=("3:247434-251178")

# Initialize or clear region-specific consensus FASTA file
for REGION in "${REGIONS[@]}"; do
    IFS=':' read -ra ADDR <<< "$REGION"
    CHROM=${ADDR[0]}
    IFS='-' read -ra POS <<< "${ADDR[1]}"
    START=${POS[0]}
    END=${POS[1]}
    REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/vcf_sequences_${CHROM}_${START}_${END}.fasta"
    > ${REGION_FASTA_FILE} # Clear or initialize the file
    echo "Initialized VCF-based sequence FASTA file: ${REGION_FASTA_FILE}"
done

# Create an empty file to track processed VCF files
PROCESSED_VCF_FILES="${OUTPUT_DIR_CONSENSUS}/processed_vcf_files.txt"
> ${PROCESSED_VCF_FILES}

# Loop through VCF files and extract variant sequences for each strain
for VCF_FILE in $(ls ${VCF_DIR}/*.vcf.gz | grep -v 'sorted'); do
    STRAIN_NAME=$(basename ${VCF_FILE} | cut -d'.' -f1)
    echo "Processing VCF file: ${VCF_FILE} (Strain: ${STRAIN_NAME})"
    
    # Add the VCF file name to the processed file list
    echo "${VCF_FILE}" >> ${PROCESSED_VCF_FILES}

    for REGION in "${REGIONS[@]}"; do
        IFS=':' read -ra ADDR <<< "$REGION"
        CHROM=${ADDR[0]}
        IFS='-' read -ra POS <<< "${ADDR[1]}"
        START=${POS[0]}
        END=${POS[1]}
        REGION_STR="${CHROM}:${START}-${END}"
        REGION_FASTA_FILE="${OUTPUT_DIR_CONSENSUS}/vcf_sequences_${CHROM}_${START}_${END}.fasta"

        # Filter VCF for the specific region and extract sequences
        echo "Extracting variant information for region: ${REGION_STR} in strain: ${STRAIN_NAME}"

        # Extract sequence based on genotype, removing "," and "<NON_REF>"
        bcftools query -r ${REGION_STR} -f '[%CHROM:%POS %REF %ALT %GT\n]' ${VCF_FILE} > temp_variants.txt

        # Check if any variants were found
        VARIANT_COUNT=$(wc -l < temp_variants.txt)
        echo "Found ${VARIANT_COUNT} variants for strain ${STRAIN_NAME} in region ${REGION_STR}"

        # If no variants, skip this strain
        if [[ "${VARIANT_COUNT}" -eq 0 ]]; then
            echo "No variants found for strain ${STRAIN_NAME} in region ${REGION_STR}. Skipping."
            continue
        fi

        # Build sequence for the strain by interpreting the genotypes
        SEQUENCE=""

        while read -r LINE; do
            CHROM_POS=$(echo "$LINE" | cut -d ' ' -f 1)
            REF=$(echo "$LINE" | cut -d ' ' -f 2)
            ALT=$(echo "$LINE" | cut -d ' ' -f 3 | sed 's/,<NON_REF>//g')  # Remove "," and "<NON_REF>"
            GENOTYPE=$(echo "$LINE" | cut -d ' ' -f 4)

            echo "Processing variant at ${CHROM_POS} with REF=${REF}, ALT=${ALT}, GENOTYPE=${GENOTYPE}"

            # Use REF for homozygous reference (0/0), ALT for homozygous alternate (1/1), and REF for heterozygous (0/1)
            if [[ "$GENOTYPE" == "0" ]]; then
                SEQUENCE+="$REF"
            elif [[ "$GENOTYPE" == "1" ]]; then
                SEQUENCE+="$ALT"
            else
                SEQUENCE+="$REF"  # For any unknown or heterozygous cases
            fi
        done < temp_variants.txt

        # Ensure that the header with the sample ID is written before the sequence
        if [[ -n "$SEQUENCE" ]]; then
            echo ">${STRAIN_NAME}" >> ${REGION_FASTA_FILE}
            echo "$SEQUENCE" >> ${REGION_FASTA_FILE}
            echo "Appended sequence for strain ${STRAIN_NAME} to ${REGION_FASTA_FILE}"
        else
            echo "No sequence generated for strain ${STRAIN_NAME}, skipping..."
        fi

        # Clean up
        rm temp_variants.txt
    done
done

# Display the list of processed VCF files
echo "The following VCF files have been processed:"
cat ${PROCESSED_VCF_FILES}

echo "Sequence extraction from VCF completed."
