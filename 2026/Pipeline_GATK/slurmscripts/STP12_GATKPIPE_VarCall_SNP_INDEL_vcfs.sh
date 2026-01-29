#!/bin/bash
#SBATCH --job-name=VariantFiltering
#SBATCH --output=VariantFiltering_%j.log
#SBATCH --error=VariantFiltering_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load GATK module
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 

# Define variables
REF_GENOME="0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa"
GVCF_DIR="4_processing/GVCF"
FINAL_VCF="$GVCF_DIR/final_joint_called.vcf"
SNP_VCF="$GVCF_DIR/final_snps.vcf"
INDEL_VCF="$GVCF_DIR/final_indels.vcf"

# Select SNPs
gatk SelectVariants \
   -R $REF_GENOME \
   -V $FINAL_VCF \
   --select-type-to-include SNP \
   -O $SNP_VCF

echo "SNP selection completed."

# Select Indels
gatk SelectVariants \
   -R $REF_GENOME \
   -V $FINAL_VCF \
   --select-type-to-include INDEL \
   -O $INDEL_VCF

echo "Indel selection completed."

# Apply filters to SNPs
gatk VariantFiltration \
   -R $REF_GENOME \
   -V $SNP_VCF \
   --filter-expression "QD < 4.0" --filter-name "QD_filter" \
   --filter-expression "FS > 20.0" --filter-name "FS_filter" \
   --filter-expression "MQ < 55.0" --filter-name "MQ_filter" \
   --filter-expression "SOR > 3.0" --filter-name "SOR_filter" \
   --filter-expression "MQRankSum < -10.0" --filter-name "MQRankSum_filter" \
   --filter-expression "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum_filter" \
   -O $GVCF_DIR/filtered_snps.vcf


echo "SNP filtering completed."

# Apply  filters to Indels
gatk VariantFiltration \
   -R $REF_GENOME \
   -V $INDEL_VCF \
   --filter-expression "QD < 2.5" --filter-name "QD_filter" \
   --filter-expression "FS > 150.0" --filter-name "FS_filter" \
   --filter-expression "SOR > 8.0" --filter-name "SOR_filter" \
   -O $GVCF_DIR/filtered_indels.vcf

echo "Indel filtering completed."

echo "Job finished"