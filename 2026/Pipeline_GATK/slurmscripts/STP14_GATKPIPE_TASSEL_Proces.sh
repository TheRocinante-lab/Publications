#!/bin/bash
#SBATCH --job-name=TASSELAnalysis
#SBATCH --output=TASSELAnalysis_%j.log
#SBATCH --error=TASSELAnalysis_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="your_email_here"

# Load required modules
module load TASSEL/5.2.86-Java-1.8.0_212
module load BCFtools/1.9-GCC-8.2.0-2.31.1

# Increase TASSEL's memory allocation
export _JAVA_OPTIONS="-Xms256m -Xmx60g"

# Set paths
ANNOTATED_SNPS_VCF="5_annotation/snpEff_snps/geno_filtered_snps.ann.vcf.gz"
ANNOTATED_INDELS_VCF="5_annotation/snpEff_indels/geno_filtered_indels.ann.vcf.gz"
OUT_DIR="6_tassel_analysis"
INDEL_OUT_DIR="6_tassel_analysis_indels"

# Create output directories
mkdir -p $OUT_DIR
mkdir -p $INDEL_OUT_DIR

# Remove samples with "sorted" in their names from SNP VCF
bcftools view -s ^$(bcftools query -l $ANNOTATED_SNPS_VCF | grep "sorted" | paste -sd,) $ANNOTATED_SNPS_VCF -Oz -o $OUT_DIR/filtered_snps.vcf.gz
bcftools index $OUT_DIR/filtered_snps.vcf.gz

# Remove samples with "sorted" in their names from indel VCF
bcftools view -s ^$(bcftools query -l $ANNOTATED_INDELS_VCF | grep "sorted" | paste -sd,) $ANNOTATED_INDELS_VCF -Oz -o $INDEL_OUT_DIR/filtered_indels.vcf.gz
bcftools index $INDEL_OUT_DIR/filtered_indels.vcf.gz

# Convert filtered SNP VCF to TASSEL format
gunzip -c $OUT_DIR/filtered_snps.vcf.gz > $OUT_DIR/temp_snps.vcf
run_pipeline.pl -Xmx60g -importGuess $OUT_DIR/temp_snps.vcf -export $OUT_DIR/tassel
rm $OUT_DIR/temp_snps.vcf

# Haplotype finding for SNPs
run_pipeline.pl -Xmx60g -FILLINFindHaplotypesPlugin -hmp $OUT_DIR/tassel.hmp.txt -o $OUT_DIR/snp_matrix_haplotypes

# Imputation using FILLIN for SNPs
run_pipeline.pl -Xmx60g -FILLINImputationPlugin -hmp $OUT_DIR/tassel.hmp.txt -d $OUT_DIR/snp_matrix_haplotypes -o $OUT_DIR/snp_matrix_imputed.hmp.txt.gz

# Filter SNPs with sliding window
FILTERED_SNP_OUT="$OUT_DIR/filtered_tassel"
run_pipeline.pl -Xmx60g -fork1 -importGuess $OUT_DIR/snp_matrix_imputed.hmp.txt.gz \
    -filterAlign -filterAlignMinFreq 0.05 -filterAlignMaxFreq 0.95 \
    -filterAlignMinCount 10 \
    -export $FILTERED_SNP_OUT

run_pipeline.pl -Xmx60g -importGuess $OUT_DIR/snp_matrix_imputed.hmp.txt.gz \
-genotypeSummary all -export $OUT_DIR/snp_genotypeSummary

# Ensure the steps below use the filtered SNPs output
# Numeric Allele Counts for filtered SNPs
run_pipeline.pl -Xmx60g -fork1 -importGuess $OUT_DIR/snp_matrix_imputed.hmp.txt.gz \
    -NumericalGenotypePlugin -endPlugin -export $OUT_DIR/numeric_allele_counts_filtered \
    -exportType ReferenceProbability

# Kinship Matrix Calculation for filtered SNPs
run_pipeline.pl -Xmx60g -fork1 -importGuess $OUT_DIR/snp_matrix_imputed.hmp.txt.gz \
    -KinshipPlugin -method Centered_IBS -endPlugin -export $OUT_DIR/kinship_matrix_filtered

# PCA for filtered SNPs
run_pipeline.pl -Xmx60g -fork1 -importGuess $OUT_DIR/snp_matrix_imputed.hmp.txt.gz \
    -PrincipalComponentsPlugin -covariance true -endPlugin -export $OUT_DIR/pca_filtered

echo "TASSEL SNP analysis completed."

# Convert filtered indel VCF to TASSEL format
gunzip -c $INDEL_OUT_DIR/filtered_indels.vcf.gz > $INDEL_OUT_DIR/temp_indels.vcf
run_pipeline.pl -Xmx60g -importGuess $INDEL_OUT_DIR/temp_indels.vcf -export $INDEL_OUT_DIR/indel_tassel
rm $INDEL_OUT_DIR/temp_indels.vcf

# Haplotype finding for indels
run_pipeline.pl -Xmx60g -FILLINFindHaplotypesPlugin -hmp $INDEL_OUT_DIR/indel_tassel.hmp.txt -o $INDEL_OUT_DIR/indel_snp_matrix_haplotypes

# Imputation using FILLIN for indels
run_pipeline.pl -Xmx60g -FILLINImputationPlugin -hmp $INDEL_OUT_DIR/indel_tassel.hmp.txt -d $INDEL_OUT_DIR/indel_snp_matrix_haplotypes -o $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz

# Filter indels
FILTERED_INDEL_OUT="$INDEL_OUT_DIR/filtered_indel_tassel"
run_pipeline.pl -Xmx60g -fork1 -importGuess $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz \
    -filterAlign -filterAlignMinFreq 0.05 -filterAlignMaxFreq 0.95 \
    -filterAlignMinCount 10 -export $FILTERED_INDEL_OUT

run_pipeline.pl -Xmx60g -importGuess $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz -genotypeSummary all -export $INDEL_OUT_DIR/indel_genotypeSummary

# Ensure the steps below use the filtered indels output
# Numeric Allele Counts for filtered indels
run_pipeline.pl -Xmx60g -fork1 -importGuess $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz \
    -NumericalGenotypePlugin -endPlugin -export $INDEL_OUT_DIR/numeric_allele_counts_filtered_indel \
    -exportType ReferenceProbability

# Kinship Matrix Calculation for filtered indels
run_pipeline.pl -Xmx60g -fork1 -importGuess $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz \
    -KinshipPlugin -method Centered_IBS -endPlugin -export $INDEL_OUT_DIR/kinship_matrix_filtered_indel

# PCA for filtered indels
run_pipeline.pl -Xmx60g -fork1 -importGuess $INDEL_OUT_DIR/indel_snp_matrix_imputed.hmp.txt.gz \
    -PrincipalComponentsPlugin -covariance true -endPlugin -export $INDEL_OUT_DIR/pca_filtered_indel

echo "Indel TASSEL analysis completed."

# Zip all the files created by TASSEL
zip -r 6_tassel_analysis.zip 6_tassel_analysis
zip -r 6_tassel_analysis_indels.zip 6_tassel_analysis_indels

echo "TASSEL analysis files zipped."

# To download the files, use the following command, replacing the username and server address with your own
#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/6_tassel_analysis.zip ~/Desktop
#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/6_tassel_analysis_indels.zip ~/Desktop
#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/5_annotation/snpEff_snps/geno_filtered_snps.ann.vcf.gz ~/Desktop
#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/5_annotation/snpEff_indels/geno_filtered_indels.ann.vcf.gz ~/Desktop

geno_filtered_snps.ann.vcf.gz
geno_filtered_indels.ann.vcf.gz
