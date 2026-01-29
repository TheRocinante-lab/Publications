#!/bin/bash
#SBATCH --job-name=VariantAnnotation
#SBATCH --output=VariantAnnotation_%j.log
#SBATCH --error=VariantAnnotation_%j.err
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user="deniz.akdemir.work@gmail.com"

# Load required modules
module load Java/1.8.0_212
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8 
module load picard/2.26.10-Java-15

# Ensure SnpEff is downloaded and unzipped
if [ ! -d "snpEff" ]; then
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip -d snpEff
fi

# Define paths
PICARD="$EBROOTPICARD/picard.jar"
SNPEFF="snpEff/snpEff/snpEff.jar"
SNPEFF_DATA_DIR="snpEff/snpEff/data"

# Update the snpEff.config to include Zymoseptoria tritici
SNPEFF_CONFIG="snpEff/snpEff/snpEff.config"
if ! grep -q "zymoseptoria_tritici.genome" "$SNPEFF_CONFIG"; then
    echo "# Zymoseptoria tritici" >> "$SNPEFF_CONFIG"
    echo "zymoseptoria_tritici.genome : Zymoseptoria_tritici.MG2" >> "$SNPEFF_CONFIG"
fi

# Define directories for VCF files and annotations
GVCF_DIR="4_processing/GVCF"
SNPEFF_DIR="5_annotation/snpEff_snps"
mkdir -p $SNPEFF_DIR

# Download genome and annotation files for Zymoseptoria tritici
GENOME_DIR="$SNPEFF_DATA_DIR/zymoseptoria_tritici"
mkdir -p $GENOME_DIR

wget -P $GENOME_DIR https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-58/fasta/zymoseptoria_tritici/dna/Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz
wget -P $GENOME_DIR https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-58/gff3/zymoseptoria_tritici/Zymoseptoria_tritici.MG2.58.gff3.gz

gunzip $GENOME_DIR/Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz
gunzip $GENOME_DIR/Zymoseptoria_tritici.MG2.58.gff3.gz

# Rename the genome sequence file for SnpEff compatibility
mv $GENOME_DIR/Zymoseptoria_tritici.MG2.dna.toplevel.fa $GENOME_DIR/sequences.fa

# Rename the GFF3 file to 'genes.gff' for SnpEff compatibility
mv $GENOME_DIR/Zymoseptoria_tritici.MG2.58.gff3 $GENOME_DIR/genes.gff

# Build the SnpEff database
java -jar $SNPEFF build -gff3 -v zymoseptoria_tritici

# Annotation of filtered SNPs variants
java -jar $PICARD GatherVcfs -I $GVCF_DIR/filtered_snps.vcf -O $SNPEFF_DIR/geno_filtered_snps.vcf 
java -Xmx32g -jar $SNPEFF Zymoseptoria_tritici $SNPEFF_DIR/geno_filtered_snps.vcf > $SNPEFF_DIR/geno_filtered_snps.ann.vcf
gzip $SNPEFF_DIR/geno_filtered_snps.ann.vcf
echo "SNP annotation completed."

# Similarly for filtered indels
SNPEFF_INDELS_DIR="5_annotation/snpEff_indels"
mkdir -p $SNPEFF_INDELS_DIR

java -jar $PICARD GatherVcfs -I $GVCF_DIR/filtered_indels.vcf -O $SNPEFF_INDELS_DIR/geno_filtered_indels.vcf 

java -Xmx32g -jar $SNPEFF Zymoseptoria_tritici $SNPEFF_INDELS_DIR/geno_filtered_indels.vcf > $SNPEFF_INDELS_DIR/geno_filtered_indels.ann.vcf
gzip $SNPEFF_INDELS_DIR/geno_filtered_indels.ann.vcf
echo "Indel annotation completed."


# On the local machine, download the annotated VCF files using scp:
#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/5_annotation/snpEff_snps/geno_filtered_snps.ann.vcf.gz  ~/Desktop


#scp -P 33322 deniz@loginhpckairos.cbgp.upm.es:GenomicData/github/GATKPIPE2/5_annotation/snpEff_indels/geno_filtered_indels.ann.vcf.gz  ~/Desktop

