# GATKPIPE2

A bioinformatics pipeline for analyzing *Zymoseptoria tritici* sequencing data, closely following the GATK variant calling best practices (Van der Auwera & O'Connor, 2020; Poplin et al., 2017).

## Pipeline Overview

1. **Quality Control and Trimming:**
   - **Decompression and Initial QC:** Raw FASTQ files were decompressed and subjected to quality control using FastQC, with results compiled via MultiQC.
   - **Trimming:** Using Trimmomatic, we trimmed adapters and low-quality sequences, followed by additional quality control.

2. **Read Alignment:**
   - **Alignment:** Reads were aligned to the *Z. tritici* reference genome (IPO323) from EnsemblGenomes.org using BWA-MEM.
   - **SAM to BAM Conversion:** SAM files were converted to sorted BAM files with SAMtools.
   - **Duplicate Marking:** Duplicate reads were marked, and metrics were collected using Picard.

3. **Variant Calling:**
   - **Raw VCF Generation:** GATK's HaplotypeCaller was employed to generate raw VCF files, applying stringent filters to identify high-confidence SNPs.
   - **Base Quality Score Recalibration (BQSR):** Systematic errors were corrected, and recalibrated BAM files were used to produce GVCF files.
   - **Joint Genotyping:** GVCF files were merged into a GenomicsDB using GATK's GenomicsDBImport tool, and joint genotyping across all samples was performed to create a consolidated VCF.
   - **SNP Annotation:** SnpEff was used for SNP annotation.

4. **Downstream Analysis:**
   - **Conversion and Haplotype Analysis:** TASSEL (Bradbury et al., 2007) was utilized to convert annotated VCFs into hapmap format, perform haplotype finding and imputation, and calculate kinship matrices and PCA for population structure assessment.
   - **Consensus Sequence Generation:** Consensus sequences for specific genomic regions were generated, providing insights into the genetic variability of *Z. tritici* strains.

This pipeline integrates robust tools and methodologies, adhering to bioinformatics best practices for high-throughput sequencing data analysis (McKenna et al., 2010; DePristo et al., 2011).

## Directory Structure

Folders other than `slurmscripts` are empty and are present to establish the directory structure. All the code is located in the `slurmscripts` folder.

### Data Placement

Please place the relevant data files in the following directories:

1. **Raw FASTQ Files:**
   - `1_data/fastq_set1`: Place the `X204SC22100973-Z01-F001.tar` file here.
   - `1_data/fastq_set2`: Place the `X204SC22100973-Z01-F002.tar` file here.

2. **Reference Genomes:**
   - `1_data/referenceIPO323`: Place the `Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz` file here.
   - `1_data/reference3D7`: Place the `Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz` file here.

You can use the `get_reference_stp3.sh` script and run it in the terminal to download the reference genome files.

### Unzipping Instructions

To unzip the files and prepare them for the pipeline, run the following script in the terminal:

```bash
# STP1_GATKPIPE_PrepData.sh is a bash script to be run in the terminal
# To run this script, type the following command in the terminal, from the directory containing the script:
# bash slurmscripts/STP1_GATKPIPE_PrepData.sh
```

The remaining steps are to be run on an HPC cluster using **sbatch**.