#!/bin/bash

# Create the reference directory if it doesn't exist
mkdir -p referenceIPO323

# Download the reference genome files to the reference directory
# Using wget to download all files from the FTP directory
wget -r -np -nH --cut-dirs=7 -P referenceIPO323 ftp://ftp.ensemblgenomes.org/pub/fungi/release-57/fasta/zymoseptoria_tritici/dna/

# Verify that the files have been copied
echo "Checking files in /referenceIPO323:"
ls -l referenceIPO323/


# Create the reference directory if it doesn't exist
mkdir -p reference3D7

# Download the reference genome files to the reference directory
# Using wget to download all files from the FTP directory
wget -r -np -nH --cut-dirs=7 -P reference3D7 ftp://ftp.ensemblgenomes.org/pub/fungi/release-57/fasta/fungi_ascomycota4_collection/zymoseptoria_tritici_st99ch_3d7_gca_900091695/dna/

# Verify that the files have been copied
echo "Checking files in /reference3D7:"
ls -l reference3D7/
