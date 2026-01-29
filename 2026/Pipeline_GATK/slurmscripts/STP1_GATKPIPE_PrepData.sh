# This code unzips the fastq files and creates a list of the fastq files to be used in the next steps of the pipeline

# This is a bash script to be run in the terminal
# to run this script type the following command in the terminal, from the directory containing the script:
# bash slurmscripts/STP1_GATKPIPE_PrepData.sh



# 1- unzip the fastq file :X204SC22100973-Z01-F001.tar under the directory 1_data/fastq_set1

mkdir -p 1_data/fastq_set1
tar -xf 1_data/X204SC22100973-Z01-F001.tar -C 1_data/fastq_set1

# 2- unzip the fastq file :X204SC22100973-Z01-F002.tar under the directory 1_data/fastq_set2

mkdir -p 1_data/fastq_set2
tar -xf 1_data/X204SC22100973-Z01-F002.tar -C 1_data/fastq_set2

# 3- Unzip the reference genome Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz under the directory 0_index/referenceIPO323

mkdir -p 0_index/referenceIPO323
gunzip -c 1_data/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz > 0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa

# 4- Unzip the reference genome Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz under the directory 0_index/reference3D7

mkdir -p 0_index/reference3D7
gunzip -c 1_data/reference3D7/Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz > 0_index/reference3D7/Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa
