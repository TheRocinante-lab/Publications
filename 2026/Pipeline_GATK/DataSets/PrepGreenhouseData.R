# Set working directory
setwd("~/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/DataSets")

# Load necessary libraries
library(readxl)
library(tidyverse)
library(data.table)
library(rrBLUP)

load("GPdata.RData")
rm(imputed_wheat_genotypes)
rm(septoria_combined_pheno)
gc()

ls()

# Load phenotypes and ensure matching with genotypes
phenotypes <- read.csv("RawData/GreenHouse/raw_phenotypes.csv")
phenotypes[1:5,]

library(magrittr)
library(tidyverse)
# Extract information from the Picture column
phenotypes <- phenotypes %>%
  mutate(
    WheatGeno = str_extract(Picture, "^[^_]+"),
    Year = str_extract(Picture, "(?<=_)\\d+")
  )

table(phenotypes$WheatGeno)
table(phenotypes$Year)

pcs2022<-pc_septoria_combined[pc_septoria_combined$Year=="2022",]
pcs2023<-pc_septoria_combined[pc_septoria_combined$Year=="2023",]

colnames(phenotypes)
phenotypes[1:5,]
phenotypes$Isolate<-NA
## -------------------------------------------------------------------------------------------
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypes)){
  for (j in 1:nrow(pcs2022)){
    if (grepl(pcs2022$Isolate[j], phenotypes$Picture[i])){
      phenotypes$Isolate[i]<-pcs2022$Isolate[j]
    }
  }
}
table(phenotypes$Isolate)



SNPSummary<-read.table("RawData/6_tassel_analysis/snp_genotypeSummary3.txt", header = TRUE, sep = "\t", as.is = TRUE, skip=0)
SNPSummary[1:5,]


# Now a map file
SNPMAP<- data.frame(
  rs = SNPSummary$Site.Name,
  chrom = SNPSummary$Chromosome,
  pos = SNPSummary$Physical.Position
)

SNPMAP<-SNPMAP[SNPMAP$rs%in%colnames(geno_septoria_with_mixes),]
dim(SNPMAP)
SNPMAP<-SNPMAP[SNPMAP$chrom<=13,]
dim(geno_septoria_with_mixes)
geno_septoria_with_mixes[1:5,1:5]
geno_septoria_with_mixes<-geno_septoria_with_mixes[, SNPMAP$rs]
dim(geno_septoria_with_mixes)


save(phenotypes, geno_septoria_with_mixes, SNPMAP,pcs2022,pcs2023, file="Pheno_geno_Map_GreenhouseTraits.RData")


dim(geno_septoria_with_mixes)

