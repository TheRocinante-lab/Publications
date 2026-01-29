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

rownames(geno_septoria_with_mixes)





# Load phenotypes and ensure matching with genotypes
phenotypesR3E3<-read.csv("RawData/PilarEndophyte/PhenoData/20240423_R3E3_R1&R2_GWAS_astringent&loose.csv", header=TRUE)
phenotypesR3E3$Zt_Strain <- gsub("-", ".", phenotypesR3E3$Zt_Strain)

print(table(phenotypesR3E3$Removal))
phenotypesR3E3<-phenotypesR3E3[phenotypesR3E3$Removal=="Astringent",]

# Load phenotypes and ensure matching with genotypes
phenotypesR3E5<-read.csv("RawData/PilarEndophyte/PhenoData/20240820_R3E5_R1&R2_GWAS_astringent_Deniz.csv", header=TRUE)
phenotypesR3E5$Zt_Strain <- gsub("-", ".", phenotypesR3E5$Zt_Strain)
phenotypesR3E5$Removal<-NA

pcs2022<-pc_septoria_combined[pc_septoria_combined$Year=="2022",]
pcs2023<-pc_septoria_combined[pc_septoria_combined$Year=="2023",]

## -------------------------------------------------------------------------------------------
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypesR3E3)){
  for (j in 1:nrow(pcs2022)){
    if (grepl(pcs2022$Isolate[j], phenotypesR3E3$Zt_Strain[i])){
      phenotypesR3E3$Isolate[i]<-pcs2022$Isolate[j]
    }
  }
}
table(phenotypesR3E3$Isolate)


## -------------------------------------------------------------------------------------------
# put that name in the phenotype file in a new column
for (i in 1:nrow(phenotypesR3E5)){
  for (j in 1:nrow(pcs2022)){
    if (grepl(pcs2022$Isolate[j], phenotypesR3E5$Zt_Strain[i])){
      phenotypesR3E5$Isolate[i]<-pcs2022$Isolate[j]
    }
  }
}
table(phenotypesR3E5$Isolate)

head(phenotypesR3E5)


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
print(colnames(phenotypesR3E3))
print(colnames(phenotypesR3E5))
colnamestokeep<-intersect(colnames(phenotypesR3E5), colnames(phenotypesR3E3))
phenotypes<-rbind(phenotypesR3E3[,colnamestokeep], phenotypesR3E5[,colnamestokeep])


save(phenotypes, geno_septoria_with_mixes, SNPMAP,pcs2022,pcs2023, file="Pheno_geno_Map_PillarTraits.RData")




