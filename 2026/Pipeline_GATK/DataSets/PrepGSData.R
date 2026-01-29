# Set working directory
setwd("~/Library/CloudStorage/GoogleDrive-deniz.akdemir.work@gmail.com/.shortcut-targets-by-id/1n5uwNs72GdTeZxmWZKLfFkgE0FIjPxVG/Akdemir_Deniz/github/GATKPIPE2/DataSets")

# Load necessary libraries
library(readxl)
library(tidyverse)
library(data.table)
library(rrBLUP)

# Function to process the data
process_data <- function(file_path, cm_value) {
  data <- read_excel(file_path) %>%
    mutate(
      Plant = sprintf("PyrSep.%03s", as.character(Plant)),
      leafArea = as.numeric(leafArea),
      necrosisArea = as.numeric(necrosisArea),
      PLACL = as.numeric(PLACL),
      pycnidiaCount = as.numeric(pycnidiaCount),
      meanPycnidiaArea = as.numeric(meanPycnidiaArea),
      pycnidiaPerCm2Leaf = as.numeric(pycnidiaPerCm2Leaf),
      pycnidiaPerCm2Lesion = as.numeric(pycnidiaPerCm2Lesion),
      pycnidiaGreyValue = as.numeric(pycnidiaGreyValue),
      Picture = as.character(Picture),
      Leaf = as.character(Leaf),
      cm = as.factor(cm_value)
    )
  return(data)
}

# Process datasets with the correct cm values
septoria_12cm <- process_data("RawData/GSData/Output_reps_unidas_12cm_mock_checks.xlsx", "12")
septoria_17cm <- process_data("RawData/GSData/Output_reps_unidas_17cm_mock_checks.xlsx", "17")
septoria_12cm$cm<-"12"
septoria_17cm$cm<-"17"

# Load wheat variety SNP matrix and rename it to wheat_genotype_matrix
load("RawData/WheatGeno/Genotyping/SNPmatrixFilt_inbreed.RData")
wheat_genotype_matrix <- snp_matrix
rm(snp_matrix)

# Intersect Plants that are present in both datasets and SNP matrix
common_plants <- intersect(union(septoria_12cm$Plant, septoria_17cm$Plant), rownames(wheat_genotype_matrix))

# Filter datasets to include only Plants in SNP matrix
septoria_12cm <- septoria_12cm %>% filter(Plant %in% common_plants)
septoria_17cm <- septoria_17cm %>% filter(Plant %in% common_plants)

# Combine both datasets
septoria_combined_pheno <- bind_rows(septoria_12cm, septoria_17cm)

# Convert all character variables to factors
septoria_combined_pheno <- septoria_combined_pheno %>%
  mutate(across(where(is.character), as.factor))

# Load septoria genotype data
septoria_genotype_data <- fread("RawData/6_tassel_analysis/numeric_allele_counts_filtered.txt",
                                skip = 1, header = TRUE, sep = "\t", na.strings = c("", "NA", 0.5, "0.5"))
septoria_genotype_data <- as.data.frame(septoria_genotype_data)

# Rename the first column to "Taxa" and remove the last 13 characters
colnames(septoria_genotype_data)[1] <- "Taxa"
septoria_genotype_data$Taxa <- substr(septoria_genotype_data$Taxa, 1, nchar(septoria_genotype_data$Taxa) - 13)

namesGeno<-septoria_genotype_data$Taxa
#see which column names are duplicated
duplicatednames<-duplicated(namesGeno)

# now for each duplicated column name, find the index of the first occurence of the column name and the sec
for (i in 1:length(duplicatednames)) {
  print(i)
  if(duplicatednames[i]) {
    firstOccurence<-which(namesGeno==namesGeno[i])[1]
    secondOccurence<-which(namesGeno==namesGeno[i])[2]
    # now take the mean of the two columns
    septoria_genotype_data[firstOccurence,2:ncol(septoria_genotype_data)]<-colMeans(septoria_genotype_data[c(firstOccurence, secondOccurence),2:ncol(septoria_genotype_data)], na.rm=TRUE)
    # remove the second column
  }
}

sum(duplicatednames)
dublicatedGenotypeNames<-septoria_genotype_data$Taxa[duplicatednames]


septoria_genotype_data<-septoria_genotype_data[!duplicatednames,]
# Convert septoria genotype data to a matrix
septoria_genotype_matrix <- as.matrix(septoria_genotype_data[-1]) * 2 - 1
rownames(septoria_genotype_matrix) <- septoria_genotype_data$Taxa

# Remove rows with "sorted" in rownames
septoria_genotype_matrix <- septoria_genotype_matrix[!grepl("sorted", rownames(septoria_genotype_matrix)), ]

# Perform genotype imputation for septoria and wheat genotypes
imputed_septoria_genotypes <- rrBLUP::A.mat(septoria_genotype_matrix, impute.method = "mean", return.imputed = TRUE, min.MAF = 0.05, max.missing = 0.3)$imputed
imputed_wheat_genotypes <- rrBLUP::A.mat(wheat_genotype_matrix, impute.method = "mean", return.imputed = TRUE, min.MAF = 0.05, max.missing = 0.3)$imputed

# Clean up
rm(septoria_genotype_matrix, wheat_genotype_matrix)
gc()

# Load isolate lists
isolateslist1 <- read_xlsx("RawData/GenotypeTables/isolated list 2021-2022 for sequencing.xlsx")
isolateslist2 <- read_xlsx("RawData/GenotypeTables/isolated list 2023 for sequencing.xlsx")

# Perform SVD and calculate principal components (PCs) for septoria genotypes
svd_septoria <- svd(scale(imputed_septoria_genotypes, scale = FALSE, center = TRUE), nu = 5, nv = 5)
pc_septoria <- scale(imputed_septoria_genotypes, scale = FALSE, center = TRUE) %*% svd_septoria$v

# Convert PCs to a data frame and add Taxa information
pc_septoria <- data.frame(Taxa = rownames(imputed_septoria_genotypes), pc_septoria)
colnames(pc_septoria)[2:6] <- paste("PC", 1:5, sep = "")

# Assign Year based on the pattern in Taxa
pc_septoria$Year <- ifelse(grepl("EKDN2300427", pc_septoria$Taxa), "2023", "2022")

# Split data by Year and process Taxa
pcs2022 <- pc_septoria[pc_septoria$Year == "2022", ]
pcs2023 <- pc_septoria[pc_septoria$Year == "2023", ]

pcs2022$TaxaShort <- gsub("combined_", "", substr(pcs2022$Taxa, 1, 4))
pcs2023$TaxaShort <- gsub("combined_", "", substr(pcs2023$Taxa, 1, 4))

# there are some "_*" at the end of TaxaShort, remove them, * is any character, sometimes the name has no "_" at the end
pcs2022$TaxaShort <- gsub("_.*$", "", pcs2022$TaxaShort)
pcs2023$TaxaShort <- gsub("_.*$", "", pcs2023$TaxaShort)





pcs2022 <- pcs2022[order(as.numeric(gsub("S", "", pcs2022$TaxaShort))), ]
pcs2023 <- pcs2023[order(as.numeric(gsub("S", "", pcs2023$TaxaShort))), ]

pcs2022$Isolate <- isolateslist1$Isolate[match(pcs2022$TaxaShort, isolateslist1$Name)]
pcs2023$Isolate <- isolateslist2$Isolate[match(pcs2023$TaxaShort, isolateslist2$Name)]

# Clean up isolate names
pcs2022$Isolate <- gsub(",", ".", pcs2022$Isolate)
pcs2022$Isolate <- recode(pcs2022$Isolate,
                          "22_Conil_Fer_L1" = "22_ConilFer_L1",
                          "22_Jerez Val_L1" = "22_JerezVal_L1",
                          "22_EcijaSecOrt _L1" = "22_EcijaSecOrt_L1",
                          "22_EcijaSecSim_L1" = "22_EcijaSecSim_L2",
                          "22_EscIca_L2" = "22_EscIca_L2",
                          "22_EcijaSecSha_L1" = "22_EcijaSecSah_L1")

pcs2023$Isolate <- gsub(",", ".", pcs2023$Isolate)

# Combine pcs2022 and pcs2023
pc_septoria_combined <- rbind(pcs2022, pcs2023)

# Replace the rownames of imputed_septoria_genotypes with the Isolate column of pc_septoria_combined
rownames_septoria_genotypes <- pc_septoria_combined$Isolate[match(pc_septoria_combined$Taxa, rownames(imputed_septoria_genotypes))]
rownames(imputed_septoria_genotypes) <- rownames_septoria_genotypes

# Remove any duplicate row names in imputed_septoria_genotypes
imputed_septoria_genotypes <- imputed_septoria_genotypes[!duplicated(rownames(imputed_septoria_genotypes)), ]

# Define the component mixes
mix1 <- c('22_EcijaSec83Ica_L2', '22_EcijaSecCris_L1', '22_EcijaRegTej_L1')
mix2 <- c('22_CorKiko_L1', '22_CorCale_L1', '22_Cor3927_L1')
mix3 <- c('22_ConilAmi_L1', '22_Conil3806_L1', '22_Jerez3927_L1')
mix4 <- c('22_CorVal_L1')

mixgeno <- rbind(
  apply(imputed_septoria_genotypes[mix1, ], 2, function(x) mean(na.omit(x))),
  apply(imputed_septoria_genotypes[mix2, ], 2, function(x) mean(na.omit(x))),
  apply(imputed_septoria_genotypes[mix3, ], 2, function(x) mean(na.omit(x))),
  imputed_septoria_genotypes[mix4, ]
)

# Set rownames for the mixes
rownames(mixgeno) <- c("mix1", "mix2", "mix3", "mix4")

# Add the mix genotypes to the imputed_septoria_genotypes matrix
geno_septoria_with_mixes <- rbind(imputed_septoria_genotypes, mixgeno)

# Remove any duplicates again just in case
geno_septoria_with_mixes <- geno_septoria_with_mixes[!duplicated(rownames(geno_septoria_with_mixes)), ]

# Check dimensions of the final matrix
dim(geno_septoria_with_mixes)

# Save all the processed data into one RData file
save(septoria_combined_pheno, geno_septoria_with_mixes, imputed_wheat_genotypes, pc_septoria_combined, file = "GPdata.RData")

# Clean up
rm(list = ls())
gc()
