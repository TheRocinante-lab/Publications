# Load necessary libraries
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(ape)
library(Biostrings)

# Read the FASTA file
fasta_file <- "/Users/denizakdemir/Desktop/consensus_3_247434_251178.fasta"
sequences <- readDNAStringSet(fasta_file, format = "fasta")
length(sequences)
sequences[2]
# Convert DNAStringSet to character vector
sequences <- as.character(sequences)
# take only odd sequences
sequences<-sequences[seq(1,length(sequences),2)]
# Pad sequences to ensure equal length by aligning them manually
max_len <- max(nchar(sequences))  # Find the maximum sequence length
padded_sequences <- sapply(sequences, function(seq) {
  pad_length <- max_len - nchar(seq)
  if (pad_length > 0) {
    return(paste0(seq, paste(rep("-", pad_length), collapse = "")))  # Add gaps to shorter sequences
  } else {
    return(seq)
  }
}, USE.NAMES = FALSE)

# Convert padded sequences to a matrix
alignment <- do.call(rbind, strsplit(padded_sequences, ""))

# Function to identify polymorphic loci
find_polymorphic_loci <- function(alignment) {
  # Identify columns (positions) where more than one unique nucleotide is present
  polymorphic_loci <- which(apply(alignment, 2, function(col) {
    length(unique(col)) > 1
  }))
  return(polymorphic_loci)
}

# Find polymorphic loci
polymorphic_loci <- find_polymorphic_loci(alignment)

# Print results
cat("Polymorphic loci positions:\n")
print(polymorphic_loci)

# calculate PCA after converting the alignment to a numeric matrix (only on polymorphic loci)
pca_data <- as.data.frame(t(as.matrix(alignment[, polymorphic_loci])))
dim(pca_data)
namestouse<-rownames(as.matrix(alignment[, polymorphic_loci]))
# only take odd rows
namestouse<-namestouse[seq(1,length(namestouse),2)]
colnames(pca_data) <-namestouse

# Perform PCA
# convert to numeric, from A, C, G, T to 1, 2, 3, 4
pca_data_num <- apply(pca_data, 2, function(col) {
  as.numeric(factor(col, levels = c("A", "C", "G", "T")))
})

# do mean imputation
pca_data_num[is.na(pca_data_num)] <- colMeans(pca_data_num, na.rm = TRUE)

# remove columns with NA values
pca_data_num <- pca_data_num[, colSums(is.na(pca_data_num)) == 0]

# Plot PCA
pca_result <- prcomp(pca_data_num, scale. = TRUE)
plot(pca_result$x[, 1], pca_result$x[, 2], xlab = "PC1", ylab = "PC2", main = "PCA of Polymorphic Loci")

round((pca_result$sdev/sum(pca_result$sdev))[1:10], 3)
