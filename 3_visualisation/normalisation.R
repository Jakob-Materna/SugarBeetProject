# normalisation.R
# 
# Description: Script used for normalisation and visualisation of correlation 
# between replicates
# 
# Input: count.matrix.tsv
# Output: PCA plot,correlatio n between replicates and a file of averaged gene counts  
# 
# Version: 1.00
# Date: 28.05.2024
# Author: Jakob Materna

library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(tibble)
library(reshape2)
library(stringr)

# Get characteristic (strain, treatment, replicate) from sample name
get_field <- function(name, field) {
  characteristic <- strsplit(name, "_")[[1]][field]
  return(characteristic)
}


count_matrix <- as.matrix(read.csv("input/count.matrix.tsv", sep="\t", row.names="gene_id"))
count_matrix <- count_matrix[, colnames(count_matrix) != "gene_name"]
# Removing outlier
count_matrix <- count_matrix[, colnames(count_matrix) != "X69_2_B"]

# Renaming and reordering the matrix
colnames(count_matrix) <- gsub("X", "", colnames(count_matrix))
count_matrix <- count_matrix[, str_sort(colnames(count_matrix), numeric = TRUE)]

# Turning values into integers
mode(count_matrix) <- "integer"

# Creating table of conditions
sample_names <- as.vector(colnames(count_matrix))
sample_conditions <- sapply(sample_names, get_field, field=2)
sample_table <- data.frame(treatment=as.factor(sample_conditions))
row.names(sample_table) <- sample_names

# Creating dds object
dds <-  DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = sample_table,
                               design = ~ treatment)

# Calculating size factors
dds <- estimateSizeFactors(dds)

# Median of Ratios normalization
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))

# Averaging mrn values of replicates
i <- 1
avg <- data.frame(gene_id=as.factor(rownames(normalized_counts)))
repeat {
  if (is.na(sample_names[i])) {
    break
  }
  strain <- get_field(sample_names[i], field=1)
  treatment <- get_field(sample_names[i], field=2)
  name <- paste0(strain, "_", treatment)
  if (treatment == get_field(sample_names[i+2], field=2) && strain == get_field(sample_names[i+2], field=1)) {
    avg[name] <- rowMeans(normalized_counts[,i:(i+2)])
    # cat("\t", name, "\t", i:(i+2), "\n")
    i <- i+3
  }
  else {
    avg[name] <- rowMeans(normalized_counts[,i:(i+1)])
    # cat("2:\t", name, "\t", i:(i+1), "\n")
    i <- i+2
  }
}

# Write Tabe of averaged gene counts in a table
write.table(avg, file="avg.matrix.tsv", row.names=FALSE, quote=FALSE, sep="\t")

# Regularized log transformation
vst <- vst(dds, blind=TRUE)

print(plotPCA(vst, intgroup=c("treatment"), ntop=500))

# visualizing pairwise correlation
count_df <- tibble::rownames_to_column(normalized_counts, "gene_id")
count_df <- melt(count_df, id.vars = "gene_id", variable.name = "mrn")
par(mfrow=c(5, 3), pty="s")
i <- 1
repeat {
  if (is.na(sample_names[i])) {
    break
  }
  treatment <- get_field(sample_names[i], field=2)
  if (treatment == get_field(sample_names[i+2], field=2)) {
    main <- paste0("R=", round(cor(1+normalized_counts[,i+1], 1+normalized_counts[,i], method = "pearson"), digits=4))
    plot(log2(1+normalized_counts[,c(i+1, i)]), col="black", pch=20, cex=0.3, main=main)
    main <- paste0("R=", round(cor(1+normalized_counts[,i+2], 1+normalized_counts[,i], method = "pearson"), digits=4))
    plot(log2(1+normalized_counts[,c(i+2, i)]), col="black", pch=20, cex=0.3, main=main)
    main <- paste0("R=", round(cor(1+normalized_counts[,i+2], 1+normalized_counts[,i+1], method = "pearson"), digits=4))
    plot(log2(1+normalized_counts[,c(i+2, i+1)]), col="black", pch=20, cex=0.3, main=main)
    i <- i+3
  }
  else {
    main <- paste0("R=", round(cor(1+normalized_counts[,i+1], 1+normalized_counts[,i], method = "pearson"), digits=4))
    plot(log2(1+normalized_counts[,c(i+1, i)]), col="black", pch=20, cex=0.3, main=main)
    plot(NULL, xlim=c(0, 15), ylim=c(0, 15), ylab="", xlab="", main="")
    plot(NULL, xlim=c(0, 15), ylim=c(0, 15), ylab="", xlab="", main="")
    i <- i+2
  }
}
i <- 1
repeat {
  if (is.na(sample_names[i])) {
    break
  }
  treatment <- get_field(sample_names[i], field=2)
  if (treatment == get_field(sample_names[i+2], field=2)) {
    main <- paste0("R=", round(cor(assay(vst)[,i+1], assay(vst)[,i], method = "pearson"), digits=4))
    plot(assay(vst)[,c(i+1, i)], col="black", pch=20, cex=0.3, main=main)
    main <- paste0("R=", round(cor(assay(vst)[,i+2], assay(vst)[,i], method = "pearson"), digits=4))
    plot(assay(vst)[,c(i+2, i)], col="black", pch=20, cex=0.3, main=main)
    main <- paste0("R=", round(cor(assay(vst)[,i+2], assay(vst)[,i+1], method = "pearson"), digits=4))
    plot(assay(vst)[,c(i+2, i+1)], col="black", pch=20, cex=0.3, main=main)
    i <- i+3
  }
  else {
    main <- paste0("R=", round(cor(assay(vst)[,i+1], assay(vst)[,i], method = "pearson"), digits=4))
    plot(assay(vst)[,c(i+1, i)], col="black", pch=20, cex=0.3, main=main)
    plot(NULL, xlim=c(0, 15), ylim=c(0, 15), ylab="", xlab="", main="")
    plot(NULL, xlim=c(0, 15), ylim=c(0, 15), ylab="", xlab="", main="")
    i <- i+2
  }
}

