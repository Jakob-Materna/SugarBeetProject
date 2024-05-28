# visualisation.R
# 
# Description: Script used for visualisation of distances between samples 
# between replicates
# 
# Input: avg.matrix.tsv
# Output: Heatmap of distances between samples
# 
# Version: 1.00
# Date: 28.05.2024
# Author: Jakob Materna

library(DESeq2)
library(gplots)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(tibble)
library(reshape2)
library(stringr)
library(paletteer)

# Get characteristic (strain, treatment, replicate) from sample name
get_field <- function(name, field) {
  characteristic <- strsplit(name, "_")[[1]][field]
  return(characteristic)
}

get_biostimulation <- function(name) {
  treat <- strsplit(name, "_")[[1]][1]
  if (treat %in% c("69", "77", "F")) {
    return("positive")
  }
  if (treat %in% c("54")) {
    return("negative")
  }
  if (treat %in% c("A", "66", "70")) {
    return("neutral")
  }
}

get_biocontrol <- function(name) {
  treat <- strsplit(name, "_")[[1]][1]
  if (treat %in% c("69", "A", "54")) {
    return("positive")
  }
  if (treat %in% c("F", "66", "70")) {
    return("negative")
  }
  if (treat %in% c("77")) {
    return("neutral")
  }
}

get_susceptibility <- function(name) {
  treat <- strsplit(name, "_")[[1]][1]
  if (treat %in% c("69", "A", "54", "F")) {
    return("susceptible")
  }
  if (treat %in% c("77", "66", "70")) {
    return("resistant")
  }
}

# import count matrix
count_matrix <- as.matrix(read.csv("avg.matrix.tsv", sep="\t", row.names="gene_id"))
count_matrix <- count_matrix[, colnames(count_matrix) != "gene_id"]

# renaming and reordering the matrix
colnames(count_matrix) <- gsub("X", "", colnames(count_matrix))
count_matrix <- count_matrix[, str_sort(colnames(count_matrix), numeric=TRUE)]

# Create a subset of samples
# samples <- c("54_1", "54_2", 
#              "66_1", "66_2", 
#              "69_1", "69_2", 
#              "70_1", "70_2", 
#              "77_1", "77_2", 
#              "A_1", "A_2", 
#              "F_1", "F_2")

# samples <- c("54_3", "54_4", "54_5", "54_6",
#              "66_3", "66_4", "66_5", "66_6",
#              "69_3", "69_4", "69_5", "69_6",
#              "70_3", "70_4", "70_5", "70_6",
#              "77_3", "77_4", "77_5", "77_6",
#              "A_3", "A_4", "A_5", "A_6",
#              "F_3", "F_4", "F_5", "F_6")

samples <- c("54_7", "54_8", "54_9", "54_10",
             "66_7", "66_8", "66_9", "66_10",
             "69_7", "69_8", "69_9", "69_10",
             "70_7", "70_8", "70_9", "70_10",
             "77_7", "77_8", "77_9", "77_10",
             "A_7", "A_8", "A_9", "A_10",
             "F_7", "F_8", "F_9", "F_10")

count_matrix <- count_matrix[, colnames(count_matrix) %in% samples]


# turning values into numeric
mode(count_matrix) <- "integer"

# creating table of conditions
sample_names <- as.vector(colnames(count_matrix))
sample_treatment <- sapply(sample_names, get_field, field=2)
sample_strain <- sapply(sample_names, get_field, field=1)
sample_biostimulation <- sapply(sample_names, get_biostimulation)
sample_biocontrol <- sapply(sample_names, get_biocontrol)
sample_susceptibility <- sapply(sample_names, get_susceptibility)
sample_table <- data.frame(treatment=as.factor(sample_treatment), 
                           strain=as.factor(sample_strain),
                           biostimulation=as.factor(sample_biostimulation),
                           biocontrol=as.factor(sample_biocontrol),
                           susceptibility=as.factor(sample_susceptibility))
row.names(sample_table) <- sample_names

# creating dds object
dds <-  DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = sample_table,
                               design = ~ treatment)

# calculating size factors
dds <- estimateSizeFactors(dds)

# regularized log transformation
vst <- vst(dds, blind=TRUE)

# calculating distances
dist <- dist(t(assay(vst)))
mat <- as.matrix(dist)
hc <- hclust(dist)

# color palette
annotation <- sample_table[, colnames(sample_table) %in% c("biostimulation", "biocontrol", "susceptibility")]
annotation_colors <- list(biostimulation=c(negative="#D73027", neutral="#C6C6C6", positive="#1A9850"),
                          biocontrol=c(negative="#D73027", neutral="#C6C6C6", positive="#1A9850"),
                          susceptibility=c(resistant="#1A9850", susceptible="#D73027"))
heatmap_colors <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(81))
genes_colors <- rev(colorRampPalette(brewer.pal(9, "Blues"))(300)[100:200])

# Plot the heatmap
breaks <- seq(min(0), max(80), by=1)
jpeg(file="heatmap.jpeg")
pheatmap(mat,
         breaks=breaks,
         Rowv=as.dendrogram(hc),
         color=heatmap_colors,
         symm=TRUE,
         trace="none",
         annotation=annotation,
         angle_col = c("45"),
         annotation_colors=annotation_colors,
         main="48h after 2nd treatment")
dev.off()

