# DE_analysis.R
# 
# Description: Script used for differential expression analysis
# 
# Input: File with Gene names and descriptions and count matrix
# Output: Heatplots, MA plots and Tables for all samples
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

# Define Phenotype vectors
biostimulation <- list(
  positive = c("69", "77", "F"),
  neutral = c("66", "70", "A"),
  negative = c("54")
)

biocontrol <- list(
  positive = c("54", "69", "A"),
  neutral = c("77"),
  negative = c("66", "70", "F")
)

susceptibility <- list(
  susceptible = c("54", "69", "A", "F"),
  resistant = c("66", "70", "77")
)

# Define Time points
time <- list(point_1 <- c("1", "2"),
             point_2 <- c("3", "4", "5", "6"),
             point_3 <- c("7", "8", "9", "10"))

# Define Replicate identifiers
replicates <- c("A", "B", "C")

# Returns field (strain, treatment, replicate) from sample name
get_field <- function(name, field) {
  return(strsplit(name, "_")[[1]][field])
}
# Returns phenotype from sample name in a given phenotype list
get_phenotype <- function(name, phenotype_list) {
  sample <- get_field(name, 1)
  for (phenotype in names(phenotype_list)) {
    if (sample %in% phenotype_list[[phenotype]]) {
      return(phenotype)
    }
  }
}

# Returns a comprehensive figure title
main_renamer <- function(expression, grouping, contrast, group_1, group_2) {
  # Create a named vector for expression mapping
  expression_map <- c(
    positive = "positive",
    neutral = "neutral",
    negative = "negative",
    susceptible = "susceptible",
    resistant = "resistant"
  )
  
  # Create a named vector for treatment mapping
  treatment_map <- c(
    "6" = "Trichoderma and Aphanomyces",
    "10" = "Trichoderma and Aphanomyces",
    "5" = "Aphanomyces",
    "9" = "Aphanomyces",
    "2" = "Trichoderma",
    "4" = "Trichoderma",
    "8" = "Trichoderma",
    "1" = "control",
    "3" = "control",
    "7" = "control"
  )
  
  # Map the expression and groups
  ex <- expression_map[[expression]]
  treat_1 <- treatment_map[[as.character(group_1)]]
  treat_2 <- treatment_map[[as.character(group_2)]]
  
  # Construct the title
  result <- paste0("Most significant DEGs in strains with ", ex, " ", grouping, 
                   " when comparing treatments with ", treat_1, " vs ", treat_2)
  
  return(result)
}

# Load input files
count_matrix <- as.matrix(read.csv("input/count.matrix.tsv", sep="\t", row.names="gene_id"))
count_matrix <- count_matrix[, colnames(count_matrix) != "gene_name"]
descriptions <- as.data.frame(read.csv("input/mapped_genes.tsv", sep="\t"))

# Renaming and reordering the matrix
colnames(count_matrix) <- gsub("X", "", colnames(count_matrix))
count_matrix <- count_matrix[, str_sort(colnames(count_matrix), numeric = TRUE)]

# Remove outlier
count_matrix <- count_matrix[, colnames(count_matrix) != "69_2_B"]

# turning values into integers
mode(count_matrix) <- "integer"

# Loop over phenotypes (biostimulation, biocontrol, susceptibility)
calculate_contrasts <- function(grouping, time_point) {
  summary <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(summary) <- c("contrast", "phenotype", "total upregulated", "total downregulated", "sig. upregulated", "sig. downregulated")
  
  # Create out directories
  dir <- getwd()
  out_dir <- paste(dir, "DEG", substitute(time_point), deparse(substitute(grouping)), sep="/")
  table_dir <- paste(out_dir, "tables", sep="/")
  ma_plots_dir <- paste(out_dir, "ma_plots", sep="/")
  heatmaps_dir <- paste(out_dir, "heatmaps", sep="/")
  dir.create(file.path(table_dir), recursive=TRUE)
  dir.create(file.path(ma_plots_dir), recursive=TRUE)
  dir.create(file.path(heatmaps_dir), recursive=TRUE)
  
  # Loop over expressions 
  for (expression in names(grouping)) {
    phenotype <- grouping[[expression]]
    combinations <- expand.grid(phenotype, time_point, replicates)
    samples <- str_sort(paste0(combinations$Var1, "_", combinations$Var2, "_", combinations$Var3), numeric=TRUE)
    # samples <- str_sort(paste0(combinations$Var1, "_", combinations$Var2), numeric=TRUE)
    sample_matrix <- count_matrix[, colnames(count_matrix) %in% samples]
    
    # create table of conditions of the phenotype phenotypes
    sample_names <- as.vector(colnames(sample_matrix))
    sample_table <- data.frame(treatment = as.factor(sapply(sample_names, get_field, 2)), 
                               replicate = as.factor(sapply(sample_names, get_field, 3)), 
                               strain = as.factor( sapply(sample_names, get_field, 1)),
                               biostimulation = as.factor(sapply(sample_names, get_phenotype, biostimulation)),
                               biocontrol = as.factor(sapply(sample_names, get_phenotype, biocontrol)),
                               susceptibility = as.factor(sapply(sample_names, get_phenotype, susceptibility)))
    row.names(sample_table) <- sample_names
    
    # Create dds object
    dds <-  DESeqDataSetFromMatrix(countData = sample_matrix,
                                   colData = sample_table,
                                   design = ~ treatment)
    
    # Calculate size factors
    dds <- estimateSizeFactors(dds)
    
    # Median of Ratios normalization
    normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
    
    # Variance stabilizing transformation 
    vst <- vst(dds, blind=TRUE)
    
    # Regularized log transformation
    # rld <- rlog(dds, blind=TRUE)
    
    # Run the analysis
    dds <- DESeq(dds)
    
    # Generate contrasts without duplicates
    all_combinations <- combn(time_point, 2)
    
    # Loop over each contrast
    for (i in (1:ncol(all_combinations))) {
      # Define contrast
      contrast <- "treatment"
      group_1 <- all_combinations[2,i]
      group_2 <- all_combinations[1,i]
      contrast_pr <- c(contrast, group_1, group_2)
      
      # Extract the results for your specified contrast
      res <- results(dds, contrast=contrast_pr)
      
      # Sort results based on significance
      res_table <- res[order(res$padj),]
      
      # Convert the rownames to the first column
      gene_id <- rownames(res_table)
      res_table <- cbind(gene_id , data.frame(res_table))
      
      # save tables of significantly DEGs
      upreg <- length(rownames(res_table)[res_table$log2FoldChange > 1])
      downreg <- length(rownames(res_table)[res_table$log2FoldChange < -1])
      resSig <- subset(res_table , padj < 0.05)
      resSig_upreg <- rownames(resSig)[resSig$log2FoldChange > 1]
      resSig_downreg <- rownames(resSig)[resSig$log2FoldChange < -1]
      file_name <- paste0(table_dir, "/", paste("diffExpr", expression, deparse(substitute(grouping)), contrast, group_1, group_2, "tab", sep="."))
      write.table(resSig, file=file_name, sep="\t", quote=FALSE, row.names=FALSE)
      # file_name_up <- paste0(table_dir, "/", paste("diffExpr.up", expression, deparse(substitute(grouping)), contrast, group_1, group_2, "tab", sep="."))
      # write.table(resSig_upreg, file=file_name_up, sep="\t", quote=FALSE, row.names=FALSE)
      # file_name_down <- paste0(table_dir, "/", paste("diffExpr.down", expression, deparse(substitute(grouping)), contrast, group_1, group_2, "tab", sep="."))
      # write.table(resSig_downreg, file=file_name_down, sep="\t", quote=FALSE, row.names=FALSE)
      
      # Initiate a summary table
      summary[nrow(summary) + 1,] = c(paste0(group_1, "vs", group_2), 
                                      paste(expression, deparse(substitute(grouping))), 
                                      upreg,
                                      downreg,
                                      length(resSig_upreg), 
                                      length(resSig_downreg))
      
      # MA plot needs the dds object
      if (length(resSig_upreg) + length(resSig_downreg) == 0) break
      main <- paste("Significant DEGs in strains with", contrast, group_1,
                    "\nrelative to strains with", contrast, group_2)
      file_name <- paste0(ma_plots_dir, "/", paste("diffExpr", expression, deparse(substitute(grouping)), contrast, group_1, group_2, "sig05", "jpg", sep="."))
      jpeg(file=file_name)
      plotMA(res, main=main, alpha=0.05)
      dev.off()
      
      # Calculate top significant genes
      top <- na.omit(res_table$gene_id[res_table$padj < 0.05 & abs(res_table$log2FoldChange) > 1])[1:25]
      top_res_table <- subset(res_table, gene_id %in% top)
      top_res_table <- top_res_table[order(-top_res_table$log2FoldChange), ]
      
      allSig <- merge(normalized_counts, top_res_table, by=0)
      sigCounts <- allSig[,2:(ncol(allSig)-7)]
      row.names(sigCounts) <- allSig$Row.names
      sigCounts <- sigCounts[, colnames(sigCounts) %in% grep(paste0("_[", group_1, group_2, "]_"), sample_names, value = TRUE)]
      sigCounts <- sigCounts[rownames(top_res_table), ]
      
      # Rename gene names
      rownames(sigCounts) <- ifelse(is.na(descriptions$title[match(rownames(sigCounts), descriptions$ensembl_id)]),
                                    rownames(sigCounts),
                                    paste(rownames(sigCounts), sub("\\(.*", "", as.character(descriptions$title[match(rownames(sigCounts), descriptions$ensembl_id)])), sep=" - "))
      
      # Color palette
      annotation <- sample_table[, as.vector(colnames(sample_table)) %in% c("biostimulation", "biocontrol", "susceptibility")]
      annotation_colors <- list(biostimulation=c(negative="#D73027", neutral="#C6C6C6", positive="#1A9850"),
                                biocontrol=c(negative="#D73027", neutral="#C6C6C6", positive="#1A9850"),
                                susceptibility=c(resistant="#1A9850", susceptible="#D73027"))
      
      # Create a heatmap 
      file_name <- paste0(heatmaps_dir, "/", paste("diffExpr", expression, deparse(substitute(grouping)), contrast, group_1, group_2, "sig05", "jpg", sep="."))
      jpeg(file=file_name, width=4000, height=900) # , height=900
      main <- main_renamer(expression, deparse(substitute(grouping)), contrast, group_1, group_2)
      steps <- seq(-2.5,2.5,length.out=100)
      print(pheatmap(log2(sigCounts+1), scale="row", breaks=steps, cluster_cols=FALSE, cluster_rows=FALSE, fontsize=20, cellwidth=100 , cellheight=30, angle_col=c("45"), main=main))
      dev.off()
    }
  }
  # Write a summary statistic
  write.table(summary, file=paste0(out_dir, "/summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

# Calculate contrasts for all time points
calculate_contrasts(grouping=biostimulation, time_point=point_1)
calculate_contrasts(grouping=biocontrol, time_point=point_1)
calculate_contrasts(grouping=susceptibility, time_point=point_1)

calculate_contrasts(grouping=biostimulation, time_point=point_2)
calculate_contrasts(grouping=biocontrol, time_point=point_2)
calculate_contrasts(grouping=susceptibility, time_point=point_2)

calculate_contrasts(grouping=biostimulation, time_point=point_3)
calculate_contrasts(grouping=biocontrol, time_point=point_3)
calculate_contrasts(grouping=susceptibility, time_point=point_3)



