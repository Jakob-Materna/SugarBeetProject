# shared_DEGs.R
# 
# Description: Script used for differential expression analysis of shared genes
# 
# Input: Files created with DE_analysis.R
# Output: Figure and table of shared DEGs between phenotypes
# 
# Version: 1.00
# Date: 28.05.2024
# Author: Jakob Materna

library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(stringr)
library(RColorBrewer)

# Create possible combinations
combinations_biocontrol <- expand.grid(expression = c("positive", "negative"),
                                       phenotype = c("biocontrol"),
                                       treatment = c("6.3", "5.3", "4.3"))

combinations_biostimulation <- expand.grid(expression = c("positive", "negative"),
                                           phenotype = c("biostimulation"),
                                           treatment = c("6.3", "5.3", "4.3"))

combinations_susceptibility <- expand.grid(expression = c("susceptible", "resistant"),
                                           phenotype = c("susceptibility"),
                                           treatment = c("6.3", "5.3", "4.3"))

# Including groups of intrest
combinations <- rbind(combinations_biocontrol, combinations_biostimulation)

# Create a list of all files 
time_point <- "point_2"
file_paths <- paste0("DEG/", time_point, "/", combinations$phenotype, "/tables/diffExpr.", combinations$expression, ".", combinations$phenotype, ".treatment.", combinations$treatment, ".tab")

# Read all files into a list of data frames
data_frames <- lapply(file_paths, read.csv, sep = "\t")

# Extract gene_id column from each data frame
gene_id_lists <- lapply(data_frames, `[[`, "gene_id")

# create a result table
all_combinations <- combn(1:length(gene_id_lists), 2)
result <- as.data.frame(matrix(ncol=4, nrow=0))
colnames(result) <- c("group_1", "group_2", "number", "genes")

# Calculate the number of shared genes between combinations
for (i in (1:ncol(all_combinations))) {
  group_1 <- all_combinations[1,i]
  group_2 <- all_combinations[2,i]
  shared <- intersect(gene_id_lists[[group_1]], gene_id_lists[[group_2]])

  if (combinations$treatment[group_1] %in% c("6.3", "10.7")) {
    treat_1 <- "+T-22+Ac"
  }
  if (combinations$treatment[group_1] %in% c("5.3", "9.7")) {
    treat_1 <- "+Ac"
  }
  if (combinations$treatment[group_1] %in% c("4.3", "8.7")) {
    treat_1 <- "+T-22"
  }  
  if (combinations$treatment[group_2] %in% c("6.3", "10.7")) {
    treat_2 <- "+T-22 +Ac"
  }
  if (combinations$treatment[group_2] %in% c("5.3", "9.7")) {
    treat_2 <- "+Ac"
  }
  if (combinations$treatment[group_2] %in% c("4.3", "8.7")) {
    treat_2 <- "+T-22"
  }
  # Create the names for shared groups
  group_1_name <- paste(combinations$expression[group_1], combinations$phenotype[group_1], treat_1, sep=" ")
  group_2_name <- paste(combinations$expression[group_2], combinations$phenotype[group_2], treat_2, sep=" ")
  
  # Add the number of shared genes to a list 
  result[nrow(result) + 1,] = c(group_1_name, group_2_name, length(shared), str_c(shared, collapse = ";"))
}

# Write result table
out_dir <- paste("DEG", time_point, "shared", sep="/")
dir.create(file.path(out_dir), recursive=TRUE)
write.table(result, file=paste0(out_dir, "/shared_DEGs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Create result table
result$number <- as.numeric(result$number)
result <- result[, colnames(result) != "genes"]

# Genereate Wheel diagram
G <- graph_from_data_frame(result) %>% as_tbl_graph()
G <- delete_edges(G, which(E(G)$number <0.1))
jpeg(file=paste0(out_dir, "/bioc_2.jpeg"), width=1000, height=1000)
G %>% ggraph(layout = "circle") +
        geom_edge_fan(aes(width=number, color=number, label=number)) + # , label=number
        scale_edge_colour_gradient(low="darkgray", high="darkgray") +
        geom_node_point(color = "black", size = 4) +
        # geom_node_text(aes(label = name), nudge_x=0.1, nudge_y=0.05) +
        geom_node_text(aes(label = name), size=4, nudge_y=0.05) +
        scale_color_viridis() +
        theme_void() +
        theme(legend.position="none")
dev.off()
