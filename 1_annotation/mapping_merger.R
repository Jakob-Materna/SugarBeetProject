# mapping_merger.R
# 
# Description: Help  function to merge gene names and descriptions from UniProt 
# with blast results.
# 
# Input: Blast results and UniProt results
# Output: Merged file of mapped gene names
# 
# Version: 1.00
# Date: 28.05.2024
# Author: Jakob Materna

# Load the input data
blast <- as.data.frame(read.csv("at_blast", sep="\t"))
uniprot <- as.data.frame(read.csv("idmapping.tsv", sep="\t"))

# Improve the column names
colnames(blast) <- c("ensembl_id", "id", "evalue", "bitscore", "title")
colnames(uniprot) <- c("id", "entry", "reviewed", "entry.name", "protein.name", "gene.name", "organism", "length", "Gene.Names..synonym", "Gene.Names..primary.")

# Combining the two files
merged <- merge(blast, uniprot, by="id")
merged <- subset(merged, reviewed == "reviewed")
merged$title <- ifelse(merged$Gene.Names..primary. == "", 
                       merged$protein.name, 
                       paste(merged$Gene.Names..primary., merged$protein.name, sep = " - "))

# Writing the merged list into a file 
write.table(merged, file="mapped_proteins.tsv", sep="\t", quote=FALSE, row.names=FALSE)
