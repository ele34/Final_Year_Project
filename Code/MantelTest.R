# 01APR25 - Comparison of reletionship between starin based on SNP and CNV
# seen in genes or intergenic regions

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(vegan)
library(proxy)

# CNV in comparable matrix format
interactive_coverage_data <- read.csv ("interactive_coverage_data_all.csv")
head(interactive_coverage_data)

new_dataset <- interactive_coverage_data %>%
  # Expand the 'Different_Columns' into multiple rows
  separate_rows(Different_Columns, sep = ",") %>% 
  # Create the new columns 'strain' and 'gene'
  mutate(strain = Different_Columns, 
         gene = case_when(
           Gene.ID == "unknown" ~ paste0("unknown-", Start.Position),
           Gene.ID == "- / unknown" ~ paste0("- / unknown-", End.Position),
           Gene.ID == "unknown / -" ~ paste0("unknown-", Start.Position),
           Gene.ID == "unknown / unknown" ~ paste0("unknown / unknown-", End.Position),
           grepl(" / unknown$", Gene.ID) ~ paste0(Gene.ID, "-", End.Position),  # Gene / unknown → Gene / unknown-End.Position
           grepl("^unknown / ", Gene.ID) ~ Gene.ID,  # unknown / gene → Keep the same
           TRUE ~ Gene.ID  # Keep other Gene.ID values unchanged
         )) %>%  
  # Select relevant columns to keep in the final dataset
  select(strain, gene)

# View the result
head(new_dataset)

# View the result
# head(new_dataset)
write.csv(new_dataset, "newdata.csv", row.names = FALSE)

# 1. Read the CSV file
data <- read.csv("newdata.csv", stringsAsFactors = FALSE) # change name to appropraite csv
data$strain <- trimws(data$strain)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(data$gene, data$strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# 3. Transpose the matrix to have strains on the y-axis
transposed_matrix_CNV <- t(presence_absence_matrix)

# 4. Clean strain names (remove extra spaces)
rownames(transposed_matrix_CNV) <- trimws(rownames(transposed_matrix_CNV))
print(transposed_matrix_CNV)

# SNP data

# Loading libraries
library(ggplot2)
library(dplyr)
library(ggtext)
library(pheatmap)
library(tidyverse)
library(vegan)
library(proxy)

# Pull
getwd()
data <- read.csv("FYP.csv")
head(data)

# Filter 1: Cancelling out for mutation that occur in the APC strain
filter1 <- unique(data$gene[grepl("APC", data$Strain)])
data2 <- data[!(data$gene %in% filter1), ]
head(data2)

# 2. Expand the data to create a binary matrix
expanded_data <- data.frame(
  gene = rep(data2$gene, lengths(data2$Strain)), 
  Strain = unlist(data2$Strain)
)

# 3. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$gene, expanded_data$Strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# Print the matrix
print(presence_absence_matrix)

# 4. Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)
print(transposed_matrix)



# COMPARING

# Making sure row order is the same
cnv_data <- transposed_matrix_CNV[order(rownames(transposed_matrix_CNV)), ]
snp_data <- transposed_matrix[order(rownames(transposed_matrix)), ]

cnv_dist <- dist(cnv_data, method = "Jaccard")
snp_dist <- dist(snp_data, method = "Jaccard")

# Checking
print(cnv_dist)

# Perform Mantel test (using Spearman correlation)
mantel_result <- mantel(cnv_dist, snp_dist, method = "spearman", permutations = 999)

# Print results
print(mantel_result)