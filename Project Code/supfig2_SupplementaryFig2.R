# 01APR25 - Hierarchical clustering based on areas of the which show CNV
# This will be a part of figure 2 to compare to areas of the genome that have SNPs

# Set working directory
setwd("/Users/evaedwards/Final-Year-Project/CNV")

# Read file
interactive_coverage_data <- read.csv ("interactive_coverage_data_all.csv")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# Differentiate between the unknown areas
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

# View the result to check if it's correct
head(new_dataset)
write.csv(new_dataset, "new.csv", row.names = FALSE)

# Read the CSV file
data <- read.csv("new.csv", stringsAsFactors = FALSE)
data$strain <- trimws(data$strain)

# Create the binary presence-absence matrix
presence_absence_matrix <- table(data$gene, data$strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)
write.csv(transposed_matrix, "transposed_matrix.csv", row.names = FALSE)

# Clean strain names (remove extra spaces)
rownames(transposed_matrix) <- trimws(rownames(transposed_matrix))

# Create strain categories based on prefixes
strain_categories <- data.frame(
  Strain = rownames(transposed_matrix),
  Category = case_when(
    grepl("^JOG", rownames(transposed_matrix)) ~ "JOG",
    grepl("^MAT", rownames(transposed_matrix)) ~ "MAT",
    grepl("^NW", rownames(transposed_matrix))  ~ "NW",
    grepl("^PP", rownames(transposed_matrix))  ~ "PP",
    grepl("^WS", rownames(transposed_matrix))  ~ "WS"
  )
)
head(strain_categories)
rownames(transposed_matrix)

# Define colors for strain categories
strain_colors <- list(
  Category = c(
    "JOG" = "red",
    "MAT" = "blue",
    "NW"  = "green",
    "PP"  = "purple",
    "WS"  = "orange"
  )
)

# Create an annotation for strains
annotation_row <- strain_categories
rownames(annotation_row) <- annotation_row$Strain
annotation_row <- annotation_row[, -1, drop = FALSE]

# Plot the heatmap
pheatmap(
  transposed_matrix,
  color = c("white", "black"),
  annotation_row = annotation_row,
  annotation_colors = strain_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("0", "1"),
  breaks = c(-0.1, 0.5, 1),
  fontsize_row = 5,
  fontsize_col = 4,
  show_colnames = TRUE,
  show_rownames = TRUE,
  main = "Hierarchical clustering of strains based on regions that have CNV"
)


