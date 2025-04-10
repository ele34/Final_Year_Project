# 07MAR25 - Figure 2a
# Create a hierarchical clustering tree to clusetr strains based on
# how similar the regions of genome are that have SNPs

# Loading libraries
library(ggplot2)
library(dplyr)
library(ggtext)
library(pheatmap)

# Check and set the correct working directory (if appropriate)
getwd()
setwd("/Users/evaedwards/Final-Year-Project/Datasets/CSV files")
# Read in dataset that has all the orginal data taken from breseq output
data <- read.csv("FYP.csv")
# Adjust data to include mutation ID in order to further categorize SNPs
data$mutation_ID <- paste(data$seq.id, data$position, data$mutation, data$gene, sep = "_")
write.csv(data, "FYP_mutid.csv", row.names = FALSE)
data2 <- read.csv("FYP_mutid.csv")
# change name to al
all <- read.csv("FYP_mutid.csv")

# Create a table of counts per exact mutation
count_table <- table(all$mutation_ID)
# Sort the table in descending order
sorted_counts <- sort(count_table, decreasing = TRUE)
# Create a data frame from the sorted counts
summary_table <- data.frame(
  annotation = names(sorted_counts),
  Count = as.numeric(sorted_counts),
  Strains = sapply(names(sorted_counts), function(mutation_ID) {
    paste(unique(all$Strain[all$mutation_ID == mutation_ID]), collapse = ", ")
  })
)
# Save the output as a new CSV
write.csv(summary_table, "mutationID_count.csv", row.names = FALSE)

# expanding these datasets
data1 <- read.csv("mutationID_count.csv")
data2 <- read.csv("FYP_mutid.csv")
# Merge the datasets based on 'mutation_ID'
colnames(data1)
colnames(data1)[colnames(data1) == "annotation"] <- "mutation_ID"
colnames(data2)
data <- merge(data1, data2, by = "mutation_ID", all = TRUE) 

# Write the combined dataset to a new CSV file
write.csv(data, "cancelling.csv", row.names = FALSE)

# Therefore, cancelling out gene that also occurin the APC strain as we are not interested
filter1 <- unique(data$gene[grepl("APC", data$Strain)])
data2 <- data[!(data$gene %in% filter1), ]
head(data2)

# 1. Expand the data to create a binary matrix 
expanded_data <- data.frame(
  gene = rep(data2$gene, lengths(data2$Strain)), # change to data2 if looking at filters, data if not
  Strain = unlist(data2$Strain)
)

# 2. Create the binary presence-absence matrix and ensure it is
presence_absence_matrix <- table(expanded_data$gene, expanded_data$Strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1
print(presence_absence_matrix)

# 3. Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)
print(transposed_matrix)

# 4. Create strain categories based on prefixes                                                                                                                               
strain_categories <- data.frame(
  Strain = rownames(transposed_matrix),
  Category = case_when(
    grepl("^JOG", rownames(transposed_matrix)) ~ "JOG", # Cocktail/JOG
    grepl("^MAT", rownames(transposed_matrix)) ~ "MAT", # Cocktail/MAT
    grepl("^NW", rownames(transposed_matrix))  ~ "NW", # Acetic/NW
    grepl("^PP", rownames(transposed_matrix))  ~ "PP", # Formic/PP
    grepl("^WS", rownames(transposed_matrix))  ~ "WS", # Acetic?WS
    TRUE ~ "Other"
  )
)

# 5. Define colors for school or background categories
strain_colors <- list(
  Category = c(
    "JOG" = "red",
    "MAT" = "blue",
    "NW"  = "green",
    "PP"  = "purple",
    "WS"  = "orange",
    "Other" = "grey"
  )
)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

# 6. Create an annotation for strains (rows)
annotation_row <- strain_categories
rownames(annotation_row) <- annotation_row$Strain  # Match rownames to row names of the matrix
annotation_row <- annotation_row[, -1, drop = FALSE]  # Remove 'Strain' column

# 7. Plot the heatmap with clustering and colored strain annotations
pheatmap(
  transposed_matrix,                 # Transposed presence-absence matrix
  color = c("white", "black"),       # White for absence (0), Black for presence (1)
  annotation_row = annotation_row,   # Add strain annotations on rows
  annotation_colors = strain_colors, # Color mapping for strain categories
  cluster_rows = TRUE,               # Hierarchical clustering of strains (y-axis)
  cluster_cols = TRUE,               # Hierarchical clustering of mutations (x-axis)
  legend_breaks = c(0, 1),              # Ensure only 0 and 1 appear in the legend
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 4,                  # Adjust font size for strain labels
  fontsize_col = 5,                  # Adjust font size for mutation labels
  show_colnames = TRUE,              # Show mutation IDs on x-axis
  show_rownames = TRUE,              # Show strain names on y-axis
  main = "Hierarchical Clustering Based on SNPs in Genes or Intergenic Regions across the genome"
)

