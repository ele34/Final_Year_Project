# 18JAN24 - Creating an interactive table

# Purpose of this interactive tbale is to make a table where we can analyse the SNP data
# But we are going to cancel out any mutation IDs that occur in the ancestor
# this is different to looking at genes that have change, these are specific chnages that occur

# Check and set the correct working directory
getwd()
# setwd("/Users/evaedwards/Final-Year-Project/Dataset/CSV files")

# Read in data
data <- read.csv("FYP.csv")

# This part delete when data is final - Adjust data to include mutation ID
data$mutation_ID <- paste(data$seq.id, data$position, data$mutation, data$gene, sep = "_")
write.csv(data, "FYP_mutid.csv", row.names = FALSE)
data2 <- read.csv("FYP_mutid.csv")

# Extract unique mutation_ID values for rows where Strain contains "APC" - so specific changes that occur
no_apc0 <- unique(data2$mutation_ID[grepl("APC", data2$Strain)])

# Subset data2 to exclude rows where mutation_ID is in no_apc0
no_apc <- data2[!(data2$mutation_ID %in% no_apc0), ]

# Checking if it has worked
head(no_apc[grepl("APC", no_apc$Strain), ])

# We are now in a position where we have the data all of strains bar the ones that also
# show up in the APC1.2, I am now going to add a column of the counts of mutations 
# so we can plot frequency on the y axis of the interactive plot

# Count the number of unique strains for each mutation_ID
strain_counts <- aggregate(Strain ~ mutation_ID, data = no_apc, FUN = function(x) length(unique(x)))

# Merge the counts back into the original dataset as a new column
no_apc$frequency <- strain_counts$Strain[match(no_apc$mutation_ID, strain_counts$mutation_ID)]

# View the updated dataset and generate csv file
head(no_apc)
write.csv(no_apc, "interactiveplot.csv", row.names = FALSE)

data <- read.csv("interactiveplot.csv") # sorta unneccessary - might delete

# additional part adding stuff for tables
# ---------------------------------------------------------------------------------

# Load necessary libraries
library(dplyr)

# Read the dataset
data <- read.csv("interactiveplot_updated.csv")

# Create a new column that concatenates strains by mutation_ID
data_with_strains <- data %>%
  group_by(mutation_ID) %>%  # Group by mutation_ID
  mutate(strains = paste(unique(Strain), collapse = ", ")) %>%  # Concatenate strains by mutation_ID
  ungroup()  # Ungroup after the operation

# View the updated data
head(data_with_strains)

# Optionally, save the updated dataset to a new CSV file
write.csv(data_with_strains, "interactiveplot_straincount.csv", row.names = FALSE)

# Read the dataset
data <- read.csv("interactiveplot.csv")

# Add the 'type' column based on 'annotation' values
data_with_type <- data %>%
  mutate(type = case_when(
    grepl("intergenic", annotation, ignore.case = TRUE) ~ "intergenic",  # Check for 'intergenic' in annotation
    grepl("noncoding", annotation, ignore.case = TRUE) ~ "noncoding",    # Check for 'noncoding' in annotation
    TRUE ~ "coding"  # Else, set as 'coding'
  ))

data_deluxe <- data_with_type %>%
  mutate(mut_type = case_when(
    grepl("->", mutation, ignore.case = TRUE) ~ "SNP",   
    grepl("\\+", mutation, ignore.case = TRUE) ~ "insertion",   # Escape "+" symbol
    grepl("diff", mutation, ignore.case = TRUE) & !grepl("\\+", mutation) ~ "deletion",  # Exclude insertion cases
    TRUE ~ "other"  # Else, set as 'coding'
  ))

# View the updated dataset with 'type' column
head(data_deluxe)

# Optionally, save the updated dataset to a new CSV file
write.csv(data_deluxe, "interactiveplot_deluxe.csv", row.names = FALSE)

# Load the dataset
data <- read.csv("interactiveplot_deluxe.csv")

# Add a column with the count of unique strains per gene
data$gene_strain_count <- ave(data$Strain, data$gene, FUN = function(x) length(unique(x)))

# Add a column with the list of strains per gene
data$gene_strains <- ave(data$Strain, data$gene, FUN = function(x) paste(unique(x), collapse = ", "))

# View the updated dataset
head(data)

# Add a new column to the dataset
data$Control <- ifelse(grepl("C", data$gene_strains), "Yes", "No")

# List of schools to check for
schools <- c("MAT", "JOG", "NW", "PP", "WS")

# Add the Schools column
data$Schools <- apply(data, 1, function(row) {
  # Extract gene_strains value for this row
  gene_strains <- row["gene_strains"]
  
  # Find matching schools
  matched_schools <- schools[sapply(schools, function(school) grepl(school, gene_strains))]
  
  # Combine matched schools with commas
  if (length(matched_schools) > 0) {
    paste(matched_schools, collapse = ", ")
  } else {
    NA  # Assign NA if no schools found
  }
})

data$Count <- sapply(data$Schools, function(schools) {
  if (is.na(schools)) {
    0  # Assign 0 if no schools are present
  } else {
    length(strsplit(schools, ", ")[[1]])  # Count the number of schools
  }
})

# View the updated dataset
head(data)

# Optionally save the updated dataset to a new file
write.csv(data, "interactiveplot_updated.csv", row.names = FALSE)

