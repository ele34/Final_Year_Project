# 08APR25 Linear Mixed Model (LMM)
# This code has been written to conduct LMMs of the various categories judging
# how similiar strains are based on: genic or intergenic regions that have SNPs, share coding genes that have SNPs,
# share teh same areas of the genome that have CNV, share the same coding genes that have CNV, share the same exact mutation ID
# (so identical SNPs), share identcial SNPs in coding genes, share mutated non-synonomous coding genes
# all categories have been filtered to not include changes that occur in the ancestor
# The parameters of the LMM include the effect of being in different school, having different background and being in different pairs

# Summary of steps
# 1) Load necessary libraries
# 2) Make an appropriate dataset for further analysis
# 3) Make distance matrices of the parameters of the LMM
# 4) Make distance matrices for the response variables
# 5) Conduct the LMM

# 1) Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggtext)
library(pheatmap)
library(tidyverse)
library(lmerTest)




# 2) Make an appropriate dataset for further analysis

# Check and set the correct working directory (if appropriate)
getwd()
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



# 3) Make distance matrices of the selective pressures

# Background
# Making a background matrix - controls mapped as their own background
strains <- data2$Strain
background <- data2$Background
unique_strains <- unique(strains)
background_similarity_matrix <- matrix(0, nrow = length(unique_strains), ncol = length(unique_strains))
rownames(background_similarity_matrix) <- unique_strains
colnames(background_similarity_matrix) <- unique_strains

# Loop through each pair of unique strains
for (i in 1:(length(unique_strains) - 1)) {  # Loop through up to the second-to-last strain
  for (j in (i + 1):length(unique_strains)) {  # Loop from i+1 to avoid double computation
    
    # Find the schools for each strain
    background_i <- background[strains == unique_strains[i]][1]  # Get the first school for strain i
    background_j <- background[strains == unique_strains[j]][1]  # Get the first school for strain j
    
    if (background_i == background_j) {
      background_similarity_matrix[i, j] <- 1
      background_similarity_matrix[j, i] <- 1  
    }
  }
}

# Print the similarity matrix
print(background_similarity_matrix)

# Converting to a distance matrix to compare later
background_matrix <- 1 - background_similarity_matrix 
print(background_matrix)


# Schools
# Extract relevant columns
strains <- data2$Strain
schools <- data2$School

# Get the unique strains
unique_strains <- unique(strains)

# Initialize the school distance matrix
school_similarity_matrix <- matrix(0, nrow = length(unique_strains), ncol = length(unique_strains))
rownames(school_similarity_matrix) <- unique_strains
colnames(school_similarity_matrix) <- unique_strains

# Loop through each pair of unique strains
for (i in 1:length(unique_strains)) {
  for (j in 1:length(unique_strains)) {
    
    # Find the school of each strain
    school_i <- schools[strains == unique_strains[i]][1]
    school_j <- schools[strains == unique_strains[j]][1]
    
    # Assign 1 if strains come from the same school, 0 if they come from different schools
    if (school_i == school_j) {
      school_similarity_matrix[i, j] <- 1
    }
  }
}

# Print the school similairity matrix
print(school_similarity_matrix)

# converting to a distance matrix to compare later
school_matrix <- 1 - school_similarity_matrix 
print(school_matrix)


# Pairs
# Read the CSV file
df <- read.csv("pairs_analysis.csv", stringsAsFactors = FALSE)
print(df)

# Split strains into a list
strain_list <- strsplit(df$Starins, ", ")

# Create a vector of unique strains
all_strains <- unique(unlist(strain_list))
print(all_strains)

# Create an empty similarity matrix
similarity_matrix <- matrix(0, nrow = length(all_strains), ncol = length(all_strains),
                            dimnames = list(all_strains, all_strains))

# Fill the matrix: if two strains share a pair, assign 1
for (pair in strain_list) {
  for (i in pair) {
    for (j in pair) {
      similarity_matrix[i, j] <- 1
    }
  }
}

# Print the resulting matrix
print(similarity_matrix)
pairs_matrix <- 1 - similarity_matrix # converting to a distance matrix to compare later
print(pairs_matrix)



# 4) Make distance matrices for the response variables

# 4a) Matrix for modelling similiarity of starins based on sharing the same areas of the genome that have CNV

# Making a distance matrix for CNV for general area change
interactive_coverage_data <- read.csv ("interactive_coverage_data_all.csv") # all changes that occur
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
write.csv(new_dataset, "new.csv", row.names = FALSE)

data <- read.csv("new.csv", stringsAsFactors = FALSE) 
data$strain <- trimws(data$strain)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(data$gene, data$strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# 3. Transpose the matrix to have strains on the y-axis
transposed_matrix_CNV <- t(presence_absence_matrix)

# 4. Clean strain names (remove extra spaces)
rownames(transposed_matrix_CNV) <- trimws(rownames(transposed_matrix_CNV))
print(transposed_matrix_CNV)
cnv_data <- transposed_matrix_CNV[order(rownames(transposed_matrix_CNV)), ]
cnv_general_dist <- dist(cnv_data, method = "Jaccard")

print(cnv_general_dist) # distance matrx of general CNV change

# 4b) Matrix for modelling similiarity of starins based on sharing the same coding genes that have CNV

interactive_coverage_datagene <- read.csv ("interactive_coverage_data.csv")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

head(interactive_coverage_datagene)

new_dataset <- interactive_coverage_datagene %>%
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
write.csv(new_dataset, "new_geneonly.csv", row.names = FALSE)

# Read the CSV file
data <- read.csv("new_geneonly.csv", stringsAsFactors = FALSE)
data$strain <- trimws(data$strain)  # Trim whitespace

# Remove rows where 'strain' contains 'C'
# data <- data[!grepl("C", data$strain), ] # keep this if using schools

# Create the binary presence-absence matrix
presence_absence_matrix <- table(data$gene, data$strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# Transpose the matrix to have strains on the y-axis
transposed_matrix_CNVgene <- t(presence_absence_matrix)

# Clean strain names (remove extra spaces)
rownames(transposed_matrix_CNVgene) <- trimws(rownames(transposed_matrix_CNVgene))
print(transposed_matrix_CNVgene)

# Making sure row order is the same
cnv_data <- transposed_matrix_CNVgene[order(rownames(transposed_matrix_CNVgene)), ]
cnv_gene_dist <- dist(cnv_data, method = "Jaccard")

# Checking
print(cnv_gene_dist) # differences between strains in terms of whether a gene has cnv or not


# 4c) Matrix for modelling similiarity of starins based on sharing genic or intergenic regions that have SNPs

data <- read.csv("cancelling.csv")

# Filter 1: Cancelling out for mutation that occur in the APC strain
filter1 <- unique(data$gene[grepl("APC", data$Strain)])
data2 <- data[!(data$gene %in% filter1), ]
filtered_data <- data[!(data$gene %in% filter1), ]
head(filtered_data)

# 2. Expand the data to create a binary matrix
expanded_data <- data.frame(
  gene = rep(filtered_data$gene, lengths(filtered_data$Strain)), 
  Strain = unlist(filtered_data$Strain)
)

# 3. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$gene, expanded_data$Strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1

# Print the matrix
print(presence_absence_matrix)

# Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)
print(transposed_matrix)

snp_data <- transposed_matrix[order(rownames(transposed_matrix)), ]
snp_general_dist <- dist(snp_data, method = "Jaccard")
print(snp_general_dist)

# 4d) Matrix for modelling similiarity of starins based on sharing share coding genes that have SNPs

data <- read.csv("cancelling.csv")
filter4 <- read.csv("interactiveplot_updated.csv") # this dataset has all the infomation be need to filter
data2 <- filter4[filter4$type == "coding", ]
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

snp_data <- transposed_matrix[order(rownames(transposed_matrix)), ]
snp_gene_dist <- dist(snp_data, method = "Jaccard")
print(snp_gene_dist)


# 4e) Matrix for modelling similiarity of starins based on sharing 
# the same areas of the genome that have CNV, share the same coding genes that have CNV, share the same exact mutation ID
# (so identical SNPs), share identcial SNPs in coding genes, share mutated non-synonomous coding genes
# different filters can be applied for each one

# Checking dataset we will use for further analysis - this has all mutation IDs that occur, including ancestor
getwd()
data <- read.csv("cancelling.csv")
head(data)

# Remove 'Strains' and 'count' columns
data <- data[, !(colnames(data) %in% c("Strains", "Count"))]

# Check the updated data - This is our data not filtered
head(data)
write.csv(data, "Question2base.csv", row.names = FALSE)

# Choose the appropriate one
# Filter 1: Cancelling out for mutation that occur in the APC strain - # just run this if general mut id for lmm
filter1 <- unique(data$mutation_ID[grepl("APC", data$Strain)])
data2 <- data[!(data$mutation_ID %in% filter1), ]
head(data2)

# Filter 2: Only include coding mutation_ID/genes
filter3 <- read.csv("interactiveplot_updated.csv") # this dataset has all the infomation be need to filter
filter3 <- filter3[filter3$type == "coding", ]
head(filter3)
data2 <- data2[data2$mutation_ID %in% filter3$mutation_ID, ] # Keep only mutation_IDs in `data` that also appear in `filter3`
head(data2)

# Filter 3: Only include mutated coding genes that are nonsynonmous
filter5 <- read.csv("interactiveplot_updated2.csv") # this dataset has all the infomation be need to filter
data2 <- filter5[filter5$SNS == "NS", ]
head(data2)

# Creating the distance matrix
# Expand the data to create a binary matrix - choose 'gene' or 'mutation_ID' depends what we are looking at
expanded_data <- data.frame(
  # mutation_ID = rep(data2$mutation_ID, lengths(data2$Strain)),
  gene = rep(data2$gene, lengths(data2$Strain)), 
  Strain = unlist(data2$Strain)
)

# Create the binary presence-absence matrix and ensure it is - choose 'gene' or 'mutation_ID' depends what we are looking at
# presence_absence_matrix <- table(expanded_data$mutation_ID , expanded_data$Strain) # again based on what we are looking at 
# presence_absence_matrix <- table(expanded_data$gene, expanded_data$Strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1
print(presence_absence_matrix)

# Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)
print(transposed_matrix)
filter_data <- transposed_matrix[order(rownames(transposed_matrix)), ]
filter_dist <- dist(filter_data, method = "Jaccard")
print(filter_dist)



# 5) Conduct the LMM 

# checking ready for lmer

dim(background_matrix)   
dim(school_matrix)       
dim(pairs_matrix)       

# Convert a distance matrix to long format (lower triangle only, no duplicates)
matrix_to_long <- function(mat, name) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("strain_1", "strain_2", name)
  df <- df %>% filter(as.character(strain_1) < as.character(strain_2))
  return(df)
}

print(background_long)
# Convert each matrix
background_long <- matrix_to_long(background_matrix, "background_dist")
school_long <- matrix_to_long(school_matrix, "school_dist")
pairs_long <- matrix_to_long(pairs_matrix, "pair_dist")


# Choose response variable - Convert dist object to a full matrix first

# mat <- as.matrix(snp_general_dist) # General SNPs in genes
# mat <- as.matrix(snp_gene_dist) # Coding SNPs in genes
# mat <- as.matrix(cnv_general_dist) # General CNVs
# mat <- as.matrix(cnv_gene_dist) # Coding CNVs
# mat <- as.matrix(filter_dist) # General Mutation_ID / Coding mutation_ID / Mutated NS genes - depends on how filter above

# Convert to long format
long <- matrix_to_long(mat, "info_dist") 


# Merge all into one dataframe
merged_df <- reduce(
  list(long, background_long, school_long, pairs_long),
  full_join,
  by = c("strain_1", "strain_2")
)

# Check data
head(merged_df)


# Fit a linear mixed model 
model <- lmer(
  info_dist ~ background_dist + school_dist + pair_dist + 
    (1 | strain_1) + (1 | strain_2),
  data = merged_df
)

# Summarize the model and compute confidence intervals
summary(model)
confint(model)

# Visualization of key results 
# change depending on what we are looking at 

# Get tidy summary of fixed effects with confidence intervals
fixed_effects <- tidy(model, effects = "fixed", conf.int = TRUE)
fixed_effects <- fixed_effects[fixed_effects$term != "(Intercept)", ]

# Plot the estimates and 95% CIs
fixed_effects$sig <- ifelse(fixed_effects$p.value < 0.05, "Significant", "Not Significant")
ggplot(fixed_effects, aes(x = term, y = estimate, color = sig)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  # geom_text(aes(label = sprintf("p = %.50f", p.value)), hjust = -0.2, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  coord_flip() +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "All SNPs",
       x = "", y = "Estimate", color = " ")



