# 18DEC24 - Noise Cancellation

getwd()
data1 <- read.csv("mutationID_count.csv")
data2 <- read.csv("FYP_Data_with_adjusted_seqid_and_mutation_ID.csv")

# Merge the datasets based on 'mutation_ID'
colnames(data1)
colnames(data1)[colnames(data1) == "annotation"] <- "mutation_ID"
colnames(data2)

combined_data <- merge(data1, data2, by = "mutation_ID", all = TRUE) 

# Write the combined dataset to a new CSV file
write.csv(combined_data, "cancelling.csv", row.names = FALSE)

# Now that dataset is made, it is time to filter

# remove anything that has a count of 97 or 1
# Step 1: Load the data
cancelling_data <- read.csv("cancelling.csv")

# Step 2: Calculate counts for each 'mutation_ID' and check if correct
mutation_counts <- table(cancelling_data$mutation_ID)
print(mutation_counts)

# Step 3: Filter out rows with counts of 1 or 97
filtered_data <- cancelling_data[!(cancelling_data$mutation_ID %in% names(mutation_counts[mutation_counts == 1 | mutation_counts == 97])), ]
head(filtered_data)

# remove any mutation_ID which have the controls or APC1.2
# Step 1: Identify mutation_IDs with 'C' in the Strain column
mutation_ids_with_C <- unique(filtered_data$mutation_ID[grepl("C", filtered_data$Strain)])
# mutation_ids_with_C <- unique(filtered_data$mutation_ID[grepl("APC", filtered_data$Strain)])


# Step 2: Remove rows with those mutation_IDs
filtered_data2 <- filtered_data[!(filtered_data$mutation_ID %in% mutation_ids_with_C), ]

# Step 3: Check the filtered data
head(filtered_data2)

write.csv(filtered_data2, "filtered_data2.csv", row.names = FALSE)

# remove strain and count column
# Remove 'Strains' and 'count' columns
filtered_data3 <- filtered_data2[, !(colnames(filtered_data2) %in% c("Strains", "Count"))]

# Check the updated data
head(filtered_data3)
write.csv(filtered_data3, "filtered_data3.csv", row.names = FALSE)


# 4. To create such a plot, ggplot2 is required so need to load this 
# also loading dplyr for colour coding
library(ggplot2)
library(dplyr)

data <- read.csv("filtered_data3.csv")

library(ggtext)

# 1. Create strain categories based on prefixes
matrix_df <- matrix_df %>%
  mutate(Strain_Category = case_when(
    grepl("^JOG", Strain) ~ "JOG",
    grepl("^MAT", Strain) ~ "MAT",
    grepl("^NW", Strain)  ~ "NW",
    grepl("^PP", Strain)  ~ "PP",
    grepl("^WS", Strain)  ~ "WS",
    TRUE ~ "Other"  # Default for unmatched strains
  ))

# 2. Define colors for strain categories
strain_colors <- c(
  "JOG" = "red",
  "MAT" = "blue",
  "NW"  = "green",
  "PP"  = "purple",
  "WS"  = "orange",
  "Other" = "grey"
)

# 3. Convert 'Strain' to a factor and reorder by Strain_Category
matrix_df$Strain <- factor(matrix_df$Strain, levels = unique(matrix_df$Strain))

# 4. Create a named color vector for strains
strain_color_mapping <- strain_colors[matrix_df$Strain_Category]
names(strain_color_mapping) <- matrix_df$Strain

# 5. Creating the heatmap with colored strain labels
ggplot(matrix_df, aes(x = mutation_ID, y = Strain, fill = Presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("0" = "white", "1" = "black"), name = "Presence") +  # Black/white for presence/absence
  scale_y_discrete(labels = function(x) {
    # Colorize labels based on strain_color_mapping
    sprintf('<span style="color:%s;">%s</span>', strain_color_mapping[x], x)
  }) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 3),  # Adjust x-axis labels
    axis.text.y = ggtext::element_markdown(size = 3)  # Use markdown for colored labels
  ) +
  labs(
    title = "Mutation Presence Across Strains",
    x = "Mutation ID",
    y = "Strain"
  )



# Install pheatmap and load it
if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)

# 1. Expand the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$Strain)),
  Strain = unlist(data$Strain)
)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)

# Ensure binary format (just in case)
presence_absence_matrix[presence_absence_matrix > 1] <- 1  

# 3. Transpose the matrix to have strains on the y-axis
transposed_matrix <- t(presence_absence_matrix)

# 4. Create strain categories based on prefixes
strain_categories <- data.frame(
  Strain = rownames(transposed_matrix),
  Category = case_when(
    grepl("^JOG", rownames(transposed_matrix)) ~ "JOG",
    grepl("^MAT", rownames(transposed_matrix)) ~ "MAT",
    grepl("^NW", rownames(transposed_matrix))  ~ "NW",
    grepl("^PP", rownames(transposed_matrix))  ~ "PP",
    grepl("^WS", rownames(transposed_matrix))  ~ "WS",
    TRUE ~ "Other"
  )
)

# 5. Define colors for strain categories
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
  fontsize_row = 3,                  # Adjust font size for strain labels
  fontsize_col = 3,                  # Adjust font size for mutation labels
  show_colnames = TRUE,              # Show mutation IDs on x-axis
  show_rownames = TRUE,              # Show strain names on y-axis
  main = "Hierarchical Clustering of Mutation Presence Across Strains"
)



# Part 3 - Statistical analysis (think this is more relevant later)
# Performing an anova

# Perform hierarchical clustering on rows (strains)
row_clustering <- hclust(dist(transposed_matrix))  # Calculate hierarchical clustering
num_clusters <- 5  # Define the number of clusters you want to cut the tree into

# Assign cluster labels to each strain
strain_clusters <- cutree(row_clustering, k = num_clusters)

# Combine the cluster labels into a data frame
cluster_data <- data.frame(
  Strain = rownames(transposed_matrix),
  Cluster = factor(strain_clusters),  # Convert cluster assignments to a factor
  Mutation_Count = rowSums(transposed_matrix)  # Calculate total mutations per strain
)

# Perform ANOVA to test for differences in mutation counts across clusters
anova_result <- aov(Mutation_Count ~ Cluster, data = cluster_data)

# View the ANOVA summary table
summary(anova_result)

# Perform Tukey's Honest Significant Difference (HSD) test
tukey_result <- TukeyHSD(anova_result)

# View the post-hoc test results
print(tukey_result)

# Boxplot to visualize mutation counts by cluster
ggplot(cluster_data, aes(x = Cluster, y = Mutation_Count, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Mutation Counts Across Clusters",
    x = "Cluster",
    y = "Total Mutation Count"
  )


# Example: A vector of strain labels matching row names
strain_labels <- rownames(transposed_matrix) # Extract strain IDs from matrix

# Create strain_groups based on your known grouping
strain_groups <- factor(ifelse(grepl("JOG", strain_labels), "JOG",
                               ifelse(grepl("MAT", strain_labels), "MAT",
                                      ifelse(grepl("NW", strain_labels), "NW",
                                             ifelse(grepl("PP", strain_labels), "PP",
                                                    ifelse(grepl("WS", strain_labels), "WS", "Other"))))))


distance_matrix <- vegdist(transposed_matrix, method = "euclidean") # Replace 'euclidean' if needed
transposed_matrix <- as.matrix(sapply(transposed_matrix, as.numeric))


length(strain_groups) == nrow(transposed_matrix)

# Load necessary libraries
library(vegan)

# Assuming your matrix rows represent strains and columns represent features
# Compute a distance matrix
distance_matrix <- vegdist(transposed_matrix, method = "euclidean") # Replace 'euclidean' with your desired metric

# Define strain groups (ensure this corresponds to the rows in your matrix)
strain_groups <- factor(c("JOG", "MAT", "NW", "PP", "WS")) # Replace with actual group labels

# Perform the dispersion analysis
dispersion <- betadisper(distance_matrix, strain_groups)

# Test for differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)

# Visualize the dispersion
plot(dispersion)



# Specifically create a heatmap to compare MAT and JOG strain mutations

data <- read.csv("filtered_data3.csv")
# MATJOG <- unique(data$mutation_ID[grepl("W|P", data$Strain)])
MATvJOG <- data[!grepl("^W|^P|^N", data$Strain), ]

# Step 2: Remove rows that include WS, NW or PP
# MATvJOG <- data[!(data$mutation_ID %in% MATJOG), ]
head(MATvJOG)


# 1. Expand the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(MATvJOG$mutation_ID, lengths(MATvJOG$Strain)),
  Strain = unlist(MATvJOG$Strain)
)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)

# Ensure binary format (just in case)
presence_absence_matrix[presence_absence_matrix > 1] <- 1  

# 3. Create strain categories based on prefixes
strain_categories <- data.frame(
  Strain = colnames(presence_absence_matrix),
  Category = case_when(
    grepl("^JOG", colnames(presence_absence_matrix)) ~ "JOG",
    grepl("^MAT", colnames(presence_absence_matrix)) ~ "MAT"
  )
)

# 4. Define colors for strain categories
strain_colors <- list(
  Category = c(
    "JOG" = "red",
    "MAT" = "blue"
  )
)

# 5. Create an annotation for strains (columns)
annotation_col <- strain_categories
rownames(annotation_col) <- annotation_col$Strain  # Match rownames to column names of the matrix
annotation_col <- annotation_col[, -1, drop = FALSE]  # Remove 'Strain' column

# 6. Plot the heatmap with clustering and colored strain annotations
pheatmap(
  presence_absence_matrix,           # Original presence-absence matrix (not transposed)
  color = c("white", "black"),       # White for absence (0), Black for presence (1)
  annotation_col = annotation_col,   # Add strain annotations on columns
  annotation_colors = strain_colors, # Color mapping for strain categories
  cluster_rows = TRUE,               # Hierarchical clustering of mutation IDs (y-axis)
  cluster_cols = TRUE,               # Hierarchical clustering of strains (x-axis)
  legend_breaks = c(0, 1),              # Ensure only 0 and 1 appear in the legend
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 5,                  # Adjust font size for mutation labels
  fontsize_col = 7,                  # Adjust font size for strain labels
  angle_col=0,
  show_colnames = TRUE,              # Show strain names on x-axis
  show_rownames = TRUE,              # Show mutation IDs on y-axis
  main = "Hierarchical Clustering of Specific Mutation Presence Comparing MAT vs JOG Strains"
)


# Specifically create a heatmap to compare NW and WS (acetic-evolved) strain mutations

data2 <- read.csv("filtered_data3.csv")
print(data2$Strain)

# Remove strains starting with MAT, JOG, or PP
# NWvWS <- data2[!grepl("^MAT|^JOG|^PP", data2$Strain), ]
NWWS <- unique(data2$mutation_ID[grepl("M|J|P", data2$Strain)])
NWvWS <- data2[!(data2$mutation_ID %in% NWWS), ]


# Print the updated Strain column
print(NWvWS$Strain)



# 1. Expand the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(NWvWS$mutation_ID, lengths(NWvWS$Strain)),
  Strain = unlist(NWvWS$Strain)
)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)

# Ensure binary format (just in case)
presence_absence_matrix[presence_absence_matrix > 1] <- 1  

# 3. Create strain categories based on prefixes
strain_categories <- data.frame(
  Strain = colnames(presence_absence_matrix),
  Category = case_when(
    grepl("^NW", colnames(presence_absence_matrix)) ~ "NW",
    grepl("^WS", colnames(presence_absence_matrix)) ~ "WS"
  )
)

# 4. Define colors for strain categories
strain_colors <- list(
  Category = c(
    "NW" = "green",
    "WS" = "orange"
  )
)

# 5. Create an annotation for strains (columns)
annotation_col <- strain_categories
rownames(annotation_col) <- annotation_col$Strain  # Match rownames to column names of the matrix
annotation_col <- annotation_col[, -1, drop = FALSE]  # Remove 'Strain' column

# 6. Plot the heatmap with clustering and colored strain annotations
pheatmap(
  presence_absence_matrix,           # Original presence-absence matrix (not transposed)
  color = c("white", "black"),       # White for absence (0), Black for presence (1)
  annotation_col = annotation_col,   # Add strain annotations on columns
  annotation_colors = strain_colors, # Color mapping for strain categories
  cluster_rows = TRUE,               # Hierarchical clustering of mutation IDs (y-axis)
  cluster_cols = TRUE,               # Hierarchical clustering of strains (x-axis)
  legend_breaks = c(0, 1),              # Ensure only 0 and 1 appear in the legend
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 3,                  # Adjust font size for mutation labels
  fontsize_col = 5,                  # Adjust font size for strain labels
  angle_col=0,
  show_colnames = TRUE,              # Show strain names on x-axis
  show_rownames = TRUE,              # Show mutation IDs on y-axis
  main = "Hierarchical Clustering of Specific Mutation Presence Comparing Acetic-evolved Strains, NW vs WS"
)

