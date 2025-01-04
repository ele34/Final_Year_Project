# 06DEC2024 - Creating a Matrix 
# Creating a matrix displaying strains and exact mutations that happen in each strain

# 1. Checking that we are in the correct working directory for what we are doing and then load the data
getwd()

data <- read.csv("FYP_data.csv")

# 2. For the heatmap, we are doing to need a column that reprsents 'exact mutations'
# Therefore, I wil now create a column in the FYP_data.csv with exact mutations so this includes 
# seqid, position, mutation and gene and will be called Mutation ID
# Create the 'mutation_ID' column by concatenating 'seqid', 'position', 'mutation', and 'gene'
data$mutation_ID <- paste(data$seq.id, data$position, data$mutation, data$gene, sep = "_")

# 3. Crate a new CSV file with the updates
write.csv(data, "FYP_Data_with_adjusted_seqid_and_mutation_ID.csv", row.names = FALSE)

# 4. To create such a plot, ggplot2 is required so need to load this 
# also loading dplyr for colour coding
library(ggplot2)
library(dplyr)

# 5. Expanding the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$Strain)),
  Strain = unlist(data$Strain)
)

# 6. Creating the binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)

# 7. Convert to a the binary matrix to a data frame for plotting
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "Strain", "Presence")

# 8. Convert 'Presence' to a factor format so it can be treated as discrete
matrix_df$Presence <- as.factor(matrix_df$Presence)

# 10. Checking for values to make sure there are no errors
unique(matrix_df$Presence)
subset(matrix_df, Presence == 2)

# 11. Ensure 'Presence' is in factor form for plotting
matrix_df$Presence <- factor(matrix_df$Presence, levels = c(0, 1))

# 12. Creating the heatmap
ggplot(matrix_df, aes(x = mutation_ID, y = Strain, fill = Presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("white", "black"), name = "Presence") + # have black and white colours
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=1, size = 1),  
    axis.text.y = element_text(size = 3)
  ) +
  labs(title = "Mutation Presence Across Strains",
       x = "Mutation ID",
       y = "Strain")

# 13. Check if there are any warnings
warnings()

# These are fine and not corrupting the data, figure can now be exported to pdf
# Important to pull the plot window across the the left when exporting
# in order to avoid the overlaps of mutation ID labels



# 17DEC24 - Colour coding and clustering

# Part 1 - Colour coding the strains
# Load necessary library 
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 1),  # Adjust x-axis labels
    axis.text.y = ggtext::element_markdown(size = 3)  # Use markdown for colored labels
  ) +
  labs(
    title = "Mutation Presence Across Strains",
    x = "Mutation ID",
    y = "Strain"
  )


# Part 2 - Clustering 

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
  fontsize_col = 1,                  # Adjust font size for mutation labels
  show_colnames = TRUE,              # Show mutation IDs on x-axis
  show_rownames = TRUE,              # Show strain names on y-axis
  main = "Hierarchical Clustering of Mutation Presence Across Strains"
)


# Part 3 - Statistical analysis (think this is more relevant later)
# Performing an anova

# Perform hierarchical clustering on rows (strains)
row_clustering <- hclust(dist(transposed_matrix))  # Calculate hierarchical clustering
num_clusters <- 6  # Define the number of clusters you want to cut the tree into

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


