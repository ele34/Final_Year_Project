# Figure 3: Compare MAT and JOG strain mutations

# Loading libraries
library(ggplot2)
library(dplyr)
library(ggtext)
library(pheatmap)

# Step 1: Filter PA_Heatmap to only include MAT and JOG - Figure 3A

# Read in dataset
data <- read.csv("filtered_data3.csv")

# Exclude WS, PP and NW schools
MATvJOG <- data[!grepl("^W|^P|^N", data$Strain), ]

# Expand the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(MATvJOG$mutation_ID, lengths(MATvJOG$Strain)),
  Strain = unlist(MATvJOG$Strain)
)

# Create the binary presence-absence matrix and ensure binary format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)
presence_absence_matrix[presence_absence_matrix > 1] <- 1  

# Create strain categories 
strain_categories <- data.frame(
  Strain = colnames(presence_absence_matrix),
  Category = case_when(
    grepl("^JOG", colnames(presence_absence_matrix)) ~ "JOG",
    grepl("^MAT", colnames(presence_absence_matrix)) ~ "MAT"
  )
)

# Define colors for MAT and JOG
strain_colors <- list(
  Category = c(
    "JOG" = "red",
    "MAT" = "blue"
  )
)

# Create an annotation for strains 
annotation_col <- strain_categories
rownames(annotation_col) <- annotation_col$Strain  # Match rownames to column names of the matrix
annotation_col <- annotation_col[, -1, drop = FALSE]  # Remove 'Strain' column

# Plot the heatmap with clustering and coloured strain annotations
pheatmap(
  presence_absence_matrix,           
  color = c("white", "black"),       
  annotation_col = annotation_col,   
  annotation_colors = strain_colors, 
  cluster_rows = TRUE,               
  cluster_cols = TRUE,              
  legend_breaks = c(0, 1),            
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 5,                  
  fontsize_col = 3,                  
  angle_col=0,
  show_colnames = TRUE,           
  show_rownames = TRUE,          
  main = "Hierarchical Clustering of Mutation Presence in MAT and JOG Strains"
)




# Step 2: Create a heatmap for filtered mutations only present in MAT and JOG

# Read in dataset
data <- read.csv("filtered_data3.csv")

# Remove mutations that are not unique to MAT or JOG
MATJOG <- unique(data$mutation_ID[grepl("W|P", data$Strain)])
MATvJOG2 <- data[!(data$mutation_ID %in% MATJOG), ]

# Expand the data to create a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(MATvJOG2$mutation_ID, lengths(MATvJOG2$Strain)),
  Strain = unlist(MATvJOG2$Strain)
)

# Create the binary presence-absence matrix and ensure binary format
presence_absence_matrix2 <- table(expanded_data$mutation_ID, expanded_data$Strain)
presence_absence_matrix2[presence_absence_matrix2 > 1] <- 1 

# Create strain categories 
strain_categories <- data.frame(
  Strain = colnames(presence_absence_matrix2),
  Category = case_when(
    grepl("^JOG", colnames(presence_absence_matrix2)) ~ "JOG",
    grepl("^MAT", colnames(presence_absence_matrix2)) ~ "MAT"
  )
)

# Define colors for strain categories
strain_colors <- list(
  Category = c(
    "JOG" = "red",
    "MAT" = "blue"
  )
)

# Create an annotation for strains
annotation_col <- strain_categories
rownames(annotation_col) <- annotation_col$Strain  # Match rownames to column names of the matrix
annotation_col <- annotation_col[, -1, drop = FALSE]  # Remove 'Strain' column

# Plot the heatmap with clustering and colored strain annotations
pheatmap(
  presence_absence_matrix2,           # Original presence-absence matrix (not transposed)
  color = c("white", "black"),       # White for absence (0), Black for presence (1)
  annotation_col = annotation_col,   # Add strain annotations on columns
  annotation_colors = strain_colors, # Color mapping for strain categories
  cluster_rows = TRUE,               # Hierarchical clustering of mutation IDs (y-axis)
  cluster_cols = TRUE,               # Hierarchical clustering of strains (x-axis)
  legend_breaks = c(0, 1),              # Ensure only 0 and 1 appear in the legend
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 7,                  # Adjust font size for mutation labels
  fontsize_col = 7,                  # Adjust font size for strain labels
  angle_col=0,
  show_colnames = TRUE,              # Show strain names on x-axis
  show_rownames = TRUE,              # Show mutation IDs on y-axis
  main = "Hierarchical Clustering of unique mutations to MAT and JOG Strains"
)

# Step 3: School Comparison for unique mutations

# For general mutations

# Create the matrix
expanded_data <- data.frame(
  mutation_ID = rep(MATvJOG$mutation_ID, lengths(MATvJOG$School)),
  School = unlist(MATvJOG$School)
)
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
matrix_df$Freq <- ifelse(matrix_df$Freq >= 1, 1, matrix_df$Freq) # When Freq is >= 1, giving it a score of 1 to represent presence
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))
matrix_df$Presence <- factor(matrix_df$Presence, levels = c(0, 1))

# PLot the matrix
ggplot(matrix_df, aes(x = School, y = mutation_ID, fill = Presence)) + # Flip x and y
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("white", "black"), name = "Presence") + # Black and white colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 360, hjust = 1, vjust = 1, size = 7),  
    axis.text.y = element_text(
      size = 5)
  ) +
  labs(
    title = "Mutation Presence Across Schools (MAT and JOG only)",
    x = "School", # Updated to reflect flipped axes
    y = "Mutation ID"
  )


# For mutations specific to MAT and JOG

# Create the matrix
expanded_data <- data.frame(
  mutation_ID = rep(MATvJOG2$mutation_ID, lengths(MATvJOG2$School)),
  School = unlist(MATvJOG2$School)
)
presence_absence_matrix2 <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix2))
matrix_df$Freq <- ifelse(matrix_df$Freq >= 1, 1, matrix_df$Freq) # When Freq is >= 1, giving it a score of 1 to represent presence
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))
matrix_df$Presence <- factor(matrix_df$Presence, levels = c(0, 1))

# Plot the matrix
ggplot(matrix_df, aes(x = School, y = mutation_ID, fill = Presence)) + # Flip x and y
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("white", "black"), name = "Presence") + # Black and white colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 360, hjust = 1, vjust = 1, size = 7),  
    axis.text.y = element_text(
      size = 5)
  ) +
  labs(
    title = "Specific mutations seen in JOG and MAT",
    x = "School", # Updated to reflect flipped axes
    y = "Mutation ID"
  )
