# 14JAN25 - Noise Cancellation Redraft
# Recreating heatmap excluding cancelling out MAT and JOG control mutations

getwd()
data1 <- read.csv("mutationID_count.csv")
data2 <- read.csv("FYP_mutid.csv")

# Merge the datasets based on 'mutation_ID'
colnames(data1)
colnames(data1)[colnames(data1) == "annotation"] <- "mutation_ID"
colnames(data2)

combined_data <- merge(data1, data2, by = "mutation_ID", all = TRUE) 

# Write the combined dataset to a new CSV file
write.csv(combined_data, "cancelling.csv", row.names = FALSE)

# Now that dataset is made, it is time to filter

# Step 1: Load the data
cancelling_data <- read.csv("cancelling.csv")

# Step 2: Calculate counts for each 'mutation_ID' and check if correct
mutation_counts <- table(cancelling_data$mutation_ID)
print(mutation_counts)

filtered_data2A <- unique(cancelling_data$mutation_ID[grepl("APC", cancelling_data$Strain)])
filtered_data2B <- cancelling_data[!(cancelling_data$mutation_ID %in% filtered_data2A), ]

head(filtered_data2B[grepl("APC", filtered_data2B$Strain), ])

# Step 3: Filter out rows with counts of 1 or 97
# filtered_data <- cancelling_data[!(cancelling_data$mutation_ID %in% names(mutation_counts[mutation_counts == 1 | mutation_counts == 96])), ]
head(filtered_data2B)

# remove any mutation_ID which have the controls or APC1.2
# Step 1: Identify mutation_IDs with 'C' in the Strain column
# mutation_ids_with_C2 <- unique(filtered_data$mutation_ID[grepl("APC", filtered_data$Strain)])
# 
# # Step 2: Remove rows with those mutation_IDs
# filtered_data2B <- filtered_data[!(filtered_data$mutation_ID %in% mutation_ids_with_C2), ]






# Step 3: Check the filtered data
head(filtered_data2B)
unique_strains <- unique(filtered_data2B$Strain)
print(unique_strains)
num_unique_strains <- length(unique_strains)
print(num_unique_strains)

write.csv(filtered_data2B, "filtered_data2B.csv", row.names = FALSE)

# remove strain and count column
# Remove 'Strains' and 'count' columns
filtered_data3A <- filtered_data2B[, !(colnames(filtered_data2B) %in% c("Strains", "Count"))]

filtered_data3B <- filtered_data3A[!grepl("^JOGC|^MATC", filtered_data3A$Strain), ]


# Check the updated data
head(filtered_data3B)
write.csv(filtered_data3B, "filtered_data3B.csv", row.names = FALSE)
unique_strains <- unique(filtered_data3B$Strain)
print(unique_strains)
num_unique_strains <- length(unique_strains)
print(num_unique_strains)
# 4. To create such a plot, ggplot2 is required so need to load this 
# also loading dplyr for colour coding
library(ggplot2)
library(dplyr)

data <- read.csv("filtered_data3B.csv")

library(ggtext)

# 5. Expanding the data to create a binary matrix
expanded_dataB <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$Strain)),
  Strain = unlist(data$Strain)
)

# 6. Creating the binary presence-absence matrix
presence_absence_matrixB <- table(expanded_dataB$mutation_ID, expanded_dataB$Strain)

# 7. Convert to a the binary matrix to a data frame for plotting
matrix_dfB <- as.data.frame(as.table(presence_absence_matrixB))
colnames(matrix_dfB) <- c("mutation_ID", "Strain", "Presence")

# 8. Convert 'Presence' to a factor format so it can be treated as discrete
matrix_dfB$Presence <- as.factor(matrix_dfB$Presence)

# 10. Checking for values to make sure there are no errors
unique(matrix_dfB$Presence)
subset(matrix_dfB, Presence == 2)

# 11. Ensure 'Presence' is in factor form for plotting
matrix_dfB$Presence <- factor(matrix_dfB$Presence, levels = c(0, 1))

# 1. Create strain categories based on prefixes
matrix_dfB <- matrix_dfB %>%
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
matrix_dfB$Strain <- factor(matrix_dfB$Strain, levels = unique(matrix_dfB$Strain))

# 4. Create a named color vector for strains
strain_color_mappingB <- strain_colors[matrix_dfB$Strain_Category]
names(strain_color_mappingB) <- matrix_dfB$Strain

# 5. Creating the heatmap with colored strain labels
ggplot(matrix_dfB, aes(x = mutation_ID, y = Strain, fill = Presence)) +
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
expanded_dataB <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$Strain)),
  Strain = unlist(data$Strain)
)

# 2. Create the binary presence-absence matrix
presence_absence_matrix <- table(expanded_dataB$mutation_ID, expanded_dataB$Strain)

# Ensure binary format (just in case)
presence_absence_matrixB[presence_absence_matrixB > 1] <- 1  

# 3. Transpose the matrix to have strains on the y-axis
transposed_matrixB <- t(presence_absence_matrixB)

# 4. Create strain categories based on prefixes
strain_categories <- data.frame(
  Strain = rownames(transposed_matrixB),
  Category = case_when(
    grepl("^JOG", rownames(transposed_matrixB)) ~ "JOG",
    grepl("^MAT", rownames(transposed_matrixB)) ~ "MAT",
    grepl("^NW", rownames(transposed_matrixB))  ~ "NW",
    grepl("^PP", rownames(transposed_matrixB))  ~ "PP",
    grepl("^WS", rownames(transposed_matrixB))  ~ "WS",
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
  transposed_matrixB,                 # Transposed presence-absence matrix
  color = c("white", "black"),       # White for absence (0), Black for presence (1)
  annotation_row = annotation_row,   # Add strain annotations on rows
  annotation_colors = strain_colors, # Color mapping for strain categories
  cluster_rows = TRUE,               # Hierarchical clustering of strains (y-axis)
  cluster_cols = TRUE,               # Hierarchical clustering of mutations (x-axis)
  legend_breaks = c(0, 1),              # Ensure only 0 and 1 appear in the legend
  legend_labels = c("0", "1"), 
  breaks = c(-0.1, 0.5, 1), 
  fontsize_row = 3,                  # Adjust font size for strain labels
  fontsize_col = 1.5,                  # Adjust font size for mutation labels
  show_colnames = TRUE,              # Show mutation IDs on x-axis
  show_rownames = TRUE,              # Show strain names on y-axis
  main = "Hierarchical Clustering of Mutation Presence Across Strains"
)








# 21DEC24 - ANOVA and post-hoc analysis of between and within variation

# Checking the working directory 
getwd()

# Read in the filtered data and load the libraries
data <- read.csv("filtered_data3B.csv")
library(tidyr)
library(ggplot2)
library(dplyr)

# Part 1 - Comparison of Schools in terms of mutation ID presence
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Perform ANOVA to check if mutation presence varies by school
anova_result <- aov(Presence ~ School, data = matrix_df)
summary(anova_result)

# Visualize the results 
# Boxplot to show mutation presence by school
library(ggplot2)
ggplot(matrix_df, aes(x = School, y = as.numeric(Presence), fill = School)) +
  geom_boxplot() +
  labs(title = "Mutation Presence Across Schools After Controlling For Noise", x = "School", y = "Mutation Presence (0=Absent, 1=Present)") +
  theme_minimal()

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
print(tukey_result$School)

# Create a data frame for visualization
tukey_result <- data.frame(
  comparison = c("MAT-JOG", "NW-JOG", "PP-JOG", "WS-JOG", "NW-MAT", "PP-MAT", "WS-MAT", "PP-NW", "WS-NW", "WS-PP"),
  diff = c(-1.7226463, 3.5470738, 4.6692112, -1.0050891, 5.2697201, 6.3918575, 0.7175573, 1.1221374, -4.5521628, -5.6743003),
  lwr = c(-3.2375322, 2.0321879, 3.1543253, -2.5199750, 3.7548342, 4.8769716, -0.7973287, -0.3927485, -6.0670488, -7.1891862),
  upr = c(-0.2077604, 5.0619597, 6.1840971, 0.5097969, 6.7846060, 7.9067434, 2.2324432, 2.6370233, -3.0372769, -4.1594143),
  p_adj = c(1.653226e-02, 2.063799e-09, 3.713674e-11, 3.670969e-01, 3.713074e-11, 3.708700e-11, 6.954322e-01, 2.554998e-01, 3.714795e-11, 3.713219e-11)
)


# Visualizing which schools have significant difference and which do not - tukey's result
ggplot(tukey_result, aes(x = comparison, y = -log10(p_adj))) +
  geom_point(color = "red", size = 3) +  # Points for the p-values
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 1) +  # Line for p = 0.05
  labs(title = "Significance of Pairwise Comparisons",
       x = "Comparison",
       y = "-log10 Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Part 2 - Adding background as a variable
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Add 'background' to the expanded data
expanded_data <- merge(
  expanded_data,
  data[, c("mutation_ID", "Background")],
  by = "mutation_ID",
  all.x = TRUE
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Add background to matrix_df
matrix_df <- merge(matrix_df, data[, c("mutation_ID", "Background")], by = "mutation_ID")

# Perform ANOVA to check if mutation presence varies by school
anova_result <- aov(Presence ~ School + Background, data = matrix_df)
summary(anova_result)

anova_result <- aov(Presence ~ Background, data = matrix_df)
summary(anova_result)

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
head(tukey_result)


# Part 3 - Understanding the reletionship between strain vs school
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Add 'background' to the expanded data
expanded_data <- merge(
  expanded_data,
  data[, c("mutation_ID", "Strain")],
  by = "mutation_ID",
  all.x = TRUE
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Add background to matrix_df
matrix_df <- merge(matrix_df, data[, c("mutation_ID", "Strain")], by = "mutation_ID")

# Perform ANOVA to check if mutation presence varies by school
model <- aov(Presence ~ School + Error(School/Strain), data = matrix_df)
summary(model)

library(lme4)

# Fit a mixed-effects model
model <- lmer(Presence ~ (1|School) + (1|School:Strain), data = matrix_df)

# Extract variance components
variance_components <- as.data.frame(VarCorr(model))

# Calculate proportion of variance
total_variance <- sum(variance_components$vcov)
between_school_variance <- variance_components$vcov[1] / total_variance
within_school_variance <- variance_components$vcov[2] / total_variance

print(c(Between_School = between_school_variance, Within_School = within_school_variance))

