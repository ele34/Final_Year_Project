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
unique_strains <- unique(filtered_data$Strain)
print(unique_strains)
num_unique_strains <- length(unique_strains)
print(num_unique_strains)

# remove any mutation_ID which have the controls or APC1.2
# Step 1: Identify mutation_IDs with 'C' in the Strain column
mutation_ids_with_C <- unique(filtered_data$mutation_ID[grepl("APC|MATC|JOGC", filtered_data$Strain)])

# Step 2: Remove rows with those mutation_IDs
filtered_data2 <- filtered_data[!(filtered_data$mutation_ID %in% mutation_ids_with_C), ]


# Step 3: Check the filtered data
head(filtered_data2)
unique_strains <- unique(filtered_data2$Strain)
print(unique_strains)
num_unique_strains <- length(unique_strains)
print(num_unique_strains)

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


# 21DEC24 - ANOVA and post-hoc analysis of between and within variation

# Checking the working directory 
getwd()

# Read in the filtered data and load the libraries
data <- read.csv("filtered_data3.csv")
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
  diff = c(-0.1942857, 2.6971429, 3.3828571, 0.5314286, 2.8914286, 3.5771429, 0.7257143, 0.6857143, -2.1657143, -2.8514286),
  lwr = c(-1.4732856, 1.4181429, 2.1038572, -0.7475714, 1.6124286, 2.2981429, -0.5532856, -0.5932856, -3.4447142, -4.1304285),
  upr = c(1.0847142, 3.9761428, 4.6618571, 1.8104285, 4.1704285, 4.8561428, 2.0047142, 1.9647142, -0.8867144, -1.5724286),
  p_adj = c(0.9937703, 1.134731e-07, 1.064981e-11, 0.7874982, 9.850990e-09, 6.187273e-13, 0.5294673, 0.5852519, 4.165012e-05, 1.649283e-08)
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

# possible delete
anova_result <- aov(Presence ~ School + Strain, data = matrix_df)
summary(anova_result)

anova_result2 <- aov(Presence ~ Strain, data = matrix_df)
summary(anova_result2)

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
head(tukey_result)

tukey_result2 <- TukeyHSD(anova_result2)
head(tukey_result2)

anova_result <- aov(Presence ~ School + Error(School/Strain), data = matrix_df)
summary(anova_result)


