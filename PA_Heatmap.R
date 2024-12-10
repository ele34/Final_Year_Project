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
library(ggplot2)

# 5. Expand the data for a binary matrix
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$Strain)),
  Strain = unlist(data$Strain)
)

# 6. Create a binary presence-absence matrix
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$Strain)

# 7. Convert to a data frame for plotting
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "Strain", "Presence")

# 8. Convert 'Presence' to a factor to ensure it is treated as discrete
matrix_df$Presence <- as.factor(matrix_df$Presence)

# 9. Checking for values to make sure no errors - might not keep in
unique(matrix_df$Presence)
subset(matrix_df, Presence == 2)

# 10. Ensure 'Presence' is a factor for plotting
matrix_df$Presence <- factor(matrix_df$Presence, levels = c(0, 1))

# 11. Creating a heatmap
ggplot(matrix_df, aes(x = mutation_ID, y = Strain, fill = Presence)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("white", "black"), name = "Presence") + # have black and white colours
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 1),  # Smaller font for x-axis
    axis.text.y = element_text(size = 3)                         # Smaller font for y-axis
  ) +
  labs(title = "Mutation Presence Across Strains",
       x = "Mutation ID",
       y = "Strain")

# 12. Check if there are any warnings
warnings()

# These are fine and not corrupting the data, figure can now be exported to pdf




