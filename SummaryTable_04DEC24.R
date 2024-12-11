#Setting the working directory
setwd("/Users/evaedwards/Strain Data")
getwd()

#Read in the dataset and call it 'all'
all <- read.csv("FYP_Data_with_adjusted_seqid_and_mutation_ID.csv")

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




# Create a table of counts per gene
# count_table <- table(all$gene)

# Sort the table in descending order
# sorted_counts <- sort(count_table, decreasing = TRUE)

# Create a data frame from the sorted counts
# summary_table <- data.frame(
#   gene = names(sorted_counts),
#   Count = as.numeric(sorted_counts),
#   StrainCount = sapply(names(sorted_counts), function(gene) {
#     length(unique(all$Strain[all$gene == gene]))  # Count of unique strains
#   }),
#   Strains = sapply(names(sorted_counts), function(gene) {
#     paste(unique(all$Strain[all$gene == gene]), collapse = ", ")  # List of strains
#   })
# )
# 
# # Save the output as a new CSV
# write.csv(summary_table, "gene_count.csv", row.names = FALSE)


# Filter out rows where StrainCount is 96, 1, or where StrainCount equals Count
filtered_table <- summary_table[!(summary_table$Count %in% c(96, 1) | 
                                    summary_table$Count == summary_table$Count), ]

# Add a new column: difference between Count and StrainCount
# filtered_table$Difference = filtered_table$Count - filtered_table$StrainCount

# Rearrange columns to make Difference the 4th column
# filtered_table <- filtered_table[, c("gene", "Count", "StrainCount", "Difference", "Strains")]

# Sort the table by the Difference column in descending order
ordered_table <- filtered_table[order(-filtered_table$Difference), ]

# Save the ordered output as a new CSV
write.csv(ordered_table, "filtered_gene_count.csv", row.names = FALSE)



