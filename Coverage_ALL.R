getwd()

library(dplyr)
library(ggplot2)

setwd("/Users/evaedwards")

file_path <- "/Users/evaedwards/BAMoutputs/TXT2"
files <- list.files(file_path, pattern = "*.txt", full.names = TRUE)
print(files)

# Combine all files into a single data frame
combined_data <- lapply(files, function(file) {
  # Read the file without headers and assign column names
  df <- read.table(file, header = FALSE, col.names = c("Chromosome_name", "Position", "Coverage"))
  
  # Extract strain name from the file name
  strain_name <- gsub(pattern = "\\.txt$", replacement = "", basename(file))
  
  # Add the strain column
  df$Strain <- strain_name
  
  return(df)
})

# Combine all individual data frames into one
combined_data <- do.call(rbind, combined_data)

# Reorder columns to match the desired order
combined_data <- combined_data[, c("Strain", "Chromosome_name", "Position", "Coverage")]
colnames(combined_data) <- gsub("^\\s+|\\s+$", "", colnames(combined_data))


output_file <- "combined_coverage_data.txt"
write.table(combined_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
head(combined_data)
str(combined_data)
# Define the window size for summarization
window_size <- 5000

# Summarize the data by strain and window
combined_data_summarized <- combined_data %>%
  mutate(Window = floor(Position / window_size)) %>%
  group_by(Strain, Window) %>%  # Group by both Strain and Window
  summarize(
    Mean_Position = mean(Position),
    Mean_Coverage = mean(Coverage),
    .groups = 'drop'
  )

# Plot the summarized data with different colors per strain
ggplot(combined_data_summarized, aes(x = Mean_Position, y = Mean_Coverage, color = Strain)) +
  geom_line(linewidth = 1) +  # Add lines for each strain
  scale_color_manual(values = c("jog1" = "blue", "jog2" = "red", "jog3" ='orange', "jog5"='green')) +  # Customize strain colors
  labs(
    x = "Position",
    y = "Mean Coverage (5 kbp window)",
    title = "JOG Strains Read Coverage Plot (Summarized)",
    color = "Strain"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    legend.position = "top"
  )

# Ideally we would adjust this code to group all the background
# and then have the last strain as the ancestor - might need
# multiple plots as a panel - but potential here
getwd()

cd /Users/evaedwards/BAMoutputs


