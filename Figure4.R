# 06JAN25 - Creating the 4th figure: Copy read coverage plot in order to
# demonstarte depth of coverage

# Figure 4A: Single strain coverage plot
# Check working directory
getwd()

# Change to appropriate directory and check
setwd("/Users/evaedwards/Final-Year-Project")
getwd()

# Load coverage plot without headers
coverage_data <- read.table("APC1.2.txt", header = FALSE) # fill with appropriate .txt files

# Assign custom column names
colnames(coverage_data) <- c("Chromosome name", "Position", "Coverage")

# Inspect the data to ensure it's loaded correctly
head(coverage_data)

# Ensure 'Position' and 'Coverage' columns are numeric
coverage_data$Position <- as.numeric(coverage_data$Position)
coverage_data$Coverage <- as.numeric(coverage_data$Coverage)

# Check column names and structure
print(colnames(coverage_data))
str(coverage_data)

# Edit format 
colnames(coverage_data) <- gsub("^\\s+|\\s+$", "", colnames(coverage_data))

# Load appropriate libraries
library(ggplot2)
library(dplyr)

# Using a 5kbp region for noise reduction and efficancy
window_size <- 5000

# Calculate the window group for each position in the dataset
coverage_data_summarized <- coverage_data %>%
  mutate(Window = floor(Position / window_size)) %>%
  group_by(Window) %>%
  summarize(
    Mean_Position = mean(Position),
    Mean_Coverage = mean(Coverage),
    .groups = 'drop'
  )

# Create the plot using the summarized data - remeber to change title each strain
ggplot(coverage_data_summarized, aes(x = Mean_Position, y = Mean_Coverage)) +
  geom_line() +
  labs(x = "Position", y = "Mean Coverage (5 kbp window)", title = "APC1.2 Read Coverage Plot (Summarized)") +
  theme_minimal()



# Figure 4B: All strains
# To the correct file path for the strains of interest
file_path <- "/Users/evaedwards/BAMoutputs/JOG_txt_files"
files <- list.files(file_path, pattern = "*.txt", full.names = TRUE)
print(files)

# Combine all files into one single data frame
combined_data <- lapply(files, function(file) {
  # Read the file without headers and assign column names
  df <- read.table(file, header = FALSE, col.names = c("Chromosome_name", "Position", "Coverage"))
  
  # Extract strain name from the file name
  strain_name <- gsub(pattern = "\\.txt$", replacement = "", basename(file))
  
  # Add the strain column 
  df$Strain <- strain_name
  
  return(df)
})

# might do this as schools, but not sure

# Combine all individual data frames into one
combined_data <- do.call(rbind, combined_data)

# Reorder columns to match order
combined_data <- combined_data[, c("Strain", "Chromosome_name", "Position", "Coverage")]
colnames(combined_data) <- gsub("^\\s+|\\s+$", "", colnames(combined_data))

# Create the output file
output_file <- "combined_coverage_data.txt"
write.table(combined_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
head(combined_data)
str(combined_data)
# Define the window size for summarization
window_size <- 5000

# Summarize the data according to strain and window
combined_data_summarized <- combined_data %>%
  mutate(Window = floor(Position / window_size)) %>%
  group_by(Strain, Window) %>%  # Group by both Strain and Window
  summarize(
    Mean_Position = mean(Position),
    Mean_Coverage = mean(Coverage),
    .groups = 'drop'
  )

# Plot the data on a coverage graph with colours reflecting different strains
ggplot(combined_data_summarized, aes(x = Mean_Position, y = Mean_Coverage, color = Strain)) +
  geom_line(linewidth = 1) +  # Add lines for each strain
  scale_color_manual(values = c("jog1" = "blue", "jog2" = "red", "jog3" ='orange', "jog5"='green')) +  # Customize strain colors or school
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

# Plot 2 - including all the jog strains with colours
file_path <- "/Users/evaedwards/BAMoutputs/TXT3"
files <- list.files(file_path, pattern = "*.txt", full.names = TRUE)
print(files)

# Combine all files into one single data frame
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

# Create the output file
output_file <- "combined_coverage_data.txt"
write.table(combined_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
head(combined_data)
str(combined_data)
# Define the window size for summarization
window_size <- 5000

# Summarize the data according to strain and window
combined_data_summarized <- combined_data %>%
  mutate(Window = floor(Position / window_size)) %>%
  group_by(Strain, Window) %>%  # Group by both Strain and Window
  summarize(
    Mean_Position = mean(Position),
    Mean_Coverage = mean(Coverage),
    .groups = 'drop'
  )

# Plot the data on a coverage graph with colours reflecting different strains
ggplot(combined_data_summarized, aes(x = Mean_Position, y = Mean_Coverage, color = Strain)) +
  geom_line(linewidth = 1) +  # Add lines for each strain
  scale_color_manual(values = c("jog1" = "blue", "jog2" = "red", "jog3" ='orange', "jog5"='green', "jog11"='purple', "jog13"='pink', "jog15"='darkgreen', "jog2"='orchid', "jog21"= )) +  # Customize strain colors
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

