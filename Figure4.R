# 06JAN25 - Creating the 4th figure: Copy read coverage plot in order to
# demonstarte depth of coverage

# Check working directory
getwd()

# Change to appropriate directory and check
setwd("/Users/evaedwards/Final-Year-Project")
getwd()

# Load coverage plot without headers
coverage_data <- read.table("NW15.txt", header = FALSE) # fill with appropriate .txt files

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
  labs(x = "Position", y = "Mean Coverage (5 kbp window)", title = "NW15 Read Coverage Plot (Summarized)") +
  theme_minimal()