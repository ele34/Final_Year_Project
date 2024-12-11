getwd()

# Load data without headers
coverage_data <- read.table("jog11.txt", header = FALSE)

# Assign custom column names
colnames(coverage_data) <- c("Chromosome name", "Position", "Coverage")

# Inspect the data to ensure it's loaded correctly
head(coverage_data)

# Ensure 'Position' and 'Coverage' columns are numeric
coverage_data$Position <- as.numeric(coverage_data$Position)
coverage_data$Coverage <- as.numeric(coverage_data$Coverage)

# Check column names to ensure they're correct
print(colnames(coverage_data))

# Check the structure of your data
str(coverage_data)

# Remove leading/trailing spaces from column names
colnames(coverage_data) <- gsub("^\\s+|\\s+$", "", colnames(coverage_data))

# Load ggplot2 library
library(ggplot2)

# Plot the data using ggplot2
# ggplot(coverage_data, aes(x = Position, y = Coverage)) +
#   geom_line() +
#   labs(x = "Position", y = "Coverage", title = "Read Coverage Plot") +
#   theme_minimal()

library(dplyr)
# 
# # Filter coverage data for the 5 kbp window
# coverage_data_filtered <- coverage_data %>%
#   filter(Position >= start_position & Position <= start_position + 5000)
# 
# # Create the plot
# ggplot(coverage_data_filtered, aes(x = Position, y = Coverage)) +
#   geom_line() +
#   labs(x = "Position", y = "Coverage", title = "Read Coverage Plot (5 kbp region)") +
#   theme_minimal()

# Define the window size (5 kbp)
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

# Create the plot using the summarized data
ggplot(coverage_data_summarized, aes(x = Mean_Position, y = Mean_Coverage)) +
  geom_line() +
  labs(x = "Position", y = "Mean Coverage (5 kbp window)", title = "JOG11 Read Coverage Plot (Summarized)") +
  theme_minimal()




