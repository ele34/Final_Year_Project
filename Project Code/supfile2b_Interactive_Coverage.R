# 25MAR25 - Interactive Coverage Table for Genes only

# Read in file

# Load original dataset (keep decimal values)
df <- read.csv("Table_CRC_CH62.csv", check.names = FALSE, stringsAsFactors = FALSE)  # Preserve column names

# Convert columns 6 to 102 to numeric (ensure correct type)
df[, 6:102] <- lapply(df[, 6:102], as.numeric)

# Apply function row-wise to count strains where value is at least +1 or -1 from the ancestor
df$Differences <- apply(df, 1, function(row) sum(abs(as.numeric(row[7:102]) - as.numeric(row[6])) >= 1, na.rm = TRUE))

# Find column names where the strain is at least +1 or -1 from the ancestor
df$Different_Columns <- apply(df, 1, function(row) {
  diff_cols <- names(df)[7:102][abs(as.numeric(row[7:102]) - as.numeric(row[6])) >= 1]  # Check ±1 rule
  paste(diff_cols, collapse = ", ")  # Convert to a comma-separated string
})

# Function to extract school prefixes and count unique ones
extract_schools <- function(diff_col_str) {
  if (diff_col_str == "") {
    return(c(Count = 0, Names = ""))  # No differences, no schools
  }
  diff_cols <- unlist(strsplit(diff_col_str, ", "))  # Split column names
  prefixes <- unique(substr(diff_cols, 1, 2))  # Get first two letters
  return(c(Count = length(prefixes), Names = paste(prefixes, collapse = ", ")))  # Return count & names
}

# Apply the function to create Schools and School_Names columns
school_info <- t(sapply(df$Different_Columns, extract_schools))  # Apply function row-wise
df$Schools <- as.numeric(school_info[, "Count"])  # Number of unique schools
df$School_Names <- as.character(school_info[, "Names"])  # Convert to character to avoid list issue

# Check if "c" appears in Different_Columns and set Control column
df$Control <- ifelse(grepl("c", df$Different_Columns, ignore.case = TRUE), "Yes", "No")

# Remove rows where Differences = 1 or Differences = 0
df <- df[df$Differences != 0, ]

# Save the modified dataset (keep decimal values)
write.csv(df, "filtered_ch62.csv", row.names = FALSE)

# Load necessary library
library(dplyr)

# Read in the original table and important genes table
original_table <- read.csv("Table_CRC_CH62.csv")  # Replace with actual file name
important_genes <- read.csv("filtered_ch62.csv")  # Replace with actual file name

# Check column names to ensure they match
head(important_genes)

# Filter the original table to only include genes from the important genes table
filtered_table <- original_table %>%
  filter(Start.Position %in% important_genes$Start.Position)

# Merge the Differences column from important_genes
filtered_table <- filtered_table %>%
  left_join(important_genes[, c("Start.Position", "Differences", "Different_Columns", "School_Names", "Control")], by = "Start.Position")

# Save the new table with the Differences column
write.csv(filtered_table, "filtered_genesch62.csv", row.names = FALSE)






library(dplyr)

folder_path <- "/Users/evaedwards/Final-Year-Project/Interactive_Coverage"


# Get the list of all CSV files in the folder
file_list <- list.files(path = folder_path, pattern = "^filtered_genesch\\d{2}\\.csv$", full.names = TRUE)

# Print the file list to check if the files are being selected correctly
print(file_list)

# Initialize an empty list to store data frames
combined_data <- list()

# Loop through each file, read it, and add the last two digits of the filename as the first column
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file)
  
  # Print the structure of the first file to check column names and types
  print(str(data))
  
  # Extract the last two digits of the filename
  file_name <- basename(file)  # Extract only the filename
  last_two_digits <- gsub("\\D", "", file_name)  # Remove all non-digits
  last_two_digits <- substr(last_two_digits, nchar(last_two_digits) - 1, nchar(last_two_digits))
  
  # Add the last two digits as a new column and rename it to 'Chromosome_Number'
  data$Chromosome_Number <- last_two_digits
  
  # Define the desired column order
  desired_order <- c("Chromosome_Number", "Gene.ID", "Gene.Function", 
                     "Start.Position", "End.Position", "Differences", 
                     "Different_Columns", "School_Names", "Control")
  
  # Reorder columns: First the desired ones, then everything else
  data <- data %>% select(all_of(intersect(desired_order, colnames(data))), everything())
  
  # Add this data frame to the list
  combined_data[[length(combined_data) + 1]] <- data
}

# Combine all data frames into one
final_combined_data <- bind_rows(combined_data)

# Save the combined data to a new CSV file
write.csv(final_combined_data, "plusminus_coverage_data.csv", row.names = FALSE)

# Load necessary libraries
library(dplyr)

# Load your data
plusminus_coverage_data <- read.csv("plusminus_coverage_data.csv")

# Check the column names of your dataset to ensure 'school_names' exists
colnames(plusminus_coverage_data)

# Define the function to assign backgrounds
assign_background <- function(School_Names) {
  backgrounds <- c()  # Empty vector to store assigned backgrounds
  
  if (grepl("PP", School_Names)) {
    backgrounds <- c(backgrounds, "Formic")
  }
  if (grepl("JO|MA", School_Names)) {
    backgrounds <- c(backgrounds, "Cocktail")
  }
  if (grepl("NW|WS", School_Names)) {
    backgrounds <- c(backgrounds, "Acetic")
  }
  
  # Return the backgrounds as a single string, joined by commas if multiple
  return(paste(backgrounds, collapse = ", "))
}

# Apply the function to the 'school_names' column and create a new 'Background' column
plusminus_coverage_data <- plusminus_coverage_data %>%
  mutate(Background = sapply(`School_Names`, assign_background))

# Reorder columns to place 'Background' after 'School_Names'
plusminus_coverage_data <- plusminus_coverage_data %>%
  select(Chromosome_Number, Gene.ID, Gene.Function, Start.Position, End.Position, School_Names, Background, everything())

# Check the result
head(plusminus_coverage_data)
write.csv(plusminus_coverage_data, "interactive_coverage_data.csv", row.names = FALSE)




  
