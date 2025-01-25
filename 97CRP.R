# This section is where going to combine all txt files into one big figure
# with 7 columns and 97 strains

# Get and set the appropriate working directory
getwd()
setwd("/Users/evaedwards/Final-Year-Project/Datasets/TXT/TXTcompressed")

# Load data - read all compressed txt files
coverage_data <- read.table("NW13.txt", header = FALSE)
head(coverage_data)
# Assign custom column names
colnames(coverage_data) <- c("Chromosome_name", "Position", "Coverage")

# Inspect the data to ensure it's loaded correctly
head(coverage_data)

# Ensure 'Position' and 'Coverage' columns are numeric
coverage_data$Position <- as.numeric(coverage_data$Position)
coverage_data$Coverage <- as.numeric(coverage_data$Coverage)
library(dplyr)
coverage_data <- coverage_data %>%
  mutate(Chromosome_name = as.numeric(sub("^CP", "", Chromosome_name)))  # Remove "CP" and convert to numeric

# Making a note of the mean coverage
median_y <- median(coverage_data$Coverage)
print(median_y)

# Check column names to ensure they're correct
print(colnames(coverage_data))

# Check the structure of your data
str(coverage_data)

# Remove leading/trailing spaces from column names
colnames(coverage_data) <- gsub("^\\s+|\\s+$", "", colnames(coverage_data))

# Load ggplot2 library
library(ggplot2)

# A plot of all mapped on one plot

# Define the window size (5 kbp)
window_size <- 5000

# Create an offset for each chromosome so that their positions are continuous across chromosomes
coverage_data <- coverage_data %>%
  arrange(Chromosome_name, Position) %>%  # Ensure sorted data
  group_by(Chromosome_name) %>%  # Group by chromosome
  mutate(
    Chromosome_max = max(Position)  # Get the max position for each chromosome
  ) %>%
  ungroup() %>%
  mutate(
    Chromosome_offset = cumsum(lag(Chromosome_max, default = 0)),  # Cumulative sum of max positions
    Adjusted_Position = Position + Chromosome_offset  # Add the offset to each position
  )
head(coverage_data)

# Summarize the data: calculate the mean coverage in 5kbp windows, grouped by Chromosome_name and Window
coverage_data_summarized <- coverage_data %>%
  mutate(Window = floor(Position / window_size)) %>%
  group_by(Chromosome_name, Window) %>%
  summarize(
    Mean_Coverage = mean(Coverage),
    Start_position = min(Position),
    .groups = 'drop'
  )
head(coverage_data_summarized)

# Save as a tab-delimited text file
write.table(coverage_data_summarized, "NW13compressed.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Combine all txt

# Plot the data
ggplot(coverage_data_summarized, aes(x = Start_position, y = Mean_Coverage)) +
  geom_line() +
  geom_hline(aes(yintercept = median_y, linetype = "Median Coverage"), color = "purple") +
  geom_hline(aes(yintercept = median_y*2, linetype = "2x Median Read Depth"), color = "red") +
  geom_hline(aes(yintercept = median_y*3, linetype = "3x Median Read Depth"), color = "blue") +
  facet_wrap(~Chromosome_name, scales = "free_x") +  # Add ~ to properly use the column
  labs(
    x = "Position",
    y = "Mean Coverage",
    title = "JOG1 Read Coverage Plot (Summarized in 5 kbp Windows)", # change for strain
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_linetype_manual(values = c("dashed", "dashed")) +  
  ylim(0, 600)



# Load necessary libraries
library(ggplot2)

# Set the directory where the TXT files are stored
data_dir <- "/Users/evaedwards/Final-Year-Project/Datasets/TXT/TXTcompressed"  # Change this to your folder path
file_list <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
file_list <- file_list[order(
  substr(basename(file_list), 1, 3),  # Sort by the prefix (e.g., "JOG")
  as.numeric(gsub("\\D", "", basename(file_list)))  # Sort by the numeric part of the filename
)]
print(file_list)

# Create a directory for individual plots
individual_plots_dir <- "individual_plots"
dir.create(individual_plots_dir, showWarnings = FALSE)

# Loop through each file and save individual plots as PNG
for (file in file_list) {
  tryCatch({
    # Extract strain name from the file name
    strain_name <- tools::file_path_sans_ext(basename(file))
    
    # Remove the 'compressed' suffix if present
    strain_name <- gsub("compressed$", "", strain_name)
    
    # Determine the prefix based on the length of the strain name
    if (nchar(strain_name) >= 4 && substr(strain_name, 1, 4) %in% c("JOGC", "MATC")) {
      prefix <- substr(strain_name, 1, 4)  # Get the first 4 characters (for "JOGC" or "MATC")
    } else {
      prefix <- substr(strain_name, 1, 2)  # Get the first 2 characters for other cases (e.g., "JOG", "NW")
    }
    
    # Set the custom title based on the prefix
    if (prefix == "JO") {
      custom_title <- paste0(strain_name, " (Cocktail-evolved)")
    } else if (prefix == "NW") {
      custom_title <- paste0(strain_name, " (Acetic-evolved)")
    } else if (prefix == "MA") {
      custom_title <- paste0(strain_name, " (Cocktail-evolved)")
    } else if (prefix == "WS") {
      custom_title <- paste0(strain_name, " (Acetic-evolved)")
    } else if (prefix == "JOGC") {
      custom_title <- paste0(strain_name, " (Cocktail-evolved - control)")
    } else if (prefix == "MATC") {
      custom_title <- paste0(strain_name, " (Cocktail-evolved - control)")
    } else if (prefix == "PP") {
      custom_title <- paste0(strain_name, " (Formic-evolved)")
    } else {
      custom_title <- paste0(strain_name, " (Ancestor)")
    }
  
    print(custom_title)
    
    # Read the data
    coverage_data <- read.table(file, header = TRUE, sep = "\t")  # Adjust separator if needed
    
    # Verify the data has the required columns
    if (!all(c("Chromosome_name", "Start_position", "Mean_Coverage") %in% names(coverage_data))) {
      stop("Required columns missing in file: ", file)
    }
    
    # Calculate the median for your plot
    median_y <- median(coverage_data$Mean_Coverage, na.rm = TRUE)
    
    # Create the ggplot
    p <- ggplot(coverage_data, aes(x = Start_position, y = Mean_Coverage)) +
      geom_line() +
      geom_hline(aes(yintercept = median_y, linetype = "Median Coverage"), color = "purple") +
      geom_hline(aes(yintercept = median_y*2, linetype = "2x Median Read Depth"), color = "red") +
      geom_hline(aes(yintercept = median_y*3, linetype = "3x Median Read Depth"), color = "blue") +      
      facet_wrap(~Chromosome_name, scales = "free_x") +
      labs(
        x = "Position",
        y = "Mean Coverage",
        title = custom_title  # Use the dynamic custom title here
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_linetype_manual(values = c("dashed", "dashed", "dashed")) +
      ylim(0, 600)
    
    # Save the plot as a PNG file
    output_file <- file.path(individual_plots_dir, paste0(strain_name, ".png"))
    ggsave(output_file, p, width = 16, height = 9, dpi = 300)  # Landscape orientation
    
    cat("Saved individual plot:", output_file, "\n")
  }, error = function(e) {
    warning(paste("Error processing file:", file, ":", e$message))
  })
}

warnings()

library(magick)

# Set the directory containing individual PNGs
individual_plots_dir <- "individual_plots"
output_pdf <- "combined_plots.pdf"

# Get the list of individual PNGs
png_files <- list.files(individual_plots_dir, pattern = "\\.png$", full.names = TRUE)

# Ensure there are files to process
if (length(png_files) == 0) {
  stop("No PNG files found in directory:", individual_plots_dir)
}

# Split PNGs into groups of 4
png_chunks <- split(png_files, ceiling(seq_along(png_files) / 4))

# Initialize an empty magick image object
pdf_pages <- image_blank(width = 1, height = 1)  # Placeholder

# Combine 4 PNGs per page
for (chunk in png_chunks) {
  # Read the images
  images <- lapply(chunk, image_read)
  
  # Combine images into a 2x2 grid
  grid <- image_montage(do.call(c, images), tile = "2x2", geometry = "+10+10")
  
  # Append to the PDF pages
  pdf_pages <- c(pdf_pages, grid)
}

# Save the combined plots to a PDF
pdf_pages <- pdf_pages[-1]  # Remove the placeholder
image_write(pdf_pages, path = output_pdf, format = "pdf")

cat("Combined plots saved to:", output_pdf)
output_file <- file.path(individual_plots_dir, paste0(strain_name, ".png"))
ggsave(output_file, p, width = 16, height = 9, dpi = 300)  # Landscape orientation
cat("Saved individual plot to:", output_file, "\n")
image_write(pdf_pages, path = output_pdf, format = "pdf")
cat("Combined plots saved to:", normalizePath(output_pdf), "\n")
