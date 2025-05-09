# Creating a figure for chromosome 58 


setwd("/Users/evaedwards/Final-Year-Project/Datasets/TXT/TXTcompressed")


# Load libraries
library(ggplot2)
library(dplyr)
library(magick)

# Read in the data
coverage_data_summarized <- read.delim("JOG2compressed.txt", header = TRUE, sep = "\t")
head(coverage_data_summarized)

# Calculate global median and mean of coverage
median_y <- median(coverage_data_summarized$Mean_Coverage, na.rm = TRUE)

# Filter for the chromosome of interest
coverage_data_filtered <- coverage_data_summarized %>%
  filter(Chromosome_name == "34458")

# Create the plot - to check if we are doing it correctly
ggplot(coverage_data_filtered, aes(x = Start_position, y = Mean_Coverage)) +
  geom_line() +
  geom_hline(aes(yintercept = median_y, linetype = "Median Coverage"), color = "purple") +
  geom_hline(aes(yintercept = median_y * 2, linetype = "2x Median Coverage"), color = "red") +
  geom_hline(aes(yintercept = median_y * 3, linetype = "3x Median Coverage"), color = "blue") +
  facet_wrap(~Chromosome_name, scales = "free_x") +
  labs(
    x = "Position",
    y = "Mean Coverage",
    title = "CP034458 JOG1 Read Coverage Plot"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "dotted")) +
  ylim(0, 1000)


# Load libraries
library(ggplot2)
library(dplyr)

# Set the directory where the TXT files are stored
data_dir <- "/Users/evaedwards/Final-Year-Project/Datasets/TXT/TXTcompressed/58"  # Adjust if needed
file_list <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
file_list <- file_list[order(
  substr(basename(file_list), 1, 3),
  as.numeric(gsub("\\D", "", basename(file_list)))
)]
print(file_list)

# Create a directory for individual plots
individual_plots_dir <- "individual_plots"
dir.create(individual_plots_dir, showWarnings = FALSE)

# Loop through each file
for (file in file_list) {
  tryCatch({
    # Extract strain name
    strain_name <- tools::file_path_sans_ext(basename(file))
    strain_name <- gsub("compressed$", "", strain_name)
    
    # Identify prefix
    if (nchar(strain_name) >= 4 && substr(strain_name, 1, 4) %in% c("JOGC", "MATC")) {
      prefix <- substr(strain_name, 1, 4)
    } else {
      prefix <- substr(strain_name, 1, 2)
    }
    
    # Set custom title
    custom_title <- switch(prefix,
                           "JO" = paste0(strain_name, " (Cocktail-evolved)"),
                           "NW" = paste0(strain_name, " (Acetic-evolved)"),
                           "MA" = paste0(strain_name, " (Cocktail-evolved)"),
                           "WS" = paste0(strain_name, " (Acetic-evolved)"),
                           "JOGC" = paste0(strain_name, " (Cocktail-evolved - control)"),
                           "MATC" = paste0(strain_name, " (Cocktail-evolved - control)"),
                           "PP" = paste0(strain_name, " (Formic-evolved)"),
                           paste0(strain_name, " (Ancestor)")
    )
    
    print(custom_title)
    
    # Read the data
    coverage_data <- read.table(file, header = TRUE, sep = "\t")
    
    # Ensure required columns exist
    if (!all(c("Chromosome_name", "Start_position", "Mean_Coverage") %in% names(coverage_data))) {
      stop("Required columns missing in file: ", file)
    }
    
    # Compute global statistics
    median_y <- median(coverage_data$Mean_Coverage, na.rm = TRUE)
    
    # Filter for CP034458 only
    filtered_data <- coverage_data %>%
      filter(Chromosome_name == "34458")
    
    # Create plot
    p <- ggplot(filtered_data, aes(x = Start_position, y = Mean_Coverage)) +
      geom_line() +
      geom_hline(aes(yintercept = median_y, linetype = "Median Coverage"), color = "purple") +
      geom_hline(aes(yintercept = median_y * 2, linetype = "2x Median Coverage"), color = "red") +
      geom_hline(aes(yintercept = median_y * 3, linetype = "3x Median Coverage"), color = "blue") +
      facet_wrap(~Chromosome_name, scales = "free_x") +
      labs(
        x = "Position",
        y = "Mean Coverage",
        title = custom_title
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_linetype_manual(values = c("dashed", "dashed", "dashed", "dotted")) +
      ylim(0, 1000) +
      theme(legend.position = "none") 
    
    # Save the plot
    output_file <- file.path(individual_plots_dir, paste0(strain_name, ".png"))
    ggsave(output_file, p, width = 16, height = 9, dpi = 300)
    
    cat("Saved individual plot:", output_file, "\n")
  }, error = function(e) {
    warning(paste("Error processing file:", file, ":", e$message))
  })
}

# 6) Combine all individual plopts into one big pdf files
# Load necessary library
library(magick)

# Set the directory containing individual PNG files
individual_plots_dir <- "individual_plots"
output_pdf <- "combined_plots.pdf"
output_png <- "combined_plots_A4_highres.png"  # High-res PNG output

# Get the list of individual PNG files
png_files <- list.files(individual_plots_dir, pattern = "\\.png$", full.names = TRUE)
print(png_files)

# Ensure there are files to process
if (length(png_files) == 0) {
  stop("No PNG files found in directory: ", individual_plots_dir)
}

# Prioritize specific files first (if needed)
priority_files <- png_files[basename(png_files) %in% c("JOGC1.png", "NW31.png")]
other_files <- png_files[!basename(png_files) %in% c("JOGC1.png", "NW31.png")]

# Final ordered list of files
ordered_pngs <- c(priority_files, other_files)

# Split into groups of 9 for 3x3 layout
png_chunks <- split(ordered_pngs, ceiling(seq_along(ordered_pngs) / 9))

# Create and collect PDF pages with 3x3 grid
pdf_pages <- lapply(png_chunks, function(chunk) {
  images <- lapply(chunk, image_read)
  image_montage(do.call(c, images), tile = "3x3", geometry = "+10+10")
})

# Combine all pages into one document
pdf_document <- image_join(pdf_pages)

# Save as a multi-page PDF
image_write(pdf_document, path = output_pdf, format = "pdf")
cat("✅ Combined PDF saved to:", normalizePath(output_pdf), "\n")

# Save high-resolution PNG version (only first page or all as needed)
final_image <- image_resize(pdf_document, "3300x2400!")  # A4 landscape
image_write(final_image, path = output_png, format = "png")
cat("✅ High-res PNG saved to:", normalizePath(output_png), "\n")











