# 12APR25 - Analysis of starin by strain
# This is a copy of SbyS analysis.csv but for background

# 1) Read in the dataset
getwd()
data <- read.csv("sbys.csv")
head(data)
library(ggplot2)
library(readr)
library(patchwork)

# 2) Plot the jitter plots

# MUT
# Keep only the two columns needed for analysis
# Keep only the two columns needed for analysis
# Keep only the two columns needed for analysis
data <- data[, 1:2]
head(data)

# Add a 'Prefix' column from the first four letters of each ID
data$Prefix <- substr(data$X, 1, 4)

# Set a dummy x value for all (e.g., "All")
data$x <- ""

# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))  # Default to "Other"

# Create the plot with jitter and boxplot overlay
plot1a <- ggplot(data, aes(x = x, y = MUT)) +
  # Jitter plot with different colors for each category
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +  
  # Boxplot summarizing all values (without distinction by prefix)
  geom_boxplot(aes(y = MUT), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  # Labels and scales
  labs(x = "MUT", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$MUT), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",    # JO or MA gets blue (without C)
                                "Acetic" = "green",     # NW or WS gets green
                                "Formic" = "purple",    # PP gets purple
                                "Other" = "gray"        # Default "Other" gets gray
  )) +   
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("MUT_plot_with_single_boxplot.png", plot1a, width = 8, height = 6, dpi = 300)





# UMUT

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 3)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot2a <- ggplot(data, aes(x = x, y = UMUT)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
# boxplot  
  geom_boxplot(aes(y = UMUT), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  # Labels and scales
  labs(x = "UMUT", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$UMUT), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
                              )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("UMUT_plot_school_box.png", plot2a, width = 8, height = 6, dpi = 300)


# CHGREG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 4)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot3a <- ggplot(data, aes(x = x, y = CHGREG)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = CHGREG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "CHGREG", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$CHGREG), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("CHGREG_plot_school.png", plot3a, width = 8, height = 6, dpi = 300)


# UCHGREG


data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 5)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot4a <- ggplot(data, aes(x = x, y = UCHGREG)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = UCHGREG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "UCHGREG", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGREG), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("UCHGREG_plot_school_box.png", plot4a, width = 8, height = 6, dpi = 300)


# CHGGEN


data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 6)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot5a <- ggplot(data, aes(x = x, y = CHGGEN)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = CHGGEN), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "CHGGEN", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$CHGGEN), by = 1), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none") 
# Save the plot
ggsave("CHGGEN_plot_schoo_boxl.png", plot5a, width = 8, height = 6, dpi = 300)


# UCHGGEN

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 7)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot6a <- ggplot(data, aes(x = x, y = UCHGGEN)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = UCHGGEN), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "UCHGGEN", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGGEN), by = 1), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none") 
# Save the plot
ggsave("UCHGGEN_plot_school_box.png", plot6a, width = 8, height = 6, dpi = 300)


# CNV

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 8)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot7a <- ggplot(data, aes(x = x, y = CNV)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = CNV), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "CNV", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$CNV), by = 10), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("CNV_plot_school_box.png", plot7a, width = 8, height = 6, dpi = 300)


# UCNV

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 9)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot8a <- ggplot(data, aes(x = x, y = UCNV)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = UCNV), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "UCNV", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$UCNV), by = 10), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("UCNV_plot_schoo_boxl.png", plot8a, width = 8, height = 6, dpi = 300)


# CNVG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 10)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot9a <- ggplot(data, aes(x = x, y = CNVG)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = CNVG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "CNVG", y = "", color = "Prefix") +
  scale_y_continuous(breaks = seq(0, max(data$CNVG), by = 10), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() +
  theme(legend.position = "none")
# Save the plot
ggsave("CNVG_plot_school_box.png", plot9a, width = 8, height = 6, dpi = 300)


# UCNVG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 11)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 4)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Assign a color based on the conditions
data$Color <- ifelse(grepl("C", data$Prefix), "Control", 
                     ifelse(grepl("JO|MA", data$Prefix) & !grepl("C", data$Prefix), "Cocktail", 
                            ifelse(grepl("NW|WS", data$Prefix), "Acetic", 
                                   ifelse(grepl("PP", data$Prefix), "Formic", "Other"))))
# Define a custom color scale for the categories
plot10a <- ggplot(data, aes(x = x, y = UCNVG)) +
  geom_jitter(aes(color = Color), width = 0.1, height = 0, alpha = 1) +
  # boxplot  
  geom_boxplot(aes(y = UCNVG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # Exclude outliers
               outlier.size = 0) +  # Make sure outliers are not displayed
  labs(x = "UCNVG", y = "", color = "Background") +
  scale_y_continuous(breaks = seq(0, max(data$UCNVG), by = 10), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  scale_color_manual(values = c("Control" = "red",       # C gets red
                                "Cocktail" = "blue",  # JO or MA gets blue
                                "Acetic" = "green", # NW or WS gets green
                                "Formic" = "purple"   # PP gets purple
  )) +  
  theme_minimal() 
# Save the plot
ggsave("UCNVG_plot_school_box.png", plot10a, width = 8, height = 6, dpi = 300)


# Making a combined plot

combined_plot <- (plot1a | plot2a | plot3a | plot4a | plot5a | plot6a | plot7a | plot8a | plot9a | plot10a) 


ggsave("strainbystrain_school.pdf", combined_plot, width = 20, height = 20)
