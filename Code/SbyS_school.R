# 12APR25 - Analysis of starin by strain 
# The table include the individual attributes of genetic chnange, for example, amount of SNPs
# that the respective strain has across every strain so we can compare between individual
# starins, pairs of strains, between backgrounds and schools
# we are going to do a jitter plot

# 1) Read in the dataset
getwd()
data <- read.csv("sbys.csv")
head(data)
library(ggplot2)
library(readr)
library(patchwork)

# 2) Plot the jitter plots
# MUT

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, 1:2]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot1 <- ggplot(data, aes(x = x, y = MUT)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = MUT), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,  # This removes the outliers from the boxplot
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "MUT", y = "", title = "Jitter Plot with Boxplot for MUT") +
  scale_y_continuous(breaks = seq(0, max(data$MUT), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("MUT_plot_with_boxplot_all.pdf", plot1, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(MUT ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$MUT, data$Prefix, method = "bonferroni") 


# UMUT

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 3)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot2 <- ggplot(data, aes(x = x, y = UMUT)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = UMUT), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA,
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "UMUT", y = "", title = "Jitter Plot with Boxplot for UMUT") +
  scale_y_continuous(breaks = seq(0, max(data$UMUT), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("UMUT_plot_with_boxplot_all.png", plot2, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(UMUT ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$UMUT, data$Prefix, method = "bonferroni") 

# CHGREG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 4)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot3 <- ggplot(data, aes(x = x, y = CHGREG)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = CHGREG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "CHGREG", y = "", title = "Jitter Plot with Boxplot for CHGREG") +
  scale_y_continuous(breaks = seq(0, max(data$CHGREG), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("CHGREG_plot_with_boxplot_all.png", plot3, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(CHGREG ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$CHGREG, data$Prefix, method = "bonferroni") 

# UCHGREG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 5)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot4 <- ggplot(data, aes(x = x, y = UCHGREG)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = UCHGREG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "UCHGREG", y = "", title = "Jitter Plot with Boxplot for UCHGREG") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGREG), by = 1), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("UCHGREG_plot_with_boxplot_all.png", plot4, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(UCHGREG ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$UCHGREG, data$Prefix, method = "bonferroni") 

# CHGGEN

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 6)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot5 <- ggplot(data, aes(x = x, y = CHGGEN)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = CHGGEN), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "CHGGEN", y = "", title = "Jitter Plot with Boxplot for CHGGEN") +
  scale_y_continuous(breaks = seq(0, max(data$CHGGEN), by = 1), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("CHGGEN_plot_with_boxplot_all.png", plot5, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(CHGGEN ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$CHGGEN, data$Prefix, method = "bonferroni") 

# UCHGGEN

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 7)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot6 <- ggplot(data, aes(x = x, y = UCHGGEN)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = UCHGGEN), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "UCHGGEN", y = "", title = "Jitter Plot with Boxplot for UCHGGEN") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGGEN), by = 1), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("UCHGGEN_plot_with_boxplot_all.png", plot6, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(UCHGGEN ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$UCHGGEN, data$Prefix, method = "bonferroni") 

# CNV

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 8)]  # Keep only the necessary columns
data$Prefix <- substr(data$X, 1, 2)  # Extract 'Prefix' from the ID
data$x <- ""  # Set a dummy x value for all (e.g., "All")

# Create the plot with jitter and boxplot overlay
plot7 <- ggplot(data, aes(x = x, y = CNV)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = CNV), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "CNV", y = "", title = "Jitter Plot with Boxplot for CNV") +
  scale_y_continuous(breaks = seq(0, max(data$CNV), by = 10), limits = c(0, NA)) +  # Set breaks at 10, 20, 30, etc.
  theme_minimal() +
  theme(legend.position = "none")  # Remove the legend for simplicity

# Save the plot
ggsave("CNV_plot_with_boxplot_all.png", plot7, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(CNV ~ Prefix, data = data)
# Display the result
print(kruskal_result)


# UCNV

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 9)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot8 <- ggplot(data, aes(x = x, y = UCNV)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = UCNV), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "UCNV", y = "", title = "Jitter Plot with Boxplot for UCNV") +
  scale_y_continuous(breaks = seq(0, max(data$UCNV), by = 10), limits = c(0, NA)) +  # Set breaks at 10, 20, 30, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("UCNV_plot_with_boxplot_all.png", plot8, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(UCNV ~ Prefix, data = data)
# Display the result
print(kruskal_result)

# CNVG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 10)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot9 <- ggplot(data, aes(x = x, y = CNVG)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = CNVG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "CNVG", y = "", title = "Jitter Plot with Boxplot for CNVG") +
  scale_y_continuous(breaks = seq(0, max(data$CNVG), by = 10), limits = c(0, NA)) +  # Set breaks at 10, 20, 30, etc.
  theme_minimal() +
  theme(legend.position = "none")
ggsave("CNVG_plot_with_boxplot_all.png", plot9, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(CNVG ~ Prefix, data = data)
# Display the result
print(kruskal_result)
dunn_result <- dunn.test(data$CNVG, data$Prefix, method = "bonferroni") 


# UCNVG

data <- read.csv("sbys.csv")
# Keep only the two columns needed for analysis
data <- data[, c(1, 11)]
head(data)
# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all (e.g., "All")
data$x <- ""
# Create the plot with jitter and boxplot overlay
plot10 <- ggplot(data, aes(x = x, y = UCNVG)) +
  geom_jitter(aes(color = Prefix), width = 0.1, height = 0, alpha = 1) +  # Jitter plot for all points, colored by Prefix
  geom_boxplot(aes(y = UCNVG), 
               width = 0.3, 
               alpha = 0.2, 
               outlier.colour = NA, 
                
               outlier.size = 0) +  # Boxplot summarizing the whole dataset
  labs(x = "UCNVG", y = "", title = "Jitter Plot with Boxplot for UCNVG") +
  scale_y_continuous(breaks = seq(0, max(data$UCNVG), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal()
ggsave("UCNVG_plot_with_boxplot_all.png", plot10, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(UCNVG ~ Prefix, data = data)
# Display the result
print(kruskal_result)
dunn_result <- dunn.test(data$UCNVG, data$Prefix, method = "bonferroni") 


# Making a combined plot

combined_plot <- (plot1 | plot2 | plot3 | plot4 | plot5 | plot6 | plot7 | plot8 | plot9 | plot10) 


ggsave("strainbystrain_wbox.pdf", combined_plot, width = 20, height = 20)
