# 12APR25 - Analysis of starin by strain 
# The table include the individual attributes of genetic chnange, for example, amount of SNPs
# that the respective strain has across every strain so we can compare between individual
# starins, pairs of strains, between backgrounds and schools
# we are going to do a jitter plot

# 1) Read in the dataset
getwd()
setwd("/Users/evaedwards/Final-Year-Project/Datasets/CSV files")
data <- read.csv("sbys.csv")
head(data)
library(ggplot2)
library(readr)
library(patchwork)

# 2) Plot the jitter plots
# MUT

library(ggplot2)

data <- read.csv("sbys.csv")
data <- data[, 1:2]
head(data)

# Add a 'Prefix' column from the first two letters of each ID
data$Prefix <- substr(data$X, 1, 2)
# Set a dummy x value for all
data$x <- ""
# Create the plot with jitter and violin overlay
plot1 <- ggplot(data, aes(x = x, y = MUT)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter plot for all points
  geom_violin(aes(y = MUT), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "MUT", y = "", title = "Jitter Plot with Violin for MUT") +
  scale_y_continuous(breaks = seq(0, max(data$MUT), by = 5), limits = c(0, NA)) +  
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Change font size here,
    axis.text.y = element_text(size = 20)
  )
ggsave("MUT_plot_with_violin_all.pdf", plot1, width = 8, height = 6, dpi = 300)

kruskal_result <- kruskal.test(MUT ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$MUT, data$Prefix, method = "bonferroni") 


# UMUT

# UMUT Violin Plot

data <- read.csv("sbys.csv")
data <- data[, c(1, 3)]  # Keep only the necessary columns
data$Prefix <- substr(data$X, 1, 2)  # Extract 'Prefix' from the ID
data$x <- ""  # Dummy x value for a single violin

plot2 <- ggplot(data, aes(x = x, y = UMUT)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +
  geom_violin(aes(y = UMUT), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +
  labs(x = "UMUT", y = "", title = "Jitter Plot with Violin Plot for UMUT") +
  scale_y_continuous(breaks = seq(0, max(data$UMUT), by = 5), limits = c(0, NA)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Change font size here,
    axis.text.y = element_text(size = 20)
  )

# Save the plot
ggsave("UMUT_plot_with_violin_all.pdf", plot2, width = 8, height = 6, dpi = 300)


kruskal_result <- kruskal.test(UMUT ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$UMUT, data$Prefix, method = "bonferroni") 

# CHGREG

data <- read.csv("sbys.csv")
data <- data[, c(1, 4)]
data$Prefix <- substr(data$X, 1, 2)
data$x <- ""

plot3 <- ggplot(data, aes(x = x, y = CHGREG)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = CHGREG), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "CHGREG", y = "", title = "Jitter Plot with Violin Plot for CHGREG") +
  scale_y_continuous(breaks = seq(0, max(data$CHGREG), by = 5), limits = c(0, NA)) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("CHGREG_plot_with_violin_all.pdf", plot3, width = 8, height = 6, dpi = 300)

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

# Create the plot with jitter and violin overlay
plot4 <- ggplot(data, aes(x = x, y = UCHGREG)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = UCHGREG), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "UCHGREG", y = "", title = "Jitter Plot with Violin Plot for UCHGREG") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGREG), by = 1), limits = c(0, NA)) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("UCHGREG_plot_with_violin_all.pdf", plot4, width = 8, height = 6, dpi = 300)


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

# Create the plot with jitter and violin overlay
plot5 <- ggplot(data, aes(x = x, y = CHGGEN)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = CHGGEN), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "CHGGEN", y = "", title = "Jitter Plot with Violin Plot for CHGGEN") +
  scale_y_continuous(breaks = seq(0, max(data$CHGGEN), by = 1), limits = c(0, NA)) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("CHGGEN_plot_with_violin_all.pdf", plot5, width = 8, height = 6, dpi = 300)


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

# Create the plot with jitter and violin overlay
plot6 <- ggplot(data, aes(x = x, y = UCHGGEN)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = UCHGGEN), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "UCHGGEN", y = "", title = "Jitter Plot with Violin Plot for UCHGGEN") +
  scale_y_continuous(breaks = seq(0, max(data$UCHGGEN), by = 1), limits = c(0, NA)) +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("UCHGGEN_plot_with_violin_all.pdf", plot6, width = 8, height = 6, dpi = 300)


kruskal_result <- kruskal.test(UCHGGEN ~ Prefix, data = data)
# Display the result
print(kruskal_result)
# Perform post-hoc test
dunn_result <- dunn.test(data$UCHGGEN, data$Prefix, method = "bonferroni") 

# CNV

# CNV Violin Plot

data <- read.csv("sbys.csv")
data <- data[, c(1, 8)]  # Keep only the necessary columns
data$Prefix <- substr(data$X, 1, 2)  # Extract 'Prefix' from the ID
data$x <- ""  # Set a dummy x value for all (e.g., "All")

# Create the violin plot with jitter overlay
plot7 <- ggplot(data, aes(x = x, y = CNV)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = CNV), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "CNV", y = "", title = "Jitter Plot with Violin Plot for CNV") +
  scale_y_continuous(breaks = seq(0, max(data$CNV), by = 10), limits = c(0, NA)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Change font size here,
    axis.text.y = element_text(size = 20)
  )

# Save the plot
ggsave("CNV_plot_with_violin_all.pdf", plot7, width = 8, height = 6, dpi = 300)


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

# Create the plot with jitter and violin overlay
plot8 <- ggplot(data, aes(x = x, y = UCNV)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = UCNV), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "UCNV", y = "", title = "Jitter Plot with Violin Plot for CNV") +
  scale_y_continuous(breaks = seq(0, max(data$UCNV), by = 10), limits = c(0, NA)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),  # Change font size here,
    axis.text.y = element_text(size = 20)
  )


# Save the plot
ggsave("UCNV_plot_with_violin_all.pdf", plot8, width = 8, height = 6, dpi = 300)

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

# Create the plot with jitter and violin overlay
plot8 <- ggplot(data, aes(x = x, y = UCNV)) +
  geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +  # Jitter points
  geom_violin(aes(y = UCNV), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "UCNV", y = "", title = "Jitter Plot with Violin Plot for UCNV") +
  scale_y_continuous(breaks = seq(0, max(data$UCNV), by = 10), limits = c(0, NA)) +
  theme_bw() +  # Change background to white
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Y-axis grid lines
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),  # Minor grid lines
    panel.background = element_rect(fill = "white")  # Ensure white background
  )


# Save the plot
ggsave("CNVG_plot_with_violin_all.pdf", plot9, width = 8, height = 6, dpi = 300)


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

# Create the plot with jitter and violin plot overlay
plot10 <- ggplot(data, aes(x = x, y = UCNVG)) +
geom_jitter(aes(color = Prefix), width = 0.5, height = 0, alpha = 1, size = 4) +
  geom_violin(aes(y = UCNVG), 
              width = 0.5, 
              alpha = 0.2, 
              fill = "gray", 
              color = NA) +  # Violin plot summarizing the whole dataset
  labs(x = "UCNVG", y = "", title = "Jitter Plot with Violin Plot for UCNVG") +
  scale_y_continuous(breaks = seq(0, max(data$UCNVG), by = 5), limits = c(0, NA)) +  # Set breaks at 5, 10, 15, etc.
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave("UCNVG_plot_with_violin_all.pdf", plot10, width = 8, height = 6, dpi = 300)


kruskal_result <- kruskal.test(UCNVG ~ Prefix, data = data)
# Display the result
print(kruskal_result)
dunn_result <- dunn.test(data$UCNVG, data$Prefix, method = "bonferroni") 


# Making a combined plot - choose one based on what you want to display
# I went for number 3 to summerise the findings
combined_plot <- (plot1 | plot2 | plot3 | plot4 | plot5 | plot6 | plot7 | plot8 | plot9 | plot10) 
combined_plot <- (plot1 | plot2 | plot7 | plot8 | plot3 | plot4 | plot5 ) 
combined_plot <- (plot1 | plot2 | plot7 | plot8  ) 


ggsave("strainbystrain_schoolviolin.pdf", combined_plot, width = 20, height = 20, dpi = 300)

