# Check and set the correct working directory (if appropriate)
getwd()
setwd("/Users/evaedwards/Final-Year-Project/Datasets/CSV files")

# Read in dataset that has all the orginal data taken from breseq output
data <- read.csv("tage.csv")
head(data)

density_values <- density(data$CNV)
values <- data$CNV

# Plot the density
plot(density_values,
     main = "Density Plot with Jittered Points",
     xlab = "Value",
     ylab = "Density",
     col = "blue",
     lwd = 2)



data$Prefix<- substr(data$Strain, 1, 2)
head(data)

stripchart(values,
           method = "jitter",
           jitter = 0.2,
           color = data$Prefix,
           pch = 16,
           main = "Jitter Plot of Values",
           xlab = "Value",
           vertical = TRUE)

library(ggplot2)

ggplot(data, aes(x = "", y = values, color = Prefix)) +
  geom_jitter(width = 0.5, size = 4) +
  labs(title = "", x = "", y = "NCNV") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),       # Axis tick labels
    axis.title = element_text(size = 16),      # Axis titles
  )
