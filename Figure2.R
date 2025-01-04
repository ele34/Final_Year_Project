# 21DEC24 - ANOVA and post-hoc analysis of between and within variation

# Checking the working directory 
getwd()

# Read in the filtered data and load the libraries
data <- read.csv("filtered_data3.csv")
library(tidyr)
library(ggplot2)
library(dplyr)

# Part 1 - Comparison of Schools in terms of mutation ID presence
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Perform ANOVA to check if mutation presence varies by school
anova_result <- aov(Presence ~ School, data = matrix_df)
summary(anova_result)

# Visualize the results 
# Boxplot to show mutation presence by school
library(ggplot2)
ggplot(matrix_df, aes(x = School, y = as.numeric(Presence), fill = School)) +
  geom_boxplot() +
  labs(title = "Mutation Presence Across Schools After Controlling For Noise", x = "School", y = "Mutation Presence (0=Absent, 1=Present)") +
  theme_minimal()

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
print(tukey_result$School)

# Create a data frame for visualization
tukey_result <- data.frame(
  comparison = c("MAT-JOG", "NW-JOG", "PP-JOG", "WS-JOG", "NW-MAT", "PP-MAT", "WS-MAT", "PP-NW", "WS-NW", "WS-PP"),
  diff = c(-0.1942857, 2.6971429, 3.3828571, 0.5314286, 2.8914286, 3.5771429, 0.7257143, 0.6857143, -2.1657143, -2.8514286),
  lwr = c(-1.4732856, 1.4181429, 2.1038572, -0.7475714, 1.6124286, 2.2981429, -0.5532856, -0.5932856, -3.4447142, -4.1304285),
  upr = c(1.0847142, 3.9761428, 4.6618571, 1.8104285, 4.1704285, 4.8561428, 2.0047142, 1.9647142, -0.8867144, -1.5724286),
  p_adj = c(0.9937703, 1.134731e-07, 1.064981e-11, 0.7874982, 9.850990e-09, 6.187273e-13, 0.5294673, 0.5852519, 4.165012e-05, 1.649283e-08)
)

# Visualizing which schools have significant difference and which do not - tukey's result
ggplot(tukey_result, aes(x = comparison, y = -log10(p_adj))) +
  geom_point(color = "red", size = 3) +  # Points for the p-values
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 1) +  # Line for p = 0.05
  labs(title = "Significance of Pairwise Comparisons",
       x = "Comparison",
       y = "-log10 Adjusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Part 2 - Adding background as a variable
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Add 'background' to the expanded data
expanded_data <- merge(
  expanded_data,
  data[, c("mutation_ID", "Background")],
  by = "mutation_ID",
  all.x = TRUE
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Add background to matrix_df
matrix_df <- merge(matrix_df, data[, c("mutation_ID", "Background")], by = "mutation_ID")

# Perform ANOVA to check if mutation presence varies by school
anova_result <- aov(Presence ~ School + Background, data = matrix_df)
summary(anova_result)

anova_result <- aov(Presence ~ Background, data = matrix_df)
summary(anova_result)

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
head(tukey_result)


# Part 3 - Understanding the reletionship between strain vs school
# Change into necessary format
expanded_data <- data.frame(
  mutation_ID = rep(data$mutation_ID, lengths(data$School)),
  School = unlist(data$School)
)

# Add 'background' to the expanded data
expanded_data <- merge(
  expanded_data,
  data[, c("mutation_ID", "Strain")],
  by = "mutation_ID",
  all.x = TRUE
)

# Transfer to matrix format
presence_absence_matrix <- table(expanded_data$mutation_ID, expanded_data$School)
matrix_df <- as.data.frame(as.table(presence_absence_matrix))
colnames(matrix_df) <- c("mutation_ID", "School", "Presence")
matrix_df$Presence <- as.factor(matrix_df$Presence)
matrix_df$Presence <- as.numeric(as.character(matrix_df$Presence))

# Add background to matrix_df
matrix_df <- merge(matrix_df, data[, c("mutation_ID", "Strain")], by = "mutation_ID")

# Perform ANOVA to check if mutation presence varies by school
model <- aov(Presence ~ School + Error(School/Strain), data = matrix_df)
summary(model)

library(lme4)

# Fit a mixed-effects model
model <- lmer(Presence ~ (1|School) + (1|School:Strain), data = matrix_df)

# Extract variance components
variance_components <- as.data.frame(VarCorr(model))

# Calculate proportion of variance
total_variance <- sum(variance_components$vcov)
between_school_variance <- variance_components$vcov[1] / total_variance
within_school_variance <- variance_components$vcov[2] / total_variance

print(c(Between_School = between_school_variance, Within_School = within_school_variance))

# possible delete
anova_result <- aov(Presence ~ School + Strain, data = matrix_df)
summary(anova_result)

anova_result2 <- aov(Presence ~ Strain, data = matrix_df)
summary(anova_result2)

# Tukey's HSD post-hoc test to identify which schools are most different
tukey_result <- TukeyHSD(anova_result)
head(tukey_result)

tukey_result2 <- TukeyHSD(anova_result2)
head(tukey_result2)

anova_result <- aov(Presence ~ School + Error(School/Strain), data = matrix_df)
summary(anova_result)


