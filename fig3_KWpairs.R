# repeating sbys - but now for pairs

# Read your main data
data <- read.csv("sbys.csv")
head(data)

# Merge by ID column — assuming 'X' is the ID column in sbys.csv and 'ID' in pairs.csv
data <- merge(data, pairs, by.x = "X", by.y = "ID", all.x = TRUE)

# Check the merged data
head(data)

library(tidyr)
library(dplyr)
library(stringr)

# Read the pairs file
pairs_raw <- read.csv("pairs_analysis.csv")

# Split the 'Starins' column by comma, unnest into long format
pairs_long <- pairs_raw %>%
  mutate(Starins = str_split(Starins, ",\\s*")) %>%  # split by ", "
  unnest(Starins) %>%
  rename(ID = Starins)

# Trim any whitespace
pairs_long$ID <- str_trim(pairs_long$ID)

# Check the new long format
head(pairs_long)

# Merge by ID column — assuming 'X' is the ID column in sbys.csv and 'ID' in pairs.csv
data <- merge(data, pairs_long, by.x = "X", by.y = "ID", all.x = TRUE)
head(data)

kruskal.test(CNVG ~ Pair, data = data) # Run this for each different way of testing


dunn_result <- dunn.test(data$CHGREG, data$Pair.x, method = "bonferroni") 
# Extract the results table
dunn_table <- as.data.frame(dunn_result$res)
# Filter for adjusted p-values < 0.05
sig_dunn <- dunn_table[dunn_table$P.adj < 0.05, ]
# Show only significant comparisons
print(sig_dunn)