# Load necessary libraries
library(DESeq2)
library(airway)
library(tidyverse)

# Load the built-in 'airway' dataset
# This data comes from an RNA-seq experiment on airway smooth muscle cells
data("airway")

# The 'airway' object is a SummarizedExperiment, a common genomics container
# Extract the count data (matrix of gene counts per sample)
count_data <- assay(airway)

# Extract the sample information (metadata)
col_data <- colData(airway) %>% 
  as.data.frame() %>%
  select(cell, dex) # 'dex' is the treatment column (dexamethasone vs control)

# View what we have
cat("Dimensions of count data:", dim(count_data), "(genes x samples)\n")
cat("Sample information:\n")
print(col_data)

# Save the processed data for the next steps
save(count_data, col_data, file = "data/processed/loaded_data.RData")
cat("Data loaded and saved successfully.\n")
