# Install required Bioconductor and CRAN packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximeta", "SummarizedExperiment", "org.Hs.eg.db"))
install.packages(c("tidyverse", "ggplot2", "pheatmap", "ggrepel", "knitr", "rmarkdown"))

cat("All required packages installed successfully!\n")