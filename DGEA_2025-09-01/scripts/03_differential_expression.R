# Load libraries and data
library(DESeq2)
library(tidyverse)
load("data/processed/loaded_data.RData")

# 1. Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ dex) # Model: differences explained by treatment

# 2. Pre-filtering: remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 3. Run the DESeq2 analysis (does normalization and statistics in one step)
dds <- DESeq(dds)

# 4. Get the results
# contrast: compare 'treated' vs 'untreated' levels of the 'dex' factor
results <- results(dds, contrast = c("dex", "trt", "untrt"))
results_df <- as.data.frame(results) %>%
  rownames_to_column("gene_id") # Move gene IDs from row names to a column

# 5. Apply a significance filter (Common thresholds: p-value < 0.05, |log2FC| > 1)
significant_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) # padj is adjusted p-value

# 6. Save results
write_csv(results_df, "output/tables/all_genes_results.csv")
write_csv(significant_genes, "output/tables/significant_genes.csv")

cat("Differential expression analysis complete!\n")
cat("Total number of significant genes (padj < 0.05, |log2FC| > 1):", nrow(significant_genes), "\n")
