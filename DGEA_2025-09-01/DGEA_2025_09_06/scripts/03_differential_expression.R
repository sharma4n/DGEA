# Load libraries and data
library(DESeq2)
library(tidyverse)
load("data/processed/breast_cancer_loaded_data.RData")

# 1. Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)

# 2. Pre-filtering: remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "\n")

# 3. Run the DESeq2 analysis
dds <- DESeq(dds)

# 4. Get the results - compare "luminal" vs "basal"
results <- results(dds, contrast = c("condition", "luminal", "basal"))
results_df <- as.data.frame(results) %>%
  rownames_to_column("gene_id")

# 5. Apply a significance filter
significant_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# 6. Save results
write_csv(results_df, "output/tables/breast_cancer_all_genes_results.csv")
write_csv(significant_genes, "output/tables/breast_cancer_significant_genes.csv")

cat("Differential expression analysis complete!\n")
cat("Total number of significant genes (padj < 0.05, |log2FC| > 1):", nrow(significant_genes), "\n")

# 7. Save the DESeq2 object for visualization
save(dds, file = "output/tables/breast_cancer_dds_object.RData")