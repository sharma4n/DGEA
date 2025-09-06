# Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)

# Load the data
load("output/tables/breast_cancer_dds_object.RData")
de_results <- read_csv("output/tables/breast_cancer_all_genes_results.csv")
sig_genes <- read_csv("output/tables/breast_cancer_significant_genes.csv")

# 1. MA Plot
png("output/figures/breast_cancer_ma_plot.png", width=800, height=600)
plotMA(results(dds), main="MA Plot: Luminal vs Basal Breast Cancer", ylim=c(-8,8))
dev.off()

# 2. Volcano Plot
volcano_plot <- de_results %>%
  mutate(significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Luminal vs Basal Breast Cancer",
       subtitle = "Differential Gene Expression Analysis",
       x = "log2(Fold Change) [Luminal/Basal]",
       y = "-log10(Adjusted p-value)",
       color = "Significant") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "darkred") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "darkblue") +
  theme(legend.position = "bottom")

ggsave("output/figures/breast_cancer_volcano_plot.png", volcano_plot, width=8, height=6)

# 3. Heatmap of top significant genes
if(nrow(sig_genes) > 0) {
  # Get normalized counts
  norm_counts <- counts(dds, normalized=TRUE)
  
  # Select top 20 most significant genes
  top_genes <- sig_genes %>%
    arrange(padj) %>%
    head(20) %>%
    pull(gene_id)
  
  top_norm_counts <- norm_counts[top_genes, ]
  
  # Create the heatmap
  png("output/figures/breast_cancer_heatmap_top_genes.png", width=800, height=1000)
  pheatmap(top_norm_counts,
           scale = "row",
           main = "Top 20 Significant Genes: Luminal vs Basal",
           show_rownames = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 8,
           annotation_col = data.frame(Condition = col_data$condition, 
                                       row.names = col_data$sample))
  dev.off()
}

cat("Breast cancer visualizations created and saved to 'output/figures/'.\n")