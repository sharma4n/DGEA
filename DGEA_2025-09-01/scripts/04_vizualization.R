# Load libraries and data
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
load("data/processed/loaded_data.RData")
read_csv("output/tables/all_genes_results.csv") # Or recreate from dds if needed

# 1. MA Plot (Visualizes log2FC vs mean expression)
plot_ma <- function(){
  plotMA(results(dds), main="MA Plot", ylim=c(-5,5))
}
png("output/figures/ma_plot.png", width=800, height=600)
plot_ma()
dev.off()

# 2. Volcano Plot (A classic in genomics)
volcano_data <- read_csv("output/tables/all_genes_results.csv")

volcano_plot <- volcano_data %>%
  mutate(significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed")

ggsave("output/figures/volcano_plot.png", volcano_plot, width=8, height=6)

# 3. Heatmap of top significant genes
# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Select top 20 most significant genes
top_genes <- significant_genes %>%
  arrange(padj) %>%
  head(20) %>%
  pull(gene_id)

top_norm_counts <- norm_counts[top_genes, ]

# Create the heatmap
png("output/figures/heatmap_top_genes.png", width=800, height=1000)
pheatmap(top_norm_counts,
         scale = "row", # Scale by row (gene) to see relative expression
         main = "Expression of Top 20 Significant Genes",
         show_rownames = TRUE,
         cluster_cols = TRUE)
dev.off()

cat("Visualizations created and saved to 'output/figures/'.\n")
