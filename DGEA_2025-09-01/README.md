# Differential Gene Expression Analysis

## Project Description
A project demonstrating a standard RNA-seq differential gene expression analysis workflow using R and the DESeq2 package.

## Analysis Workflow
1.  **Data Loading:** Used the built-in `airway` dataset from Bioconductor.
2.  **Quality Control:** Pre-filtered genes with low counts.
3.  **Differential Expression:** Used DESeq2 to identify genes significantly altered by dexamethasone treatment.
4.  **Visualization:** Created standard plots (MA Plot, Volcano Plot, Heatmap).

## Key Findings
- Identified significant differentially expressed genes (adj. p-value < 0.05, |log2FC| > 1).
- The results are consistent with the known effects of dexamethasone.

## How to Run
1.  Clone this repository.
2.  Open the RProject file.
3.  Run `source("scripts/run_all.R")` to execute the entire analysis.
4.  Open `final_report.Rmd` and knit to generate the HTML report.
