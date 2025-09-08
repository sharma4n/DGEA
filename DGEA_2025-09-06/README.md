Differential Gene Expression Analysis

Project Description
A project demonstrating a standard RNA-seq differential gene expression analysis workflow using R and the DESeq2 package.

Analysis Workflow

01. Data Loading: Used synthetic Breast cancer data built with rnbinom function.

02. Quality Control: Pre-filtered genes with low counts.

03. Differential Expression: Used DESeq2 to identify genes that are significantly altered.

04. Visualization: Created standard plots (MA Plot, Volcano Plot, Heatmap).

Key Findings
Identified significantly differentially expressed genes (adj. p-value < 0.05, |log2FC| > 1).

How to Run
* Clone this repository.
* Open the RProject file.
* Run source("scripts/run_all.R") to execute the entire analysis.
* Open final_report.Rmd and knit to generate the HTML report. ####Working on it.
