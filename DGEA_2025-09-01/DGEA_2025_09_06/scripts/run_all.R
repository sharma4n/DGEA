# Master Script to run the entire Pasilla analysis
source("scripts/01_install_packages.R")
source("scripts/02_data_loading.R")
source("scripts/03_differential_expression.R")
source("scripts/04_visualization.R")

print("All scripts executed successfully! Pasilla knockdown analysis complete.")