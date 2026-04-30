# ============================================================
# SCRIPT 00: Install All Required Packages
# Project: Microglial Heterogeneity in Alzheimer's Disease
# Author:  Eiman Meer
# Run this script ONCE before starting the project
# ============================================================

# ---- CRAN Packages ----
cran_packages <- c(
  "Seurat",        # Core scRNA-seq analysis
  "ggplot2",       # Plotting
  "patchwork",     # Combine plots
  "dplyr",         # Data manipulation
  "tidyr",         # Data tidying
  "readr",         # Read files
  "Matrix",        # Sparse matrix support
  "scales",        # Color scales
  "RColorBrewer",  # Color palettes
  "viridis",       # Viridis color palette
  "cowplot",       # Plot themes
  "ggrepel",       # Non-overlapping labels
  "pheatmap",      # Heatmaps
  "openxlsx",      # Export Excel tables
  "GEOquery"       # Download GEO datasets
)

installed <- rownames(installed.packages())
to_install <- cran_packages[!cran_packages %in% installed]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
} else {
  message("All CRAN packages already installed.")
}

# ---- Bioconductor Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "SingleR",          # Automated cell type annotation
  "celldex",          # Reference datasets for SingleR
  "clusterProfiler",  # Pathway enrichment analysis
  "org.Hs.eg.db",     # Human gene annotation database
  "enrichplot",       # Enrichment visualization
  "dittoSeq",         # Publication-quality scRNA-seq plots
  "limma",            # Differential expression support
  "edgeR",            # Pseudobulk DE analysis
  "MAST"              # Hurdle model for scRNA-seq DE
)

bioc_to_install <- bioc_packages[!bioc_packages %in% installed]
if (length(bioc_to_install) > 0) {
  BiocManager::install(bioc_to_install, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

message("
====================================================
 Package installation complete!
 Next step: Run 01_data_download_QC.R
====================================================
")
