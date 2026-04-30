# ============================================================
# SCRIPT 01: Data Download & Quality Control
# Project: Microglial Heterogeneity in Alzheimer's Disease
# Dataset: Mathys et al. 2019, Nature (GSE138852)
#          ~80,000 single-nucleus RNA-seq profiles
#          Human prefrontal cortex — AD vs Control
# Author:  Eiman Meer
# ============================================================

# ---- 0. Setup: Load Libraries ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(readr)
library(patchwork)
library(GEOquery)

# ---- 1. Set Project Paths ----
base_dir   <- "E:/sc_RNA_Project"
data_dir   <- file.path(base_dir, "data", "raw")
output_dir <- file.path(base_dir, "outputs")
fig_dir    <- file.path(base_dir, "figures", "01_QC")

# Create all directories if they don't exist
for (d in c(data_dir, output_dir, fig_dir,
            file.path(base_dir, "figures", "02_clustering"),
            file.path(base_dir, "figures", "03_annotation"),
            file.path(base_dir, "figures", "04_DE"),
            file.path(base_dir, "figures", "05_pathways"),
            file.path(base_dir, "figures", "06_ncrna"),
            file.path(base_dir, "rds"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("Project directories created.\n")

# ---- 2. Download Dataset from GEO ----
# Dataset: Mathys et al. 2019, Nature
# GEO Accession: GSE138852
# This downloads two files:
#   - GSE138852_counts.csv.gz  (raw count matrix)
#   - GSE138852_meta.csv.gz    (cell metadata with cell type labels)
#
# NOTE: Files are ~500MB total. Allow 10-20 min depending on connection.
# If download fails, download manually from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852

cat("Downloading dataset from GEO (GSE138852)...\n")
cat("This may take 10-20 minutes. Please be patient.\n")

counts_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_counts.csv.gz"
meta_url   <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/GSE138852_meta.csv.gz"

counts_file <- file.path(data_dir, "GSE138852_counts.csv.gz")
meta_file   <- file.path(data_dir, "GSE138852_meta.csv.gz")

# Only download if files don't already exist
if (!file.exists(counts_file)) {
  tryCatch({
    download.file(counts_url, destfile = counts_file, mode = "wb", timeout = 600)
    cat("Count matrix downloaded successfully.\n")
  }, error = function(e) {
    cat("Auto-download failed. Please download manually:\n")
    cat("URL:", counts_url, "\n")
    cat("Save to:", counts_file, "\n")
  })
} else {
  cat("Count matrix already exists. Skipping download.\n")
}

if (!file.exists(meta_file)) {
  tryCatch({
    download.file(meta_url, destfile = meta_file, mode = "wb", timeout = 600)
    cat("Metadata downloaded successfully.\n")
  }, error = function(e) {
    cat("Auto-download failed. Please download manually:\n")
    cat("URL:", meta_url, "\n")
    cat("Save to:", meta_file, "\n")
  })
} else {
  cat("Metadata already exists. Skipping download.\n")
}

# ---- 3. Load Data ----
cat("\nLoading count matrix (this may take a few minutes)...\n")
counts <- read.csv(gzfile(counts_file), row.names = 1, check.names = FALSE)
counts <- as(as.matrix(counts), "dgCMatrix")   # Convert to sparse matrix (saves RAM)
cat("Count matrix loaded:", nrow(counts), "genes x", ncol(counts), "cells\n")

cat("Loading metadata...\n")
meta <- read.csv(gzfile(meta_file), row.names = 1)
cat("Metadata loaded:", nrow(meta), "cells\n")

# Preview metadata columns
cat("\nMetadata columns available:\n")
print(colnames(meta))
cat("\nCell types in dataset:\n")
print(table(meta$broad.cell.type))

# ---- 4. Create Seurat Object (Full Dataset) ----
cat("\nCreating full Seurat object...\n")
seurat_full <- CreateSeuratObject(
  counts   = counts,
  meta.data = meta,
  project  = "AD_snRNAseq",
  min.cells = 3,    # Keep genes expressed in at least 3 cells
  min.features = 200  # Keep cells with at least 200 genes
)

cat("Full Seurat object created:", ncol(seurat_full), "cells\n")

# ---- 5. Subset to Microglia Only ----
# We focus on microglia for targeted analysis
# The metadata column for broad cell type is 'broad.cell.type'
# Microglia label in this dataset: "Mic" 

cat("\nSubsetting to microglia...\n")

# Check exact label used
print(table(seurat_full$broad.cell.type))

# Subset — adjust label if your metadata uses a different name
microglia <- subset(seurat_full, subset = broad.cell.type == "Mic")
cat("Microglia subset:", ncol(microglia), "cells\n")

# Clean up full object to free RAM
rm(seurat_full, counts)
gc()

# ---- 6. Add Key Metadata Columns ----
# Ensure we have AD vs Control status clearly labeled
# In Mathys 2019: 'diagnosis' column contains "AD" or "Control" (or similar)
cat("\nDisease status distribution:\n")
print(table(microglia$diagnosis))

# Rename for clarity if needed (adjust column name to match your metadata)
microglia$disease_status <- ifelse(
  grepl("AD|Alzheimer", microglia$diagnosis, ignore.case = TRUE),
  "AD", "Control"
)
cat("\nDisease status after recoding:\n")
print(table(microglia$disease_status))

# ---- 7. Quality Control ----
cat("\nCalculating QC metrics...\n")

# Calculate percent mitochondrial reads
# Mitochondrial genes start with "MT-" in human
microglia[["percent.mt"]] <- PercentageFeatureSet(microglia, pattern = "^MT-")

# Calculate percent ribosomal reads
microglia[["percent.ribo"]] <- PercentageFeatureSet(microglia, pattern = "^RP[SL]")

# View QC summary
cat("\nQC Summary:\n")
summary(microglia@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])

# ---- 8. Visualize QC Metrics (Before Filtering) ----
cat("\nGenerating QC plots...\n")

p_vln_before <- VlnPlot(
  microglia,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0,
  group.by = "disease_status"
) +
  plot_annotation(title = "QC Metrics Before Filtering — Microglia",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

ggsave(file.path(fig_dir, "QC_violin_before_filtering.pdf"),
       p_vln_before, width = 14, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "QC_violin_before_filtering.png"),
       p_vln_before, width = 14, height = 6, dpi = 300)

# Scatter plot: nCount vs nFeature (detect doublets — unusually high both)
p_scatter <- FeatureScatter(microglia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                             group.by = "disease_status") +
  ggtitle("nCount vs nFeature (Microglia)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(fig_dir, "QC_scatter_nCount_nFeature.pdf"),
       p_scatter, width = 8, height = 6, dpi = 300)

cat("QC plots saved.\n")

# ---- 9. Apply QC Filters ----
# Thresholds based on Mathys 2019 and standard practice:
#   - Minimum 200 genes per cell (remove empty droplets)
#   - Maximum 4000 genes per cell (remove likely doublets)
#   - Maximum 10% mitochondrial reads (remove dying cells)

cat("\nApplying QC filters...\n")
cat("Cells before filtering:", ncol(microglia), "\n")

microglia <- subset(
  microglia,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 4000 &
           percent.mt < 10
)

cat("Cells after filtering:", ncol(microglia), "\n")
cat("Cells removed:", (ncol(microglia) - ncol(microglia)), "\n")

# ---- 10. Visualize QC After Filtering ----
p_vln_after <- VlnPlot(
  microglia,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0,
  group.by = "disease_status"
) +
  plot_annotation(title = "QC Metrics After Filtering — Microglia",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

ggsave(file.path(fig_dir, "QC_violin_after_filtering.pdf"),
       p_vln_after, width = 14, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "QC_violin_after_filtering.png"),
       p_vln_after, width = 14, height = 6, dpi = 300)

# ---- 11. Save Filtered Object ----
cat("\nSaving filtered Seurat object...\n")
saveRDS(microglia, file.path(base_dir, "rds", "01_microglia_filtered.rds"))

cat("
====================================================
 Script 01 Complete!
 Saved: rds/01_microglia_filtered.rds
 Figures saved to: figures/01_QC/
 Next step: Run 02_normalization_clustering.R
====================================================
")
