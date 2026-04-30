# ============================================================
# SCRIPT 02: Normalization, Dimensionality Reduction & Clustering
# Project: Microglial Heterogeneity in Alzheimer's Disease
# Author:  Eiman Meer
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)

# ---- 0. Paths ----
base_dir <- "E:/sc_RNA_Project"
fig_dir  <- file.path(base_dir, "figures", "02_clustering")
rds_dir  <- file.path(base_dir, "rds")

# ---- 1. Load Filtered Object ----
cat("Loading filtered microglia object...\n")
microglia <- readRDS(file.path(rds_dir, "01_microglia_filtered.rds"))
cat("Loaded:", ncol(microglia), "cells,", nrow(microglia), "genes\n")

# ---- 2. Normalization (SCTransform) ----
# SCTransform is superior to basic log-normalization
# It accounts for sequencing depth differences between cells
# vars.to.regress: removes technical variation from mitochondrial reads
cat("\nRunning SCTransform normalization...\n")
cat("(This may take 5-10 minutes)\n")

microglia <- SCTransform(
  microglia,
  vars.to.regress = "percent.mt",   # Regress out mitochondrial variation
  verbose = TRUE,
  return.only.var.genes = FALSE
)

cat("SCTransform complete.\n")

# ---- 3. PCA ----
cat("\nRunning PCA...\n")
microglia <- RunPCA(
  microglia,
  npcs = 50,           # Compute 50 PCs
  verbose = FALSE
)

# Elbow plot to determine how many PCs to use
p_elbow <- ElbowPlot(microglia, ndims = 50) +
  ggtitle("Elbow Plot — Selecting Number of PCs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  annotate("text", x = 21, y = max(ElbowPlot(microglia, ndims=50)$data$stdev)*0.9,
           label = "Selected cutoff", color = "red", hjust = 0, size = 3.5)

ggsave(file.path(fig_dir, "PCA_elbow_plot.pdf"), p_elbow, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "PCA_elbow_plot.png"), p_elbow, width = 8, height = 5, dpi = 300)
cat("Elbow plot saved — examine it to confirm PC cutoff (we use 20)\n")

# PCA plot colored by disease status
p_pca <- DimPlot(microglia, reduction = "pca", group.by = "disease_status",
                  cols = c("AD" = "#E63946", "Control" = "#457B9D")) +
  ggtitle("PCA — AD vs Control Microglia") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "PCA_disease_status.pdf"), p_pca, width = 8, height = 6, dpi = 300)

# Heatmap of top genes driving PC1 and PC2
pdf(file.path(fig_dir, "PCA_heatmap_PC1_PC2.pdf"), width = 12, height = 8)
DimHeatmap(microglia, dims = 1:2, cells = 500, balanced = TRUE)
dev.off()

# ---- 4. Nearest Neighbor Graph ----
# Using first 20 PCs (adjust based on your elbow plot)
n_pcs <- 20

cat("\nComputing nearest neighbor graph (", n_pcs, "PCs)...\n")
microglia <- FindNeighbors(
  microglia,
  dims = 1:n_pcs,
  verbose = FALSE
)

# ---- 5. Clustering at Multiple Resolutions ----
# Resolution controls granularity: lower = fewer, larger clusters
# We test multiple resolutions and choose the most biologically meaningful
cat("Testing clustering resolutions: 0.2, 0.4, 0.6, 0.8...\n")

microglia <- FindClusters(
  microglia,
  resolution = c(0.2, 0.4, 0.6, 0.8),
  verbose = FALSE
)

# ---- 6. UMAP Embedding ----
cat("Running UMAP...\n")
microglia <- RunUMAP(
  microglia,
  dims = 1:n_pcs,
  seed.use = 42,     # Set seed for reproducibility
  verbose = FALSE
)

# ---- 7. UMAP Plots at Different Resolutions ----
cat("Generating UMAP plots...\n")

# Choose a working resolution for downstream analysis
# Resolution 0.4 typically gives ~5-8 clusters for microglia — biologically meaningful
Idents(microglia) <- "SCT_snn_res.0.4"

# Set a nice color palette
cluster_colors <- brewer.pal(n = max(as.integer(levels(Idents(microglia)))) + 1, "Set2")

# Main UMAP by cluster
p_umap_clusters <- DimPlot(
  microglia,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  cols = cluster_colors
) +
  ggtitle("Microglial Clusters — Resolution 0.4") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

ggsave(file.path(fig_dir, "UMAP_clusters_res0.4.pdf"),
       p_umap_clusters, width = 9, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_clusters_res0.4.png"),
       p_umap_clusters, width = 9, height = 7, dpi = 300)

# UMAP colored by disease status
p_umap_disease <- DimPlot(
  microglia,
  reduction = "umap",
  group.by  = "disease_status",
  cols = c("AD" = "#E63946", "Control" = "#457B9D"),
  pt.size = 0.5
) +
  ggtitle("UMAP — AD vs Control Microglia") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_disease_status.pdf"),
       p_umap_disease, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_disease_status.png"),
       p_umap_disease, width = 8, height = 6, dpi = 300)

# Side-by-side comparison
p_combined <- p_umap_clusters | p_umap_disease
ggsave(file.path(fig_dir, "UMAP_combined.pdf"),
       p_combined, width = 16, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_combined.png"),
       p_combined, width = 16, height = 6, dpi = 300)

# Split by disease status
p_split <- DimPlot(
  microglia,
  reduction = "umap",
  split.by  = "disease_status",
  label = TRUE,
  label.size = 4,
  cols = cluster_colors
) +
  ggtitle("Microglial Clusters Split by Disease Status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_split_by_disease.pdf"),
       p_split, width = 14, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_split_by_disease.png"),
       p_split, width = 14, height = 6, dpi = 300)

# ---- 8. Known Microglial Marker Genes ----
# Homeostatic microglia markers: P2RY12, CX3CR1, TMEM119, CSF1R, HEXB
# Disease-associated microglia (DAM) markers: TREM2, APOE, LPL, CD9, SPP1
# Activation markers: CD68, AIF1 (IBA1)

homeostatic_markers <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "HEXB", "SALL1")
dam_markers         <- c("TREM2", "APOE", "LPL", "CD9", "SPP1", "CLEC7A", "CTSD")
activation_markers  <- c("CD68", "AIF1", "MHC2", "IL1B", "TNF", "CCL2")

# Feature plots for key markers
p_homeostatic <- FeaturePlot(
  microglia,
  features  = homeostatic_markers[homeostatic_markers %in% rownames(microglia)],
  reduction = "umap",
  ncol = 3,
  order = TRUE,
  cols = c("lightgrey", "#2166AC")
) +
  plot_annotation(title = "Homeostatic Microglial Markers",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

ggsave(file.path(fig_dir, "FeaturePlot_homeostatic_markers.pdf"),
       p_homeostatic, width = 14, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "FeaturePlot_homeostatic_markers.png"),
       p_homeostatic, width = 14, height = 9, dpi = 300)

p_dam <- FeaturePlot(
  microglia,
  features  = dam_markers[dam_markers %in% rownames(microglia)],
  reduction = "umap",
  ncol = 3,
  order = TRUE,
  cols = c("lightgrey", "#B2182B")
) +
  plot_annotation(title = "Disease-Associated Microglia (DAM) Markers",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

ggsave(file.path(fig_dir, "FeaturePlot_DAM_markers.pdf"),
       p_dam, width = 14, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "FeaturePlot_DAM_markers.png"),
       p_dam, width = 14, height = 9, dpi = 300)

# Dot plot of all key markers across clusters
all_key_markers <- c(homeostatic_markers, dam_markers, activation_markers)
all_key_markers <- all_key_markers[all_key_markers %in% rownames(microglia)]

p_dot <- DotPlot(
  microglia,
  features = all_key_markers,
  group.by = "SCT_snn_res.0.4",
  cols = c("lightgrey", "#B2182B")
) +
  coord_flip() +
  ggtitle("Key Microglial Markers Across Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(fig_dir, "DotPlot_key_markers.pdf"),
       p_dot, width = 10, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "DotPlot_key_markers.png"),
       p_dot, width = 10, height = 10, dpi = 300)

# ---- 9. Cell Proportions per Cluster ----
prop_table <- microglia@meta.data %>%
  group_by(SCT_snn_res.0.4, disease_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SCT_snn_res.0.4) %>%
  mutate(proportion = n / sum(n))

p_prop <- ggplot(prop_table, aes(x = SCT_snn_res.0.4, y = proportion, fill = disease_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("AD" = "#E63946", "Control" = "#457B9D"),
                    name = "Condition") +
  labs(x = "Cluster", y = "Proportion of Cells",
       title = "Cell Composition per Cluster by Disease Status") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 11))

ggsave(file.path(fig_dir, "Barplot_cell_proportions.pdf"),
       p_prop, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Barplot_cell_proportions.png"),
       p_prop, width = 9, height = 6, dpi = 300)

# ---- 10. Save ----
saveRDS(microglia, file.path(rds_dir, "02_microglia_clustered.rds"))

cat("
====================================================
 Script 02 Complete!
 Saved: rds/02_microglia_clustered.rds
 Figures saved to: figures/02_clustering/
 Next step: Run 03_annotation.R
====================================================
")
