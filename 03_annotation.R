# ============================================================
# SCRIPT 03: Microglial Subtype Annotation
# Project: Microglial Heterogeneity in Alzheimer's Disease
# Author:  Eiman Meer
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)
library(pheatmap)
library(RColorBrewer)

# ---- 0. Paths ----
base_dir <- "E:/sc_RNA_Project"
fig_dir  <- file.path(base_dir, "figures", "03_annotation")
rds_dir  <- file.path(base_dir, "rds")

# ---- 1. Load Clustered Object ----
cat("Loading clustered microglia object...\n")
microglia <- readRDS(file.path(rds_dir, "02_microglia_clustered.rds"))
cat("Loaded:", ncol(microglia), "cells\n")

Idents(microglia) <- "SCT_snn_res.0.4"

# ---- 2. Find Cluster Marker Genes ----
# FindAllMarkers identifies genes that are elevated in each cluster
# compared to all other clusters combined
cat("\nFinding marker genes for each cluster...\n")
cat("(This may take 10-15 minutes)\n")

cluster_markers <- FindAllMarkers(
  microglia,
  only.pos  = TRUE,       # Only upregulated markers
  min.pct   = 0.25,       # Gene must be expressed in ≥25% of cells in the cluster
  logfc.threshold = 0.25, # Minimum log fold-change
  test.use  = "wilcox",   # Wilcoxon rank-sum test (default, robust)
  verbose   = FALSE
)

# Filter for significant markers
cluster_markers <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# Save full marker table
write.csv(cluster_markers,
          file.path(base_dir, "outputs", "cluster_markers_all.csv"),
          row.names = FALSE)

# Top 10 markers per cluster
top10 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10,
          file.path(base_dir, "outputs", "cluster_markers_top10.csv"),
          row.names = FALSE)

cat("Marker genes saved.\n")

# ---- 3. Heatmap of Top Marker Genes ----
cat("Generating cluster marker heatmap...\n")

top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

p_heatmap <- DoHeatmap(
  microglia,
  features = top5$gene,
  group.by = "SCT_snn_res.0.4",
  angle    = 45,
  size     = 3
) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Scaled\nExpression") +
  ggtitle("Top 5 Marker Genes per Cluster") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "Heatmap_cluster_markers.pdf"),
       p_heatmap, width = 14, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "Heatmap_cluster_markers.png"),
       p_heatmap, width = 14, height = 10, dpi = 300)

# ---- 4. Manual Annotation Based on Known Biology ----
# Microglial states in AD (based on Mathys 2019 + Keren-Shaul 2017 + Hammond 2019):
#
# Homeostatic microglia:
#   High: P2RY12, CX3CR1, TMEM119, HEXB, CSF1R, SALL1
#   These are the "resting" surveillance microglia
#
# Disease-Associated Microglia (DAM):
#   High: TREM2, APOE, LPL, SPP1, CLEC7A, CD9, CTSD
#   These are activated in response to amyloid plaques
#   DAM have two stages: DAM-1 (TREM2-independent) and DAM-2 (TREM2-dependent)
#
# Inflammatory microglia:
#   High: IL1B, TNF, CCL2, NLRP3, CD86
#   Pro-inflammatory, elevated in AD
#
# Proliferating microglia:
#   High: MKI67, TOP2A, PCNA
#   Dividing cells

# Based on marker genes, assign identities
# NOTE: You will need to adjust these mappings based on YOUR actual cluster markers
# Look at cluster_markers_top10.csv to guide your decisions

# Example annotation (adjust cluster numbers based on your results):
new_cluster_ids <- c(
  "0" = "Homeostatic",
  "1" = "Homeostatic",
  "2" = "DAM (TREM2-high)",
  "3" = "Transitional",
  "4" = "Inflammatory",
  "5" = "DAM (APOE-high)",
  "6" = "Proliferating",
  "7" = "Unassigned"
)

# Apply only clusters that exist
existing_clusters <- levels(Idents(microglia))
new_cluster_ids   <- new_cluster_ids[names(new_cluster_ids) %in% existing_clusters]

microglia <- RenameIdents(microglia, new_cluster_ids)
microglia$microglial_state <- Idents(microglia)

cat("\nMicroglial state distribution:\n")
print(table(microglia$microglial_state))

# ---- 5. UMAP with Annotated States ----
state_colors <- c(
  "Homeostatic"      = "#4DAF4A",
  "DAM (TREM2-high)" = "#E41A1C",
  "DAM (APOE-high)"  = "#FF7F00",
  "Transitional"     = "#984EA3",
  "Inflammatory"     = "#F781BF",
  "Proliferating"    = "#A65628",
  "Unassigned"       = "#999999"
)

p_umap_annotated <- DimPlot(
  microglia,
  reduction = "umap",
  group.by  = "microglial_state",
  label     = TRUE,
  label.size = 4,
  repel     = TRUE,
  cols      = state_colors
) +
  ggtitle("Annotated Microglial States — Alzheimer's Disease") +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right")

ggsave(file.path(fig_dir, "UMAP_annotated_states.pdf"),
       p_umap_annotated, width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_annotated_states.png"),
       p_umap_annotated, width = 10, height = 7, dpi = 300)

# Split by disease
p_split_annotated <- DimPlot(
  microglia,
  reduction = "umap",
  group.by  = "microglial_state",
  split.by  = "disease_status",
  cols      = state_colors,
  label     = FALSE
) +
  ggtitle("Microglial States: AD vs Control") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(fig_dir, "UMAP_annotated_split_disease.pdf"),
       p_split_annotated, width = 14, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "UMAP_annotated_split_disease.png"),
       p_split_annotated, width = 14, height = 6, dpi = 300)

# ---- 6. State Proportions by Disease Status ----
prop_state <- microglia@meta.data %>%
  group_by(microglial_state, disease_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(disease_status) %>%
  mutate(proportion = n / sum(n))

p_prop_state <- ggplot(prop_state,
                        aes(x = disease_status, y = proportion, fill = microglial_state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = state_colors, name = "Microglial State") +
  labs(x = "Condition", y = "Proportion of Microglia",
       title = "Microglial State Composition: AD vs Control") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text  = element_text(size = 12))

ggsave(file.path(fig_dir, "Barplot_state_proportions.pdf"),
       p_prop_state, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Barplot_state_proportions.png"),
       p_prop_state, width = 7, height = 6, dpi = 300)

# ---- 7. Violin Plots for Key State Markers ----
key_features <- c("P2RY12", "TREM2", "APOE", "SPP1", "IL1B", "MKI67")
key_features <- key_features[key_features %in% rownames(microglia)]

p_vln_states <- VlnPlot(
  microglia,
  features  = key_features,
  group.by  = "microglial_state",
  ncol      = 3,
  pt.size   = 0,
  cols      = state_colors
) +
  plot_annotation(title = "Key Marker Expression Across Microglial States",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")))

ggsave(file.path(fig_dir, "Violin_state_markers.pdf"),
       p_vln_states, width = 14, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "Violin_state_markers.png"),
       p_vln_states, width = 14, height = 9, dpi = 300)

# ---- 8. Save Annotated Object ----
saveRDS(microglia, file.path(rds_dir, "03_microglia_annotated.rds"))

cat("
====================================================
 Script 03 Complete!
 Saved: rds/03_microglia_annotated.rds
 Figures saved to: figures/03_annotation/
 Next step: Run 04_differential_expression.R
====================================================
")
