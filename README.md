# Microglia-snRNAseq-Analysis
# Microglial Heterogeneity in Alzheimer's Disease — snRNA-seq Analysis

**Author:** Eiman Meer  
**Tools:** R, Seurat, clusterProfiler, SingleR, ggplot2  
**Dataset:** Mathys et al. 2019, *Nature* ([GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852))

## Project Overview

This project performs a focused single-nucleus RNA-sequencing (snRNA-seq) analysis of **microglial heterogeneity in Alzheimer's disease (AD)** using a landmark human postmortem prefrontal cortex dataset. The analysis characterizes distinct microglial transcriptional states, identifies disease-associated gene expression changes, and integrates findings with ncRNA target networks, connecting to my published review on ncRNA https://pubmed.ncbi.nlm.nih.gov/39258255/.

## Research Question

> *Do distinct microglial transcriptional states in human Alzheimer's disease prefrontal cortex correlate with disease-associated gene expression programs, and can ncRNA-regulated hub genes be identified within these states?*

## Dataset
| Property | Detail |
|---|---|
| Publication | Mathys et al. 2019, *Nature* |
| GEO Accession | GSE138852 |
| Tissue | Human prefrontal cortex (postmortem) |
| Technology | Single-nucleus RNA-seq (snRNA-seq) |
| Full dataset | ~80,000 nuclei (6 cell types) |
| This analysis | Microglia subset (~3,000–5,000 cells) |
| Conditions | AD patients vs age-matched controls |
## Methods
### Analysis Pipeline
Raw count matrix (GEO)
    ↓
Quality Control (Seurat)
    — Filter: nFeature 200–4000, percent.mt < 10%
    ↓
Normalization (SCTransform)
    — Regressed out: percent mitochondrial reads
    ↓
Dimensionality Reduction
    — PCA (50 PCs) → selected 20 PCs (ElbowPlot)
    — UMAP (seed = 42, for reproducibility)
    ↓
Clustering (Louvain, resolution 0.4)
    ↓
Cell Annotation
    — Manual annotation based on canonical microglial markers
    — Homeostatic: P2RY12, CX3CR1, TMEM119, HEXB
    — DAM: TREM2, APOE, LPL, SPP1, CLEC7A
    — Inflammatory: IL1B, TNF, CCL2
    — Proliferating: MKI67, TOP2A
    ↓
Differential Expression (MAST test, Seurat)
    — AD vs Control: global and per microglial state
    — Threshold: adj. p < 0.05, |Log2FC| > 0.25
    ↓
Pathway Enrichment (clusterProfiler)
    — GO Biological Process
    — KEGG Pathways
    ↓
ncRNA Integration
    — Cross-referenced DE genes with miR-124-3p,
      miR-155, and miR-146a predicted targets
      (miRDB score > 80, TargetScan validated)
## Repository Structure
sc_RNA_Project/
├── 00_install_packages.R       # Install all dependencies
├── 01_data_download_QC.R       # Download GSE138852 + QC filtering
├── 02_normalization_clustering.R  # SCTransform + PCA + UMAP + clustering
├── 03_annotation.R             # Marker gene identification + annotation
├── 04_differential_expression.R   # MAST-based DE: AD vs Control
├── 05_pathway_enrichment.R     # GO + KEGG enrichment analysis
├── 06_ncrna_integration.R      # miRNA target cross-referencing
├── 07_final_figures.R          # Publication-quality combined figure
├── data/raw/                   # Raw downloaded data (not tracked by git)
├── rds/                        # Intermediate Seurat objects
├── figures/                    # All output figures (PDF + PNG)
│   ├── 01_QC/
│   ├── 02_clustering/
│   ├── 03_annotation/
│   ├── 04_DE/
│   ├── 05_pathways/
│   ├── 06_ncrna/
│   └── final_panel/
└── outputs/                    # Tables, DE results, enrichment CSVs
## Dependencies
### R Packages
| Package | Version | Purpose |
|---|---|---|
| Seurat | ≥5.0 | Core scRNA-seq analysis |
| SCTransform | latest | Normalization |
| ggplot2 | ≥3.4 | Visualization |
| patchwork | latest | Plot composition |
| clusterProfiler | ≥4.0 | Pathway enrichment |
| org.Hs.eg.db | latest | Human gene annotation |
| SingleR | latest | Cell type annotation |
| MAST | latest | DE testing |
| enrichplot | latest | Enrichment visualization |
| dittoSeq | latest | scRNA-seq plotting |
| ggrepel | latest | Non-overlapping labels |
Install all: `source("00_install_packages.R")`
## Related Work

- **Published review:** Meer E. (2024). Role of Non-coding RNAs in Modulating Microglial Phenotype. *Global Medical Genetics*, 11:304–311. — The conceptual foundation for the ncRNA integration in Script 06.

- **miRNA target prediction pipeline:** [miRNA-Target-Predictor](https://github.com/EimanMeer/miRNA-Target-Predictor) — Python/Jupyter pipeline for miR-124-3p binding site prediction, directly informing the target lists used in this project.
## How to Run
1. Clone this repository
2. Open RStudio and set working directory
3. Run scripts in order: `00` → `01` → `02` → ... → `07`
4. Adjust file paths in each script if needed
5. Examine `figures/` and `outputs/` for results
> **Note:** Script 01 downloads ~500MB of data. Allow 20–30 minutes on first run.

eimanmeer720@gmail.com | [GitHub](https://github.com/EimanMeer)
