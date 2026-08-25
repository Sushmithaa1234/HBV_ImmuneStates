# Comparative Single-cell Profiling of Intrahepatic Immune States across Clinical States of Chronic Hepatitis B

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Status](https://img.shields.io/badge/Status-In%20Progress-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.5.0-blue)

## 📋 Overview

A **fully reproducible single-cell RNA-seq analysis** characterizing **CD8 T-cell** and **macrophage transcriptional states** across 23 individuals representing five distinct clinical states of chronic hepatitis B infection. This project integrates differential expression, pathway enrichment, co-expression network analysis, and predicted cell-cell communication across clinically-defined HBV disease stages.

---

## 🔬 Clinical Context

This study investigates intrahepatic immune cell populations across the HBV disease spectrum:

| Clinical State | Abbr. | N | Description |
|---|---|---|---|
| Healthy/Seronegative | **NL** | 6 | Anti-HBc antibody negative |
| Immunotolerant | **IT** | 6 | High viral load, minimal inflammation |
| Immune Active | **IA** | 5 | Active viral replication, inflammation |
| Anti-HBe Seroconversion | **AR** | 3 | Transitional immune phase |
| Anti-HBc Seroconversion | **AC** | 3 | Resolved/chronic inactive state |

---

## 📊 Dataset

| Property | Value |
|---|---|
| **Source** | Gene Expression Omnibus (GEO) |
| **Accession** | [GSE182159](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182159) |
| **Technology** | 10x Genomics Chromium (v2/v3) |
| **Samples** | 23 liver biopsies (1 per donor) |
| **Post-QC Cells** | **106,592 cells** |
| **Genes Detected** | **18,925 genes** |

**Data Access:** Raw expression matrices are downloaded automatically via `TCGAbiolinks` during pipeline execution. No raw data are stored in this repository.

---

## 📈 Analysis Workflow

| Phase | Objective | Status |
|---|---|---|
| **0** | Pre-analysis dataset audit | ✅ **Complete** |
| **1** | Project infrastructure setup | ✅ **Complete** |
| **2** | Data loading + cell-level QC | ✅ **Complete** |
| **3** | Global atlas + cell-type annotation | 🟢 **In Progress** |
| **4** | CD8 T-cell subset analysis | ⏳ Planned |
| **5** | Macrophage/myeloid subset analysis | ⏳ Planned |
| **6** | Differential expression analysis | ⏳ Planned |
| **7** | Pathway enrichment (GO/KEGG/GSEA) | ⏳ Planned |
| **8** | Cell-cell communication (CellChat) | ⏳ Planned |
| **9** | Figure generation + integration | ⏳ Planned |

**Important:** Workflow and analysis decisions may evolve as data are processed. All updates will be reflected in this README and corresponding scripts.

---

## 📁 Repository Structure

```
HBV_ImmuneStates/
│
├── README.md                          # This file
├── .gitignore                         # Git exclusions
│
├── scripts/
│   ├── 00_setup.R                     # R environment + package setup
│   ├── 01_load_and_merge.R            # Load + merge 23 samples
│   ├── 02_qc_filtering.R              # Cell-level QC + metrics
│   ├── 03_metadata_annotation.R       # Global atlas + annotation
│   ├── 04_subset_cd8.R                # CD8 T-cell subset
│   ├── 05_subset_myeloid.R            # Macrophage/myeloid subset
│   ├── 06_de_analysis.R               # Differential expression
│   ├── 07_pathway_enrichment.R        # Functional enrichment
│   ├── 08_cellchat_analysis.R         # Cell-cell communication
│   └── run_full_pipeline.R            # Master control script
│
├── data/
│   ├── raw/                           # [Auto-downloaded at runtime]
│   └── processed/                     # Intermediate objects, checkpoints
│
└── results/
    ├── rds_objects/                   # Seurat checkpoints (.rds)
    ├── tables/                        # CSV/TSV result tables
    └── figures/                       # PNG/PDF figures
```

---

## 🛠️ Software Requirements

### **R Version**
```
R >= 4.5.0
```

### **Core Packages**

#### Bioconductor
```r
BiocManager::install(c(
  "TCGAbiolinks",        # GEO data retrieval
  "DESeq2",              # Differential expression
  "clusterProfiler",     # Gene enrichment
  "org.Hs.eg.db",        # Gene annotations
  "WGCNA",               # Co-expression networks
  "SingleR",             # Reference-based annotation
  "scDblFinder"          # Doublet detection
))
```

#### CRAN
```r
install.packages(c(
  "Seurat",              # scRNA-seq toolkit
  "ggplot2",             # Visualization
  "dplyr",               # Data manipulation
  "pheatmap",            # Heatmaps
  "survival",            # Survival analysis
  "scales",              # Plotting utilities
  "tidyr"                # Data tidying
))
```

**See `sessionInfo.txt` for exact package versions** (generated after reproducibility run).

---

## 🚀 Quick Start

### **Step 1: Clone Repository**
```bash
git clone https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA.git
cd HBV-ImmuneStates
```

### **Step 2: Set Up R Environment**
```r
# Open R in repository root
source("scripts/00_setup.R")
```

This will:
- ✅ Verify R version
- ✅ Install/load required packages
- ✅ Record package versions
- ✅ Confirm reproducibility environment

### **Step 3: Run Full Pipeline**

**Option A — Run all phases sequentially:**
```r
source("scripts/run_full_pipeline.R")
```

**Option B — Run individual phases:**
```r
source("scripts/01_load_and_merge.R")
source("scripts/02_qc_filtering.R")
source("scripts/03_metadata_annotation.R")
# ... etc
```

### **Step 4: Check Results**
All outputs are written to `results/` with subdirectories for figures, tables, and RDS objects.

---

## ⏱️ Runtime

**Expected:** ~4–6 hours (on 32GB RAM system)  
**Bottleneck:** Phase 3 annotation + projection (most compute-intensive)

---

## 📊 Results

### **Phase 3 (Current — In Progress)**

Upon completion of Phase 3, the following deliverables will be available:

- ✅ **Global UMAP:** All 106,592 cells in low-dimensional space
- ✅ **17 Major Cell Clusters:** Identified via graph-based clustering on leverage-score sketch
- ✅ **Cell-Type Annotations:** Multi-evidence framework (cluster markers + canonical genes + SingleR)
- ✅ **Transfer Confidence Scores:** Projection validation (99.79% sketch cluster retention)
- ✅ **Sample/Donor Metadata:** Linked to all cells
- ✅ **Canonical Marker Validation:** T cells, CD8 T cells, macrophages, dendritic cells, B cells, NK cells, etc.

### **Phases 4–9**

Results tables, figures, and biological interpretations will be added incrementally as each phase completes. See `results/` directory for current outputs.

---

## 🔑 Key Methodological Decisions

### **Leverage-Score Sketch for Scalability**
- Reference sketch: **21,756 representative cells** (up to 1,000 cells per sample)
- Smallest sample (66 cells) fully retained
- Full dataset projected into learned sketch-derived PCA/UMAP space
- **Benefit:** Maintains traceability while preserving computational tractability

### **Multi-Evidence Cell-Type Annotation**
Cell type assignments integrate:
1. **Cluster-specific markers** (FindAllMarkers)
2. **Canonical lineage gene panels** (T cells, macrophages, etc.)
3. **Module scoring** (curated marker programs)
4. **SingleR validation** (reference-based annotation)
5. **Manual biological interpretation**

This avoids reliance on automated classifiers alone.

### **Transparent QC Decisions**
- Phase 0 flags samples; **no automatic removal**
- Flagged samples/cells carried forward into Phase 2
- Final filtering decisions made **explicitly with justification**
- Sensitivity analyses included

### **Reproducibility-First**
- All thresholds documented and justified
- No post-hoc result removal for convenience
- Session information recorded for full transparency

---

## ⚠️ Important Interpretive Notes

### **This is Observational Single-Cell Analysis**

1. **CellChat Predictions:** Identify predicted ligand-receptor interactions, not experimentally proven functional signaling
2. **Differential Expression:** Identifies transcriptional associations, not causal mechanisms
3. **Network Topology:** Hub genes are highly connected within co-expression networks; centrality ≠ biological importance

### **Recommended Language**
Use: *associated with*, *suggests*, *characterized by*, *predicted*  
Avoid: *proves*, *causes*, *demonstrates mechanism*

---

## 📚 Citation

If you use this analysis or repository, please cite:

```bibtex
@article{Chandrasekar2026,
  author = {Sushmithaa Chandrasekar},
  title = {Single-cell immunotranscriptomic profiling of intrahepatic immune states in chronic hepatitis B},
  year = {2026},
  url = {https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA}
}
```

Or as text:

```
Sushmithaa Chandrasekar (2026). Single-cell immunotranscriptomic profiling 
of intrahepatic immune states in chronic hepatitis B. 
GitHub: https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA
```

---

## 📄 License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## 💬 Contact & Support

**Author:** Sushmithaa Chandrasekar  
**GitHub:** [@Sushmithaa1234](https://github.com/Sushmithaa1234)  
**GitHub Issues:** [Report bugs or request features](https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA/issues)

For methodology questions or analysis collaboration inquiries, please open an issue or contact directly.

---

## 🙏 Acknowledgments

- GEO repository and GSE182159 contributors for public data access
- Seurat, CellChat, and R/Bioconductor communities for open-source tools
- Inspiration from single-cell immunology literature

---

**Last Updated:** August 2026  
**Status:** Analysis in progress (Phase 3)
