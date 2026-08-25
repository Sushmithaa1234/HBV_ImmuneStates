# Comparative Single-cell Profiling of Intrahepatic Immune States across Clinical States of Chronic Hepatitis B

Overview

A reproducible single-cell RNA-seq analysis characterizing CD8 T-cell and macrophage transcriptional states across 23 individuals representing five distinct clinical states of chronic hepatitis B infection.

Clinical Context

This study investigates intrahepatic immune populations across HBV disease stages:

NL — Healthy/Seronegative (n=6)
IT — Immunotolerant (n=6)
IA — Immune Active (n=5)
AR — Anti-HBe Seroconversion (n=3)
AC — Anti-HBc Seroconversion (n=3)
Dataset

Source: Gene Expression Omnibus (GEO)
Accession: GSE182159
Platform: 10x Genomics Chromium
Samples: 23 liver biopsies (one per donor)
Post-QC Cells: 106,592 cells
Genes: 18,925

Raw data are downloaded automatically from GEO via TCGAbiolinks. No raw matrices are included in this repository.

Analysis Workflow
Phase	Goal	Status
0	Pre-analysis audit	✅ Complete
1	Project infrastructure	✅ Complete
2	Data loading + QC	✅ Complete
3	Global atlas + cell-type annotation	🟢 In Progress
4	CD8 T-cell characterization	⏳ Planned
5	Macrophage/myeloid characterization	⏳ Planned
6	Differential expression analysis	⏳ Planned
7	Pathway enrichment	⏳ Planned
8	Cell-cell communication (CellChat)	⏳ Planned
9	Figure generation + integration	⏳ Planned

Note: Workflow and specific analyses may evolve as data are processed. Updates will be reflected in this README and individual scripts.

Repository Structure
HBV_ImmuneStates/
├── README.md
├── .gitignore
│
├── scripts/
│   ├── 00_setup.R                     # Environment + packages
│   ├── 01_load_and_merge.R            # Load 23 samples, merge
│   ├── 02_qc_filtering.R              # Cell QC
│   ├── 03_metadata_annotation.R       # Global atlas + annotation
│   ├── 04_subset_cd8.R                # CD8 subset
│   ├── 05_subset_myeloid.R            # Myeloid subset
│   ├── 06_de_analysis.R               # Differential expression
│   ├── 07_pathway_enrichment.R        # GO, KEGG, GSEA
│   ├── 08_cellchat_analysis.R         # Cell-cell communication
│   └── run_full_pipeline.R            # Master script
│
├── data/
│   ├── raw/                           # [Downloaded at runtime]
│   └── processed/                     # Checkpoints, intermediate objects
│
└── results/
    ├── rds_objects/                   # Seurat checkpoints
    ├── tables/                        # CSV/TSV tables
    └── figures/                       # PNG/PDF figures
Software Requirements

R: ≥ 4.2.0

Key Packages:

r
# Bioconductor
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("WGCNA")
BiocManager::install("SingleR")

# CRAN
install.packages("Seurat")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("survival")

See sessionInfo.txt for exact versions after reproducibility run.

Quick Start
1. Clone repository
bash
git clone https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA.git
cd HBV-ImmuneStates
2. Set up environment
r
source("scripts/00_setup.R")
3. Run pipeline
r
# Full pipeline
source("scripts/run_full_pipeline.R")

# Or run individual phases
source("scripts/01_load_and_merge.R")
source("scripts/02_qc_filtering.R")
# ... etc

Runtime: ~4–6 hours (32GB RAM recommended)

Results

Results from completed phases will be added here as analysis progresses.

Phase 3 (Current)
Global UMAP of all 106,592 cells
17 major cell clusters
Cell-type annotations with confidence scores
[Full results pending completion]
Phase 4–9

Results tables, figures, and interpretations will be added as each phase completes.

Key Analysis Decisions
Leverage-score sketch: 21,756 representative cells used to establish global structure; full dataset projected into learned space for scalability and traceability.
Multi-evidence annotation: Cell types assigned using cluster markers + canonical lineage genes + reference validation (SingleR) + module scoring, not automated classifiers alone.
No blindside removals: Flagged samples/cells from Phase 0 carried forward; final filtering decisions made explicitly with justification.
Reproducibility first: All thresholds documented; sensitivity analyses included.
Important Notes

This is observational scRNA-seq analysis:

Findings describe transcriptional associations, not causal mechanisms.
CellChat predicts communication; experimental validation required.
Results should use language like "associated with," "suggests," "characterized by" rather than "proves" or "causes."
Citation
Sushmithaa Chandrasekar.
Single-cell immunotranscriptomic profiling of intrahepatic immune states in chronic hepatitis B.
[Preprint in preparation]
GitHub: https://github.com/Sushmithaa1234/HBV-ImmuneStates-scRNA
License

MIT License — See LICENSE for details.

Contact

GitHub: @Sushmithaa1234

For questions or issues, please open a GitHub issue.
