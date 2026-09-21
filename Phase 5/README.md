# Phase 5 — Myeloid State Analysis

## Overview

Phase 5 performs focused transcriptional state discovery and biological adjudication within the **Monocyte_Myeloid compartment** identified during Phase 3 of the HBV liver single-cell RNA-seq analysis.

The analysis asks:

> **What transcriptionally distinct myeloid-associated populations are present within the liver Monocyte_Myeloid compartment, and what biological programs characterize them?**

This phase is deliberately exploratory and descriptive. It does **not** infer clinical-state differences, disease progression, causality, or functional phenotypes.

Clinical-state inference is reserved for **Phase 6**, where donor identity is treated as the biological replicate.

---

## Dataset Context

The project uses liver-only single-cell RNA-seq data from:

**GEO accession:** `GSE182159`

The dataset contains liver samples from **23 individuals** representing five clinically defined HBV states:

| Clinical state | Description |
|---|---|
| AC | Acute carrier |
| AR | Active replication |
| IA | Immune active |
| IT | Immune tolerant |
| NL | Normal liver |

> The five groups are treated as clinically defined comparison groups. They are **not modeled as a progression trajectory**.

Phase 5 operates exclusively on the cells classified as:

```
Final_Cell_Type == "Monocyte_Myeloid"
```

from the Phase 3 broad immune atlas.

---

## Input

Phase 5 uses the final Phase 3 full-dataset object:

```
results/rds_objects/phase3_final_full_dataset.rds
```

**Phase 3 dataset contains:**

- 106,592 total cells
- 24,452 genes
- 23 donors
- 5 clinical states

**Phase 5 Monocyte_Myeloid compartment:**

- 2,489 cells
- 24,452 genes
- 23 donors
- 5 clinical states

All 23 donors are represented in the input compartment.

---

## Analytical Principles

Several principles were locked before biological interpretation.

### 1. Donors are the biological replicates

Individual cells are used for:

- transcriptional state discovery
- clustering
- marker identification
- program characterization

Individual donors are the biological replicates for:

- clinical-state inference
- differential abundance
- differential expression/state analysis

Therefore, Phase 5 does **not** treat thousands of cells from the same individual as thousands of independent biological replicates.

---

### 2. Processed expression values are preserved

The GEO supplementary matrices contain processed expression values corresponding to log-transformed CP10K expression.

Phase 5 therefore uses the existing processed expression data.

No attempt is made to reconstruct raw UMI counts.

In particular:

- no reverse transformation to pseudo-counts
- no fabricated count matrix
- no re-normalization from reconstructed counts
- no claim that stored expression values are raw counts

---

### 3. Biological interpretation is evidence-driven

Cluster identities are adjudicated using multiple sources of evidence:

- differential marker genes
- lineage-associated marker panels
- transcriptional program scores
- antigen-presentation programs
- inflammatory programs
- macrophage-associated programs
- donor representation
- cluster size
- potential lineage ambiguity

No biological identity is assigned solely because a particular marker appears among the top-ranked genes.

---

## Computational Workflow

Phase 5 follows the same analysis architecture established for the CD8 compartment in Phase 4.

### 1. Subset the Phase 3 myeloid compartment

The Phase 3 full dataset is subset to:

```
Final_Cell_Type == "Monocyte_Myeloid"
```

**Expected result:**

```
2,489 cells
24,452 genes
23 donors
5 clinical states
```

The source expression layers are preserved before joining the processed expression layers for downstream analysis.

---

### 2. Preserve the existing expression representation

The Phase 5 analysis uses the processed log-CP10K expression supplied with the GEO dataset.

The analysis does **not** call `NormalizeData()` and does not attempt to reconstruct raw counts.

---

### 3. Highly variable feature selection

A set of 2,000 highly variable genes is selected for dimensionality reduction.

---

### 4. Scaling and PCA

The selected features are scaled and principal component analysis is performed.

The analysis retains:

```
50 principal components
```

for downstream structure discovery.

---

### 5. Neighbor graph construction

A nearest-neighbor graph is constructed using:

```
PC1–PC20
k = 20
```

The Phase 5 graph names are:

```
Myeloid_nn
Myeloid_snn
```

---

### 6. Graph-based clustering

Clusters are identified using the shared-nearest-neighbor graph with:

```
Resolution = 0.4
Algorithm = Louvain
```

This produced **five computational clusters**.

Importantly, these are computational clusters first. They are not automatically interpreted as five distinct biological cell states.

---

### 7. UMAP visualization

A UMAP embedding is generated from the PCA representation to visualize the transcriptional structure of the myeloid compartment.

The UMAP is used as a visualization of transcriptional relationships and is **not** treated as an independent source of biological identity.

---

### 8. Marker identification

Cluster-level differential markers are identified using the Wilcoxon rank-sum test.

Marker lists are subsequently interpreted alongside:

- lineage panels
- transcriptional programs
- donor representation
- biological context

---

## Phase 5 Computational Clusters

The final clustering produced five populations:

| Cluster | Cells | Fraction |
|---:|---:|---:|
| 0 | 1,153 | 46.3% |
| 1 | 649 | 26.1% |
| 2 | 336 | 13.5% |
| 3 | 283 | 11.4% |
| 4 | 68 | 2.7% |

These cluster labels are computational identifiers and should not be confused with final biological annotations.

---

## Biological Adjudication

The final biological annotations were assigned after reviewing:

- top-ranked markers
- lineage-associated panels
- transcriptional program scores
- macrophage-associated signal
- antigen-presentation signal
- inflammatory signal
- donor representation
- evidence of non-myeloid populations

### Final annotations

| Cluster | Biological annotation | Confidence | Interpretation |
|---:|---|---|---|
| 0 | Inflammatory monocyte-associated population | High | Interpretable, but strongly donor concentrated |
| 1 | Activated/antigen-presenting monocyte-associated population | High | Broad donor representation |
| 2 | Antigen-presenting myeloid population | Moderate–High | Broad representation; macrophage identity not forced |
| 3 | FCGR3A-associated myeloid population | Moderate | Strong donor concentration and mixed granulocyte signal |
| 4 | Platelet-associated population | High | Non-myeloid-associated; entirely donor concentrated |

---

## Cluster 0 — Inflammatory Monocyte-Associated Population

Cluster 0 is characterized by a strong inflammatory monocyte-associated transcriptional program.

**Notable markers:**

`CCR2`, `S100A8`, `S100A12`, `FCN1`, `CD14`, `VCAN`, `NFE2`, `ALOX5AP`

The cluster also shows strong core myeloid and antigen-presentation programs.

**Biological annotation:**

> **Inflammatory monocyte-associated population**

### Important limitation

Cluster 0 is strongly donor concentrated:

```
1,086 / 1,153 cells (94.2%)
```

originate from the largest contributing donor.

Therefore, the existence of this transcriptional population can be described within the dataset, but its prevalence or relationship to clinical state should **not** be generalized from the cell count alone.

Donor-aware analysis is deferred to Phase 6.

---

## Cluster 1 — Activated/Antigen-Presenting Monocyte-Associated Population

Cluster 1 displays a coherent activated and antigen-presenting monocyte-associated profile.

**Notable markers:**

`IL1A`, `SERPINB2`, `C15orf48`, `GPR183`, `IL1R2`, `TRAF1`, `MIR155HG`, `THBS1`

The population also demonstrates strong HLA class II and antigen-presentation expression:

`HLA-DRA`, `HLA-DRB1`, `HLA-DPA1`, `HLA-DPB1`, `CD74`

**Biological annotation:**

> **Activated/antigen-presenting monocyte-associated population**

The population is represented across multiple donors, making the annotation less dependent on a single donor than clusters 0 and 3.

However, the annotation remains descriptive. It does not establish a causal or functional activation mechanism.

---

## Cluster 2 — Antigen-Presenting Myeloid Population

Cluster 2 shows strong core myeloid identity together with prominent antigen-presentation programs.

**Representative evidence:**

`LYZ`, `TYROBP`, `FCER1G`, `CTSS`, `AIF1`, `CST3`, `HLA-DRA`, `HLA-DRB1`, `HLA-DPA1`, `HLA-DPB1`, `CD74`

The population also shows some macrophage-associated expression. However, the overall macrophage-associated program is not strong enough to justify forcing a specific macrophage identity.

**Biological annotation:**

> **Antigen-presenting myeloid population**

This intentionally avoids unsupported labels such as:

- macrophage
- Kupffer cell
- disease-associated macrophage
- TREM2+ macrophage
- FOLR2+ macrophage

---

## Cluster 3 — FCGR3A-Associated Myeloid Population

Cluster 3 shows a strong FCGR3A-associated myeloid signature.

**Notable markers:**

`FCGR3A`, `CX3CR1`, `MS4A4A`, `LYZ`, `TYROBP`, `FCER1G`, `AIF1`, `CST3`

Strong HLA class II expression is also present.

However, the cluster contains some granulocyte-associated signal:

`FCGR3B`, `CSF3R`, `FPR1`

and is strongly donor concentrated:

```
259 / 283 cells (91.5%)
```

originate from one donor.

**Biological annotation:**

> **FCGR3A-associated myeloid population**

Rather than being forced into a neutrophil, dendritic-cell, or macrophage identity.

The donor concentration is an important limitation for downstream cohort-level interpretation.

---

## Cluster 4 — Platelet-Associated Population

Cluster 4 is qualitatively different from the other four populations.

**Strongest markers:**

`PF4`, `PPBP`, `TUBB1`, `ITGA2B`, `GP9`, `TREML1`, `MPIG6B`, `PF4V1`, `CLEC1B`, `CALD1`

This is a strong platelet-associated transcriptional signature.

The population consists of:

```
68 cells
```

and all 68 cells originate from a single donor.

**Biological annotation:**

> **Platelet-associated population**

This population is **not interpreted as a biological myeloid state**.

### Why it is not simply called a doublet

The available expression data are processed log-CP10K values rather than raw counts.

Therefore, this analysis does not have the raw-count information necessary to make a definitive doublet classification.

The conservative interpretation is:

> **Platelet-associated / non-myeloid-associated population**

rather than a definitive claim of doublet contamination.

---

## Macrophage Program Assessment

One of the important findings of Phase 5 is that a strong macrophage-associated transcriptional program was **not observed** across the five computational clusters at this resolution.

Macrophage-associated markers and programs were comparatively weak.

This means the analysis does **not** force the Monocyte_Myeloid compartment into conventional macrophage subtypes simply because macrophages were part of the original biological question.

This is an intentional methodological decision.

### Interpretation

The absence of a strong macrophage-associated program at this clustering resolution is itself a result of the analysis.

It does **not** demonstrate that macrophages are absent from the liver.

It means only that the current dataset representation and analytical resolution do not provide sufficiently strong evidence to assign the observed populations to specific macrophage states.

---

## Donor Representation

Donor representation was explicitly evaluated for every computational cluster.

| Cluster | Donors represented | Largest-donor fraction | Interpretation |
|---:|---:|---:|---|
| 0 | 16 | 94.2% | Strong donor concentration |
| 1 | 18 | 26.2% | Broad representation |
| 2 | 19 | 30.7% | Broad representation |
| 3 | 8 | 91.5% | Strong donor concentration |
| 4 | 1 | 100% | Extremely donor concentrated |

These observations are important because a transcriptional population can be biologically coherent while still being poorly suited to cohort-level inference if it is dominated by one individual.

Accordingly:

> **Cell abundance within a cluster is not interpreted as evidence of clinical-state enrichment in Phase 5.**

Clinical-state comparisons require donor-aware analysis and are reserved for Phase 6.

---

## What Phase 5 Can Conclude

Within the liver Monocyte_Myeloid compartment, the analysis identifies:

1. an inflammatory monocyte-associated population
2. an activated/antigen-presenting monocyte-associated population
3. an antigen-presenting myeloid population
4. an FCGR3A-associated myeloid population
5. a small platelet-associated population (not interpreted as a myeloid state)

The analysis also indicates that macrophage-associated transcriptional programs are relatively weak at this resolution.

These findings provide a structured set of candidate myeloid populations for downstream donor-aware analysis.

---

## What Phase 5 Cannot Conclude

Phase 5 does **not** establish:

- clinical-state-specific abundance differences
- clinical-state-specific expression differences
- disease progression
- causal relationships
- pathogenic or protective functions
- HBV-specificity of any transcriptional state
- macrophage polarization
- M1/M2 identity
- Kupffer-cell identity
- disease-associated macrophage identity
- functional cell-cell interactions
- experimental validation
- generalization beyond the analyzed cohort

These questions require additional evidence and, where appropriate, donor-aware statistical analysis.

---

## Relationship to Phase 6

Phase 5 is primarily a **state-discovery and biological-adjudication phase**.

The resulting populations provide candidate states for Phase 6.

Phase 6 will address questions such as:

- Do myeloid population abundances differ across the five clinical groups?
- Do transcriptional programs differ between clinical groups?
- Which genes show donor-aware differential expression/state?
- Are observed differences robust across individuals rather than driven by one donor?

These analyses will treat:

> **donor = biological replicate**

and will not use individual cells as independent biological replicates for clinical-state inference.

---

## Final Outputs

### RDS

**Final annotated Seurat object:**

```
results/rds_objects/phase5_final_myeloid_analysis_adjudicated.rds
```

The object contains the computational clustering together with the final biological annotations.

**Final annotation metadata:**

```
Myeloid_Biological_Annotation
Myeloid_Lineage_Category
Myeloid_Annotation_Confidence
Myeloid_Donor_Representation
Myeloid_Interpretation_Eligible
```

**Analysis-level metadata:**

```
Phase5_Analysis
Phase5_Expression
Phase5_Donor_Unit
Phase5_Inference
```

---

### Final Biological Annotation Table

```
results/tables/phase5_myeloid_final_biological_annotations.csv
```

This table contains:

- computational cluster
- cell count
- biological annotation
- lineage category
- confidence
- donor representation interpretation
- myeloid interpretation eligibility
- evidence summary
- interpretation boundary

---

### Additional outputs

| Output | Purpose |
|---|---|
| `phase5_myeloid_cluster_sizes.csv` | Cluster population sizes |
| `phase5_myeloid_all_positive_markers.csv` | All significant positive markers per cluster |
| `phase5_myeloid_top20_markers.csv` | Top 20 markers per cluster |
| `phase5_myeloid_top50_markers.csv` | Top 50 markers per cluster |
| `phase5_myeloid_cluster_program_scores.csv` | Biological program scores by cluster |
| `phase5_myeloid_marker_panel_cluster_means.csv` | Lineage-panel gene expression by cluster |
| `phase5_myeloid_cluster_donor_counts.csv` | Donor and GSM composition by cluster |
| `phase5_myeloid_cluster_donor_representation.csv` | Donor representation summary by cluster |
| `phase5_myeloid_donor_cluster_proportions.csv` | Cluster proportions within each donor |
| `phase5_myeloid_umap_clusters.png` | UMAP colored by cluster |
| `phase5_myeloid_umap_clinical_state.png` | UMAP colored by clinical state |
| `phase5_myeloid_marker_validation_dotplot.png` | DotPlot of marker genes by cluster |

---

## Reproducibility

The Phase 5 workflow is designed to be reproducible from the Phase 3 final object.

The analysis preserves:

- the original cell identities
- donor identities
- clinical-state metadata
- processed expression representation
- computational cluster assignments
- marker results
- biological adjudication

No cells were removed specifically because their biological identity was inconvenient or ambiguous.

Instead, populations with weaker or conflicting evidence are explicitly flagged and interpreted conservatively.

---

## Interpretation Philosophy

This phase deliberately follows a simple principle:

> **Do not force the biology to look cleaner than the data.**

A computational cluster is not automatically a biological cell state.

A marker is not automatically a cell identity.

A cell count is not automatically a biological replicate.

A transcriptional association is not automatically a mechanism.

And a predicted or inferred identity is not automatically a demonstrated function.

The goal of Phase 5 is therefore not to maximize the number of named cell types.

The goal is to produce a **defensible representation of the transcriptional structure actually supported by the dataset**, while preserving uncertainty for downstream analysis.

---

## Software

**Core analysis performed using:**

- R
- Seurat 5.5.1
- SeuratObject 5.4.0
- dplyr
- tidyr
- tibble
- stringr
- Matrix
- ggplot2

Exact package versions are recorded in the project environment/session information.

---

## Project Position

Phase 5 is part of a broader multi-phase analysis of intrahepatic immune transcriptional states in chronic HBV.

```
Phase 0 — Pre-analysis audit
        ↓
Phase 1 — Input/data understanding
        ↓
Phase 2 — QC object construction
        ↓
Phase 3 — Broad immune atlas
        ↓
Phase 4 — CD8 T-cell state analysis
        ↓
Phase 5 — Myeloid state analysis
        ↓
Phase 6 — Donor-aware differential abundance/state analysis
        ↓
Phase 7 — Pathway enrichment
        ↓
Phase 8 — Predicted cell-cell communication
        ↓
Phase 9 — Integrated interpretation
```

Phase 5 therefore provides the **myeloid state framework** used by subsequent donor-aware analyses.

---

## Status

### ✅ Phase 5 — COMPLETE

**Final dataset:**

```
2,489 myeloid-compartment cells
24,452 genes
23 donors
5 clinical states
5 computational clusters
5 cluster-level biological annotations
4 clusters eligible for interpretation as myeloid-associated states
1 platelet-associated/non-myeloid-associated population
```

The final annotations are intentionally conservative and are designed to support, rather than pre-empt, the donor-aware analyses performed in Phase 6.

---

**Phase 5 README updated:** September 2026  
**Analysis scope:** Myeloid compartment transcriptional state discovery and biological adjudication  
**Inference unit:** Donor for downstream clinical-state inference  
**Biological interpretation:** Exploratory and observational
