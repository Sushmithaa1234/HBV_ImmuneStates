# PHASE 3 — GLOBAL CELL ATLAS CONSTRUCTION & FULL-DATASET LABEL TRANSFER

**Status:** ✅ COMPLETE — LOCKED
**Analysis date:** September 2026
**Dataset:** GEO GSE182159 — liver samples from 23 individuals with chronic HBV infection
**Input representation:** GEO processed log-CP10K expression values
**Software:** R / Seurat 5.5.1 / SeuratObject 5.4.0

---

## Overview

Phase 3 constructs a global liver immune-cell atlas from a representative **20,000-cell sketch** of the complete dataset, performs multi-evidence cell-type annotation of the resulting clusters, and transfers those annotations to all **106,592 cells** in the full dataset.

The analysis was designed for a processed-expression dataset in which the available matrices contain **log-counts-per-ten-thousand (log-CP10K)** values rather than raw UMI counts.

Accordingly, Phase 3 does **not**:

* reconstruct raw counts;
* apply a second normalization step;
* perform conventional count-based QC;
* perform additional cell filtering;
* perform batch correction;
* remove cells based on label-transfer confidence.

The final workflow is:

```text
106,592 liver cells
        │
        ▼
Existing GEO processed log-CP10K expression
        │
        ▼
Sample-aware 20,000-cell sketch
        │
        ▼
Global PCA + clustering
        │
        ▼
Multi-evidence cluster annotation
        │
        ▼
Locked 20,000-cell reference atlas
        │
        ▼
Validated PCA projection of all 106,592 cells
        │
        ▼
KNN label transfer (PC1–20, k = 30)
        │
        ▼
Transfer confidence + targeted QC
        │
        ▼
New full-dataset UMAP
        │
        ▼
Final annotated dataset
```

---

# 1. Dataset and Input Representation

The analysis uses the **liver-only portion of GEO GSE182159**.

### Dataset structure

| Property         |                         Value |
| ---------------- | ----------------------------: |
| Liver samples    |                        **23** |
| Donors           |                        **23** |
| Clinical states  |                         **5** |
| Genes            |                    **24,452** |
| Cells            |                   **106,592** |
| RNA data layers  | **23 sample-specific layers** |
| Raw count layers |                         **0** |

Clinical states:

* AC
* AR
* IA
* IT
* NL

The 23 sample-specific expression matrices were retained throughout the analysis.

### Expression representation

The GEO supplementary matrices contain processed expression values corresponding to **log-CP10K** values.

Therefore:

> **The values used in Phase 3 are not raw UMI counts.**

This distinction is critical. Count-based quantities such as raw `nCount_RNA`, conventional count-depth normalization, and count-based doublet assumptions are not reconstructed from the processed matrices.

No artificial raw counts were generated.

---

# 2. SCRIPT 01 — GLOBAL SKETCH ATLAS CONSTRUCTION

## Objective

Construct a computationally tractable but representative **20,000-cell global atlas** while preserving representation from all 23 samples.

The complete dataset contains 106,592 cells, so the sketch represents approximately **18.8%** of all cells.

## Sampling strategy

The sketch was constructed using a **sample-aware allocation strategy**:

1. Every sample was required to contribute cells.
2. Approximately proportional allocation was used across samples.
3. A minimum representation safeguard was applied.
4. Within each sample, diversity/leverage-informed sampling was performed without replacement.
5. The final sketch contained exactly **20,000 cells**.

This avoids allowing the largest samples to completely dominate the reference atlas while also avoiding artificial equalization of clinical states.

### Sketch validation

* Reference cells: **20,000**
* Samples represented: **23/23**
* Clinical states represented: **5/5**
* Reference genes: **24,452**
* No cells were removed from the full dataset as part of sketch construction.

---

# 3. Highly Variable Features

Because the merged Seurat object contains 23 sample-specific expression layers, conventional merged-object HVG calculation was not used directly.

Instead, variable features were calculated **within each sample layer** using the existing processed expression values.

The workflow was:

```text
sample-specific processed expression
        ↓
per-sample HVG identification
        ↓
HVG frequency/rank aggregation
        ↓
consensus HVG ranking
        ↓
2,000 consensus variable features
```

No additional normalization was performed for this step.

The final atlas PCA contained loadings for **1,999 features**, because one of the selected variable features did not have a usable PCA loading.

The projection workflow therefore uses the **actual 1,999 PCA-loading features**, rather than assuming that all 2,000 selected HVGs contributed to PCA.

---

# 4. PCA and Global Atlas Construction

PCA was calculated on the 20,000-cell sketch.

### PCA

| Parameter                   |        Value |
| --------------------------- | -----------: |
| Reference cells             |       20,000 |
| PCA dimensions              |       **50** |
| PCA loading features        |    **1,999** |
| Downstream graph dimensions | **PC1–PC20** |

The 50-dimensional PCA representation was retained for projection.

The first 20 PCs were used for:

* nearest-neighbor graph construction;
* clustering;
* label transfer;
* full-dataset UMAP.

No additional batch correction was applied.

---

# 5. Global Clustering

The sketch produced **20 global clusters**.

The cluster identities were determined from the complete multi-evidence annotation workflow described below.

Final atlas cluster mapping:

| Atlas cluster | Final annotation    | Confidence |
| ------------: | ------------------- | ---------- |
|             0 | CD8_T               | High       |
|             1 | NK                  | High       |
|             2 | CD4_T               | High       |
|             3 | MAIT                | High       |
|             4 | CD8_T               | High       |
|             5 | NK                  | High       |
|             6 | MAIT                | Moderate   |
|             7 | CD8_T               | Moderate   |
|             8 | CD8_T               | High       |
|             9 | B_Cell              | High       |
|            10 | Unresolved_Lymphoid | Low        |
|            11 | Treg                | High       |
|            12 | NK                  | High       |
|            13 | NK                  | Moderate   |
|            14 | NK                  | High       |
|            15 | Monocyte_Myeloid    | High       |
|            16 | Monocyte_Myeloid    | High       |
|            17 | Dendritic_Cell      | High       |
|            18 | Plasma_Cell         | High       |
|            19 | pDC                 | High       |

### Important annotation principle

The global atlas deliberately uses **broad, evidence-supported labels**.

In particular:

* `Monocyte_Myeloid` is intentionally broad.
* Macrophage-associated markers were retained as a diagnostic program.
* Detailed macrophage-state resolution is deferred to **Phase 5**.
* Cluster 10 remains `Unresolved_Lymphoid` rather than being forcibly assigned to a more specific lineage because its evidence was conflicting and its cells were highly sample-concentrated.

---

# 6. SCRIPT 02 — CLUSTER MARKERS AND CELL-TYPE ANNOTATION

## Objective

Assign biological identities to the 20 global clusters using multiple independent evidence streams rather than relying on a single automated classifier.

Four principal evidence streams were used:

1. Cluster marker discovery
2. Canonical lineage programs
3. Module scores
4. Reference-based SingleR annotation

These were followed by manual evidence adjudication.

---

# 7. Evidence Stream 1 — Cluster Marker Discovery

Cluster markers were identified using Seurat's `FindAllMarkers()` after joining the sketch layers for this cluster-level analysis.

Parameters:

```text
test: Wilcoxon rank-sum
only.pos: TRUE
min.pct: 0.25
logfc.threshold: 0.25
```

The resulting markers were used as **descriptive cluster markers**, not as independent proof of biological identity.

Marker statistics were retained for downstream inspection.

---

# 8. Evidence Stream 2 — Canonical Immune Marker Programs

Thirteen curated marker programs were evaluated:

| Program          |
| ---------------- |
| T_Cell           |
| CD8_T            |
| CD4_T            |
| Treg             |
| NK               |
| MAIT             |
| GammaDelta_T     |
| B_Cell           |
| Plasma_Cell      |
| Monocyte_Myeloid |
| Macrophage       |
| Dendritic_Cell   |
| pDC              |

These programs were used to assess lineage concordance at the cluster level.

The `Macrophage` program was retained as a **diagnostic program**, even though the global atlas uses the broader `Monocyte_Myeloid` label.

---

# 9. Evidence Stream 3 — Module Scores

Module scores were calculated for the curated marker programs.

The scores provide an additional quantitative comparison between competing lineage programs.

They were treated as supporting evidence rather than as an automatic labeling rule.

This was particularly important for closely related lymphoid populations such as:

* CD4 T cells
* CD8 T cells
* MAIT cells
* NK cells
* Treg cells
* γδ T cells

---

# 10. Evidence Stream 4 — SingleR

Cluster-level reference annotation was performed using:

**SingleR + Monaco Immune reference**

The SingleR output was used as an independent reference-based evidence stream.

SingleR predictions were not accepted blindly. Discordance between SingleR and marker/module evidence was explicitly inspected during manual adjudication.

---

# 11. Manual Evidence Adjudication

Final labels were assigned by integrating:

```text
Cluster markers
      +
Canonical marker programs
      +
Module scores
      +
SingleR
      ↓
Manual evidence adjudication
      ↓
Final cluster label
```

### Final annotation principles

* Concordant evidence → high-confidence assignment.
* Strong but incomplete evidence → moderate-confidence assignment.
* Conflicting evidence → conservative label or unresolved status.
* No cluster was forced into a biologically specific identity solely because an automated classifier produced a label.

### Final atlas categories

The final global atlas contains:

* CD8_T
* NK
* CD4_T
* MAIT
* B_Cell
* Treg
* Monocyte_Myeloid
* Dendritic_Cell
* Plasma_Cell
* pDC
* Unresolved_Lymphoid

`Unresolved_Lymphoid` is an intentional category rather than a failed annotation.

---

# 12. SCRIPT 03 — FULL-DATASET PCA PROJECTION

## Objective

Transfer the validated sketch-atlas representation to **all 106,592 cells** without rerunning full-dataset clustering.

The full dataset was projected into the PCA coordinate system learned from the 20,000-cell reference atlas.

### Full dataset

| Property            |        Value |
| ------------------- | -----------: |
| Cells               |  **106,592** |
| Genes               |   **24,452** |
| Reference cells     |   **20,000** |
| PCA dimensions      |       **50** |
| Transfer dimensions | **PC1–PC20** |
| KNN k               |       **30** |
| Cells removed       |        **0** |

---

# 13. Projection Method

The full dataset was projected using the actual PCA parameters learned from the atlas:

1. Atlas PCA loading features were extracted.
2. Atlas feature means and variances were calculated.
3. Full-dataset expression was extracted from each sample-specific layer.
4. Query values were scaled using the **atlas reference scaling parameters**.
5. The scaled expression was multiplied by the locked PCA loading matrix.

Conceptually:

```text
Full-data expression
        ↓
Atlas reference scaling
        ↓
Locked atlas PCA loadings
        ↓
50-dimensional PCA coordinates
```

No new normalization was performed.

No new PCA was fitted to the full dataset.

No batch correction was introduced during projection.

---

# 14. Projection Validation

The projection implementation was independently validated by projecting the atlas itself through the same mathematical procedure.

The resulting coordinates reproduced the stored atlas PCA coordinates essentially exactly:

* PC correlations: **1.000**
* Sign-aligned correlations: **1.000**
* RMSE: approximately **10⁻¹⁴**
* All 50 PCs successfully reproduced

This provides an internal validation that the projection mathematics correctly reproduces the reference PCA coordinate system.

### Full-data projection result

All 23 sample-specific layers were projected successfully.

Final projected PCA matrix:

**106,592 × 50**

Cell identities and ordering were checked against the Phase 2 object.

Result:

**106,592 / 106,592 cells successfully projected.**

---

# 15. KNN Label Transfer

The annotated 20,000-cell sketch was used as the reference population.

The full 106,592-cell dataset was treated as the query population.

Label transfer was performed using nearest-neighbor voting in:

**PC1–PC20**

with:

**k = 30**

For every query cell:

1. Identify the 30 nearest reference cells.
2. Retrieve their atlas cell-type labels.
3. Count label votes.
4. Assign the label with the highest vote count.
5. Calculate neighbor agreement.
6. Calculate the vote margin relative to the second-most-supported label.

### Important distinction

The transfer was performed in **PC1–PC20**, matching the dimensionality used to construct the atlas neighbor graph.

The remaining PCA dimensions (21–50) were retained in the final object but were not used for label-transfer voting.

---

# 16. Label Transfer Results

All **106,592 cells** received a transferred label.

### Overall transferred composition

| Cell type           |       Cells | Percentage |
| ------------------- | ----------: | ---------: |
| CD8_T               |      41,281 |     38.73% |
| NK                  |      25,316 |     23.75% |
| MAIT                |      16,945 |     15.90% |
| CD4_T               |      10,635 |      9.98% |
| B_Cell              |       3,432 |      3.22% |
| Unresolved_Lymphoid |       2,831 |      2.66% |
| Monocyte_Myeloid    |       2,489 |      2.34% |
| Treg                |       1,534 |      1.44% |
| Dendritic_Cell      |       1,042 |      0.98% |
| Plasma_Cell         |         812 |      0.76% |
| pDC                 |         275 |      0.26% |
| **Total**           | **106,592** |   **100%** |

These are **descriptive cell counts from the analyzed dataset**. They should not be interpreted as population-level estimates of immune-cell frequencies in the broader HBV patient population.

---

# 17. Transfer Confidence

Confidence was evaluated from the agreement of the 30 nearest reference neighbors.

### Definitions

**High confidence**

```text
Neighbor agreement ≥ 0.80
```

**Low confidence**

```text
Neighbor agreement < 0.80
```

Low-confidence cells were **flagged but retained**.

No cells were removed on the basis of transfer confidence.

### Overall result

| Transfer confidence |       Cells | Percentage |
| ------------------- | ----------: | ---------: |
| High                |      93,176 |     87.41% |
| Low                 |      13,416 |     12.59% |
| **Total**           | **106,592** |   **100%** |

At least 95% neighbor agreement was observed for:

**79,879 cells.**

---

# 18. Confidence by Cell Type

Agreement was not uniform across cell types.

| Cell type           |  Cells | Mean agreement | High-confidence |
| ------------------- | -----: | -------------: | --------------: |
| Plasma_Cell         |    812 |          0.999 |           99.9% |
| B_Cell              |  3,432 |          0.999 |           99.8% |
| pDC                 |    275 |          0.995 |           98.9% |
| Monocyte_Myeloid    |  2,489 |          0.991 |           98.4% |
| Dendritic_Cell      |  1,042 |          0.966 |           93.9% |
| NK                  | 25,316 |          0.964 |           92.9% |
| CD8_T               | 41,281 |          0.937 |           87.9% |
| MAIT                | 16,945 |          0.929 |           86.6% |
| Unresolved_Lymphoid |  2,831 |          0.874 |           74.9% |
| CD4_T               | 10,635 |          0.861 |           72.1% |
| Treg                |  1,534 |          0.829 |           64.0% |

The lower agreement observed for Treg, CD4 T-cell and unresolved lymphoid populations reflects greater overlap between neighboring lymphoid transcriptional states.

Confidence is therefore treated as a **quality-control variable**, not as proof that high-confidence cells are biologically correct or that low-confidence cells are biologically incorrect.

---

# 19. Targeted CD8 and Myeloid QC

Because the downstream scientific objectives focus particularly on CD8 T cells and myeloid populations, these compartments received targeted transfer diagnostics.

## CD8 T cells

Reference atlas clusters:

```text
0, 4, 7, 8
```

Transferred CD8 population:

* **41,281 cells**
* Mean agreement: **0.937**
* Median agreement: **1.000**
* High confidence: **87.9%**
* ≥95% agreement: **75.5%**

## Monocyte/Myeloid

Reference atlas clusters:

```text
15, 16
```

Transferred population:

* **2,489 cells**
* Mean agreement: **0.991**
* Median agreement: **1.000**
* High confidence: **98.4%**
* ≥95% agreement: **96.6%**

The broad `Monocyte_Myeloid` label is retained intentionally. Detailed macrophage/myeloid-state characterization is deferred to Phase 5.

---

# 20. Full-Dataset UMAP

A new UMAP was fitted to the **full dataset's projected PCA coordinates**.

Parameters:

* Input: `pca.full`
* Dimensions: **PC1–PC20**
* Neighbors: **30**
* Metric: cosine
* UMAP embedding: **106,592 × 2**

### Important methodological note

`umap.full` is a **new full-dataset UMAP** fitted from the projected PCA coordinates.

It is **not** the original 20,000-cell atlas UMAP projected onto the full dataset.

This distinction is preserved in the analysis documentation.

---

# 21. Clinical-State Transfer Diagnostics

Mean transfer agreement by clinical state:

| Clinical state |  Cells | Mean agreement | High confidence | Low confidence |
| -------------- | -----: | -------------: | --------------: | -------------: |
| AC             | 12,837 |          0.943 |           89.0% |          11.0% |
| AR             | 17,746 |          0.943 |           89.0% |          11.0% |
| IA             | 32,289 |          0.937 |           87.8% |          12.2% |
| IT             | 19,501 |          0.928 |           85.6% |          14.4% |
| NL             | 24,219 |          0.929 |           86.3% |          13.7% |

These values demonstrate that transfer confidence was evaluated across all clinical states rather than only in the largest or most disease-associated groups.

---

# 22. Final Integrity Check

The final Phase 3 object contains:

* **106,592 cells**
* **24,452 genes**
* **23 samples**
* **23 donors**
* **5 clinical states**
* **20 atlas clusters**
* **11 final cell-type categories including Unresolved_Lymphoid**
* `pca.full`: **106,592 × 50**
* `umap.full`: **106,592 × 2**
* all original 23 processed-expression layers
* transferred cell-type labels
* neighbor agreement scores
* vote margins
* confidence classifications

### Final validation

```text
Cells retained:             106,592 / 106,592
Genes retained:               24,452 / 24,452
Samples represented:               23 / 23
Donors represented:                23 / 23
Clinical states represented:        5 / 5
Cells successfully projected: 106,592 / 106,592
Cells receiving labels:       106,592 / 106,592
Cells filtered in Phase 3:              0
```

---

# 23. Final Output

The canonical final Phase 3 object is:

```text
results/rds_objects/phase3_final_full_dataset.rds
```

This is the object used as the starting point for downstream Phase 4 and Phase 5 analyses.

---

# 24. Key Methodological Decisions

### No additional normalization

The input data were already processed log-CP10K expression values.

`NormalizeData()` was therefore not applied again.

### No raw-count reconstruction

Raw UMI counts were not available and were not inferred from processed values.

### No additional batch correction

The Phase 3 atlas was constructed without introducing a batch-correction transformation.

### No global cell filtering

All 106,592 cells were retained.

### No automatic doublet removal in Phase 3

Phase 3 does not claim that the final atlas is doublet-free.

### No confidence-based cell removal

Low-confidence transferred cells remain in the dataset and are explicitly flagged.

### Broad global annotation

The global atlas uses conservative lineage-level labels.

Fine-grained CD8 and myeloid states are addressed in downstream phases.

---

# 25. Reproducibility Notes

The final projection implementation uses Seurat internal computational functions for sparse matrix statistics and scaling.

The analysis was therefore version-locked to:

```text
Seurat       5.5.1
SeuratObject 5.4.0
```

Changes to Seurat versions should be treated as a reproducibility event and the projection validation should be rerun before relying on the resulting coordinates.

The most important computational validation is the atlas self-projection test, which reproduced the stored atlas PCA coordinates with essentially zero numerical error.

---

# 26. Outputs

The exact output inventory should be treated as the files generated by the finalized Phase 3 scripts.

Key checkpoint objects include:

```text
results/rds_objects/phase3_atlas_sketch_umap_clustered.rds
results/rds_objects/phase3_atlas_sketch_final_annotation.rds
results/rds_objects/phase3_final_full_dataset.rds
```

The analysis also generates tables documenting:

* dataset and sketch composition;
* HVG selection;
* PCA configuration;
* cluster composition;
* marker statistics;
* marker-program overlap;
* module scores;
* SingleR results;
* final cluster annotations;
* full-dataset projection;
* label-transfer parameters;
* neighbor agreement;
* transfer confidence;
* cell-type composition;
* clinical-state composition;
* donor/sample composition.

Figures document:

* sketch representation;
* PCA diagnostics;
* atlas clustering;
* cluster annotation;
* marker programs;
* module scores;
* full-dataset UMAP;
* cell-type distribution;
* clinical-state distribution;
* donor distribution;
* transfer confidence;
* neighbor-agreement diagnostics.

---

# 27. Interpretation Boundaries

The Phase 3 cell counts describe the **composition of the analyzed liver-cell dataset**.

They do not establish:

* population prevalence in all people with chronic HBV;
* causal relationships between clinical states and cell frequencies;
* differential expression between clinical states;
* functional activation states;
* macrophage subtypes;
* disease mechanisms.

Those questions require the downstream analyses in Phases 4–9.

In particular, the observed transferred cell-type composition should not be interpreted as an independent biological finding until donor-level structure and appropriate statistical analysis are considered.

---

# 28. Downstream Analysis

## Phase 4 — CD8 T-Cell Analysis

The global atlas identifies **41,281 CD8_T cells**.

Phase 4 will subset this population for detailed characterization of CD8 transcriptional states, including:

* CD8 subclustering;
* differentiation/state programs;
* activation and cytotoxicity;
* exhaustion-associated states;
* clinical-state comparisons;
* donor-aware downstream analyses.

The global Phase 3 label is therefore a **starting population definition**, not the final CD8 biological interpretation.

---

## Phase 5 — Myeloid Analysis

The global atlas identifies **2,489 Monocyte_Myeloid cells**.

Phase 5 will perform dedicated myeloid analysis and determine whether biologically meaningful macrophage/monocyte states can be resolved within this broader global population.

This is why Phase 3 intentionally does not overclaim macrophage identities.

---

# 29. Phase 3 Summary

| Parameter                 |           Final value |
| ------------------------- | --------------------: |
| Input cells               |           **106,592** |
| Genes                     |            **24,452** |
| Samples                   |                **23** |
| Donors                    |                **23** |
| Clinical states           |                 **5** |
| Sketch size               |            **20,000** |
| Sketch fraction           |             **18.8%** |
| Atlas clusters            |                **20** |
| PCA dimensions            |                **50** |
| PCA loading features      |             **1,999** |
| Transfer dimensions       |            **PC1–20** |
| KNN k                     |                **30** |
| Cells projected           | **106,592 / 106,592** |
| Cells labeled             | **106,592 / 106,592** |
| High-confidence transfers |             **87.4%** |
| Low-confidence transfers  |             **12.6%** |
| Cells removed in Phase 3  |                 **0** |
| Final UMAP                |       **106,592 × 2** |

---

# 30. Final Status

## ✅ PHASE 3 COMPLETE

The Phase 3 global atlas and full-dataset label-transfer workflow has passed its computational integrity checks.

The final analysis preserves:

* all 23 samples;
* all 23 donors;
* all 5 clinical states;
* all 106,592 cells;
* the original processed-expression representation;
* the locked atlas PCA;
* transferred labels;
* transfer-confidence metrics;
* targeted CD8 and myeloid diagnostics.

The final Phase 3 object is:

```text
results/rds_objects/phase3_final_full_dataset.rds
```

**Phase 3 is locked.**

The next analysis begins from this final object rather than modifying the global atlas.

---

## Data Source

**GEO:** GSE182159

The source dataset should be obtained from the public GEO repository rather than redistributed through this repository.

## Reference Annotation

**Monaco Immune Cell reference** via SingleR.

## Software

* R
* Seurat 5.5.1
* SeuratObject 5.4.0
* SingleR
* Monaco Immune reference
* tidyverse
* ggplot2
* pheatmap
* RANN / nearest-neighbor utilities as applicable to the finalized implementation

---

**Phase 3 finalized:** September 2026
**Status:** 🔒 **LOCKED / COMPLETE**
