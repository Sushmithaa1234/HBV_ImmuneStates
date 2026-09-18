# Phase 2: Data Loading & Quality Control

## Overview

**Goal:** Construct and quality-assess a production-ready Seurat object from the 23 liver expression matrices while preserving the processed expression representation supplied by GEO and avoiding inappropriate raw-count assumptions.

**Input:** 23 liver samples from GSE182159
**Input cells:** 106,592
**Input genes:** 24,452
**Expression representation:** Processed log-counts-per-10,000 (log-CP10K), as supplied in the GEO supplementary expression matrices
**Output:** 23 samples, 106,592 cells, 24,452 genes
**Cells removed by Phase 2 QC:** 0
**Overall retention:** 100%

> **Important:** The supplied matrices contain processed log-CP10K expression values rather than raw UMI counts. Therefore, conventional raw-count QC metrics such as `nCount_RNA`, count-based mitochondrial percentage, and RNA complexity defined as `nFeature_RNA / nCount_RNA` were **not calculated or interpreted as raw-count metrics**.

---

## Workflow

Phase 2 consisted of two scripts:

1. Load the 23 liver expression matrices and construct a merged Seurat v5 object while preserving sample and clinical metadata.
2. Calculate QC metrics appropriate for the processed log-CP10K representation.
3. Examine QC distributions globally and by sample/clinical state.
4. Apply a conservative minimum detected-gene threshold.
5. Perform sample-level diagnostic assessment using MAD statistics.
6. Flag unusually high-complexity cells for downstream inspection rather than automatically removing them.
7. Preserve all cells because no cells met the predefined minimum-gene criterion.
8. Save the validated processed-expression object for Phase 3.

**No normalization, batch correction, clustering, biological annotation, or doublet removal was performed in Phase 2.**

---

# Methods

## Script 01: Load and Merge

### Tasks

* Load all 23 liver expression matrices.
* Read the space-delimited processed expression files.
* Create one Seurat object per sample using the supplied processed expression values.
* Preserve sample identity and Phase 0 metadata:

  * GSM
  * Donor
  * Clinical state
  * Sample/tissue information
  * Phase 0 QC flags
* Merge all 23 individual objects into a unified Seurat v5 object.
* Preserve sample-specific assay layers.
* Verify metadata integrity and complete sample representation after merging.

### Seurat v5 Representation

Each sample was constructed using a Seurat v5 assay containing the supplied processed expression matrix as a `data` layer.

After merging, the RNA assay contains 23 sample-specific processed-expression layers:

```text
data.GSM5519469
data.GSM5519471
data.GSM5519472
...
data.GSM5519512
```

No raw `counts.*` layers are present.

### Output

`seurat_merged_processed_expression.rds`

**Merged dataset:**

* 23 samples
* 23 donors
* 5 clinical states
* 106,592 cells
* 24,452 genes

---

# Script 02: Processed-Expression QC

Because the input consists of log-CP10K expression rather than raw UMI counts, QC metrics were adapted accordingly.

## A. Cell-Level QC Metrics

| Metric                      | Description                                                                       |
| --------------------------- | --------------------------------------------------------------------------------- |
| **Genes_Detected**          | Number of genes with supplied processed expression > 0 per cell                   |
| **Mito_Genes_Detected**     | Number of detected mitochondrial genes per cell                                   |
| **Mito_Detection_Fraction** | Fraction of detected genes that are mitochondrial                                 |
| **Ribo_Genes_Detected**     | Number of detected genes matching the applied ribosomal-gene identifier pattern   |
| **Ribo_Detection_Fraction** | Fraction of detected genes matching the applied ribosomal-gene identifier pattern |
| **Total_LogCP10K**          | Sum of supplied log-CP10K values per cell; descriptive only                       |

### Interpretation of these metrics

`Genes_Detected` is a detection-based metric derived from the supplied processed expression matrix. It is **not equivalent to raw-count `nFeature_RNA`**.

`Total_LogCP10K` is **not a UMI/library-size metric**. It is retained only as a descriptive measure of the supplied processed expression values.

`Mito_Detection_Fraction` measures the fraction of detected genes that are mitochondrial. It is **not equivalent to conventional mitochondrial percentage calculated from raw UMI counts**.

---

## B. Mitochondrial and Ribosomal Features

Mitochondrial genes were identified using the gene-symbol pattern:

```r
^MT-
```

This identified **13 mitochondrial genes** in the supplied feature set.

Ribosomal genes were assessed using the applied `^RP[SL]` identifier pattern. This yielded **0 matching genes**.

This does **not** imply that the biological samples contain no ribosomal transcripts. Rather, ribosomal detection was **not informative in this feature representation/identifier set** and was not used as a QC exclusion criterion.

---

## C. Conservative QC Threshold

A single conservative cell-level criterion was applied:

| Criterion              | Threshold | Purpose                                 |
| ---------------------- | --------: | --------------------------------------- |
| Minimum detected genes |  **≥200** | Exclude extremely low-information cells |

No upper detected-gene cutoff was imposed.

No mitochondrial hard cutoff was imposed.

No donor or sample was automatically excluded based on cell number.

### Rationale

The matrices are already processed log-CP10K expression data, so raw-count QC thresholds cannot be transferred directly to this representation.

The Phase 2 strategy therefore prioritizes:

* preservation of biological populations,
* diagnostic assessment before exclusion,
* avoidance of arbitrary filtering,
* explicit documentation of the processed-expression limitations,
* and downstream inspection of unusual cells during atlas construction.

---

## D. Sample-Level MAD Diagnostics

Sample-specific Median Absolute Deviation (MAD) diagnostics were calculated to identify unusual QC distributions.

These diagnostics were **informational only** and were not used for automatic filtering.

This preserves potentially meaningful populations that may naturally display unusual expression complexity while still providing an audit trail for downstream inspection.

---

## E. High-Complexity Cell Flagging

Cells with:

```text
Genes_Detected > 5,000
```

were flagged as:

```text
Potential_HighComplexity
```

These cells were **retained**.

The flag is intended for downstream inspection during Phase 3 rather than being treated as proof of doublet status.

---

## F. Doublet Detection

**No automatic doublet detection or doublet removal was performed in Phase 2.**

In particular:

* `scDblFinder` was **not applied**.
* No predicted doublet/singlet classification was generated.
* No cells were removed as predicted doublets.
* No doublet rate was inferred from the processed expression representation.

### Rationale

The supplied matrices contain processed log-CP10K expression rather than raw UMI counts. Phase 2 therefore does not force a count-dependent doublet workflow onto the data.

Potential doublets and mixed-lineage cells will instead be evaluated in the context of the global atlas during Phase 3, where clustering, marker expression, and cell-type coherence can provide biological evidence for suspicious populations.

---

# Results

## Dataset Summary

| Metric                    |   Value |
| ------------------------- | ------: |
| **Samples**               |      23 |
| **Donors**                |      23 |
| **Clinical states**       |       5 |
| **Genes**                 |  24,452 |
| **Cells**                 | 106,592 |
| **Processed data layers** |      23 |
| **Raw counts layers**     |       0 |

Clinical-state cell totals:

| Clinical State |       Cells |
| -------------- | ----------: |
| **AC**         |      12,837 |
| **AR**         |      17,746 |
| **IA**         |      32,289 |
| **IT**         |      19,501 |
| **NL**         |      24,219 |
| **TOTAL**      | **106,592** |

---

## QC Metrics

### Global Distribution

| Metric                      |   Min |      Q1 |      Median |    Mean |      Q3 |     Max |
| --------------------------- | ----: | ------: | ----------: | ------: | ------: | ------: |
| **Genes_Detected**          |   524 |     886 |   **1,195** |   1,230 |   1,469 |   5,825 |
| **Mito_Genes_Detected**     |     0 |      12 |      **12** |   11.93 |      13 |      13 |
| **Mito_Detection_Fraction** |     0 | 0.00832 | **0.01008** | 0.01085 | 0.01286 | 0.02481 |
| **Ribo_Genes_Detected**     |     0 |       0 |       **0** |       0 |       0 |       0 |
| **Ribo_Detection_Fraction** |     0 |       0 |       **0** |       0 |       0 |       0 |
| **Total_LogCP10K**          | 509.9 | 1,988.2 | **2,322.1** | 2,295.5 | 2,581.3 | 4,132.6 |

### Interpretation

The processed-expression dataset shows a median of **1,195 detected genes per cell**.

The mitochondrial detection fraction is low across the dataset, with a median of approximately **1.01% of detected genes**.

The ribosomal metric was not informative because no features matched the applied ribosomal identifier pattern.

The minimum detected-gene value was **524**, meaning that every cell already exceeded the conservative 200-gene threshold.

---

# Fixed QC Filtering Results

| Criterion                      | Cells Removed | % of Input |
| ------------------------------ | ------------: | ---------: |
| **Genes_Detected < 200**       |             0 |         0% |
| **Mitochondrial hard cutoff**  |   Not applied |          — |
| **Upper detected-gene cutoff** |   Not applied |          — |
| **Doublet removal**            |   Not applied |          — |
| **Total cells removed**        |         **0** |     **0%** |

### Result

**106,592 / 106,592 cells passed the Phase 2 QC procedure.**

**Retention = 100%.**

The 200-gene threshold was therefore non-selective in this dataset.

---

# Clinical State Retention

| Clinical State |      Before |       After |  Lost | % Retained |
| -------------- | ----------: | ----------: | ----: | ---------: |
| **AC**         |      12,837 |      12,837 |     0 |       100% |
| **AR**         |      17,746 |      17,746 |     0 |       100% |
| **IA**         |      32,289 |      32,289 |     0 |       100% |
| **IT**         |      19,501 |      19,501 |     0 |       100% |
| **NL**         |      24,219 |      24,219 |     0 |       100% |
| **TOTAL**      | **106,592** | **106,592** | **0** |   **100%** |

Because no cells were removed, the clinical-state composition of the dataset was unchanged.

---

# High-Complexity Cells

Five cells had:

```text
Genes_Detected > 5,000
```

These cells were flagged with:

```text
Potential_HighComplexity = TRUE
```

They were retained in the dataset.

Their presence will be considered during Phase 3 atlas construction and biological annotation rather than being automatically classified as doublets.

---

# Flagged Sample: GSM5519491

**GSM5519491 — P190902 — IT**

* Pre-analysis cell count: **66**
* Post-Phase-2 cell count: **66**
* Retention: **100%**
* Status: retained

This sample was identified during Phase 0 as a major low-cell-recovery outlier.

It was **not excluded** in Phase 2 because:

1. sample size alone is not a valid per-cell QC criterion;
2. no cells failed the minimum detected-gene threshold;
3. automatic donor/sample exclusion was not part of the Phase 2 protocol.

The sample should remain visible during Phase 3 clustering and atlas interpretation because its small cell count limits the amount of biological information it can contribute.

---

# Key Findings

## 1. No Cells Failed the Conservative QC Criterion

All 106,592 cells had at least 524 detected genes.

Therefore:

* minimum threshold = 200 genes;
* minimum observed = 524 genes;
* cells removed = 0;
* retention = 100%.

The QC procedure did not alter the biological cell population.

---

## 2. Mitochondrial Detection Was Low

Thirteen mitochondrial genes were present in the feature set.

The median mitochondrial detection fraction was approximately:

**1.01% of detected genes per cell.**

No mitochondrial hard cutoff was applied because the available metric is detection-based rather than raw-count-based.

---

## 3. Ribosomal QC Was Not Informative

The applied ribosomal-gene identifier pattern detected zero matching features.

This is documented as **unavailable/not informative**, rather than interpreted as absence of ribosomal transcripts.

No filtering decision was based on this metric.

---

## 4. No Doublet Filtering Was Performed

The Phase 2 object does **not** represent a doublet-free dataset.

There is:

* no scDblFinder classification,
* no estimated Phase 2 doublet rate,
* no automatic doublet removal.

Potential doublets or mixed-lineage cells will be assessed in Phase 3 using atlas-level biological evidence.

---

## 5. All Clinical States Were Preserved

All five clinical states retained their complete cell populations:

* AC: 12,837
* AR: 17,746
* IA: 32,289
* IT: 19,501
* NL: 24,219

No clinical-state representation was altered by Phase 2 filtering.

---

# QC Decision Summary

| Decision                                | Outcome            | Rationale                                                    |
| --------------------------------------- | ------------------ | ------------------------------------------------------------ |
| **Processed-expression representation** | Accepted           | GEO matrices contain processed log-CP10K values              |
| **Raw UMI QC metrics**                  | Not calculated     | Raw UMI counts are unavailable                               |
| **Minimum detected-gene filter**        | Applied            | Conservative threshold of ≥200 genes                         |
| **Cells removed by gene threshold**     | **0**              | Minimum observed = 524                                       |
| **Mitochondrial hard cutoff**           | Not applied        | Detection fraction is not raw-count mitochondrial percentage |
| **Upper gene cutoff**                   | Not applied        | High-complexity cells retained for downstream inspection     |
| **Ribosomal QC**                        | Informational only | Identifier pattern yielded zero matches                      |
| **MAD diagnostics**                     | Informational only | Used for sample-level transparency                           |
| **Doublet detection**                   | Not performed      | Avoid count-dependent inference from processed expression    |
| **Doublet removal**                     | Not performed      | No Phase 2 singlet/doublet classification                    |
| **Clinical-state exclusion**            | None               | All states retained 100%                                     |
| **Sample exclusion**                    | None               | All 23 samples retained                                      |
| **Final Phase 2 object**                | **106,592 cells**  | No cells removed                                             |

---

# Output & Checkpoints

## Primary Deliverable

### `seurat_qc_processed_expression.rds`

Final Phase 2 Seurat object containing:

* **106,592 cells**
* **24,452 genes**
* **23 samples**
* **23 donors**
* **5 clinical states**
* Processed log-CP10K expression
* Sample-specific Seurat v5 data layers
* QC metadata
* QC flags
* `Potential_HighComplexity` flags
* Original sample/donor/clinical metadata

This object is the input for Phase 3.

---

## Supporting Data Files

Phase 2 generates supporting QC outputs including:

* `sample_qc_before_filtering.csv`
* `qc_filtering_decision_summary.csv`
* `phase_cell_retention.csv`
* `sample_cell_retention.csv`
* `qc_summary.csv`
* sample-level QC summaries
* sample-specific MAD diagnostics
* high-complexity cell summary
* final QC status
* QC figures

Exact filenames should follow the files generated by the finalized Phase 2 scripts.

---

# Visualizations

The Phase 2 QC figures focus on the metrics that are valid for the supplied processed-expression representation.

Key visualizations include:

1. **Genes detected by clinical state**
2. **Genes detected by sample**
3. **Mitochondrial detection fraction by clinical state**
4. **Mitochondrial detection fraction by sample**
5. **Total log-CP10K distribution**
6. **Genes detected versus total log-CP10K**
7. **Sample-level QC distributions**
8. **Sample-specific MAD diagnostics**
9. **High-complexity cell distribution**
10. **Cell-retention summaries by clinical state**
11. **Cell-retention summaries by sample**

Raw-count-specific plots such as UMI-vs-gene complexity and raw-count mitochondrial percentage are intentionally **not included**, because raw UMI counts are unavailable.

---

# Technical Notes

## Processed Expression Versus Raw Counts

The supplied GEO matrices represent processed expression values generated after normalization and log transformation.

Consequently, the following conventional raw-count metrics should **not** be interpreted as available:

```text
nCount_RNA
raw UMI library size
raw-count RNA complexity
percent.mt based on UMI counts
```

Instead, Phase 2 uses detection-based and processed-expression metrics that can be calculated without reconstructing unavailable raw counts.

**Raw counts must not be reconstructed from log-CP10K values.**

---

## Why No Doublet Removal?

Doublet detection is most informative when supported by an appropriate expression/count representation and should not be treated as a mandatory checkbox independent of the available data.

Because this dataset is supplied as processed log-CP10K expression:

* no raw UMI counts are available;
* no count-based doublet rate can be established;
* Phase 2 does not manufacture a doublet classification from unsuitable input;
* suspicious cells can instead be evaluated during global atlas construction using expression profiles, cluster coherence, and lineage-marker relationships.

This keeps Phase 2 focused on **data integrity and conservative QC**, rather than introducing an unsupported exclusion step.

---

## Why Was No Upper Gene Threshold Applied?

Five cells exceeded 5,000 detected genes.

Rather than automatically excluding these cells, they were flagged:

```text
Potential_HighComplexity = TRUE
```

High detected-gene counts can arise from genuine biological complexity as well as technical multiplets. Without an appropriate raw-count representation and independent biological evidence, automatic removal would risk discarding legitimate cells.

Phase 3 provides a more appropriate context for evaluating these cells.

---

## MAD Diagnostics

MAD-based diagnostics are retained as a transparent diagnostic layer.

They identify samples with unusual QC distributions without automatically declaring those samples or their cells to be invalid.

This is particularly important for this cohort because cell recovery varies substantially between samples, including the 66-cell GSM5519491 sample.

---

# Phase 2 Gate

## ✅ GATE: PASS

### Criteria

* ✅ All 23 samples retained
* ✅ All 23 donors retained
* ✅ All 5 clinical states retained
* ✅ 106,592 cells retained
* ✅ 24,452 genes retained
* ✅ No cells failed the ≥200 detected-gene criterion
* ✅ Mitochondrial detection assessed
* ✅ Ribosomal metric documented as non-informative
* ✅ No inappropriate raw-count metrics inferred from processed data
* ✅ No arbitrary donor/sample exclusion
* ✅ High-complexity cells flagged rather than automatically removed
* ✅ QC metadata attached and validated
* ✅ Final object saved successfully
* ✅ Dataset ready for Phase 3

---

# Workflow Summary

```text
23 GEO liver expression matrices
        │
        │  Processed log-CP10K expression
        ↓
Load individual sample matrices
        │
        ↓
Create Seurat v5 objects
        │
        ↓
Merge 23 samples
        │
        ├── 24,452 genes
        ├── 106,592 cells
        ├── 23 donors
        └── 5 clinical states
        │
        ↓
Verify processed-expression layers
        │
        ├── 23 data layers
        └── 0 counts layers
        │
        ↓
Calculate processed-expression QC metrics
        │
        ├── Genes detected
        ├── Mitochondrial detection
        ├── Ribosomal detection
        └── Total log-CP10K
        │
        ↓
Apply conservative minimum
Genes_Detected ≥ 200
        │
        ↓
0 cells removed
        │
        ↓
Flag 5 high-complexity cells
        │
        ↓
No doublet filtering
        │
        ↓
Final Phase 2 object
106,592 cells × 24,452 genes
        │
        ↓
READY FOR PHASE 3
```

---

# Next Steps

## Phase 3: Global Atlas Construction

The Phase 3 workflow should use:

```text
seurat_qc_processed_expression.rds
```

as its starting object.

Planned Phase 3 activities:

1. Construct the global reference atlas.
2. Use leverage-score sketching to create a representative reference subset.
3. Perform dimensional reduction and clustering on the reference.
4. Identify major immune and hepatic cell populations.
5. Establish marker-based cell-type annotations.
6. Project annotations onto the full dataset.
7. Inspect high-complexity cells and potential mixed-lineage populations.
8. Perform secondary biological QC.
9. Produce the annotated global atlas for downstream CD8 T-cell and myeloid analyses.

Particular attention should be given to:

* **GSM5519491 / P190902 / IT**, which contains only 66 cells;
* the five `Potential_HighComplexity` cells;
* ambiguous or mixed-lineage clusters;
* concordance between computational clusters and canonical marker expression.

---

# Phase 2 Complete ✅

**Input:** 106,592 cells × 24,452 genes
**Expression:** Processed log-CP10K
**Cells removed:** 0
**Retention:** **100%**
**Samples retained:** 23/23
**Donors retained:** 23/23
**Clinical states retained:** 5/5
**High-complexity cells flagged:** 5
**Doublet removal:** Not performed
**Final object:** `seurat_qc_processed_expression.rds`
**Gate Status:** **READY FOR PHASE 3**
