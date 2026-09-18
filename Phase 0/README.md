# Phase 0: Pre-Analysis Audit

## Overview

**Goal:** Establish a verified, complete sample inventory before committing to downstream analysis. This phase confirms data provenance, validates cohort structure, verifies metadata integrity, and documents technical characteristics without performing cell filtering, normalization, batch correction, or biological analysis.

**Dataset:** GEO Accession **GSE182159** — Single-cell RNA-seq study of liver and peripheral immune cells associated with hepatitis B virus (HBV) infection.

**Scope:** **23 liver biopsy samples from 23 unique donors**, representing 5 clinical-state labels present in the GEO metadata.

**Total liver cells represented in the downloaded matrices:** **106,592**

---

## Phase 0 Principles

Phase 0 is a **provenance and data-integrity audit**, not a biological analysis or quality-control filtering stage.

The following principles were applied:

* GEO sample metadata were treated as the authoritative source for sample identifiers, donor identifiers, tissue, and clinical-state labels.
* Clinical-state labels were preserved exactly as provided by GEO.
* No cells or samples were removed.
* No QC thresholds were used to determine inclusion/exclusion.
* No normalization was performed because the downloaded matrices are already processed expression data.
* No batch correction was performed.
* No clustering, cell-type annotation, differential expression, or other biological analysis was performed.
* Raw count matrices were not assumed to be available.
* Descriptive matrix statistics were recorded without using them to make downstream filtering decisions.

---

# Workflow

## 1. Identify and Verify Downloaded Files

* Identified liver expression matrices within `Data/Raw/`.
* Confirmed that **23 liver expression files** were present.
* Extracted GSM accession identifiers from filenames.
* Confirmed that all 23 GSM identifiers were unique.

**Result:** 23/23 expected liver samples detected.

---

## 2. Retrieve Official GEO Metadata

GEO metadata were retrieved programmatically using the **GEOquery** R package from the GSE182159 Series Matrix.

The following metadata were extracted:

* GSM accession
* Donor identifier
* Clinical-state label (`Stage:ch1`)
* Tissue (`tissue:ch1`)
* Platform
* Instrument
* Library selection
* Library source
* Library strategy

Downloaded files were matched against the GEO metadata using their GSM accession identifiers.

---

## 3. Validate Sample and Donor Structure

The downloaded matrices and GEO metadata were cross-checked in both directions:

* Every downloaded GSM had a corresponding GEO liver sample.
* Every GEO liver sample had a corresponding downloaded matrix.
* GSM identifiers were unique.
* Donor identifiers were unique.
* All 23 samples were annotated as **Liver**.
* No donor, tissue, or clinical-state metadata were missing.

**Result:**

* **23 liver samples**
* **23 unique donors**
* **1 liver sample per donor**
* **100% bidirectional GSM matching**
* **100% liver-tissue concordance**

---

## 4. Clinical-State Composition

The clinical-state labels were taken directly from the GEO `Stage:ch1` metadata.

| GEO Label | Clinical Interpretation         | N Samples | N Donors |
| --------- | ------------------------------- | --------: | -------: |
| **AC**    | Asymptomatic Carrier*           |         3 |        3 |
| **AR**    | Acute Resolved / Acute Recovery |         3 |        3 |
| **IA**    | Immune Active                   |         5 |        5 |
| **IT**    | Immune Tolerant                 |         6 |        6 |
| **NL**    | Normal Liver                    |         6 |        6 |
| **TOTAL** |                                 |    **23** |   **23** |

*The GEO metadata use the label **AC**. Subsequent literature may harmonize the corresponding individuals with a chronic-resolved (`CR`) clinical group. **Phase 0 preserves the original GEO label and does not perform this harmonization.**

The three GEO `AC` samples are:

* GSM5519506 — P191126
* GSM5519508 — P191127
* GSM5519510 — P191210

Any future derived clinical grouping will be explicitly documented as a derived variable rather than overwriting the source metadata.

---

# Clinical Context

The five labels do not represent a simple linear disease-progression trajectory.

They represent distinct clinical contexts associated with HBV infection, resolution, or a normal liver reference population.

Broadly:

* **IT — Immune Tolerant:** chronic HBV infection characterized in classical clinical frameworks by high viral replication with relatively limited immune-mediated liver inflammation.
* **IA — Immune Active:** HBV-associated immune-active disease with evidence of active hepatic inflammation and increased ALT and/or viral replication depending on the clinical context.
* **AR — Acute Resolved/Acute Recovery:** individuals associated with acute HBV infection that subsequently resolved.
* **AC — Asymptomatic Carrier:** the clinical label used in the GEO metadata for three individuals. The terminology has subsequently been represented differently in some literature, including classification as chronic resolved (`CR`).
* **NL — Normal Liver:** HBV-free/normal liver reference samples.

Clinical-state interpretation should therefore be based on the original study's clinical metadata and definitions rather than treating the five labels as sequential stages.

Important clinical variables underlying HBV disease-state classification can include combinations of:

* HBsAg status
* HBeAg status
* HBV DNA level
* ALT
* Evidence of hepatic inflammation
* Liver histology
* Disease/resolution history

Phase 0 does not reconstruct or alter these clinical classifications. It records the source labels for subsequent analysis.

---

# Expression Matrix Characteristics

Each downloaded liver matrix contains:

* **24,452 genes**
* Sample-specific numbers of cells
* Numeric expression values
* No duplicated gene identifiers
* No duplicated cell identifiers
* No negative expression values detected
* No missing/non-finite expression values detected

The matrices were successfully parsed as space-delimited expression matrices with genes as rows and cells as columns.

---

## Data Representation

The downloaded matrices are **processed log-normalized expression data**, not raw UMI count matrices.

The data are documented by GEO as:

> **log-counts-per-10,000**

Therefore:

* `colSums(expression_matrix)` was **not interpreted as total UMI counts**.
* Per-cell nonzero gene counts were used only as a descriptive measure of **genes detected per cell**.
* Raw-count-dependent assumptions were not made during Phase 0.
* No additional normalization was performed.

**Raw count matrices:** Not available in the downloaded GEO supplementary data.

This distinction is important for downstream methodological choices because several standard scRNA-seq workflows assume access to raw count data.

---

# Cell Recovery

The downloaded matrices contain a total of:

**106,592 liver cells**

This total was obtained directly by summing the cell counts across the 23 local expression matrices:

```text
sum(final_metadata$cells_in_matrix)
[1] 106592
```

This independently reproduces the liver-cell total reported in subsequent analyses of GSE182159.

### Cells per clinical state

| GEO State | N Samples | Total Cells | Median Cells/Sample |
| --------- | --------: | ----------: | ------------------: |
| **NL**    |         6 |      24,219 |               4,488 |
| **IT**    |         6 |      19,501 |               2,889 |
| **IA**    |         5 |      32,289 |               6,070 |
| **AR**    |         3 |      17,746 |               4,122 |
| **AC**    |         3 |      12,837 |               3,293 |
| **TOTAL** |    **23** | **106,592** |                   — |

### Overall cell recovery

| Metric              |           Value |
| ------------------- | --------------: |
| Total cells         |     **106,592** |
| Mean cells/sample   |     **4,634.4** |
| Median cells/sample |       **4,122** |
| Range               | **66 – 10,164** |

No samples were excluded based on these values during Phase 0.

---

# Transcriptomic Complexity

Genes detected per cell were calculated descriptively as the number of genes with expression greater than zero in each cell.

Because the matrices are processed log-normalized expression values, this metric is referred to as:

> **genes detected per cell**

rather than `nFeature_RNA` or raw-count-derived transcriptomic complexity.

### Sample-level summary

| Metric                                            |             Value |
| ------------------------------------------------- | ----------------: |
| Median of sample-level median genes detected/cell |         **1,193** |
| Mean of sample-level median genes detected/cell   |       **1,123.7** |
| Range of sample-level medians                     | **643.5 – 1,546** |

These values are descriptive and were **not used to exclude samples in Phase 0**.

---

# Sample-Level Audit

| GSM        | Donor   | GEO Phase |  Cells | Median Genes Detected/Cell | Mean Genes Detected/Cell |
| ---------- | ------- | --------- | -----: | -------------------------: | -----------------------: |
| GSM5519469 | P190604 | IT        |  3,085 |                      1,278 |                  1,438.6 |
| GSM5519471 | P190326 | IT        |  1,379 |                        843 |                  1,157.4 |
| GSM5519472 | P190402 | IT        |  7,110 |                      1,464 |                  1,540.9 |
| GSM5519475 | P190716 | AR        |  3,460 |                      1,091 |                  1,105.3 |
| GSM5519477 | P190719 | IA        |  3,006 |                        920 |                  1,027.1 |
| GSM5519483 | Dhc570  | NL        |  4,243 |                      1,194 |                  1,207.8 |
| GSM5519484 | D528848 | NL        |  4,733 |                      1,171 |                  1,194.0 |
| GSM5519485 | D529074 | NL        |  5,956 |                      1,193 |                  1,268.1 |
| GSM5519486 | D529351 | NL        |    690 |                      674.5 |                    791.0 |
| GSM5519487 | D529354 | NL        |  5,233 |                        850 |                    945.0 |
| GSM5519488 | D529409 | NL        |  3,364 |                        729 |                    813.2 |
| GSM5519491 | P190902 | IT        |     66 |                      643.5 |                    656.5 |
| GSM5519494 | P190910 | IT        |  2,693 |                      1,546 |                  1,625.3 |
| GSM5519495 | P190808 | IT        |  5,168 |                      1,527 |                  1,547.7 |
| GSM5519496 | P190801 | IA        |  6,069 |                      1,277 |                  1,318.9 |
| GSM5519497 | P190911 | IA        |  9,082 |                      1,307 |                  1,352.8 |
| GSM5519499 | P191028 | IA        |  8,062 |                        685 |                    745.1 |
| GSM5519502 | P191112 | IA        |  6,070 |                      1,430 |                  1,451.4 |
| GSM5519504 | P191008 | AR        | 10,164 |                      1,008 |                  1,058.4 |
| GSM5519506 | P191126 | AC        |  3,293 |                      1,426 |                  1,449.4 |
| GSM5519508 | P191127 | AC        |  6,625 |                      1,423 |                  1,438.6 |
| GSM5519510 | P191210 | AC        |  2,919 |                        826 |                    893.4 |
| GSM5519512 | P191217 | AR        |  4,122 |                      1,339 |                  1,369.2 |

---

# Technical Metadata

All 23 liver samples share the following GEO technical metadata:

| Parameter                     | Value                               |
| ----------------------------- | ----------------------------------- |
| **Platform**                  | GPL20301                            |
| **Instrument**                | Illumina HiSeq 4000                 |
| **Library Strategy**          | RNA-Seq                             |
| **Library Selection**         | cDNA                                |
| **Library Source**            | transcriptomic                      |
| **Genes per matrix**          | 24,452                              |
| **Expression representation** | Processed log-normalized expression |
| **Normalization**             | log-counts-per-10,000               |
| **Raw counts available**      | No                                  |

No platform, instrument, library strategy, library selection, or library source variation was observed among the 23 downloaded liver samples.

**Important:** technical metadata consistency does not demonstrate the absence of all batch effects. It only establishes that these particular recorded technical variables do not vary across the cohort.

---

# Samples Requiring Attention

Phase 0 does **not** assign PASS/FAIL QC categories and does **not** exclude samples.

However, descriptive inspection identifies samples that should receive attention during the formal downstream QC stage.

### GSM5519491 — P190902 — IT

* **66 cells**
* Median genes detected/cell: **643.5**
* This is substantially smaller than the other liver samples.
* The sample should be explicitly evaluated during the formal QC stage.

**Phase 0 decision:** Retain in the dataset.

**Downstream decision:** Not predetermined here. Phase 1/2 QC must establish whether the sample should be excluded and document the rationale.

### Lower-complexity samples

Several samples have comparatively low sample-level median genes detected per cell, including:

* GSM5519486 — NL — 674.5
* GSM5519487 — NL — 850
* GSM5519488 — NL — 729
* GSM5519491 — IT — 643.5
* GSM5519499 — IA — 685
* GSM5519510 — AC — 826
* GSM5519471 — IT — 843

These observations are **not sufficient by themselves to justify sample exclusion**.

Formal QC will consider cell-level distributions and other appropriate metrics before any filtering decisions are made.

---

# Key Findings

## 1. Cohort Completeness

* **23/23** expected liver samples were identified.
* **23/23** downloaded samples matched GEO metadata.
* **23/23** GEO liver samples had corresponding local matrices.
* **23/23** donors were unique.
* **23/23** samples were annotated as liver tissue.

## 2. Matrix Integrity

* Every matrix contains **24,452 genes**.
* No duplicated gene identifiers were detected.
* No duplicated cell identifiers were detected.
* Expression matrices were successfully parsed as numeric data.
* Expression values ranged from zero to positive values.
* No missing or non-finite expression values were detected.

## 3. Cell Count Validation

The 23 local matrices contain exactly:

**106,592 liver cells**

This provides an independent consistency check against the published liver-cell cohort associated with GSE182159.

## 4. Clinical-State Metadata

The source GEO distribution is:

```text
AC = 3
AR = 3
IA = 5
IT = 6
NL = 6
```

No clinical-state labels were manually altered during Phase 0.

## 5. Data Representation

The downloaded matrices are **processed log-normalized expression data** rather than raw count matrices.

Consequently, raw UMI-based metrics were not reconstructed or reported.

## 6. Technical Consistency

All samples share the same recorded:

* Platform
* Instrument
* Library strategy
* Library source
* Library selection

This reduces obvious technical heterogeneity in the recorded metadata, but does **not** establish that biological and technical variation are completely separable.

---

# What Phase 0 Did NOT Do

Phase 0 deliberately did not perform:

* Cell filtering
* Sample exclusion
* QC threshold-based removal
* Normalization
* Batch correction
* Feature selection
* Dimensionality reduction
* Clustering
* Cell-type annotation
* Differential expression
* Pathway analysis
* Cell–cell communication analysis
* Clinical-state differential testing
* Biological interpretation

These decisions belong to subsequent phases and must be justified according to the relevant analytical objective.

---

# Deliverables

### Data Files

* `results/tables/phase0_final_liver_sample_metadata.csv`

  * GEO-verified sample metadata
  * Donor identifiers
  * Clinical-state labels
  * Tissue and technical metadata
  * Matrix-level descriptive statistics

* `results/tables/phase0_expression_matrix_audit.csv`

  * Per-sample matrix integrity and expression statistics

* `results/tables/phase0_geo_phase_summary.csv`

  * Sample and cell summaries by GEO clinical-state label

### Figures

* `results/figures/phase0_cells_per_sample.png`
* `results/figures/phase0_median_genes_detected.png`

### Code

* `code/phase_0_pre-analysis audit.R`

---

# Quality-Control Checkpoint

### Phase 0 Status: **PASS**

The dataset passed the Phase 0 provenance and integrity audit:

* **23 liver samples verified**
* **23 unique donors verified**
* **106,592 liver cells verified**
* **24,452 genes per matrix verified**
* **GEO metadata matched**
* **Clinical-state labels preserved**
* **Matrix integrity verified**
* **Processed-data representation documented**
* **No cells or samples filtered**

### Important downstream consideration

**GSM5519491 (P190902, IT) contains only 66 cells and requires formal evaluation during downstream QC.**

Phase 0 does **not** pre-commit to its exclusion.

---

# Technical Notes

## Data Source

**GEO Accession:** GSE182159

**Tissue analyzed:** Liver

**Number of liver donors:** 23

**Number of liver samples:** 23

**Total liver cells:** 106,592

**Genes per matrix:** 24,452

**Metadata source:** GEO Series Matrix retrieved programmatically using GEOquery.

## Expression Data

The downloaded supplementary matrices represent processed expression values normalized as log-counts-per-10,000.

Raw sequencing/count data were not available for reconstruction of raw UMI counts.

Therefore, Phase 0 reports **genes detected per cell** rather than raw-count-derived UMI metrics.

## Clinical Metadata Provenance

Clinical-state labels in the Phase 0 metadata correspond directly to the GEO `Stage:ch1` field.

No manual relabeling was performed.

The `AC` label is retained exactly as supplied by GEO. Its terminology and relationship to the chronic-resolved (`CR`) terminology used in subsequent literature will be documented separately and, if required for downstream comparisons, represented as an explicitly derived clinical variable.

---

# File Structure

```text
Phase_0_PreAnalysis/
├── README.md
├── code/
│   └── phase_0_pre-analysis audit.R
├── results/
│   ├── tables/
│   │   ├── phase0_final_liver_sample_metadata.csv
│   │   ├── phase0_expression_matrix_audit.csv
│   │   └── phase0_geo_phase_summary.csv
│   └── figures/
│       ├── phase0_cells_per_sample.png
│       └── phase0_median_genes_detected.png
```

---

# Next Step

**Phase 1: Data Loading & Formal QC**

Phase 1 will establish the appropriate downstream QC framework using the characteristics of the processed dataset documented here.

Any cell- or sample-level exclusions will be determined and documented in that phase rather than retroactively imposed by Phase 0.

---

## Phase 0 Complete ✅

**Samples verified:** 23
**Unique donors:** 23
**Liver cells:** 106,592
**Genes/matrix:** 24,452
**GEO metadata match:** 23/23
**Clinical labels altered:** No
**Cells filtered:** No
**Raw counts assumed:** No
**Gate status:** **READY FOR PHASE 1**
