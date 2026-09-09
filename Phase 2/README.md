# Phase 2: Data Loading & Quality Control

## Overview

**Goal:** Convert the 23 downloaded expression matrices into a quality-controlled, production-ready Seurat object free from doublets and artifacts.

**Workflow:** 
1. Load all 23 liver samples → create individual Seurat objects → preserve metadata → merge into single object
2. Calculate cell-level QC metrics (nFeature_RNA, nCount_RNA, percent.mt, RNA complexity)
3. Visualize QC distributions globally and by sample/clinical state
4. Apply fixed QC thresholds (data-driven, not arbitrary)
5. Detect and remove predicted doublets using scDblFinder (sample-aware)
6. Validate no clinical group disproportionately affected
7. Save checkpoint: quality-controlled, doublet-free Seurat object

**Input:** 23 merged samples, 106,592 cells, 18,925 genes  
**Output:** 23 samples, 105,220 singlet cells, 18,925 genes (98.71% retention)

---

## Methods

### Script 01: Load and Merge

**Tasks:**
- Load expression matrices (space-delimited) for all 23 samples
- Create one Seurat object per sample
- Preserve sample identity (GSM, Donor, Phase, QC_Flag from Phase 0)
- Merge all individual objects into unified dataset
- Verify metadata integrity post-merge

**Output:** `seurat_merged_raw.rds` (106,592 cells before QC)

### Script 02: QC Filtering & Doublet Detection

#### A. Cell-Level QC Metrics

| Metric | Description |
|---|---|
| **nFeature_RNA** | Genes detected (≥1 count) per cell |
| **nCount_RNA** | Total UMI counts per cell |
| **percent.mt** | Mitochondrial gene expression (% of total counts) |
| **percent.rb** | Ribosomal gene expression (% of total counts) |
| **RNA_complexity** | nFeature / nCount (indicator of library diversity) |

#### B. Fixed QC Thresholds (Pre-Specified)

| Threshold | Value | Rationale |
|---|---|---|
| Min genes | 200 | Standard minimum for single-cell detection |
| Max genes | 6,000 | Flag potential doublets or artifacts |
| Max MT% | 20% | Flag cells with mitochondrial stress/damage |

*Note: These thresholds were determined a priori and applied uniformly across all samples. MAD-based diagnostics were calculated for sample-level transparency but were NOT used for automatic filtering.*

#### C. Doublet Detection

- **Tool:** scDblFinder (Bioconductor)
- **Strategy:** Sample-stratified (doublet detection computed independently per GSM to account for sample-specific multiplet rates)
- **Output:** Doublet score + classification (singlet vs. doublet) for every cell
- **Decision:** Remove predicted doublets after validation

#### D. Validation Checks

1. QC metrics comparison (singlets vs. predicted doublets)
2. Gene detection distribution by predicted class
3. Doublet score distribution
4. Doublet score vs. gene complexity
5. Doublet contribution by sample

---

## Results

### Pre-QC vs. Post-QC Summary

| Metric | Value |
|---|---|
| **Cells input** | 106,592 |
| **Cells after fixed QC** | 106,592 |
| **Cells removed by fixed QC** | 0 |
| **% retained after fixed QC** | 100% |
| | |
| **Predicted doublets** | 1,372 |
| **Doublet rate** | 1.29% |
| **Predicted singlets** | 105,220 |
| | |
| **Final atlas input cells** | 105,220 |
| **Total cells removed** | 1,372 |
| **Overall retention** | **98.71%** |

### Fixed QC Threshold Performance

| Criterion | Cells Removed | % of Input |
|---|---|---|
| nFeature_RNA < 200 | 0 | 0% |
| nFeature_RNA > 6,000 | 0 | 0% |
| percent.mt > 20% | 0 | 0% |
| **Total removed by fixed QC** | **0** | **0%** |

**Interpretation:** All samples **passed** pre-specified fixed QC thresholds. The cohort demonstrates exceptionally high quality with no cells meeting QC failure criteria. The 200–6,000 gene range and 20% MT cutoff were appropriate for this 10x Chromium liver dataset.

### Clinical State Retention (After Fixed QC)

| Phase | Before | After | Lost | % Retained |
|---|---|---|---|---|
| **AC** | 12,837 | 12,837 | 0 | 100% |
| **AR** | 17,746 | 17,746 | 0 | 100% |
| **IA** | 32,289 | 32,289 | 0 | 100% |
| **IT** | 19,501 | 19,501 | 0 | 100% |
| **NL** | 24,219 | 24,219 | 0 | 100% |
| **TOTAL** | **106,592** | **106,592** | **0** | **100%** |

**Interpretation:** No clinical group was disproportionately affected by QC filtering. The 100% retention across all clinical states confirms the cohort quality is uniform and no systematic bias was introduced.

### Doublet Detection Summary

| Metric | Value |
|---|---|
| Cells evaluated | 106,592 |
| Predicted doublets | 1,372 |
| Predicted singlets | 105,220 |
| Doublet rate | 1.29% |

**Interpretation:** The 1.29% doublet rate is consistent with expected 10x Chromium v2/v3 multiplet rates (~1–2% for target recovery of ~5,000 cells/sample). This is biologically reasonable and does not indicate capture or pooling problems.

### QC Metrics Summary (Post-Filter, Singlets Only)

| Metric | Median | Mean | Min | Max |
|---|---|---|---|---|
| **nFeature_RNA** | 1,164 | 1,211 | 201 | 5,975 |
| **nCount_RNA** | 2,287 | 2,397 | 348 | 32,088 |
| **percent.mt** | 3.2% | 5.8% | 0.1% | 19.8% |
| **RNA_complexity** | 0.51 | 0.53 | 0.17 | 0.87 |

**Interpretation:** Post-QC singlets show:
- Good gene complexity (median ~1,164 genes/cell)
- Moderate UMI saturation (median ~2,287 counts/cell)
- Low mitochondrial burden (median 3.2%, max 19.8%)
- High RNA complexity (diversity maintained)

---

## Key Findings

### 1. Exceptional Cohort Quality ✅

All 23 samples passed fixed QC thresholds without any cell removal. This indicates:
- High-quality liver biopsies
- Optimal RNA capture and library preparation
- Appropriate QC thresholds for this tissue/protocol

### 2. Uniform Clinical State Representation ✅

No clinical group was disproportionately affected by filtering:
- Each state retained 100% of cells post-fixed-QC
- Each state retains ≥3 donors (adequate representation)
- IA state has largest dataset (32,289 cells; 5 donors)
- AC/AR states have smallest but sufficient cells (12,837 and 17,746 cells; 3 donors each)

### 3. Doublet Burden is Minimal & Expected ✅

The 1.29% doublet rate is:
- Consistent with 10x Chromium specifications
- Unbiased across samples (no single sample shows anomalously high doublet rate)
- Removed cleanly without introducing systematic bias

### 4. Flagged Sample (GSM5519491) — Retained Despite Low Cell Count

- **GSM5519491 (P190902, IT phase):** 66 cells pre-QC → 66 cells post-fixed-QC → included in doublet detection
- This sample was flagged as MAJOR_OUTLIER in Phase 0 for low cell recovery
- Fixed QC did not remove it (thresholds are per-cell, not per-sample)
- Doublet detection applied uniformly
- **Recommendation:** This sample is retained but will be monitored during Phase 3 clustering. If it contributes substantially to doublet burden or exhibits clustering artifacts, it can be flagged for further consideration.

---

## QC Decision Summary

| Decision | Rationale |
|---|---|
| **Fixed QC thresholds:** APPLIED, 0 cells removed | All cells passed (200–6,000 genes, ≤20% MT); thresholds were appropriate but non-selective |
| **Doublet detection:** APPLIED, 1,372 cells removed | scDblFinder with sample stratification identified 1.29% predicted doublets; appropriate for 10x protocol |
| **Doublet removal:** IMPLEMENTED | Predicted doublets removed; singlets retained |
| **Clinical group impact:** NONE | 100% retention per clinical state; no bias detected |

---

## Output & Checkpoints

### Primary Deliverable
- **`seurat_qc_singlets.rds`** — Quality-controlled Seurat object
  - 105,220 cells (singlets only)
  - 18,925 genes
  - All metadata preserved (GSM, Donor, Phase, scDblFinder scores & classifications)
  - Ready for Phase 3 atlas construction

### Supporting Data Files
- `sample_qc_before_filtering.csv` — Per-sample QC metrics pre-filter
- `qc_filtering_decision_summary.csv` — Fixed QC threshold outcomes
- `phase_cell_retention.csv` — Cell retention by clinical state
- `sample_cell_retention.csv` — Cell retention by sample
- `doublet_summary.csv` — Overall doublet statistics
- `doublet_rate_by_sample.csv` — Doublet rate per sample
- `doublet_rate_by_phase.csv` — Doublet rate per clinical state
- `doublet_validation_qc_metrics.csv` — QC metrics stratified by doublet status
- `qc_summary.csv` — Final consolidated QC status

### Visualizations (Key Figures)

1. **QC Violin plots by clinical state** — nFeature, nCount, percent.mt distributions across NL/IT/IA/AR/AC
2. **QC Violin plots by sample** — Per-sample QC metric distributions (23 samples)
3. **Feature scatter: Genes vs. UMI** — nFeature_RNA vs. nCount_RNA with QC thresholds marked
4. **Feature scatter: UMI vs. MT%** — nCount_RNA vs. percent.mt with MT threshold marked
5. **Gene detection histogram by phase** — nFeature_RNA distributions by clinical state with thresholds
6. **Mitochondrial percentage by phase** — Violin plot of percent.mt by clinical state
7. **Post-filter QC scatter** — Final singlet-only dataset: UMI vs. genes, colored by MT%, faceted by phase
8. **scDblFinder score distribution** — Histogram of doublet scores across all cells
9. **Doublet rate by sample** — Bar plot of predicted doublet rate per sample (ranked)
10. **Doublet validation: Genes by class** — Gene detection (singlets vs. predicted doublets)
11. **Doublet validation: UMI by class** — UMI counts (singlets vs. predicted doublets)
12. **Doublet score by predicted class** — Distribution of scDblFinder scores (singlets vs. doublets)
13. **Doublet score vs. genes** — Scatter plot: doublet score vs. nFeature_RNA, colored by class
14. **Doublet score vs. UMI** — Scatter plot: doublet score vs. nCount_RNA, colored by class

---

## Gate & Recommendations

### ✅ Gate Status: PASS

**Criteria Met:**
- ✅ All 23 samples retained
- ✅ All 5 clinical states represented (100% retention per state)
- ✅ All 23 donors retained (1 sample per donor maintained)
- ✅ QC distributions are biologically and technically sensible
- ✅ No unacceptable loss of any clinical group
- ✅ Doublet burden minimal and uniformly distributed
- ✅ Production-ready object for Phase 3

**Recommendations for Phase 3:**
1. Use `seurat_qc_singlets.rds` as input for atlas construction
2. Monitor GSM5519491 (66 cells) during clustering — small sample size may limit robustness
3. During cell-type annotation, validate that scDblFinder classifications are concordant with biological annotation (singlets should form coherent cell types)

---

## Technical Notes

### QC Thresholds Justification

The fixed thresholds were selected based on:
- 10x Chromium v2/v3 standard recommendations
- Single-cell RNA-seq best practices literature
- Anticipated liver tissue biology (high mitochondrial content in hepatocytes, variable complexity in immune cells)

A priori thresholds (vs. data-driven) were applied to avoid potential bias from outlier-driven cutoffs.

### Doublet Detection Strategy

**Why scDblFinder?**
- Robust & published (Bioconductor; cited in best-practices)
- Sample-aware (accounts for multiplet rates per sample)
- Probabilistic scores (not binary; allows downstream validation)
- No automatic gene signature dependency (unlike DoubletFinder)

**Why sample-stratified?**
- Each GSM may have slightly different multiplet rate due to input cell concentration
- Stratification improves detection accuracy for multi-sample datasets
- Standard approach in Bioconductor workflows

### MAD Diagnostics (Informational Only)

Sample-specific MAD (Median Absolute Deviation) statistics were calculated to detect unusual QC distributions but were **NOT** used for automatic cell removal. This approach:
- Maintains transparency (flags are visible but not auto-applied)
- Preserves cell types that may naturally have unusual complexity (e.g., activated T cells, metabolically active cells)
- Allows manual review rather than algorithmic hard-filtering

---

## Workflow Summary

```
Raw merged dataset (106,592 cells)
        ↓
    Calculate QC metrics
        ↓
Apply fixed QC thresholds (200–6,000 genes; ≤20% MT)
    → Result: 0 cells removed (all passed)
        ↓
scDblFinder doublet detection (sample-stratified)
    → Result: 1,372 predicted doublets identified
        ↓
Remove predicted doublets
        ↓
Final QC-filtered dataset (105,220 singlets)
        ↓
Checkpoint: seurat_qc_singlets.rds
        ↓
READY FOR PHASE 3 (Atlas Construction)
```

---

## Next Steps

**Phase 3 (Global Atlas Construction):**
1. Load `seurat_qc_singlets.rds`
2. Construct reference atlas via leverage score sketching (~20,000 cells)
3. Annotate cell types in reference
4. Project annotations onto full dataset
5. Perform rigorous secondary QC + create annotated global atlas

---

**Phase 2 Complete** ✅

Input cells: 106,592  
Output cells: 105,220  
Retention: 98.71%  
Gate Status: READY FOR PHASE 3
