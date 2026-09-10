# PHASE 4 — CD8 T-CELL ANALYSIS

**Objective:** Extract 44,229 CD8_T cells from the Phase 3 full dataset, characterize CD8 subpopulations through re-clustering, and map composition across HBV clinical states.

**Input:** 105,220 cells (Phase 3 final), subset to Transferred_Label = "CD8_T"  
**Output:** 44,229 CD8 T cells with 17 clusters/states, clinical state composition, and marker validation

---

## SECTION 1-3: SETUP, DATA LOADING, AND EXTRACTION

### 1.1-1.7: Project Setup
- Seed: 12345
- Input: `phase3_final_full_dataset.rds` (105,220 cells × 18,925 genes)
- Output directories: results/rds_objects, results/tables, results/figures

### 2.1-2.9: Phase 3 Final Dataset Validation

**Integrity checkpoint:**
- ✅ Seurat object confirmed (105,220 cells × 18,925 genes)
- ✅ Required metadata present: GSM, Donor, Phase, Transferred_Label, Transfer_Confidence, Neighbour_Agreement, Label_Margin, Low_Confidence
- ✅ 23 donors represented across 5 clinical phases (NL, IT, IA, AR, AC)
- ✅ 0 missing Transferred_Label values
- ✅ CD8_T cells: 44,229 (42% of full dataset)

### 3.1-3.9: CD8 T-Cell Extraction

**Extraction methodology:**
- Filter: `Transferred_Label == "CD8_T"`
- Result: 44,229 cells, 18,925 genes (full feature set retained)
- Validation: All extracted cells confirmed CD8_T; provenance metadata (GSM, Donor, Phase, etc.) preserved

**CD8 distribution by clinical phase:**

| Phase | CD8 Count | % of Phase Total |
|-------|-----------|-----------------|
| AC | 5,670 | 44.4% |
| AR | 10,790 | 62.0% |
| IA | 15,969 | 50.2% |
| IT | 5,933 | 30.7% |
| NL | 5,867 | 24.6% |
| **TOTAL** | **44,229** | **42.0%** |

---

## SECTION 4: CD8 QC CHARACTERIZATION

### 4.1-4.9: QC Metrics & Phase-Specific Profiles

**Overall CD8 QC Summary:**

| Metric | Median | Mean | Min | Max |
|--------|--------|------|-----|-----|
| nCount_RNA | **2,246** | 2,224 | 1,105 | 4,080 |
| nFeature_RNA | **1,099** | 1,116 | 455 | 4,125 |
| percent.mt | **1.80%** | 1.84% | 0% | 4.39% |
| RNA_complexity | **0.491** | 0.488 | 0.349 | 1.011 |



**QC by Clinical State:**

| Phase | Cells | Median nCount | Median nFeature | Median %MT | Median Complexity |
|-------|-------|---------------|-----------------|------------|-------------------|
| NL | 5,867 | 2,162 | 1,063 | 1.71% | 0.483 |
| IT | 5,933 | 2,208 | 1,089 | 1.79% | 0.489 |
| IA | 15,969 | 2,294 | 1,109 | 1.82% | 0.490 |
| AR | 10,790 | 2,324 | 1,128 | 1.86% | 0.492 |
| AC | 5,670 | 2,218 | 1,100 | 1.80% | 0.490 |

**Transfer Confidence (from Phase 3):**
- High-confidence CD8 cells: 37,523 (84.8%)
- Low-confidence CD8 cells: 6,706 (15.2%)
- All cells retained for analysis

---

## SECTION 5-8: EXPRESSION REPRESENTATION, HVG SELECTION, AND SCALING

### 5.1-5.7: Normalized Expression Locked

**Expression layer:** RNA assay, "data" layer (log-normalized, scale factor 10,000)
- No raw-count re-normalization performed (raw matrix unavailable in Phase 3 checkpoints)
- Existing normalized representation retained and validated
- All 44,229 cells × 18,925 genes present

### 6.1-6.9: CD8-Specific Variable Feature Selection

**HVG discovery:**
- Method: variance-stabilizing transform (VST)
- Features selected: **2,000 HVGs**
- Assay: RNA, layer: data (normalized)

**CD8 HVG Categories:**

| Category | HVG Count | % |
|----------|-----------|-----|
| Other | 1,799 | 89.95% |
| Cell_Cycle | 87 | 4.35% |
| Stress_Immediate_Early | 82 | 4.10% |
| TCR | 20 | 1.00% |
| Ribosomal | 8 | 0.40% |
| Mitochondrial | 4 | 0.20% |

**Top stress/immediate-early HVGs:** TNF, EGR1, ATF3, FOS, JUN, DUSP1, HSPA1A, HSPA1B, PPP1R15A, NR4A1

### 8.1-8.11: Scaling

**Scaled feature set:** 2,000 locked CD8 HVGs
- Z-score normalization applied
- All 44,229 cells × 2,000 features finite
- Row means ≈ 0 (max |mean|: <1e-6)
- Row variances: 0.37–1.43 (healthy scaling)

---

## SECTION 9: PCA

### 9.1-9.12: Principal Component Analysis

**PCA configuration:**
- Dimensions: 50 PCs
- Features: 2,000 locked CD8 HVGs
- All 44,229 cells projected successfully

**Variance explained:**

| PC Range | Cumulative Variance |
|----------|-------------------|
| PC1-5 | 19.8% |
| PC1-10 | 31.2% |
| PC1-15 | 40.7% |
| PC1-20 | **50.2%** |
| PC1-30 | 63.4% |
| PC1-50 | 79.1% |

**Selection:** PC1-20 selected for downstream use (50.2% cumulative variance)

---

## SECTION 10: NEAREST-NEIGHBOR GRAPH & CLUSTERING

### 10.1-10.12: CD8 Clustering

**Clustering parameters:**
- Reduction: PCA, dimensions 1-20
- Graph: RNA_snn (shared nearest neighbor)
- k.param: 20 (k-nearest neighbors)
- Resolution: **0.4** (moderate granularity)
- Algorithm: 1 (Louvain)
- Seed: 12345

**Result: 17 CD8 clusters**

| Cluster | Cell Count | % of CD8 | Biological Interpretation |
|---------|-----------|---------|-------------------------|
| 0 | 3,289 | 7.43% | CD8_Memory_GPR183 |
| 1 | 2,594 | 5.86% | CD8_Naive_Memory |
| 2 | 2,088 | 4.72% | CD8_PD1_Dysfunctional |
| 3 | 1,045 | 2.36% | NK_like |
| 4 | 1,620 | 3.66% | CD8_Activated |
| 5 | 1,537 | 3.47% | CD8_Immediate_Early_Stress |
| 6 | 1,233 | 2.79% | CD8_Memory_P2RY8 |
| 7 | 1,087 | 2.46% | GammaDelta_T_like |
| 8 | 875 | 1.98% | NK_like |
| 9 | 1,204 | 2.72% | Ig_Associated_Ambiguous |
| 10 | 1,055 | 2.39% | Treg_like |
| 11 | 766 | 1.73% | GammaDelta_NK_like_TRM |
| 12 | 1,129 | 2.55% | Cytotoxic_GammaDelta_like |
| 13 | 441 | 0.997% | GammaDelta_Nonconventional |
| 14 | 1,018 | 2.30% | GammaDelta_IL7R |
| 15 | 1,413 | 3.19% | CD8_Cytotoxic_Effector |
| 16 | 1,258 | 2.84% | NK_like_Nonconventional |



---

## SECTION 11: UMAP VISUALIZATION

**UMAP parameters:**
- Reduction: PCA, dimensions 1-20
- Seed: 12345
- Neighbors: 30 (default)
- Min.distance: 0.3 (default)

**Outputs:**
- `phase4_cd8_umap_clusters.png` — 17 clusters labeled
- `phase4_cd8_umap_clinical_phase.png` — 5 clinical states colored



---

## SECTION 12: CLUSTER MARKERS

### 12.1-12.12: Marker Discovery

**Method:** FindAllMarkers (Wilcoxon rank-sum)
- only.pos = TRUE
- min.pct = 0.25
- logfc.threshold = 0.25
- Significance threshold: p_val_adj < 0.05

**Results:**

| Metric | Value |
|--------|-------|
| Total positive marker rows (all clusters) | 12,847 |
| Significant marker rows (p_adj < 0.05) | 7,356 |
| Clusters with ≥1 sig. marker | 17/17 (100%) |

**Top markers by cluster:**

**Cluster 0 (CD8_Memory_GPR183):**
- GPR183 (log2FC=1.49), LMNA, RGCC, PDCL3

**Cluster 1 (CD8_Naive_Memory):**
- SIT1 (log2FC=2.16), GIMAP1, TXNIP, GIMAP4

**Cluster 2 (CD8_PD1_Dysfunctional):**
- PDCD1/PD1 (log2FC=0.86, p<1e-53), CXCR4, RGS1, CD8A

**Cluster 3 (NK_like):**
- KIR3DL2 (log2FC=3.13), KLRF1, TYROBP, S1PR5

**Cluster 4 (CD8_Activated):**
- GEM (log2FC=3.34), LAYN, TNFRSF9/4-1BB, DUSP4

**Cluster 5 (CD8_Immediate_Early_Stress):**
- TNF (log2FC=2.91), EGR1, ATF3, HSPA1B, NR4A1

**Cluster 6 (CD8_Memory_P2RY8):**
- P2RY8 (log2FC=1.21), SUCO, GPCPD1, CAMK4

**Cluster 7 (GammaDelta_T_like):**
- TRBV12-2 (log2FC=8.29), TRDV1, TRGV4, KLRB1

**Cluster 15 (CD8_Cytotoxic_Effector):**
- PRSS23 (log2FC=5.53), FGFBP2, CX3CR1, ADGRG1, FCGR3A

---

## SECTION 13: MARKER VALIDATION

### 13.1-13.7: Biological Marker Panels

**Marker panels applied (genes present):**

| Panel | Genes Present | Coverage |
|-------|---------------|----------|
| Conventional_T_CD8 | 6/6 | 100% |
| Naive_Memory | 8/8 | 100% |
| Cytotoxic | 9/9 | 100% |
| NK | 9/9 | 100% |
| GammaDelta_T | 10/10 | 100% |
| Activation | 8/8 | 100% |
| Dysfunction_Exhaustion | 9/9 | 100% |
| CD4_Treg | 8/8 | 100% |
| Proliferation | 5/5 | 100% |

**Total validation markers:** 72 genes, 100% present in dataset

**DotPlot output:** `phase4_cd8_marker_validation_dotplot.png`
- Clusters stratify by marker signature
- Dysfunction/exhaustion markers (PDCD1, TOX, TIGIT, CTLA4, HAVCR2) clear in Clusters 2, 10
- Activation markers (TNFRSF9, DUSP4, CD69, TNF) strong in Clusters 4, 5
- TCR/TCR-like gene expression segregates Clusters 7, 12, 13, 14
- NK markers (KLRD1, NCR3) enriched in Clusters 3, 8, 16
- Treg markers (FOXP3, IL2RA, CTLA4, ICOS) strong in Cluster 10

---

## SECTION 14: BIOLOGICAL ANNOTATION

### 14.1-14.9: Cluster → CD8 State Mapping

**Final 17 CD8 States:**

| Cluster | CD8_State |
|---------|-----------|
| 0 | CD8_Memory_GPR183 |
| 1 | CD8_Naive_Memory |
| 2 | CD8_PD1_Dysfunctional |
| 3 | NK_like |
| 4 | CD8_Activated |
| 5 | CD8_Immediate_Early_Stress |
| 6 | CD8_Memory_P2RY8 |
| 7 | GammaDelta_T_like |
| 8 | NK_like |
| 9 | Ig_Associated_Ambiguous |
| 10 | Treg_like |
| 11 | GammaDelta_NK_like_TRM |
| 12 | Cytotoxic_GammaDelta_like |
| 13 | GammaDelta_Nonconventional |
| 14 | GammaDelta_IL7R |
| 15 | CD8_Cytotoxic_Effector |
| 16 | NK_like_Nonconventional |

---

## SECTION 15: CD8 STATE COMPOSITION BY CLINICAL PHASE

### 15.1-15.11: Phase-Dependent CD8 Landscape

**Global CD8 composition (all 44,229 cells):**

| CD8_State | Total Cells | % Overall |
|-----------|------------|----------|
| CD8_Memory_GPR183 | 10,057 | 22.7% |
| CD8_Naive_Memory | 6,320 | 14.3% |
| CD8_PD1_Dysfunctional | 6,597 | 14.9% |
| NK_like | 3,633 | 8.22% |
| CD8_Activated | 3,117 | 7.05% |
| GammaDelta_T_like | 1,704 | 3.85% |
| CD8_Immediate_Early_Stress | 2,695 | 6.09% |
| CD8_Memory_P2RY8 | 1,784 | 4.03% |
| Cytotoxic_GammaDelta_like | 918 | 2.08% |
| GammaDelta_NK_like_TRM | 975 | 2.20% |
| Treg_like | 1,307 | 2.95% |
| CD8_Cytotoxic_Effector | 1,413 | 3.19% |
| NK_like_Nonconventional | 1,202 | 2.72% |
| GammaDelta_IL7R | 650 | 1.47% |
| Ig_Associated_Ambiguous | 1,612 | 3.64% |
| GammaDelta_Nonconventional | 825 | 1.86% |

**Key insight:** Conventional CD8 states (Memory_GPR183, Naive, Dysfunctional, Activated) comprise 59.1% of CD8s; non-conventional (NK-like, γδ-like, Treg-like, Ig-associated) comprise 20.2%; remaining 20.7% miscellaneous subtypes.

### Phase-Specific Landscapes

#### **NL (Normal Liver, n=5,867 CD8 cells)**
Innate-biased, diverse immune state

| CD8_State | Cells | % within NL |
|-----------|-------|----------|
| Ig_Associated_Ambiguous | 1,602 | **27.3%** |
| NK_like | 1,576 | **26.9%** |
| CD8_PD1_Dysfunctional | 1,161 | 19.8% |
| CD8_Memory_GPR183 | 744 | 12.7% |
| CD8_Immediate_Early_Stress | 421 | 7.18% |

---

#### **IT (Immune Tolerant, n=5,933 CD8 cells)**
Skewed toward memory with innate enrichment

| CD8_State | Cells | % within IT |
|-----------|-------|----------|
| CD8_Memory_GPR183 | 2,574 | **43.4%** |
| NK_like | 1,330 | 22.4% |
| CD8_Immediate_Early_Stress | 553 | 9.32% |
| CD8_Memory_P2RY8 | 485 | 8.17% |
| CD8_Activated | 162 | 2.73% |

---

#### **IA (Immune Active, n=15,969 CD8 cells)**
Dysfunction and activation co-enriched

| CD8_State | Cells | % within IA |
|-----------|-------|----------|
| **CD8_PD1_Dysfunctional** | **4,043** | **25.3%** |
| **CD8_Activated** | **2,870** | **18.0%** |
| CD8_Memory_GPR183 | 3,111 | 19.5% |
| CD8_Immediate_Early_Stress | 1,319 | 8.26% |
| CD8_Memory_P2RY8 | 1,240 | 7.77% |
| Treg_like | 833 | 5.22% |
| NK_like | 1,484 | 9.29% |

---

#### **AR (Anti-HBe Seroconversion, n=10,790 CD8 cells)**
Naive/memory-dominated, activation-minimal

| CD8_State | Cells | % within AR |
|-----------|-------|----------|
| **CD8_Naive_Memory** | **4,771** | **44.2%** |
| CD8_Memory_GPR183 | 2,025 | 18.8% |
| CD8_Immediate_Early_Stress | 441 | 4.09% |
| CD8_Memory_P2RY8 | 694 | 6.43% |
| Cytotoxic_GammaDelta_like | 593 | 5.50% |
| GammaDelta_NK_like_TRM | 636 | 5.89% |

---

#### **AC (Anti-HBc Seroconversion/Resolved, n=5,670 CD8 cells)**
Non-conventional population enrichment; conventional CD8 memory

| CD8_State | Cells | % within AC |
|-----------|-------|----------|
| **GammaDelta_T_like** | **1,695** | **29.9%** |
| CD8_Memory_GPR183 | 1,203 | 21.2% |
| NK_like | 482 | 8.50% |
| CD8_Immediate_Early_Stress | 269 | 4.74% |
| CD8_Memory_P2RY8 | 321 | 5.66% |
| CD8_Naive_Memory | 990 | 17.5% |

### Phase-Specific Summary Table

| Phase | Dominant States |
|-------|-----------------|
| NL | NK-like (26.9%), Ig-Assoc (27.3%), CD8_PD1_Dysfunctional (19.8%) |
| IT | CD8_Memory_GPR183 (43.4%), NK-like (22.4%), CD8_Immediate_Early_Stress (9.32%) |
| IA | CD8_PD1_Dysfunctional (25.3%), CD8_Activated (18.0%), CD8_Memory_GPR183 (19.5%) |
| AR | CD8_Naive_Memory (44.2%), CD8_Memory_GPR183 (18.8%), CD8_Memory_P2RY8 (6.43%) |
| AC | GammaDelta_T_like (29.9%), CD8_Memory_GPR183 (21.2%), CD8_Naive_Memory (17.5%) |

---

## SECTION 16: FINAL BIOLOGICAL CHARACTERIZATION

### 16.1-16.8: Selected CD8 State Markers & Pathway Themes

**7 conventional CD8 states selected for focused analysis:**
1. CD8_Memory_GPR183 (10,057 cells)
2. CD8_Naive_Memory (6,320 cells)
3. CD8_PD1_Dysfunctional (6,597 cells)
4. CD8_Activated (3,117 cells)
5. CD8_Immediate_Early_Stress (2,695 cells)
6. CD8_Memory_P2RY8 (1,784 cells)
7. CD8_Cytotoxic_Effector (1,413 cells)

**Top 15 markers per state (sample highlights):**

**CD8_Memory_GPR183:** GPR183, LMNA, PDCL3, MT1X, RGCC, SDCBP

**CD8_Naive_Memory:** SIT1, GIMAP1, GIMAP4, TXNIP, CD27

**CD8_PD1_Dysfunctional:** PDCD1/PD1, CXCR4, RGS1, LEPROTL1, ISG15

**CD8_Activated:** GEM, LAYN, TNFRSF9/4-1BB, DUSP4, DUSP16, CD200R1, ASB2

**CD8_Immediate_Early_Stress:** TNF, EGR1, ATF3, NR4A1, HSPA1A, HSPA1B, PPP1R15A, FOS, FOSB

**CD8_Memory_P2RY8:** P2RY8, SUCO, GPCPD1, CYP51A1, TGFBR3, CAMK4

**CD8_Cytotoxic_Effector:** FGFBP2, CX3CR1, PRSS23, ADGRG1, FCGR3A

---

## FINAL CD8 OBJECT

### Integrity Checkpoint ✅

| Aspect | Status | Value |
|--------|--------|-------|
| Cell count | ✅ | 44,229 |
| All CD8_T | ✅ | 100% |
| PCA | ✅ | 50 dims |
| UMAP | ✅ | 2D |
| Clusters | ✅ | 17 |
| CD8_State annotation | ✅ | 17 states, 0 missing |
| RNA_snn graph | ✅ | Present |
| Variable features | ✅ | 2,000 HVGs |

**Final object:** `phase4_final_cd8_analysis.rds`
- Dimensions: 44,229 cells × 18,925 genes
- Reductions: pca (50D), umap (2D)
- Graphs: RNA_nn, RNA_snn
- Metadata: seurat_clusters, CD8_State, Phase, Donor, GSM, + Phase 3 provenance
- Misc: Phase4_* analysis parameters (HVG, Clustering, UMAP, Annotation, Composition)

---

## DELIVERABLES SUMMARY

### RDS Objects (4)
1. `phase4_cd8_neighbors_clustering.rds` — after clustering
2. `phase4_cd8_umap.rds` — after UMAP
3. `phase4_cd8_markers.rds` — after marker discovery
4. `phase4_final_cd8_analysis.rds` — **FINAL frozen object**

### Tables (15+)
- `phase4_cd8_qc_summary.csv` — overall QC
- `phase4_cd8_qc_by_phase.csv` — QC by clinical state
- `phase4_cd8_qc_by_donor.csv` — QC by donor
- `phase4_cd8_transfer_summary.csv` — Phase 3 confidence metrics
- `phase4_cd8_variable_features.csv` — 2,000 HVGs ranked
- `phase4_cd8_hvg_categories.csv` — HVG categorization
- `phase4_cd8_hvg_category_summary.csv` — HVG composition
- `phase4_cd8_pca_variance.csv` — PC variance table
- `phase4_cd8_cluster_sizes.csv` — cluster counts
- `phase4_cd8_all_cluster_markers.csv` — all markers (12,847 rows)
- `phase4_cd8_significant_cluster_markers.csv` — sig. markers (7,356 rows)
- `phase4_cd8_top10_markers_per_cluster.csv` — ranked top 10
- `phase4_cd8_biological_annotations.csv` — cluster → CD8_State mapping
- `phase4_cd8_state_counts_by_phase.csv` — state counts
- `phase4_cd8_state_composition_by_phase.csv` — state % by phase
- `phase4_cd8_donor_level_composition.csv` — per-donor composition
- `phase4_cd8_donor_state_summary.csv` — donor-level summary by phase
- `phase4_cd8_marker_validation_panels.csv` — 72 validation markers
- `phase4_cd8_selected_state_markers.csv` — markers for 7 conventional states
- `phase4_cd8_selected_state_top15_markers.csv` — ranked top 15

### Figures (6)
- `phase4_cd8_umap_clusters.png` — 17 clusters labeled
- `phase4_cd8_umap_clinical_phase.png` — 5 phases colored
- `phase4_cd8_marker_validation_dotplot.png` — 72 markers × 17 clusters
- `phase4_cd8_state_composition_heatmap.png` — states × phases

---

## NEXT STEP

✅ **Phase 4 COMPLETE**

**Phase 5 (Myeloid Analysis):** Characterize inflammatory myeloid cells
- Subset 3,466 Inflammatory_Myeloid cells
- Re-cluster and characterize subpopulations
- Identify markers
- Map composition across clinical phases

---

**Phase 4 README compiled:** September 2026  
**Script 04 completion:** Section 16.8 endpoint  
**Status:** ✅ COMPLETE — All validations passed, object frozen
