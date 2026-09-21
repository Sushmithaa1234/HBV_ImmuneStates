# PHASE 4 — CD8-LABELLED COMPARTMENT ANALYSIS

**Objective:** Extract the Phase 3 CD8-labelled compartment from the full liver dataset, characterize its transcriptional heterogeneity through CD8-specific re-clustering and marker analysis, and perform biologically conservative annotation of the resulting populations.

**Input:** `results/rds_objects/phase3_final_full_dataset.rds`

**Output:** 41,281 CD8-labelled cells, 16 computational clusters, and 16 cluster-level biological annotations across 23 donors and 5 GEO-defined clinical states.

---

## SCIENTIFIC SCOPE

Phase 4 characterizes transcriptional heterogeneity within the cells labelled CD8_T by the Phase 3 global atlas.

The analysis is intentionally framed as **CD8-labelled compartment characterization**, rather than assuming that every transferred CD8_T cell represents a conventional CD8 T-cell state.

### Phase 4 is designed to answer

- What transcriptional populations are present within the Phase 3 CD8-labelled compartment?
- Which marker programs characterize those populations?
- Which populations show conventional CD8-associated versus nonconventional lineage-associated signatures?
- How broadly are the identified populations represented across donors?
- What biological annotations can be supported by the available transcriptional evidence?

### Phase 4 does NOT attempt to

- Infer clinical-state effects
- Establish disease progression or trajectory
- Make causal claims
- Identify HBV-specific T-cell states
- Claim definitive cellular lineage from a single marker
- Classify populations as pathogenic or protective
- Establish exhaustion as a functional state
- Perform donor-level statistical inference between clinical states
- Replace the Phase 3 global cell-type annotation

**Note:** Clinical-state inference is reserved for the donor-aware analyses of later phases.

---

## 1. INPUT DATA AND INTEGRITY VALIDATION

### 1.1 Project Setup

| Parameter | Value |
|-----------|-------|
| Input object | `results/rds_objects/phase3_final_full_dataset.rds` |
| Total cells | 106,592 |
| Total genes | 24,452 |
| Donors/samples | 23 |
| Seurat version | v5 |
| Includes Phase 3 annotations | Yes |

**Clinical-state labels:**
- AC
- AR
- IA
- IT
- NL

The Phase 4 analysis preserves these GEO-defined labels and does not reinterpret them as a progression trajectory.

### 1.2 Phase 3 Integrity Check

- ✅ Seurat object valid
- ✅ 106,592 cells present
- ✅ 24,452 genes present
- ✅ 23 donors represented
- ✅ 5 GEO-defined clinical states represented
- ✅ Phase 3 cell-type annotations present
- ✅ CD8_T population present
- ✅ Provenance metadata retained

---

## 2. CD8-LABELLED COMPARTMENT EXTRACTION

### 2.1 Extraction Rule

Cells were extracted using the Phase 3 final global annotation:

```
Final_Cell_Type == "CD8_T"
```

No additional biological filtering was applied at extraction.

### 2.2 Extraction Result

| Metric | Value |
|--------|-------|
| CD8-labelled cells | 41,281 |
| Genes | 24,452 |
| Donors | 23 |
| Clinical states | 5 |
| Computational clusters | 16 |

**Extraction integrity:**
- ✅ Exact cell-ID correspondence
- ✅ No unexpected cells introduced
- ✅ No extracted cells lost through ID mismatch
- ✅ Donor metadata preserved
- ✅ Clinical-state metadata preserved
- ✅ Phase 3 provenance preserved

---

## 3. EXPRESSION REPRESENTATION

### 3.1 Expression Layer

The Phase 3 dataset contains processed normalized expression rather than the original raw UMI count matrix. The Phase 4 analysis therefore uses the existing normalized expression representation.

| Aspect | Value |
|--------|-------|
| Assay | RNA |
| Data layers | Processed |
| Expression format | log-CP10K |
| Reverse transformation | No |
| Re-normalization | No |

The existing normalized representation was retained as the starting point for CD8 analysis.

### 3.2 Source-Layer Preservation

Before downstream manipulation, the sample-specific expression layers were preserved.

**Checkpoint:** `results/rds_objects/phase4_cd8_source_layers_preserved.rds`

---

## 4. CD8-SPECIFIC DIMENSION REDUCTION AND CLUSTERING

### 4.1 Highly Variable Features

| Parameter | Value |
|-----------|-------|
| HVGs selected | 2,000 |
| Purpose | Scaling and PCA |

### 4.2 Scaling

The 2,000 selected CD8 HVGs were scaled for dimensionality reduction. No biological filtering was performed on the basis of the scaled values.

### 4.3 Principal Component Analysis

| Parameter | Value |
|-----------|-------|
| Input features | 2,000 CD8 HVGs |
| PCs calculated | 50 |
| PCs for graph construction | PC1–PC20 |

PCA was performed on the CD8-specific expression space.

---

## 5. CD8-SPECIFIC NEIGHBOR GRAPH AND CLUSTERING

### 5.1 Graph Construction

| Parameter | Value |
|-----------|-------|
| Reduction | PCA |
| Dimensions | PC1–PC20 |
| k-value | 20 |
| Nearest-neighbor graph | CD8_nn |
| Shared-nearest-neighbor graph | CD8_snn |

### 5.2 Clustering

| Parameter | Value |
|-----------|-------|
| Graph | CD8_snn |
| Resolution | 0.4 |
| Algorithm | 1 (Louvain) |
| Resulting clusters | 16 |

**Cluster sizes:**

| Cluster | Cells |
|---------|-------|
| 0 | 8,765 |
| 1 | 6,100 |
| 2 | 5,632 |
| 3 | 3,545 |
| 4 | 3,362 |
| 5 | 3,106 |
| 6 | 1,942 |
| 7 | 1,703 |
| 8 | 1,530 |
| 9 | 1,098 |
| 10 | 995 |
| 11 | 958 |
| 12 | 906 |
| 13 | 778 |
| 14 | 544 |
| 15 | 317 |
| **Total** | **41,281** |

**Important:** These clusters represent computationally defined transcriptional populations. They are not automatically equivalent to stable biological cell states.

---

## 6. UMAP

A CD8-specific UMAP was generated from the CD8 PCA representation.

| Parameter | Value |
|-----------|-------|
| Reduction | PCA |
| Dimensions | PC1–PC20 |
| Distance metric | Cosine (UWOT) |
| Output file | `phase4_cd8_umap.png` |

UMAP is used for visualization of the transcriptional structure identified by PCA and graph-based clustering. UMAP coordinates are not treated as quantitative biological measurements.

---

## 7. MARKER DISCOVERY

### 7.1 Cluster-Level Marker Analysis

| Parameter | Value |
|-----------|-------|
| Statistical method | Wilcoxon rank-sum test |
| Purpose | Gene enrichment per cluster |

Marker analysis was exploratory and was not treated as donor-level clinical-state inference.

### 7.2 Marker Interpretation

Marker evidence was evaluated together with:

- Lineage-associated gene panels
- Transcriptional programs
- Cluster-level expression patterns
- Donor representation
- Competing lineage evidence

**Critical principle:** No single gene was treated as sufficient evidence for a definitive biological identity.

---

## 8. TRANSCRIPTIONAL PROGRAMS

Programs examined across the CD8-labelled compartment:

- Naive/Memory
- Cytotoxic
- Activation
- Dysfunction-associated
- Immediate-early response
- NK-like
- Nonconventional T-cell
- Tissue-associated
- Proliferation

**Important:** These are transcriptional programs, not automatically mutually exclusive cell states. An immediate-early response program may represent a transient transcriptional response rather than a stable cellular identity. Dysfunction-associated expression is not interpreted as proof of functional exhaustion.

---

## 9. LINEAGE-ASSOCIATED EVIDENCE

The following panels were used as supporting evidence.

### Conventional T/CD8-associated
`CD3D`, `CD3E`, `TRBC1`, `TRBC2`, `CD8A`, `CD8B`

### NK-associated
`NKG7`, `GNLY`, `KLRD1`, `KLRF1`, `KLRC1`, `KLRC2`, `TYROBP`, `FCER1G`, `NCR3`, `S1PR5`

### Gamma-delta T-associated
`TRDC`, `TRGC1`, `TRGC2`, `TRDV1`, `TRDV2`, `TRGV2`, `TRGV4`, `TRGV5`, `TRGV8`, `TRGV9`

### Treg-associated
`CD4`, `IL7R`, `FOXP3`, `IL2RA`, `CTLA4`, `TNFRSF4`, `TNFRSF18`, `ICOS`

### B-cell-associated
`CD79A`, `MS4A1`, `CD37`, `CD74`, `HLA-DRA`

### Myeloid-associated
`LYZ`, `S100A8`, `S100A9`, `FCGR3A`, `CTSD`

**Note:** These panels are used as evidence for annotation, not as automated cell-type classifiers.

---

## 10. BIOLOGICAL ADJUDICATION

Each computational cluster was evaluated using multiple evidence streams:

- Cluster-level marker expression
- Lineage-associated panels
- Transcriptional programs
- Competing lineage evidence
- Donor representation
- Degree of biological ambiguity

**Principle:** Annotations were deliberately conservative. Terms such as "associated," "like," and "population" are used where the available transcriptomic evidence does not justify a definitive lineage or functional claim.

---

## 11. FINAL CLUSTER-LEVEL BIOLOGICAL ANNOTATIONS

The 16 computational clusters were assigned 16 cluster-level annotations.

**Important:** 16 clusters do not imply 16 unique biological categories.

| Cluster | Final Biological Annotation | Annotation Category | Confidence |
|---------|----------------------------|----------------------|------------|
| 0 | GZMK-associated memory-like CD8 T-cell | Conventional_CD8_associated | Moderate |
| 1 | Memory-associated conventional T-cell | Conventional_T_associated | Moderate |
| 2 | Activated/dysfunction-associated CD8 T-cell | Conventional_CD8_associated | High |
| 3 | Immediate-early-response CD8-associated population | Transcriptional_program | High |
| 4 | Activated/dysfunction-associated CD8 T-cell | Conventional_CD8_associated | High |
| 5 | Activation/tissue-associated CD8-associated population | Conventional_CD8_associated | Moderate |
| 6 | NK-like cytotoxic nonconventional T-cell population | Noncanonical_T_associated | High |
| 7 | Gamma-delta-like cytotoxic T-cell population | GammaDelta_associated | High |
| 8 | NK-like/nonconventional cytotoxic T-cell population | Noncanonical_T_associated | High |
| 9 | NK-like cytotoxic population | NK_associated | High |
| 10 | NK/gamma-delta-like nonconventional T-cell population | Noncanonical_T_associated | High |
| 11 | Gamma-delta-like T-cell population | GammaDelta_associated | High |
| 12 | Activated/dysfunction-associated tissue-associated CD8 T-cell | Conventional_CD8_associated | High |
| 13 | Gamma-delta-like cytotoxic T-cell population | GammaDelta_associated | High |
| 14 | NK-like cytotoxic population | NK_associated | High |
| 15 | Gamma-delta-like nonconventional T-cell population | GammaDelta_associated | High |

**Annotation summary:**
- 16/16 computational clusters annotated
- 13 unique biological annotation categories
- No cluster left without a final annotation
- No claim that every population represents a canonical CD8 T-cell state

---

## 12. DONOR REPRESENTATION

Donor representation was explicitly evaluated because the 23 donors, rather than individual cells, represent the biological replicates.

### Examples of donor concentration

| Cluster | Largest donor representation |
|---------|------------------------------|
| 7 | ~97.6% |
| 13 | ~99.6% |
| 14 | ~96.1% |
| 9 | ~91.6% |
| 1 | ~77% |
| 12 | ~79.8% |

**Interpretation:** These clusters are retained because donor concentration does not invalidate their existence as observed transcriptional populations. However, donor concentration limits the strength of cohort-level biological generalization. Consequently, the presence of a cluster is not interpreted as evidence that the corresponding population is uniformly present across all clinical states or donors.

---

## 13. CD8 COMPARTMENT INTERPRETATION

The Phase 4 analysis demonstrates that the Phase 3 CD8-labelled compartment is transcriptionally heterogeneous.

### Conventional CD8/T-associated populations

Including populations characterized by:
- Memory-associated programs
- Activation-associated programs
- Dysfunction-associated programs
- Tissue-associated programs
- Immediate-early response programs

### Nonconventional T-cell-associated populations

Including:
- NK-like cytotoxic populations
- Gamma-delta-like populations
- Mixed NK/gamma-delta-like populations

**Why this matters:** The Phase 3 transfer label describes the cells' position within the global atlas, whereas Phase 4 provides a more detailed characterization of their transcriptional structure.

---

## 14. IMPORTANT INTERPRETATION BOUNDARIES

### No exhaustion claim
"Dysfunction-associated" transcriptional evidence is not treated as proof of functional exhaustion.

### No tissue-resident memory claim
"Tissue-associated" expression is not treated as definitive evidence of TRM identity.

### No HBV-specificity claim
The analysis does not establish antigen specificity or HBV-specific T-cell receptor activity.

### No pathogenic/protective claim
The analysis does not determine whether any population is pathogenic, protective, beneficial, or harmful.

### No clinical-state inference
Differences in cell composition across clinical states are not interpreted as statistically supported clinical-state effects in Phase 4. Those questions require donor-aware analyses in later phases.

### No trajectory/progression model
The five GEO-defined clinical states are treated as distinct clinical groups. Phase 4 does not infer a trajectory among them.

### No causal inference
All observations are descriptive and exploratory.

---

## 15. FINAL OBJECT

### Integrity Checkpoint

| Aspect | Status | Value |
|--------|--------|-------|
| Cell count | ✅ | 41,281 |
| Feature count | ✅ | 24,452 |
| Donors | ✅ | 23 |
| Clinical states | ✅ | 5 |
| Computational clusters | ✅ | 16 |
| Biological annotations | ✅ | 16/16 clusters |
| Unique biological categories | ✅ | 13 |
| PCA | ✅ | 50 dimensions |
| CD8 UMAP | ✅ | Present |
| CD8 neighbor graph | ✅ | Present |
| CD8 SNN graph | ✅ | Present |
| Variable features | ✅ | 2,000 HVGs |
| Missing biological annotations | ✅ | 0 |
| Seurat object validity | ✅ | Passed |
| Cell-ID integrity | ✅ | Passed |

### FINAL FROZEN OBJECT

**Primary final object:**
```
results/rds_objects/phase4_final_cd8_analysis_adjudicated.rds
```

**Dimensions:** 41,281 cells × 24,452 genes

**Reductions:**
- `pca`
- `umap_cd8`

**Graphs:**
- `CD8_nn`
- `CD8_snn`

**Key metadata:**
- `GSM`
- `Donor`
- `Phase`
- `Final_Cell_Type`
- `Transfer_Confidence`
- `Neighbour_Agreement`
- `Low_Confidence`
- `CD8_Cluster`
- `CD8_Biological_Annotation`

---

## 16. PHASE 4 OUTPUTS

### RDS Objects

| Object | Purpose |
|--------|---------|
| `phase4_cd8_source_layers_preserved.rds` | Source-layer checkpoint created before joining CD8 normalized data layers |
| `phase4_final_cd8_analysis_adjudicated.rds` | Frozen final Phase 4 object |

### Tables

- `phase4_cd8_cluster_sizes.csv`
- `phase4_cd8_all_cluster_markers.csv`
- `phase4_cd8_significant_cluster_markers.csv`
- `phase4_cd8_top20_markers_per_cluster.csv`
- `phase4_cd8_top50_markers_per_cluster.csv`
- `phase4_cd8_program_definitions.csv`
- `phase4_cd8_program_scores.csv`
- `phase4_cd8_lineage_panel_summary.csv`
- `phase4_cd8_annotation_adjudication.csv`
- `phase4_cd8_final_biological_annotations.csv`
- `phase4_cd8_donor_annotation_representation.csv`
- `phase4_cd8_donor_cluster_proportions.csv`
- `phase4_cd8_donor_cluster_summary.csv`
- `phase4_cd8_summary.txt`

**Note:** Exact filenames may vary according to the final script output directory.

### Figures

Key visualization outputs include:
- CD8 cluster UMAP
- CD8 biological annotation UMAP
- CD8 marker heatmap
- CD8 marker/program visualizations
- Donor representation summaries

**Final biological annotation figure:** `phase4_cd8_final_biological_annotations.png`

---

## 17. REPRODUCIBILITY

The Phase 4 analysis is implemented in:
```
04_phase4_cd8_analysis.R
```

The script contains:
- Input validation
- CD8 extraction
- Source-layer preservation
- Normalized-expression handling
- HVG selection
- Scaling
- PCA
- Neighbor graph construction
- Clustering
- UMAP
- Marker discovery
- Transcriptional program scoring
- Lineage evidence assessment
- Annotation adjudication
- Donor representation analysis
- Final integrity validation
- Final RDS export

---

## 18. BIOLOGICAL ROLE OF PHASE 4 IN THE OVERALL PROJECT

Phase 4 provides the detailed transcriptional characterization of the CD8-labelled compartment identified in the global Phase 3 atlas.

It establishes the populations and programs that will later be examined using donor-aware analyses.

**Key methodological principle carried forward:**

> Cells are used to discover and characterize transcriptional states; donors are the biological replicates for clinical-state inference.

Therefore, Phase 4 does not treat thousands of individual cells as independent biological replicates.

---

## 19. NEXT STEP

### Phase 5 — MYELOID ANALYSIS

Phase 5 will characterize the myeloid compartment identified in the Phase 3 global atlas.

**Planned objectives:**
- Extract the Phase 3 Monocyte_Myeloid compartment
- Preserve donor and clinical-state provenance
- Perform myeloid-specific dimensionality reduction and clustering
- Characterize transcriptional heterogeneity
- Identify inflammatory/macrophage-associated populations
- Evaluate lineage and state-associated marker programs
- Assess donor representation
- Produce conservative biological annotations

As with Phase 4, Phase 5 will distinguish:
- Computational clusters → Transcriptional programs → Biological interpretation

Rather than treating every computational cluster as an automatically defined biological cell state.

---

## PHASE 4 STATUS

### ✅ PHASE 4 COMPLETE

**Final checkpoint:**
- 41,281 CD8-labelled cells
- 24,452 genes
- 23 donors
- 5 GEO-defined clinical states
- 16 computational clusters
- 16/16 clusters biologically annotated
- 13 unique biological annotation categories
- 2,000 CD8 HVGs
- 50-dimensional PCA
- CD8-specific UMAP
- Marker discovery completed
- Lineage evidence evaluated
- Donor representation evaluated
- Final Seurat validity check passed
- Cell-ID integrity check passed
- Final adjudicated RDS frozen

**Final object:**
```
results/rds_objects/phase4_final_cd8_analysis_adjudicated.rds
```

**Status:** ✅ **COMPLETE — FINAL OBJECT FROZEN**

---

**Phase 4 README updated:** September 2026  
**Analysis scope:** CD8-labelled compartment characterization  
**Inference unit:** Donor for downstream clinical-state inference  
**Biological interpretation:** Exploratory and observational
