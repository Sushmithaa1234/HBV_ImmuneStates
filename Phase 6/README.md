# Phase 6 — Donor-Aware Differential-State Analysis

## Overview

Phase 6 performs donor-aware differential-state analysis of intrahepatic CD8 T-cell populations across five clinically defined HBV-related groups using publicly available single-cell RNA-sequencing data from GEO accession **GSE182159**.

The objective is to identify genes exhibiting differences in modeled expression between clinical states within individual CD8 T-cell populations and establish a statistical foundation for subsequent pathway enrichment and biological interpretation.

Unlike conventional cell-level differential expression analysis, this phase treats the **donor as the biological unit of inference**, reducing the risk of pseudoreplication associated with treating individual cells as independent biological replicates.

### Clinical groups

| Abbreviation | Clinical state |
|---|---|
| NL | Normal Liver |
| IT | Immune Tolerant |
| IA | Immune Active |
| AC | Asymptomatic Carrier |
| AR | Acute Resolved |

These groups are treated as distinct clinical categories, not as a chronological disease progression.

## 1. Analysis scope

- **Dataset:** GSE182159
- **Tissue:** Liver only
- **Biological replicates:** 23 donors
- **Inferential compartment:** CD8 T-cell populations
- **Primary populations:** 8 computationally identified CD8 clusters
- **Clinical groups:** 5
- **Pairwise contrasts:** 10
- **Total population–contrast combinations:** 80

Myeloid populations were characterized descriptively in Phase 5 but excluded from downstream inferential analysis because donor/state representation was insufficient for robust five-state comparisons.

## 2. Input data

The analysis uses the following previously generated objects:

| Input | Description |
|---|---|
| `phase3_final_full_dataset.rds` | Full-dataset cell atlas and metadata |
| `phase4_final_cd8_analysis_adjudicated.rds` | Adjudicated CD8 subclustering object |

The expression input consists of processed **log-CP10K values**.

These are non-count expression values. They are not reconstructed into raw counts, reverse-transformed, or treated as conventional count data.

## 3. CD8 population selection

Primary populations were selected using predefined donor-representation criteria.

A donor–population combination was retained when it contained at least 10 cells.

A CD8 population was considered eligible for primary inference when at least two retained donors were represented in every clinical state.

Eight populations met these criteria:

- CD8_Cluster_0
- CD8_Cluster_1
- CD8_Cluster_2
- CD8_Cluster_3
- CD8_Cluster_5
- CD8_Cluster_6
- CD8_Cluster_8
- CD8_Cluster_11

The remaining clusters were excluded from primary five-state inference because of insufficient donor/state representation.

The two-donor minimum is an eligibility threshold, not a guarantee of adequate statistical power. Donor representation and concentration are retained for interpretation.

## 4. Donor-level aggregation

Processed expression values were aggregated using:

`dreamlet::aggregateNonCountSignal()`

Aggregation was performed at the:

**Donor × CD8 population level**

Parameters:

| Parameter | Value |
|---|---:|
| `min.cells` | 10 |
| `min.signal` | 0.01 |
| `min.samples` | 4 |
| `min.prop` | 0.4 |

The resulting donor-level expression summaries were used for downstream differential-state modelling. Precision weights were retained through the dreamlet workflow.

No pseudo-count reconstruction was performed.

## 5. Statistical modelling

Differential-state analysis was performed using **dreamlet** with a cell-means model:

`~ 0 + Clinical_State`

The model estimates expression means for each clinical state and evaluates ten planned pairwise contrasts.

### Planned contrasts

| Reference group | Comparison group |
|---|---|
| NL | IT |
| NL | IA |
| NL | AR |
| NL | AC |
| IT | IA |
| IT | AR |
| IT | AC |
| IA | AR |
| IA | AC |
| AR | AC |

Each contrast was evaluated across all eight primary CD8 populations.

### Effect interpretation

The reported `logFC` represents the limma/dreamlet coefficient on the modeled expression scale.

It should not automatically be interpreted as a conventional log2 fold change.

For a contrast such as `NL_vs_IA`:

- Positive logFC: higher modeled expression in NL.
- Negative logFC: higher modeled expression in IA.

The sign interpretation follows the order of the specified contrast.

### Statistical significance

Gene-level results include:

- `logFC`: estimated expression difference.
- `AveExpr`: average modeled expression.
- `t`: moderated t-statistic.
- `P.Value`: unadjusted p-value.
- `adj.P.Val`: multiple-testing-adjusted p-value.
- `B`: log-odds statistic.
- `z.std`: standardized statistic.

Genes with `adj.P.Val < 0.05` are classified as statistically significant for descriptive summaries.

The adjusted p-values are extracted from the complete dreamlet result for each contrast. The intended correction scope is across the eight primary populations within each contrast. This correction scope should be verified against the installed dreamlet/limma implementation and documented alongside the final results. The values should not be assumed to represent correction across all ten contrasts.

## 6. Outputs

### RDS objects

| File | Contents |
|---|---|
| `phase6_cd8_donor_population_noncount_aggregated.rds` | Donor–population aggregated expression |
| `phase6_cd8_dreamlet_contrasts.rds` | Fitted dreamlet result |
| `phase6_cd8_final_differential_state_analysis.rds` | Consolidated Phase 6 analysis object |

### Tables

| File | Description |
|---|---|
| `phase6_cd8_cell_metadata.csv` | CD8 cell-level analysis metadata |
| `phase6_cd8_donor_population_counts.csv` | Donor–population cell counts |
| `phase6_cd8_state_coverage.csv` | Donor coverage across clinical states |
| `phase6_cd8_population_eligibility.csv` | Primary, secondary and excluded population adjudication |
| `phase6_cd8_sce_donor_population_counts.csv` | SCE-level donor–population counts |
| `phase6_cd8_sce_threshold_audit.csv` | Aggregation threshold audit |
| `phase6_cd8_dreamlet_model_details.csv` | Model-fitting diagnostics |
| `phase6_cd8_all_gene_level_results.csv` | Complete gene-level results across all contrasts and populations |
| `phase6_cd8_significant_gene_summary.csv` | Significant-gene counts and effect summaries |

### Figures

| File | Description |
|---|---|
| `phase6_cd8_donor_state_coverage.png` | Donor representation across populations and clinical states |
| `phase6_cd8_fdr05_gene_counts.png` | Number of FDR-significant gene results across comparisons |

## 7. Quality control

The script includes checks for:

- Expected clinical-state representation.
- One clinical state per donor.
- Unique cell identifiers.
- Expected donor and population counts.
- Valid processed expression values.
- Successful construction of the SingleCellExperiment object.
- Expected primary assay representation.
- Successful dreamlet model fitting.
- Finite gene-level statistics.
- Valid adjusted p-values.
- Expected result-row counts.

The model-fitting diagnostics are retained for subsequent review.

## 8. Interpretation and limitations

Phase 6 identifies **differential-state expression**, not differential population abundance.

The results describe differences in modeled expression within CD8 populations across clinical groups. They do not establish:

- Causal effects of HBV infection.
- HBV antigen specificity.
- Functional exhaustion or activation.
- Protective or pathogenic cellular roles.
- Changes in CD8 population abundance.
- Direct clinical utility.

Cluster annotations remain provisional transcriptional classifications rather than definitive functional identities.

The analysis is restricted to liver-derived samples and the donor representation available in the public dataset. The number of donors, particularly in individual clinical groups, constrains statistical power and generalizability.

**Phase status:** Differential-state modelling completed successfully across eight primary CD8 populations and ten pairwise clinical-state contrasts. Downstream pathway enrichment and biological interpretation remain to be completed.
