# Phase 6 — Donor-Aware Differential Abundance and Expression Analysis

## Overview

Phase 6 performs **statistical inference on myeloid population abundance and transcriptional state** across the five clinically defined HBV groups.

The key analytical principle is:

> **Donors are the biological replicates. Cells are units for state discovery.**

This phase asks:

* Do myeloid population abundances differ across clinical groups when properly accounting for donor variation?
* Do transcriptional programs differ across clinical groups?
* Which genes show statistically significant differential expression across clinical groups, controlling for donor-level heterogeneity?
* Are observed differences robust to donor structure, or are they artifacts of one or two outlier individuals?

Phase 6 is the **first phase** in which clinical-state inference is performed.

All results explicitly account for donor identity as a random effect or blocking factor.

---

## Input Data

Phase 6 uses the final Phase 5 annotated myeloid Seurat object:

```text
results/rds_objects/phase5_final_myeloid_analysis_adjudicated.rds
```

The object contains:

* **2,489 myeloid cells**
* **24,452 genes**
* **23 donors (biological replicates)**
* **5 clinical states: AC, AR, IA, IT, NL**
* **5 computational clusters with biological annotations**
* **Processed log-CP10K expression**

---

## Critical Issue: Expression Representation

### The Problem

The Phase 5 input data contain **processed log-CP10K expression** rather than raw UMI counts.

This creates a fundamental challenge for standard count-based differential expression methods:

* Standard DE tools (edgeR, DESeq2) assume raw counts as input
* Reversing log-transformation introduces artifacts and information loss
* Reconstructing pseudo-counts from log-CP10K is not mathematically defensible
* Standard library-size normalization is not applicable to pre-normalized data

### Why This Matters

Differential expression analysis typically follows this pipeline:

```text
raw counts
    ↓
[library-size normalization]
    ↓
[log transformation]
    ↓
[differential expression testing]
```

Because the input data already exist at the log-transformed stage, attempting to reverse this creates:

* artificial variance structures
* spurious count distributions
* incorrect dispersion estimation
* invalid p-value calibration

### Solution Approach

Phase 6 uses **log-space linear regression** with proper modeling of the pseudobulk structure:

1. **Aggregate cells to pseudobulk** (per donor × cluster combination)
2. **Normalize the pseudobulk matrix** using appropriate log-scale methods
3. **Model in log-space** using empirical-Bayes linear regression
4. **Apply mixed-model framework** to account for donor as a random effect

This avoids reconstructing counts while preserving statistical validity.

---

## Analytical Framework

### Unit of Analysis: Pseudobulk Aggregation

Instead of testing 2,489 individual cells, Phase 6 aggregates to:

```text
One pseudobulk sample per donor per cluster
```

This produces:

| Cluster | Donors represented | Pseudobulk samples |
| ------: | -----------------: | -----------------: |
|       0 |                 16 |                 16 |
|       1 |                 18 |                 18 |
|       2 |                 19 |                 19 |
|       3 |                  8 |                  8 |
|       4 |                  1 |                  1 |

Each pseudobulk sample represents:

```text
summed expression across all cells from [Donor X] in [Cluster Y]
```

This gives each pseudobulk sample a biologically coherent origin — it is a true replicate from one individual.

### Normalization Strategy

**Pseudobulk construction:**

1. For each cluster, sum the log-CP10K expression values across all cells from each donor
2. This yields a donor × gene matrix where each entry is the sum of log-transformed values
3. (Note: summing log values is not standard; the alternative is to work in the linear space by exponentiating first, summing, then re-log-transforming — this is more principled but computationally intensive and requires stored raw values)

**Recommended approach for your data:**

Because you have access to the individual cell-level log-CP10K values and donor labels, prefer:

1. Exponentiate the cell-level log-CP10K values to recover per-cell linear-space approximations
2. Aggregate these linear-space values (sum) to pseudobulk per donor × cluster
3. Re-apply log transformation to the pseudobulk aggregates
4. Normalize using quantile normalization or VST on the donor × gene pseudobulk matrix

This ensures that pseudobulk aggregates represent true summed expression and are appropriately normalized.

Alternatively, if you retain the original per-sample library-size information from the GEO matrices:

1. Aggregate cell-level log-CP10K sums per donor × cluster
2. Rescale by the per-sample library-size factor
3. Apply standard count-based normalization (TMM, edgeR, or DESeq2) to the rescaled pseudobulk matrix
4. Proceed with log-space testing

---

### Statistical Model

Phase 6 uses a **mixed-effects linear model** on log-scale pseudobulk data:

```text
log(Pseudobulk expression) ~ Clinical State + (1 | Donor) + Batch/Technical Variables (if applicable)
```

Breaking this down:

| Component | Interpretation |
| --- | --- |
| **log(Pseudobulk expression)** | Outcome: log-scaled aggregated expression per donor × cluster |
| **Clinical State** | Fixed effect: AC, AR, IA, IT, NL (the biological comparison) |
| **(1 \| Donor)** | Random intercept: allows each donor to have its own baseline expression level |
| **Batch/Technical** | Optional: if samples were processed in different batches, include batch as a fixed effect or random effect depending on design |

This model structure ensures:

* **Clinical-state differences are detected after accounting for baseline donor-to-donor variation**
* **Donor-specific outliers (like P190719's strong inflammatory signal in cluster 0) do not dominate the test statistic**
* **Statistical power is appropriate for N=16–19 pseudobulk samples per cluster**

---

## Cluster Eligibility

Not all Phase 5 clusters are equally suitable for Phase 6 analysis.

| Cluster | Sample size (donors) | Donor concentration | Eligibility | Recommendation |
| ------: | -------------------: | -------------------: | --- | --- |
|       0 |                   16 |              94.2% | Conditional | **Include with donor blocking; flag results as P190719-influenced** |
|       1 |                   18 |              26.2% | Primary | **Include as primary analysis** |
|       2 |                   19 |              30.7% | Primary | **Include as primary analysis** |
|       3 |                    8 |              91.5% | Weak | **Include only as sensitivity; use caution** |
|       4 |                    1 |             100% | Exclude | **Exclude from clinical-state inference** |

### Cluster 0 Caveat

Cluster 0 contains 1,086 cells from donor P190719 (IA state) and only 67 cells from other donors.

When donor P190719 is used as a random intercept in the model, its elevated baseline myeloid expression is absorbed into the random effect.

**However:**

* Any IA-state effect in cluster 0 may be confounded with P190719-specific biology
* Sensitivity analyses should exclude P190719 to test whether observed effects persist
* Results should be reported with the caveat that cluster 0 is largely donor-driven

### Cluster 3 Caveat

Cluster 3 has weak sample size (N=8 pseudobulk samples) and high donor concentration (91.5% from one donor).

Including it in Phase 6 may be uninformative.

Options:

* **Exclude from primary analysis** and treat as exploratory only
* **Include as a sensitivity** with explicit sample-size caveat in methods/results
* **Combine with cluster 1 or 2** if biological interpretation supports aggregation

### Cluster 4 Exclusion

Cluster 4 is platelet-associated and not myeloid-interpretable.

**Exclude from Phase 6.**

---

## Workflow Overview

Phase 6 proceeds through the following steps:

### Step 1: Pseudobulk Aggregation

For each myeloid cluster separately:

1. Subset the Phase 5 Seurat object to that cluster
2. Extract the expression matrix (cells × genes)
3. Create a donor-level index
4. Sum expression values (or aggregate through exponentiation/re-transformation) by donor
5. Output: One donor × gene matrix per cluster

Result: A normalized pseudobulk expression matrix per cluster where:

```text
rows = donors (N = 8 to 19 depending on cluster)
columns = genes (N = ~24,000)
values = log-scale normalized pseudobulk expression
```

---

### Step 2: Metadata Construction

Create a per-pseudobulk-sample metadata table containing:

```text
Pseudobulk_ID (unique identifier)
Donor (D528848, D529074, ..., P191217)
Clinical_State (AC, AR, IA, IT, NL)
Cluster (0, 1, 2, 3, or 4)
Sample_N_Cells (how many cells went into this pseudobulk aggregate)
```

This metadata will be used to fit the statistical model.

---

### Step 3: Statistical Framework Setup

Use a **mixed-effects linear regression framework** that can:

1. Model on log-scale data (not count data)
2. Fit random intercepts per donor
3. Test fixed effects (clinical state)
4. Provide empirical-Bayes shrinkage for small-sample inference (recommended)

**Recommended tools:**

* **limma + dream()** (R): Empirical-Bayes linear regression with mixed-effects extensions
* **lme4 + lmerTest** (R): Classical mixed-effects regression with Satterthwaite or Kenward-Roger degrees-of-freedom approximation
* **statsmodels** (Python): Mixed-effects regression via REML

The choice between them depends on whether you want:

* **limma/dream**: Faster, borrows strength across genes, designed for genomics
* **lme4/lmerTest**: More conservative, classical statistical inference, better for small sample sizes

For your scenario (N=16–19 pseudobulk samples per cluster, ~24,000 genes), **limma with dream()** is ideal because it applies empirical-Bayes shrinkage across the large gene space.

---

### Step 4: Model Fitting (Per Cluster)

For each eligible cluster:

1. Create the design matrix:

   ```text
   Formula: ~ Clinical_State + (1 | Donor)
   ```

2. Fit the mixed model to the pseudobulk matrix:

   ```text
   For each gene:
       log(Pseudobulk[donor, gene]) ~ Clinical_State_Fixed + Donor_RandomIntercept
   ```

3. Extract:

   * Gene-level effect sizes (log2-fold-change) for each clinical-state contrast
   * Gene-level p-values and adjusted p-values (FDR)
   * Model diagnostics (residuals, variance components)

---

### Step 5: Hypothesis Testing

Define biological contrasts of interest, such as:

```text
IA vs NL
AC vs NL
IT vs NL
AR vs NL
IA vs AR
(and any others of biological interest)
```

For each contrast and each cluster, the model produces:

| Gene | log2FC | p-value | FDR-adjusted p-value | Mean expression (IA) | Mean expression (comparison group) |
| --- | ---: | ---: | ---: | ---: | ---: |
| S100A8 | 0.85 | 0.0023 | 0.031 | 2.34 | 1.49 |
| ... | ... | ... | ... | ... | ... |

---

### Step 6: Validation and Sensitivity Analysis

#### Assumption Checking

1. **Linearity:** Examine residual-versus-fitted plots for each contrast
2. **Normality:** Q-Q plots of residuals per gene
3. **Homogeneity of variance:** Levene's test or visual inspection of variance across groups
4. **No hidden batch effects:** Check whether pseudobulk samples cluster by processing batch (if known)

#### Sensitivity Analyses

1. **Exclude donor P190719:** Re-fit models for clusters 0, 3 to test whether effects persist when the most extreme donor is removed

2. **Exclude cluster 0:** Re-fit to test whether analyses in clusters 1 & 2 are robust (they should not change)

3. **Apply stricter significance thresholds:** Use FDR < 0.01 instead of FDR < 0.05 to identify only the most robust signals

4. **Effect-size filtering:** Report only genes with |log2FC| > 1.0 alongside p-values to avoid overinterpreting tiny effects

---

### Step 7: Interpretation and Visualization

#### Differential Expression Tables

Produce final DE result tables per cluster per contrast, including:

```text
Gene
log2FC
p-value
FDR-adjusted p-value
Mean pseudobulk expression in group A
Mean pseudobulk expression in group B
Interpretation (significant yes/no at FDR < 0.05)
```

#### Volcano Plots

For each cluster and major contrast (e.g., IA vs NL):

```text
x-axis: log2 fold-change
y-axis: -log10(FDR-adjusted p-value)
points: one per gene
highlighted: genes with |log2FC| > 0.5 and FDR < 0.05
```

#### Heatmaps

For each cluster:

```text
Rows: top 30–50 DE genes (sorted by p-value or log2FC)
Columns: pseudobulk samples (donors), colored by clinical state
Values: normalized log-scale expression
```

This shows which DE genes are truly cluster-coherent and which are noisy.

#### Donor-level Visualization

Create plots showing:

```text
For each significant gene in cluster 1:
    y-axis: pseudobulk expression (log scale)
    x-axis: clinical state
    points: individual donors
    colored: by donor
```

This directly visualizes whether observed group differences are driven by one or many donors.

---

## What Phase 6 Addresses

Within the constraint of the processed log-CP10K data and the donor structure:

✅ **Differential abundance:** Do clusters 0–3 differ in cell count across clinical groups?
✅ **Differential expression:** Which genes are up/down-regulated in specific clusters across clinical groups?
✅ **Program-level differences:** Do transcriptional programs (inflammatory, antigen-presentation, etc.) differ across groups?
✅ **Donor robustness:** Are observed differences driven by one outlier donor, or seen across independent individuals?
✅ **Effect sizes:** What is the magnitude of expression change, not just statistical significance?

---

## What Phase 6 Does NOT Address

❌ **Causality:** Differential expression does not establish that expression changes are causal to disease state
❌ **Mechanism:** Expression differences do not reveal mechanism without additional functional data
❌ **Function:** Transcriptional state does not prove functional phenotype
❌ **Cell-cell interaction:** This phase does not model interaction between cell types; that is Phase 8
❌ **Pathway activation:** Gene-level DE does not automatically indicate pathway activation; that is Phase 7
❌ **HBV-specific responses:** The analysis cannot distinguish HBV-driven effects from other chronic viral responses without external validation

---

## Important Caveats

### 1. Sample Size

Even with 23 donors total, cluster-specific sample sizes are modest:

* Clusters 1 & 2: N = 18–19 ✅
* Cluster 0: N = 16 (⚠️ P190719-dominated)
* Cluster 3: N = 8 (❌ underpowered)

At these sample sizes:

* Large, consistent effects are reliably detected
* Small effects (|log2FC| < 0.3) may be missed
* Single-donor artifacts can appear significant

**Mitigation:** Report effect sizes alongside p-values and prioritize biological validation.

### 2. Donor Concentration in Clusters 0 & 3

The Phase 5 analysis revealed that:

* Cluster 0 is 94% donor P190719 (IA state)
* Cluster 3 is 92% donor P190719 (IA state)

While the random-intercept model absorbs baseline donor differences, P190719's dominance means:

* Any strong IA-enriched signal in these clusters may reflect P190719 biology
* Statistical significance does not prove the signal is IA-state-specific

**Mitigation:** Always present results with and without P190719 as a sensitivity check.

### 3. Processed Expression Data

Standard count-based DE tests are not applicable to processed log-CP10K data.

The workflows recommended in Phase 6 use log-space linear regression, which is:

* Statistically valid
* Appropriate for normalized data
* Well-established for microarray-type analysis

However, it does **not** recover information lost during the initial log transformation and normalization.

**Implication:** p-values and fold-changes reflect signal in the pseudobulk log-scale space, not in raw counts.

---

## Statistical Rigor Checklist

Before reporting Phase 6 results, verify:

- [ ] Pseudobulk aggregation preserves cell counts and donor identities correctly
- [ ] Normalization method is appropriate for the log-scale pseudobulk matrix
- [ ] Model formula explicitly includes (1 | Donor) or equivalent random effect
- [ ] Model assumptions (linearity, normality, homogeneity) are checked visually and statistically
- [ ] Multiple testing correction is applied (FDR or Benjamini-Hochberg)
- [ ] Effect sizes are reported alongside p-values
- [ ] Sensitivity analyses (remove P190719, increase FDR threshold) are performed and reported
- [ ] Cluster-specific caveats (sample size, donor concentration) are documented
- [ ] Results are not over-interpreted beyond what the data support

---

## Outputs

### Tables

Per-cluster differential expression results:

```text
results/tables/phase6_[CLUSTER]_DE_results.csv
```

Containing columns:

```text
Gene
log2FC
p_value
FDR_adjusted_p_value
Mean_PseudoBulk_Expression_GroupA
Mean_PseudoBulk_Expression_GroupB
Significant_FDR_0.05
```

### Figures

Per-cluster, per-contrast visualizations:

```text
results/figures/phase6_[CLUSTER]_[CONTRAST]_volcano.png
results/figures/phase6_[CLUSTER]_top_DE_genes_heatmap.png
results/figures/phase6_[CLUSTER]_[CONTRAST]_donor_expression_plot.png
```

### Model Diagnostics

Per-cluster:

```text
results/figures/phase6_[CLUSTER]_residual_diagnostics.png
results/tables/phase6_[CLUSTER]_model_fit_summary.csv
```

---

## Relationship to Subsequent Phases

Phase 6 identifies genes and programs that differ across clinical groups **at the level of transcriptional state and abundance.**

### Phase 7 — Pathway Enrichment

Phase 7 takes the Phase 6 DE gene lists and:

* Maps genes to biological pathways (KEGG, Reactome, GO)
* Identifies enriched pathways per cluster per contrast
* Contextualizes individual gene changes within broader biological processes

### Phase 8 — Cell-Cell Communication

Phase 8 uses all phases (including Phase 6 DE results) to:

* Predict myeloid-lymphoid interactions via ligand-receptor pairs
* Model communication networks across clusters
* Suggest functional consequences of cluster-state changes

### Phase 9 — Integrated Interpretation

Phase 9 synthesizes results from all phases to:

* Construct a coherent narrative of intrahepatic myeloid biology across HBV states
* Distinguish robust, multi-donor signals from outlier-driven artifacts
* Propose testable hypotheses for experimental validation

---

## Interpretation Philosophy

Phase 6 uses statistics to estimate, not dictate, biological reality.

A gene with FDR < 0.05 and |log2FC| > 0.5 is **likely to be true in this cohort** — but not necessarily in other cohorts or in functional assays.

A gene with FDR > 0.05 is **not detected in this cohort** — but may still be biologically important if the effect is small or donor variance is high.

Outlier donors (like P190719) are not artifacts to be hidden; they are part of the biological reality of human disease.

Therefore, Phase 6 results are presented with explicit acknowledgment of:

* Which donors drive observed effects
* Which effects disappear when outliers are removed
* Which signals are robust across independent individuals
* Which findings are cluster-specific versus shared

---

## Final Status

Phase 6 is the **first statistical-inference phase** in the project.

It addresses **whether observed transcriptional differences between clinical groups are statistically supported when donor structure is properly modeled.**

The results feed into Phase 7 (pathway), Phase 8 (interaction), and Phase 9 (integration).

---

**Phase 6 — In development**

**Expected outputs:** Differential expression tables and visualizations per cluster per clinical-state contrast, with explicit donor-level caveats and sensitivity analyses.
