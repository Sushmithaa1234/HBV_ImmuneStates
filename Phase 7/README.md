# Phase 7 — CD8 T-cell Differential-State Interpretation and Pathway Enrichment

## Overview

Phase 7 investigates transcriptional differences and associated biological pathways within intrahepatic CD8 T-cell populations across five clinical states of chronic hepatitis B infection:

- **NL** — Normal liver
- **IA** — Immune active
- **IT** — Immune tolerant
- **AR** — Acute resolved
- **AC** — Asymptomatic carrier

The objective is to identify biological programmes associated with differences in CD8 T-cell transcriptional states across clinical groups and to interpret these findings alongside the population representation results from Phase 4.

This phase builds upon the donor-aware differential expression analysis performed in Phase 6.

The primary questions are:

1. Which genes exhibit differential expression between clinical states within each eligible CD8 T-cell population?
2. Which biological pathways and transcriptional programmes are associated with these differences?
3. How do these programmes vary across the five clinical states?
4. How can these findings contribute to hypotheses about intrahepatic CD8 T-cell biology during chronic HBV infection?

---

## 1. Analytical scope

### Dataset

The analysis uses the liver-only subset of GEO dataset **GSE182159**.

| Characteristic                                               | Description                           |
| ------------------------------------------------------------ | ------------------------------------- |
| Individuals                                                  | 23                                    |
| Clinical states                                              | 5                                     |
| CD8 T-cell populations identified in Phase 4                 | 16                                    |
| Populations eligible for Phase 6 differential-state analysis | 8                                     |
| Expression data                                              | Processed log-CP10K expression values |

### CD8 T-cell populations

Eight populations met the donor-representation criteria established for differential-state analysis.

| Cluster | Annotation                                             | Annotation confidence |
| ------- | ------------------------------------------------------ | --------------------- |
| 0       | GZMK-associated memory-like CD8 T-cell                 | Moderate              |
| 1       | Memory-associated conventional T-cell                  | Moderate              |
| 2       | Activated/dysfunction-associated CD8 T-cell            | High                  |
| 3       | Immediate-early-response CD8-associated population     | High                  |
| 5       | Activation/tissue-associated CD8-associated population | Moderate              |
| 6       | NK-like cytotoxic nonconventional T-cell population    | High                  |
| 8       | NK-like/nonconventional cytotoxic T-cell population    | High                  |
| 11      | Gamma-delta-like T-cell population                     | High                  |

**Annotation note:** Cluster 3 represents a transcriptional programme rather than a conventional biological cell class. Its annotation should therefore be interpreted at the programme level.

The remaining eight Phase 4 clusters were not included in Phase 6 differential-state inference because they did not satisfy the established donor-representation eligibility criteria.

---

## 2. Workflow

Phase 7 consists of two principal analytical scripts.

### Script 01 — Differential expression results audit

**Script:** `phase_07_01_de_results_audit.R`

This script audits the differential expression results generated in Phase 6.

The audit includes:

- Structural integrity of differential expression result tables.
- Verification of population and contrast coverage.
- Duplicate-key detection.
- Missing statistical values and invalid adjusted p-values.
- Reconciliation of reported result summaries.
- Donor representation across populations and clinical states.
- Examination of differential expression effect directions.
- A targeted donor-level diagnostic for Cluster 2 in the NL vs IA comparison.

#### Key audit findings

- 783,280 differential expression result rows were audited.
- Eight CD8 populations and ten pairwise clinical-state contrasts were represented.
- All 80 population × contrast combinations were accounted for.
- No duplicate result keys, missing statistical values or invalid adjusted p-values were identified.
- All eight differential-state models completed without errors.
- 37,076 results had adjusted p-values below 0.05.
- 19,857 results had adjusted p-values below 0.01.

The targeted Cluster 2 NL vs IA diagnostic examined the donor-level directions of the ten leading genes. All ten showed directionally concordant donor mean differences relative to the modelled effect direction. Seven of ten showed complete group separation.

These findings supported the reliability of the leading Cluster 2 signal in this comparison. They do not establish that every significant gene across the full analysis is equally robust.

**Audit outcome: PASS.**

### Script 02 — Pathway enrichment analysis

**Script:** `phase_07_02_pathway_enrichment.R`

This script performs pathway enrichment analysis using the differential expression results from Phase 6.

Two complementary approaches were used:

**A. Gene Set Enrichment Analysis (GSEA)**

GSEA evaluates whether genes belonging to predefined biological pathways are concentrated towards either end of a ranked differential expression gene list.

This allows the analysis to identify coordinated transcriptional programmes, including programmes that may involve genes with individually modest differential expression signals.

The analysis includes:

- Hallmark gene sets.
- Reactome pathways.
- Gene Ontology Biological Process (GO BP).

**B. Over-Representation Analysis (ORA)**

ORA evaluates whether genes meeting the differential expression significance criteria are over-represented in predefined biological pathways.

ORA provides complementary pathway-level evidence to support interpretation of the GSEA findings.

All 80 population × contrast combinations were processed for both GSEA and ORA.

#### Interpretation of enrichment direction

For GSEA:

- Positive NES indicates enrichment towards the first clinical state in the contrast.
- Negative NES indicates enrichment towards the second clinical state.

For ORA:

- Positive logFC indicates higher expression towards the first clinical state.
- Negative logFC indicates higher expression towards the second clinical state.

#### Multiple-testing correction

GSEA adjusted p-values were calculated within each population × contrast × library.

ORA adjusted p-values were calculated within each population × contrast × direction × library.

Adjusted p-values should therefore be interpreted within their respective correction scopes, rather than as globally adjusted significance across the entire Phase 7 analysis.

**Enrichment analysis outcome: COMPLETED.**

---

## 3. Quality control and analytical limitations

The enrichment analysis completed successfully, with the following computational observations documented.

### GSEA calculation warnings

A total of 2,194 result rows had missing NES, p-value, adjusted p-value and log2err values because of unbalanced positive and negative ranked statistics.

These rows could not be interpreted as enrichment results and were excluded from biological interpretation.

An additional 206 rows had an available NES but missing log2err values. All had adjusted p-values of at least 0.6743 and none were significant.

Ties in the preranked statistics were approximately 0.01%, representing a negligible proportion of the ranked genes.

These observations did not require rerunning the enrichment analysis.

### Statistical limitations

- AC and AR each contain only three donors.
- Several population-state combinations had only two qualifying donors in the differential-state analysis.
- Donors, rather than individual cells, were treated as the biological replicates in the primary differential-state framework.
- Some cluster-state combinations have limited donor representation in the abundance analysis.
- Descriptive donor-level population proportions should not be interpreted as statistically established differential abundance.

### Biological interpretation limitations

- Clinical states are treated as categorical groups, not as a chronological disease-progression sequence.
- Pathway enrichment indicates coordinated expression patterns, not direct measurement of pathway activity or cellular function.
- GSEA and ORA results are complementary but are not independent biological discoveries.
- Overlapping pathways across Hallmark, Reactome and GO BP should not be counted as separate independent mechanisms.
- Pathway enrichment does not establish causality, HBV specificity, protective or pathogenic function, exhaustion, tissue residency or antiviral activity.
- Moderate-confidence annotations for Clusters 0, 1 and 5 must be considered when interpreting their biological identities.
- The results are observational and hypothesis-generating.

---

## 4. Major biological findings

### 4.1 Recurring transcriptional themes

The significant GSEA results identified several recurring biological themes across the eight CD8 populations and ten clinical-state contrasts.

**Mitochondrial metabolism and oxidative phosphorylation**

Oxidative phosphorylation, respiratory electron transport, aerobic respiration and ATP-related pathways were among the most recurrent enriched programmes.

These findings indicate differences in mitochondrial energy-metabolism-associated gene expression across multiple CD8 populations and clinical comparisons.

**Inflammatory signalling**

TNF/NF-κB-associated pathways appeared in multiple populations and comparisons, indicating differences in inflammatory signalling-related transcription.

**Transcriptional and cellular regulation**

MYC-associated transcription, mTORC1, hypoxia and apoptosis-related programmes were identified in selected populations and comparisons.

**RNA processing**

RNA splicing, mRNA processing and ribonucleoprotein biogenesis were prominent among pathways enriched towards AC, particularly within Cluster 0 in comparisons involving AR and IA.

**Antigen processing and presentation**

Proteasome, antigen-presentation and protein-processing pathways were identified in Cluster 2 in the NL vs IA comparison, enriched towards IA.

These findings suggest that differences between clinical states involve multiple biological dimensions rather than a single shared transcriptional programme.

### 4.2 Interpretation across clinical states

The following summary integrates the descriptive population representation results from Phase 4 with the transcriptional programmes identified in Phase 7.

**Normal liver (NL)**

- Cluster 2 showed a particularly high median representation (31.5%).
- Cluster 3 exhibited oxidative phosphorylation, inflammatory, MYC-associated, translation and hypoxia-related programmes enriched towards NL in NL vs AR.
- In NL vs IA, Cluster 2 showed antigen-processing, proteasome and protein-processing programmes enriched towards IA.

These findings indicate a distinct CD8 population composition in NL, alongside population-specific transcriptional differences.

**Immune active (IA)**

- Cluster 5 exhibited oxidative phosphorylation, aerobic respiration and ATP-related programmes enriched towards IA in IA vs AC.
- MYC-associated transcription was enriched towards IA in Cluster 5 in IA vs AR.
- Cluster 3 showed TNF/NF-κB-associated programmes enriched towards IA in IA vs AR.
- Cluster 2 showed antigen-processing and protein-processing programmes enriched towards IA in NL vs IA.

These findings suggest metabolic, inflammatory and antigen-processing-related transcriptional differences across multiple CD8 populations in IA.

**Immune tolerant (IT)**

- Cluster 0 had the highest median representation in IT (44.9%).
- Cluster 3 exhibited TNF/NF-κB, mTORC1, hypoxia, MYC and apoptosis-related programmes enriched towards IT in IT vs AR.
- Cluster 6 showed oxidative phosphorylation-related programmes enriched towards IT in IT vs IA.

These findings indicate that IT is associated with a distinct CD8 population composition and measurable transcriptional differences involving cellular regulation, inflammatory signalling and metabolism.

**Acute resolved (AR)**

- Cluster 5 had its highest median representation in AR (13.0%).
- Cluster 5 exhibited mitochondrial respiration and ATP-synthesis-related programmes enriched towards AR in AR vs AC.
- Cluster 6 showed oxidative phosphorylation-related programmes enriched towards AR in IA vs AR.

These findings suggest differences in CD8 population representation and mitochondrial metabolism-associated transcription in AR.

**Asymptomatic carrier (AC)**

- Cluster 0 had a median representation of 25.9%.
- RNA splicing, mRNA processing and ribonucleoprotein biogenesis were prominent among pathways enriched towards AC in Cluster 0, particularly in comparisons involving AR and IA.

These findings suggest RNA-processing-related transcriptional differences in the GZMK-associated memory-like CD8 population in comparisons involving AC.

**Important:** These state-level summaries describe findings from specific pairwise comparisons. They do not imply that every pathway listed is universally characteristic of a clinical state.

---

## 5. Overall biological interpretation

Phase 7 demonstrates that intrahepatic CD8 T-cell populations exhibit heterogeneous transcriptional programmes across the five clinical states of chronic hepatitis B infection.

The findings suggest that differences between clinical states involve both:

1. **Population representation:** different CD8 populations show distinct descriptive representation patterns across clinical groups.
2. **Differential-state transcription:** individual CD8 populations exhibit differences in mitochondrial metabolism, inflammatory signalling, transcriptional regulation, RNA processing and antigen-processing-related programmes.

The combination of these findings supports the hypothesis that the intrahepatic CD8 compartment is differently organised across HBV clinical states, both in population composition and in transcriptional characteristics.

These observations provide biological hypotheses for further investigation. They do not independently establish the mechanisms responsible for viral persistence, immune-mediated liver injury, viral clearance or clinical disease progression.

---

### Next phase

**Phase 8 — Predicted cell-cell communication analysis**

Investigate predicted cell-cell communication patterns within the intrahepatic immune environment across the five clinical states of chronic hepatitis B infection.

The next phase will extend the analysis beyond transcriptional programmes within individual populations to examine potential signalling interactions between immune populations.

Phase 9 will integrate the findings from the preceding phases into a coherent biological interpretation of the intrahepatic immune environment during chronic HBV infection.

