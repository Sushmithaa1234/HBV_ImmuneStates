# ============================================================================
# HBV Immune States scRNA-seq Analysis
# Script 02: QC Characterization, Conservative Filtering & Diagnostics
#
# PURPOSE
#   1. Characterize cell-level QC metrics recoverable from processed
#      log-counts-per-10,000 (log-CP10K) expression data
#   2. Characterize QC distributions globally, by sample, donor and state
#   3. Apply a conservative minimum detected-gene filter
#   4. Perform sample-aware robust MAD diagnostics
#   5. Audit cell retention by GSM, donor and clinical state
#   6. Preserve QC flags and diagnostic metrics for downstream auditing
#   7. Explicitly document limitations caused by unavailable raw UMI counts
#
# DATA LIMITATION
#   The GEO supplementary matrices contain processed log-counts-per-10,000
#   expression values. Raw UMI count matrices are not available.
#
#   Therefore this script DOES NOT calculate:
#     - raw UMI counts per cell
#     - conventional mitochondrial UMI percentage
#     - UMI-based RNA complexity
#     - count-based doublet predictions
#
#   Instead, it calculates metrics directly recoverable from the processed
#   expression representation.
#
# FILTERING PHILOSOPHY
#   - Conservative
#   - Diagnostic-first
#   - No automatic mitochondrial filtering
#   - No automatic high-complexity filtering
#   - No donor/sample exclusion based solely on cell number
#
# Starting cell-level filter:
#   Genes_Detected >= 200
#
# IMPORTANT
#   This filter is a conservative starting criterion, not a claim that every
#   cell below 200 detected genes is biologically meaningless.
#
# ============================================================================


# ============================================================================
# 0. SETUP
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(here)
  library(patchwork)
})

setwd(here())

cat("\n")
cat("============================================================\n")
cat("PHASE 2: SCRIPT 02 — QC & CONSERVATIVE FILTERING\n")
cat("============================================================\n\n")


# ============================================================================
# 1. LOAD MERGED PROCESSED EXPRESSION OBJECT
# ============================================================================

cat("=== STEP 1: LOADING MERGED OBJECT ===\n\n")

seurat_merged <- readRDS(
  "results/rds_objects/seurat_merged_processed_expression.rds"
)

cat(
  "Cells:",
  ncol(seurat_merged),
  "\n"
)

cat(
  "Genes:",
  nrow(seurat_merged),
  "\n"
)

cat(
  "Samples:",
  dplyr::n_distinct(seurat_merged$GSM),
  "\n"
)

cat(
  "Donors:",
  dplyr::n_distinct(seurat_merged$Donor),
  "\n"
)

cat(
  "Clinical states:",
  dplyr::n_distinct(seurat_merged$Phase),
  "\n\n"
)

# ============================================================================
# 2. VERIFY DATA REPRESENTATION
# ============================================================================

cat("=== STEP 2: VERIFYING DATA REPRESENTATION ===\n\n")

DefaultAssay(seurat_merged) <- "RNA"

rna_layers <- Layers(
  seurat_merged[["RNA"]]
)

cat(
  "RNA layers:\n"
)

print(
  rna_layers
)

# --------------------------------------------------------------------------
# Identify processed expression layers
# --------------------------------------------------------------------------

data_layers <- grep(
  "^data(\\.|$)",
  rna_layers,
  value = TRUE
)

counts_layers <- grep(
  "^counts(\\.|$)",
  rna_layers,
  value = TRUE
)

# --------------------------------------------------------------------------
# Validate processed data layers
# --------------------------------------------------------------------------

if (length(data_layers) == 0) {
  
  stop(
    "ERROR: No RNA data layers were detected. ",
    "Expected one or more processed 'data' layers ",
    "(e.g. data.GSM5519469, data.GSM5519471, etc.)."
  )
  
}

cat(
  "\n✓ Detected ",
  length(data_layers),
  " processed RNA data layer(s).\n",
  sep = ""
)

print(
  data_layers
)

# --------------------------------------------------------------------------
# Confirm that raw counts are absent
# --------------------------------------------------------------------------

if (length(counts_layers) > 0) {
  
  stop(
    "ERROR: Raw counts layer(s) were detected:\n",
    paste(counts_layers, collapse = ", "),
    "\nThis script is designed for the processed log-CP10K dataset ",
    "with no raw UMI counts."
  )
  
}

cat(
  "\n✓ No raw UMI counts layers are present.\n"
)

cat(
  "✓ RNA expression is represented by processed data layers.\n"
)

cat(
  "✓ QC will therefore use metrics recoverable from the supplied ",
  "log-CP10K expression values.\n\n"
)

# ============================================================================
# 3. VERIFY REQUIRED METADATA
# ============================================================================

cat("=== STEP 3: VERIFYING METADATA ===\n\n")

required_metadata <- c(
  "GSM",
  "Donor",
  "Phase"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(seurat_merged@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "ERROR: Required metadata columns missing: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
  
}

cat(
  "✓ Required metadata present:",
  paste(
    required_metadata,
    collapse = ", "
  ),
  "\n\n"
)

# ============================================================================
# 4. EXTRACT PROCESSED EXPRESSION AND CALCULATE QC METRICS
#
# IMPORTANT:
#   The Seurat v5 RNA assay contains 23 sample-specific processed data layers:
#   data.GSM5519469, data.GSM5519471, etc.
#
#   Raw UMI counts are unavailable.
#   Therefore QC metrics are calculated separately within each processed
#   log-CP10K layer and then combined at the cell-metadata level.
# ============================================================================

cat("=== STEP 4: CALCULATING PROCESSED-DATA QC METRICS ===\n\n")


# --------------------------------------------------------------------------
# 4A. Identify mitochondrial and ribosomal genes
# --------------------------------------------------------------------------

gene_names <- rownames(
  seurat_merged[["RNA"]]
)

mito_genes <- grep(
  "^MT-",
  gene_names,
  value = TRUE
)

ribosomal_genes <- grep(
  "^RP[SL]",
  gene_names,
  value = TRUE
)

if (length(mito_genes) == 0) {
  
  stop(
    "ERROR: No mitochondrial genes matching '^MT-' were found."
  )
}

cat(
  "Total genes in RNA assay: ",
  length(gene_names),
  "\n",
  sep = ""
)

cat(
  "Mitochondrial genes detected: ",
  length(mito_genes),
  "\n",
  sep = ""
)

cat(
  "Ribosomal genes detected: ",
  length(ribosomal_genes),
  "\n\n",
  sep = ""
)


# --------------------------------------------------------------------------
# 4B. Calculate QC metrics separately for each processed data layer
# --------------------------------------------------------------------------

qc_metrics_list <- lapply(
  
  data_layers,
  
  function(layer_name) {
    
    cat(
      "Processing: ",
      layer_name,
      "\n",
      sep = ""
    )
    
    
    # ----------------------------------------------------------------------
    # Extract this sample's processed expression layer
    # ----------------------------------------------------------------------
    
    expression_matrix <- LayerData(
      seurat_merged[["RNA"]],
      layer = layer_name
    )
    
    
    # ----------------------------------------------------------------------
    # Validate layer
    # ----------------------------------------------------------------------
    
    if (nrow(expression_matrix) != length(gene_names)) {
      
      stop(
        "ERROR: Unexpected gene count in layer ",
        layer_name,
        "."
      )
    }
    
    
    # ----------------------------------------------------------------------
    # Genes detected
    #
    # A gene is considered detected when its supplied processed expression
    # value is > 0.
    #
    # This is NOT a raw UMI-based nFeature_RNA metric.
    # ----------------------------------------------------------------------
    
    genes_detected <- Matrix::colSums(
      expression_matrix > 0
    )
    
    
    # ----------------------------------------------------------------------
    # Mitochondrial genes detected
    # ----------------------------------------------------------------------
    
    mito_present <- intersect(
      mito_genes,
      rownames(expression_matrix)
    )
    
    if (length(mito_present) > 0) {
      
      mito_genes_detected <- Matrix::colSums(
        expression_matrix[
          mito_present,
          ,
          drop = FALSE
        ] > 0
      )
      
    } else {
      
      mito_genes_detected <- rep(
        0,
        ncol(expression_matrix)
      )
    }
    
    
    # ----------------------------------------------------------------------
    # Ribosomal genes detected
    # ----------------------------------------------------------------------
    
    ribo_present <- intersect(
      ribosomal_genes,
      rownames(expression_matrix)
    )
    
    if (length(ribo_present) > 0) {
      
      ribo_genes_detected <- Matrix::colSums(
        expression_matrix[
          ribo_present,
          ,
          drop = FALSE
        ] > 0
      )
      
    } else {
      
      ribo_genes_detected <- rep(
        0,
        ncol(expression_matrix)
      )
    }
    
    
    # ----------------------------------------------------------------------
    # Detection fractions
    #
    # These are fractions of DETECTED GENES.
    # They are NOT fractions of UMIs or transcripts.
    # ----------------------------------------------------------------------
    
    mito_detection_fraction <- (
      mito_genes_detected /
        pmax(
          genes_detected,
          1
        )
    )
    
    ribo_detection_fraction <- (
      ribo_genes_detected /
        pmax(
          genes_detected,
          1
        )
    )
    
    
    # ----------------------------------------------------------------------
    # Total log-CP10K signal
    #
    # Descriptive only.
    # This is NOT library size or UMI depth.
    # ----------------------------------------------------------------------
    
    total_logCP10K <- Matrix::colSums(
      expression_matrix
    )
    
    
    # ----------------------------------------------------------------------
    # Return cell-level metrics
    # ----------------------------------------------------------------------
    
    data.frame(
      
      Cell = colnames(
        expression_matrix
      ),
      
      Genes_Detected = as.numeric(
        genes_detected
      ),
      
      Mito_Genes_Detected = as.numeric(
        mito_genes_detected
      ),
      
      Mito_Detection_Fraction = as.numeric(
        mito_detection_fraction
      ),
      
      Ribo_Genes_Detected = as.numeric(
        ribo_genes_detected
      ),
      
      Ribo_Detection_Fraction = as.numeric(
        ribo_detection_fraction
      ),
      
      Total_LogCP10K = as.numeric(
        total_logCP10K
      ),
      
      stringsAsFactors = FALSE
    )
  }
)


# ============================================================================
# 5. COMBINE AND VALIDATE QC METRICS
# ============================================================================

cat("\n=== STEP 5: COMBINING AND VALIDATING QC METRICS ===\n\n")


qc_metrics <- do.call(
  rbind,
  qc_metrics_list
)

rownames(qc_metrics) <- qc_metrics$Cell

qc_metrics$Cell <- NULL


# --------------------------------------------------------------------------
# Validate total number of cells
# --------------------------------------------------------------------------

if (
  nrow(qc_metrics) != ncol(seurat_merged)
) {
  
  stop(
    "ERROR: Number of QC metric rows does not match number of cells.\n",
    "QC rows: ",
    nrow(qc_metrics),
    "\n",
    "Seurat cells: ",
    ncol(seurat_merged)
  )
}


# --------------------------------------------------------------------------
# Validate cell identities
# --------------------------------------------------------------------------

if (
  !identical(
    sort(rownames(qc_metrics)),
    sort(colnames(seurat_merged))
  )
) {
  
  stop(
    "ERROR: QC metric cell IDs do not match Seurat cell IDs."
  )
}


# --------------------------------------------------------------------------
# Validate missing values
# --------------------------------------------------------------------------

if (
  any(
    !is.finite(
      as.matrix(qc_metrics)
    )
  )
) {
  
  stop(
    "ERROR: Non-finite QC metric values detected."
  )
}


cat(
  "✓ QC metrics calculated for ",
  nrow(qc_metrics),
  " cells.\n",
  sep = ""
)

cat(
  "✓ Processed layers analyzed: ",
  length(data_layers),
  "\n",
  sep = ""
)

cat(
  "✓ Gene features analyzed: ",
  length(gene_names),
  "\n",
  sep = ""
)

cat(
  "✓ Cell identities validated.\n"
)

cat(
  "✓ No raw UMI metrics calculated.\n"
)

cat(
  "✓ No conventional mitochondrial UMI percentage calculated.\n\n"
)


# ============================================================================
# 6. ATTACH QC METRICS TO SEURAT METADATA
# ============================================================================

cat("=== STEP 6: ATTACHING QC METRICS TO SEURAT METADATA ===\n\n")


# --------------------------------------------------------------------------
# Match explicitly by cell ID.
#
# This avoids relying on the ordering of cells across layers.
# --------------------------------------------------------------------------

cell_match <- match(
  colnames(seurat_merged),
  rownames(qc_metrics)
)

if (
  anyNA(cell_match)
) {
  
  stop(
    "ERROR: Some Seurat cells could not be matched to QC metrics."
  )
}


seurat_merged$Genes_Detected <- (
  qc_metrics$Genes_Detected[cell_match]
)

seurat_merged$Mito_Genes_Detected <- (
  qc_metrics$Mito_Genes_Detected[cell_match]
)

seurat_merged$Mito_Detection_Fraction <- (
  qc_metrics$Mito_Detection_Fraction[cell_match]
)

seurat_merged$Ribo_Genes_Detected <- (
  qc_metrics$Ribo_Genes_Detected[cell_match]
)

seurat_merged$Ribo_Detection_Fraction <- (
  qc_metrics$Ribo_Detection_Fraction[cell_match]
)

seurat_merged$Total_LogCP10K <- (
  qc_metrics$Total_LogCP10K[cell_match]
)


cat(
  "✓ QC metrics attached to Seurat metadata.\n\n"
)


# ============================================================================
# 7. GLOBAL QC SUMMARY
# ============================================================================

cat("=== STEP 7: GLOBAL QC SUMMARY ===\n\n")


qc_metric_names <- c(
  "Genes_Detected",
  "Mito_Genes_Detected",
  "Mito_Detection_Fraction",
  "Ribo_Genes_Detected",
  "Ribo_Detection_Fraction",
  "Total_LogCP10K"
)


for (
  metric in qc_metric_names
) {
  
  cat(
    "\n",
    metric,
    ":\n",
    sep = ""
  )
  
  print(
    summary(
      seurat_merged@meta.data[[metric]]
    )
  )
}


cat("\n")


# ============================================================================
# 8. SAVE PRE-QC CELL METADATA
# ============================================================================

cat("=== STEP 8: SAVING PRE-QC CELL METADATA ===\n\n")


pre_qc_metadata <- seurat_merged@meta.data %>%
  
  rownames_to_column(
    "Cell"
  )


write_csv(
  pre_qc_metadata,
  "results/tables/pre_qc_cell_metadata.csv"
)


cat(
  "✓ Pre-QC cell metadata saved.\n\n"
)


# ============================================================================
# 9. SAMPLE-LEVEL QC SUMMARY
# ============================================================================

cat("=== STEP 9: SAMPLE-LEVEL QC SUMMARY ===\n\n")


sample_qc_before <- seurat_merged@meta.data %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Median_Genes_Detected = median(
      Genes_Detected,
      na.rm = TRUE
    ),
    
    Mean_Genes_Detected = mean(
      Genes_Detected,
      na.rm = TRUE
    ),
    
    Median_Mito_Detection_Fraction = median(
      Mito_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Mean_Mito_Detection_Fraction = mean(
      Mito_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Median_Ribo_Detection_Fraction = median(
      Ribo_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Mean_Ribo_Detection_Fraction = mean(
      Ribo_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Median_Total_LogCP10K = median(
      Total_LogCP10K,
      na.rm = TRUE
    ),
    
    Mean_Total_LogCP10K = mean(
      Total_LogCP10K,
      na.rm = TRUE
    ),
    
    Pct_Below_200_Genes = mean(
      Genes_Detected < 200,
      na.rm = TRUE
    ) * 100,
    
    Pct_Above_5000_Genes = mean(
      Genes_Detected > 5000,
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Phase,
    GSM
  )


write_csv(
  sample_qc_before,
  "results/tables/sample_qc_before_filtering.csv"
)


print(
  sample_qc_before,
  n = Inf
)


cat(
  "\n✓ Sample-level QC summary saved.\n\n"
)


# ============================================================================
# 10. SAMPLE-SPECIFIC MAD DIAGNOSTICS
#
# IMPORTANT:
#   MAD is diagnostic only.
#   It does NOT automatically remove cells.
# ============================================================================

cat("=== STEP 10: SAMPLE-SPECIFIC MAD DIAGNOSTICS ===\n\n")


calculate_mad_summary <- function(
    data,
    value_column,
    metric_name
) {
  
  data %>%
    
    group_by(
      GSM,
      Donor,
      Phase
    ) %>%
    
    summarise(
      
      Metric = metric_name,
      
      N = sum(
        is.finite(
          .data[[value_column]]
        )
      ),
      
      Median = median(
        .data[[value_column]],
        na.rm = TRUE
      ),
      
      MAD = mad(
        .data[[value_column]],
        na.rm = TRUE
      ),
      
      Lower_MAD_3 = Median - 3 * MAD,
      
      Upper_MAD_3 = Median + 3 * MAD,
      
      .groups = "drop"
    )
}


mad_genes <- calculate_mad_summary(
  seurat_merged@meta.data,
  "Genes_Detected",
  "Genes_Detected"
)


mad_mito <- calculate_mad_summary(
  seurat_merged@meta.data,
  "Mito_Detection_Fraction",
  "Mito_Detection_Fraction"
)


mad_total_expression <- calculate_mad_summary(
  seurat_merged@meta.data,
  "Total_LogCP10K",
  "Total_LogCP10K"
)


mad_diagnostics <- bind_rows(
  mad_genes,
  mad_mito,
  mad_total_expression
)


write_csv(
  mad_diagnostics,
  "results/tables/sample_specific_MAD_diagnostics.csv"
)


cat(
  "✓ Sample-specific MAD diagnostics calculated.\n"
)

cat(
  "✓ MAD diagnostics saved.\n"
)

cat(
  "IMPORTANT: MAD statistics are diagnostic only.\n\n"
)


# ============================================================================
# 11. DEFINE CONSERVATIVE QC FILTER
# ============================================================================

cat("=== STEP 11: DEFINING CONSERVATIVE QC FILTER ===\n\n")


min_genes_detected <- 200L


cat(
  "Minimum detected genes:",
  min_genes_detected,
  "\n"
)

cat(
  "No mitochondrial hard cutoff is applied.\n"
)

cat(
  "No upper detected-gene cutoff is applied.\n"
)

cat(
  "No donor/sample is excluded based solely on cell number.\n\n"
)


# ============================================================================
# 12. QC FLAGS
# ============================================================================

cat("=== STEP 12: CALCULATING QC FLAGS ===\n\n")


seurat_merged$QC_LowGenes <- (
  seurat_merged$Genes_Detected <
    min_genes_detected
)


seurat_merged$QC_Pass <- (
  !seurat_merged$QC_LowGenes
)


cells_before <- ncol(
  seurat_merged
)


cells_low_genes <- sum(
  seurat_merged$QC_LowGenes,
  na.rm = TRUE
)


cells_passing_qc <- sum(
  seurat_merged$QC_Pass,
  na.rm = TRUE
)


cat(
  "Cells before QC:",
  cells_before,
  "\n"
)

cat(
  "Cells below 200 detected genes:",
  cells_low_genes,
  "\n"
)

cat(
  "Cells passing conservative QC:",
  cells_passing_qc,
  "\n\n"
)


# ============================================================================
# 13. QC FILTERING DECISION SUMMARY
# ============================================================================

qc_decision_summary <- tibble(
  
  Criterion = c(
    "Below 200 detected genes",
    "Passing conservative QC"
  ),
  
  Cells = c(
    cells_low_genes,
    cells_passing_qc
  ),
  
  Percent_of_input = round(
    
    100 *
      c(
        cells_low_genes,
        cells_passing_qc
      ) /
      cells_before,
    
    3
  )
)


write_csv(
  qc_decision_summary,
  "results/tables/qc_filtering_decision_summary.csv"
)


print(
  qc_decision_summary
)


cat("\n")


# ============================================================================
# 14. PRE-QC FIGURES
# ============================================================================

cat("=== STEP 14: CREATING PRE-QC FIGURES ===\n\n")


# --------------------------------------------------------------------------
# 14A. Genes detected by clinical state
# --------------------------------------------------------------------------

p_genes_phase <- ggplot(
  
  seurat_merged@meta.data,
  
  aes(
    x = Phase,
    y = Genes_Detected
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  geom_hline(
    yintercept = min_genes_detected,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Genes detected per cell by clinical state",
    x = "Clinical state",
    y = "Detected genes"
  )


ggsave(
  "results/figures/01_qc_genes_detected_by_phase.png",
  p_genes_phase,
  width = 9,
  height = 6,
  dpi = 300
)


# --------------------------------------------------------------------------
# 14B. Genes detected by sample
# --------------------------------------------------------------------------

p_genes_sample <- ggplot(
  
  seurat_merged@meta.data,
  
  aes(
    x = reorder(
      GSM,
      Genes_Detected,
      FUN = median
    ),
    y = Genes_Detected
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  geom_hline(
    yintercept = min_genes_detected,
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Genes detected per cell by sample",
    x = "GSM",
    y = "Detected genes"
  )


ggsave(
  "results/figures/02_qc_genes_detected_by_sample.png",
  p_genes_sample,
  width = 10,
  height = 12,
  dpi = 300
)


# --------------------------------------------------------------------------
# 14C. Mitochondrial detection fraction by clinical state
# --------------------------------------------------------------------------

p_mito_phase <- ggplot(
  
  seurat_merged@meta.data,
  
  aes(
    x = Phase,
    y = Mito_Detection_Fraction
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Mitochondrial gene detection fraction",
    x = "Clinical state",
    y = "MT genes / detected genes"
  )


ggsave(
  "results/figures/03_qc_mito_detection_by_phase.png",
  p_mito_phase,
  width = 9,
  height = 6,
  dpi = 300
)


# --------------------------------------------------------------------------
# 14D. Total log-CP10K signal by clinical state
# --------------------------------------------------------------------------

p_log_signal <- ggplot(
  
  seurat_merged@meta.data,
  
  aes(
    x = Phase,
    y = Total_LogCP10K
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Total log-CP10K signal by clinical state",
    x = "Clinical state",
    y = "Sum of log-CP10K values"
  )


ggsave(
  "results/figures/04_qc_total_logCP10K_by_phase.png",
  p_log_signal,
  width = 9,
  height = 6,
  dpi = 300
)


# --------------------------------------------------------------------------
# 14E. Genes detected vs total log-expression
# --------------------------------------------------------------------------

p_genes_signal <- ggplot(
  
  seurat_merged@meta.data,
  
  aes(
    x = Total_LogCP10K,
    y = Genes_Detected
  )
) +
  
  geom_point(
    alpha = 0.2,
    size = 0.3
  ) +
  
  geom_hline(
    yintercept = min_genes_detected,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Detected genes vs total log-CP10K signal",
    x = "Total log-CP10K signal",
    y = "Detected genes"
  )


ggsave(
  "results/figures/05_qc_genes_vs_log_signal.png",
  p_genes_signal,
  width = 8,
  height = 6,
  dpi = 300
)


cat(
  "✓ Pre-QC figures created.\n\n"
)


# ============================================================================
# 15. APPLY CONSERVATIVE FILTER
# ============================================================================

cat("=== STEP 15: APPLYING CONSERVATIVE QC FILTER ===\n\n")


seurat_filtered <- subset(
  seurat_merged,
  subset = QC_Pass
)


cells_after_qc <- ncol(
  seurat_filtered
)


cells_lost_qc <- (
  cells_before -
    cells_after_qc
)


pct_retained_qc <- (
  100 *
    cells_after_qc /
    cells_before
)


cat(
  "Cells before QC:",
  cells_before,
  "\n"
)

cat(
  "Cells after QC:",
  cells_after_qc,
  "\n"
)

cat(
  "Cells removed:",
  cells_lost_qc,
  "\n"
)

cat(
  "Percent retained:",
  round(
    pct_retained_qc,
    2
  ),
  "%\n\n"
)


# ============================================================================
# 16. RETENTION BY CLINICAL STATE
# ============================================================================

cat("=== STEP 16: RETENTION BY CLINICAL STATE ===\n\n")


retained_cells <- colnames(
  seurat_filtered
)


retention_metadata <- seurat_merged@meta.data %>%
  
  rownames_to_column(
    "Cell"
  ) %>%
  
  mutate(
    Retained = Cell %in% retained_cells
  )


phase_loss <- retention_metadata %>%
  
  group_by(
    Phase
  ) %>%
  
  summarise(
    
    Before = n(),
    
    After = sum(
      Retained
    ),
    
    Lost = Before - After,
    
    Pct_Retained = round(
      100 * After / Before,
      2
    ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Phase
  )


print(
  phase_loss
)


write_csv(
  phase_loss,
  "results/tables/phase_cell_retention.csv"
)


# ============================================================================
# 17. RETENTION BY SAMPLE / DONOR
# ============================================================================

cat("\n=== STEP 17: RETENTION BY SAMPLE / DONOR ===\n\n")


sample_loss <- retention_metadata %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Before = n(),
    
    After = sum(
      Retained
    ),
    
    Lost = Before - After,
    
    Pct_Retained = round(
      100 * After / Before,
      2
    ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Pct_Retained
  )


print(
  sample_loss,
  n = Inf
)


write_csv(
  sample_loss,
  "results/tables/sample_cell_retention.csv"
)


# ============================================================================
# 18. SMALL-SAMPLE AUDIT
# ============================================================================

cat("\n=== STEP 18: SMALL-SAMPLE AUDIT ===\n\n")


small_sample_audit <- sample_loss %>%
  
  mutate(
    Small_Sample_Flag = Before < 100
  ) %>%
  
  filter(
    Small_Sample_Flag
  )


if (
  nrow(small_sample_audit) == 0
) {
  
  cat(
    "No samples contain fewer than 100 cells before QC.\n"
  )
  
} else {
  
  cat(
    "Samples containing fewer than 100 cells before QC:\n\n"
  )
  
  print(
    small_sample_audit
  )
}


write_csv(
  small_sample_audit,
  "results/tables/small_sample_audit.csv"
)


# ============================================================================
# 19. POST-QC SUMMARY
# ============================================================================

cat("\n=== STEP 19: POST-QC SUMMARY ===\n\n")


post_qc_summary <- seurat_filtered@meta.data %>%
  
  summarise(
    
    Cells = n(),
    
    Median_Genes_Detected = median(
      Genes_Detected,
      na.rm = TRUE
    ),
    
    Median_Mito_Detection_Fraction = median(
      Mito_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Median_Ribo_Detection_Fraction = median(
      Ribo_Detection_Fraction,
      na.rm = TRUE
    ),
    
    Median_Total_LogCP10K = median(
      Total_LogCP10K,
      na.rm = TRUE
    )
  )


print(
  post_qc_summary
)


# ============================================================================
# 20. POST-QC FIGURE
# ============================================================================

cat("\n=== STEP 20: CREATING POST-QC FIGURE ===\n\n")


p_post_qc <- ggplot(
  
  seurat_filtered@meta.data,
  
  aes(
    x = Total_LogCP10K,
    y = Genes_Detected,
    color = Mito_Detection_Fraction
  )
) +
  
  geom_point(
    size = 0.4,
    alpha = 0.4
  ) +
  
  facet_wrap(
    ~ Phase
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Post-QC processed-expression metrics",
    x = "Total log-CP10K signal",
    y = "Detected genes",
    color = "MT detection fraction"
  )


ggsave(
  "results/figures/06_qc_post_filter.png",
  p_post_qc,
  width = 14,
  height = 8,
  dpi = 300
)


cat(
  "✓ Post-QC figure created.\n\n"
)


# ============================================================================
# 21. DOUBLETS / HIGH-COMPLEXITY CELLS
#
# Raw UMI counts are unavailable.
#
# Therefore no automatic count-based scDblFinder filtering is performed
# in this script.
#
# High-complexity cells are flagged for downstream atlas-level inspection
# rather than automatically removed.
# ============================================================================

cat("=== STEP 21: DOUBLET / HIGH-COMPLEXITY HANDLING ===\n\n")


seurat_filtered$Potential_HighComplexity <- (
  seurat_filtered$Genes_Detected > 5000
)


high_complexity_summary <- seurat_filtered@meta.data %>%
  
  summarise(
    
    Cells = n(),
    
    High_Complexity_Cells = sum(
      Potential_HighComplexity,
      na.rm = TRUE
    ),
    
    High_Complexity_Percent = round(
      
      100 *
        High_Complexity_Cells /
        Cells,
      
      2
    )
  )


print(
  high_complexity_summary
)


write_csv(
  high_complexity_summary,
  "results/tables/high_complexity_cell_summary.csv"
)


cat(
  "\nNo cells were removed because of high detected-gene counts.\n"
)

cat(
  "Potential high-complexity cells will be inspected during atlas construction.\n"
)

cat(
  "No automatic doublet removal was performed.\n\n"
)


# ============================================================================
# 22. FINAL QC STATUS
# ============================================================================

cat("=== STEP 22: FINAL QC STATUS ===\n\n")


final_qc_status <- tibble(
  
  Metric = c(
    
    "Input cells",
    
    "Input genes",
    
    "Cells below 200 detected genes",
    
    "Post-QC cells",
    
    "Percent retained",
    
    "Samples",
    
    "Donors",
    
    "Clinical states",
    
    "Smallest sample before QC",
    
    "High-complexity cells flagged",
    
    "Mitochondrial UMI percentage",
    
    "Raw UMI counts",
    
    "Automatic doublet removal"
  ),
  
  
  Value = c(
    
    cells_before,
    
    length(gene_names),
    
    cells_low_genes,
    
    cells_after_qc,
    
    round(
      pct_retained_qc,
      2
    ),
    
    dplyr::n_distinct(
      seurat_filtered$GSM
    ),
    
    dplyr::n_distinct(
      seurat_filtered$Donor
    ),
    
    dplyr::n_distinct(
      seurat_filtered$Phase
    ),
    
    min(
      sample_qc_before$Cells
    ),
    
    sum(
      seurat_filtered$Potential_HighComplexity,
      na.rm = TRUE
    ),
    
    "Not available",
    
    "Not available",
    
    "Not performed"
  )
)


print(
  final_qc_status
)


write_csv(
  final_qc_status,
  "results/tables/final_qc_status.csv"
)


# ============================================================================
# 23. SAVE QC-FILTERED OBJECT
# ============================================================================

cat("\n=== STEP 23: SAVING QC-FILTERED OBJECT ===\n\n")


saveRDS(
  seurat_filtered,
  "results/rds_objects/seurat_qc_processed_expression.rds"
)


cat(
  "✓ Saved:\n",
  "results/rds_objects/seurat_qc_processed_expression.rds\n\n",
  sep = ""
)


# ============================================================================
# 24. SAVE QC SUMMARY
# ============================================================================

qc_summary <- tibble(
  
  metric = c(
    
    "Cells before QC",
    
    "Genes in input assay",
    
    "Cells below 200 detected genes",
    
    "Cells after conservative QC",
    
    "Percent retained",
    
    "Minimum detected genes",
    
    "Samples",
    
    "Donors",
    
    "Clinical states",
    
    "Raw UMI counts available",
    
    "Conventional mitochondrial percentage available",
    
    "Automatic doublet removal performed"
  ),
  
  
  value = c(
    
    cells_before,
    
    length(gene_names),
    
    cells_low_genes,
    
    cells_after_qc,
    
    round(
      pct_retained_qc,
      2
    ),
    
    min_genes_detected,
    
    dplyr::n_distinct(
      seurat_filtered$GSM
    ),
    
    dplyr::n_distinct(
      seurat_filtered$Donor
    ),
    
    dplyr::n_distinct(
      seurat_filtered$Phase
    ),
    
    "FALSE",
    
    "FALSE",
    
    "FALSE"
  )
)


write_csv(
  qc_summary,
  "results/tables/qc_summary.csv"
)


# ============================================================================
# 25. FINAL GATE
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PHASE 2 — SCRIPT 02 COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Input genes:",
  length(gene_names),
  "\n"
)

cat(
  "Input cells:",
  cells_before,
  "\n"
)

cat(
  "Cells below 200 detected genes:",
  cells_low_genes,
  "\n"
)

cat(
  "Final QC cells:",
  cells_after_qc,
  "\n"
)

cat(
  "Retention:",
  round(
    pct_retained_qc,
    2
  ),
  "%\n\n"
)


cat(
  "Clinical-state retention:\n"
)


print(
  phase_loss
)


cat("\n")


cat(
  "Doublet filtering: NOT performed\n"
)

cat(
  "High-complexity cells: FLAGGED, NOT REMOVED\n"
)

cat(
  "Raw UMI counts: unavailable\n"
)

cat(
  "Conventional mitochondrial percentage: unavailable\n\n"
)


cat(
  "QC figures: results/figures/\n"
)

cat(
  "QC tables: results/tables/\n"
)

cat(
  "Final QC object:\n",
  "results/rds_objects/seurat_qc_processed_expression.rds\n\n",
  sep = ""
)


cat(
  "NEXT GATE:\n"
)

cat(
  "Proceed to Phase 3 atlas construction using the QC-filtered\n",
  "processed-expression object.\n\n"
)


cat("============================================================\n")
cat("END OF SCRIPT 02\n")
cat("============================================================\n")

