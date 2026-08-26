# ============================================================================
# HBV Immune States scRNA-seq Analysis
# Script 02: QC Filtering, Outlier Assessment & Doublet Detection
#
# Purpose:
#   1. Calculate cell-level QC metrics
#   2. Characterize QC distributions globally and by sample/state
#   3. Apply pre-specified starting QC thresholds
#   4. Perform sample-aware QC diagnostics using robust MAD statistics
#   5. Assess cell retention by GSM, donor and clinical state
#   6. Evaluate potential doublets using scDblFinder when available
#   7. Preserve QC/doublet flags for downstream auditing
#   8. Save the QC-filtered Seurat object
#
# IMPORTANT:
#   - Fixed thresholds are the primary starting filter.
#   - MAD-based statistics are diagnostic and are NOT used as an
#     automatic second filtering criterion.
#   - Doublet predictions are recorded separately from QC filtering.
#   - Biological populations must not be removed solely because they
#     have unusual RNA complexity.
#
# Starting QC thresholds:
#   nFeature_RNA: 200–6000
#   percent.mt: <= 20%
#
# Cohort:
#   23 donors / samples
#   106,592 cells before QC
#
# Clinical states:
#   NL = 6 donors
#   IT = 6 donors
#   IA = 5 donors
#   AR = 3 donors
#   AC = 3 donors
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
cat("PHASE 2: SCRIPT 02 — QC FILTERING & QUALITY CONTROL\n")
cat("============================================================\n\n")


# ============================================================================
# 1. LOAD MERGED RAW OBJECT
# ============================================================================

cat("=== STEP 1: LOADING MERGED OBJECT ===\n\n")

seurat_merged <- readRDS(
  "results/rds_objects/seurat_merged_raw.rds"
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
# 2. CREATE OUTPUT DIRECTORIES
# ============================================================================

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/rds_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================================
# 3. VERIFY REQUIRED METADATA
# ============================================================================

cat("=== STEP 2: VERIFYING METADATA ===\n\n")

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
# 4. CALCULATE CELL-LEVEL QC METRICS
# ============================================================================

cat("=== STEP 3: CALCULATING QC METRICS ===\n\n")


# --------------------------------------------------------------------------
# 4A. Mitochondrial genes
# --------------------------------------------------------------------------

mito_genes <- grep(
  "^MT-",
  rownames(seurat_merged),
  value = TRUE
)

if (length(mito_genes) == 0) {
  
  stop(
    "ERROR: No mitochondrial genes matching '^MT-' were found.\n",
    "Check the gene naming convention before proceeding."
  )
  
}

cat(
  "Mitochondrial genes detected:",
  length(mito_genes),
  "\n"
)

seurat_merged[["percent.mt"]] <- PercentageFeatureSet(
  seurat_merged,
  features = mito_genes
)


# --------------------------------------------------------------------------
# 4B. Ribosomal genes
# --------------------------------------------------------------------------

ribosomal_genes <- grep(
  "^RP[SL]",
  rownames(seurat_merged),
  value = TRUE
)

if (length(ribosomal_genes) > 0) {
  
  seurat_merged[["percent.rb"]] <- PercentageFeatureSet(
    seurat_merged,
    features = ribosomal_genes
  )
  
  cat(
    "Ribosomal genes detected:",
    length(ribosomal_genes),
    "\n"
  )
  
} else {
  
  cat(
    "No ribosomal genes detected using '^RP[SL]'.\n"
  )
  
}


# --------------------------------------------------------------------------
# 4C. Complexity metric
# --------------------------------------------------------------------------

seurat_merged$RNA_complexity <- (
  seurat_merged$nFeature_RNA /
    seurat_merged$nCount_RNA
)

cat(
  "✓ Mitochondrial percentage calculated\n"
)

cat(
  "✓ Ribosomal percentage calculated\n"
)

cat(
  "✓ RNA complexity calculated\n\n"
)


# ============================================================================
# 5. GLOBAL QC SUMMARY
# ============================================================================

cat("=== STEP 4: GLOBAL QC SUMMARY ===\n\n")

qc_metrics <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "RNA_complexity"
)

if ("percent.rb" %in% colnames(seurat_merged@meta.data)) {
  
  qc_metrics <- c(
    qc_metrics,
    "percent.rb"
  )
  
}

for (metric in qc_metrics) {
  
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
# 6. SAVE PRE-QC CELL METADATA
# ============================================================================

pre_qc_metadata <- seurat_merged@meta.data %>%
  rownames_to_column("Cell")

write_csv(
  pre_qc_metadata,
  "results/tables/pre_qc_cell_metadata.csv"
)

cat(
  "✓ Pre-QC metadata saved\n\n"
)


# ============================================================================
# 7. SAMPLE-LEVEL QC SUMMARY
# ============================================================================

cat("=== STEP 5: SAMPLE-LEVEL QC SUMMARY ===\n\n")

sample_qc_before <- seurat_merged@meta.data %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Median_Genes = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Mean_Genes = mean(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_UMIs = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Mean_UMIs = mean(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_MT = median(
      percent.mt,
      na.rm = TRUE
    ),
    
    Mean_MT = mean(
      percent.mt,
      na.rm = TRUE
    ),
    
    Median_RNA_Complexity = median(
      RNA_complexity,
      na.rm = TRUE
    ),
    
    Pct_MT_Above_20 = mean(
      percent.mt > 20,
      na.rm = TRUE
    ) * 100,
    
    Pct_Below_200_Genes = mean(
      nFeature_RNA < 200,
      na.rm = TRUE
    ) * 100,
    
    Pct_Above_6000_Genes = mean(
      nFeature_RNA > 6000,
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

cat("\n✓ Sample-level QC summary saved\n\n")


# ============================================================================
# 8. ROBUST SAMPLE-SPECIFIC QC DIAGNOSTICS
#
#    MAD is used here ONLY to identify unusual distributions.
#    It is NOT used as an automatic filtering rule.
# ============================================================================

cat("=== STEP 6: SAMPLE-SPECIFIC MAD DIAGNOSTICS ===\n\n")

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
        !is.na(.data[[value_column]])
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


mad_nfeature <- calculate_mad_summary(
  seurat_merged@meta.data,
  "nFeature_RNA",
  "nFeature_RNA"
)

mad_ncount <- calculate_mad_summary(
  seurat_merged@meta.data,
  "nCount_RNA",
  "nCount_RNA"
)

mad_mt <- calculate_mad_summary(
  seurat_merged@meta.data,
  "percent.mt",
  "percent.mt"
)

mad_diagnostics <- bind_rows(
  mad_nfeature,
  mad_ncount,
  mad_mt
)

write_csv(
  mad_diagnostics,
  "results/tables/sample_specific_MAD_diagnostics.csv"
)

cat(
  "✓ Sample-specific MAD diagnostics calculated\n"
)

cat(
  "✓ MAD diagnostics saved\n"
)

cat(
  "IMPORTANT: MAD thresholds are diagnostic only and are not\n",
  "being used as automatic cell-removal criteria.\n\n"
)


# ============================================================================
# 9. DEFINE PRE-SPECIFIED QC THRESHOLDS
# ============================================================================

cat("=== STEP 7: STARTING QC THRESHOLDS ===\n\n")

min_features <- 200L
max_features <- 6000L
max_mt_pct <- 20

cat(
  "nFeature_RNA:",
  min_features,
  "–",
  max_features,
  "\n"
)

cat(
  "percent.mt <= ",
  max_mt_pct,
  "%\n\n",
  sep = ""
)


# ============================================================================
# 10. CALCULATE INDIVIDUAL QC FLAGS
# ============================================================================

cat("=== STEP 8: CALCULATING QC FLAGS ===\n\n")

seurat_merged$QC_LowGenes <- (
  seurat_merged$nFeature_RNA <
    min_features
)

seurat_merged$QC_HighGenes <- (
  seurat_merged$nFeature_RNA >
    max_features
)

seurat_merged$QC_HighMT <- (
  seurat_merged$percent.mt >
    max_mt_pct
)

seurat_merged$QC_Pass <- (
  !seurat_merged$QC_LowGenes &
    !seurat_merged$QC_HighGenes &
    !seurat_merged$QC_HighMT
)


# Count each failure category

cells_before <- ncol(seurat_merged)

cells_low_genes <- sum(
  seurat_merged$QC_LowGenes,
  na.rm = TRUE
)

cells_high_genes <- sum(
  seurat_merged$QC_HighGenes,
  na.rm = TRUE
)

cells_high_mt <- sum(
  seurat_merged$QC_HighMT,
  na.rm = TRUE
)

cells_passing_qc <- sum(
  seurat_merged$QC_Pass,
  na.rm = TRUE
)

cat(
  "Below minimum genes:",
  cells_low_genes,
  "\n"
)

cat(
  "Above maximum genes:",
  cells_high_genes,
  "\n"
)

cat(
  "Above mitochondrial threshold:",
  cells_high_mt,
  "\n"
)

cat(
  "Passing all fixed QC criteria:",
  cells_passing_qc,
  "\n\n"
)


# ============================================================================
# 11. QC FAILURE OVERLAP
# ============================================================================

cat("=== STEP 9: QC FAILURE OVERLAP ===\n\n")

qc_failure_table <- table(
  LowGenes = seurat_merged$QC_LowGenes,
  HighGenes = seurat_merged$QC_HighGenes,
  HighMT = seurat_merged$QC_HighMT
)

print(
  qc_failure_table
)

cat("\n")


# ============================================================================
# 12. QC DECISION SUMMARY
# ============================================================================

qc_decision_summary <- tibble(
  
  Criterion = c(
    "Below minimum genes",
    "Above maximum genes",
    "Above mitochondrial threshold",
    "Passing all fixed QC criteria"
  ),
  
  Cells = c(
    cells_low_genes,
    cells_high_genes,
    cells_high_mt,
    cells_passing_qc
  ),
  
  Percent_of_input = round(
    100 *
      c(
        cells_low_genes,
        cells_high_genes,
        cells_high_mt,
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
# 13. FIGURE 1 — QC VIOLINS BY CLINICAL STATE
# ============================================================================

cat("=== STEP 10: PRE-QC FIGURES ===\n\n")

p_phase <- VlnPlot(
  seurat_merged,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "Phase",
  ncol = 3,
  pt.size = 0
)

ggsave(
  "results/figures/01_qc_violin_by_phase.png",
  p_phase,
  width = 15,
  height = 5,
  dpi = 300
)


# ============================================================================
# 14. FIGURE 2 — QC VIOLINS BY SAMPLE
# ============================================================================

p_sample <- VlnPlot(
  seurat_merged,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "GSM",
  ncol = 3,
  pt.size = 0
) &
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    )
  )

ggsave(
  "results/figures/02_qc_violin_by_sample.png",
  p_sample,
  width = 20,
  height = 10,
  dpi = 300
)


# ============================================================================
# 15. FIGURE 3 — GENE DETECTION VS UMI COUNT
# ============================================================================

p_features <- FeatureScatter(
  seurat_merged,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  
  geom_hline(
    yintercept = min_features,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = max_features,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Gene detection vs UMI counts"
  )

ggsave(
  "results/figures/03_qc_scatter_features.png",
  p_features,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================================
# 16. FIGURE 4 — UMI VS MITOCHONDRIAL PERCENTAGE
# ============================================================================

p_mt <- FeatureScatter(
  seurat_merged,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
) +
  
  geom_hline(
    yintercept = max_mt_pct,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "UMI counts vs mitochondrial percentage"
  )

ggsave(
  "results/figures/04_qc_scatter_mt.png",
  p_mt,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================================
# 17. FIGURE 5 — NFEATURE DISTRIBUTION BY PHASE
# ============================================================================

p_nfeature <- seurat_merged@meta.data %>%
  
  ggplot(
    aes(
      x = nFeature_RNA
    )
  ) +
  
  geom_histogram(
    bins = 60
  ) +
  
  geom_vline(
    xintercept = min_features,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = max_features,
    linetype = "dashed"
  ) +
  
  facet_wrap(
    ~ Phase,
    scales = "free_y"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Detected genes per cell by clinical state",
    x = "Genes detected",
    y = "Cells"
  )

ggsave(
  "results/figures/05_qc_hist_nfeature.png",
  p_nfeature,
  width = 12,
  height = 8,
  dpi = 300
)


# ============================================================================
# 18. FIGURE 6 — MITOCHONDRIAL PERCENTAGE BY PHASE
# ============================================================================

p_mt_phase <- seurat_merged@meta.data %>%
  
  ggplot(
    aes(
      x = Phase,
      y = percent.mt
    )
  ) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 0.3,
    alpha = 0.2
  ) +
  
  geom_hline(
    yintercept = max_mt_pct,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Mitochondrial percentage by clinical state",
    y = "Mitochondrial reads (%)"
  )

ggsave(
  "results/figures/06_qc_mt_by_phase.png",
  p_mt_phase,
  width = 8,
  height = 6,
  dpi = 300
)


cat(
  "✓ Pre-QC figures created\n\n"
)


# ============================================================================
# 19. APPLY FIXED QC FILTER
# ============================================================================

cat("=== STEP 11: APPLYING FIXED QC FILTER ===\n\n")

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
# 20. RETENTION BY CLINICAL STATE
# ============================================================================

cat("=== STEP 12: RETENTION BY CLINICAL STATE ===\n\n")

retained_cells <- colnames(
  seurat_filtered
)

retention_metadata <- seurat_merged@meta.data %>%
  
  rownames_to_column("Cell") %>%
  
  mutate(
    Retained = Cell %in% retained_cells
  )

phase_loss <- retention_metadata %>%
  
  group_by(
    Phase
  ) %>%
  
  summarise(
    
    Before = n(),
    
    After = sum(Retained),
    
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
# 21. RETENTION BY SAMPLE / DONOR
# ============================================================================

cat("\n=== STEP 13: RETENTION BY SAMPLE ===\n\n")

sample_loss <- retention_metadata %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Before = n(),
    
    After = sum(Retained),
    
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
# 22. IDENTIFY SAMPLES WITH LOW RETENTION
# ============================================================================

low_retention_samples <- sample_loss %>%
  
  filter(
    Pct_Retained < 50
  )

cat("\nSamples retaining <50%:\n")

if (nrow(low_retention_samples) == 0) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    low_retention_samples
  )
  
}


# ============================================================================
# 23. EXPLICITLY TRACK SMALL SAMPLES
# ============================================================================

cat("\n=== STEP 14: SMALL-SAMPLE AUDIT ===\n\n")

small_sample_audit <- sample_loss %>%
  
  mutate(
    Small_Sample_Flag = Before < 100
  ) %>%
  
  filter(
    Small_Sample_Flag
  )

if (nrow(small_sample_audit) == 0) {
  
  cat(
    "No samples contain <100 cells.\n"
  )
  
} else {
  
  cat(
    "Samples containing <100 cells:\n\n"
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
# 24. POST-QC SUMMARY
# ============================================================================

cat("\n=== STEP 15: POST-QC SUMMARY ===\n\n")

post_qc_summary <- seurat_filtered@meta.data %>%
  
  summarise(
    
    Cells = n(),
    
    Median_Genes = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_UMIs = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_MT = median(
      percent.mt,
      na.rm = TRUE
    ),
    
    Mean_MT = mean(
      percent.mt,
      na.rm = TRUE
    ),
    
    Median_Complexity = median(
      RNA_complexity,
      na.rm = TRUE
    )
  )

print(
  post_qc_summary
)


# ============================================================================
# 25. POST-QC FIGURE
# ============================================================================

p_post_qc <- seurat_filtered@meta.data %>%
  
  ggplot(
    aes(
      x = nCount_RNA,
      y = nFeature_RNA,
      color = percent.mt
    )
  ) +
  
  geom_point(
    size = 0.4,
    alpha = 0.4
  ) +
  
  facet_wrap(
    ~ Phase
  ) +
  
  scale_color_viridis_c() +
  
  theme_minimal() +
  
  labs(
    title = "Post-filter QC",
    x = "UMI counts",
    y = "Genes detected",
    color = "MT %"
  )

ggsave(
  "results/figures/07_qc_post_filter_scatter.png",
  p_post_qc,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================================
# 26. DOUBLEt DETECTION
#
# IMPORTANT:
#   Doublet prediction is deliberately kept separate from the fixed QC
#   filter. Predictions are recorded first; removal is a separate decision.
# ============================================================================
cat("\n=== STEP 16: DOUBLET DETECTION ===\n\n")

scdblfinder_available <- requireNamespace(
  "scDblFinder",
  quietly = TRUE
)

if (scdblfinder_available) {
  
  cat("✓ scDblFinder is installed.\n")
  cat("Preparing temporary object for scDblFinder...\n\n")
  
  # ------------------------------------------------------------
  # Work on a copy
  # ------------------------------------------------------------
  
  seurat_dbl <- seurat_filtered
  
  DefaultAssay(seurat_dbl) <- "RNA"
  
  # ------------------------------------------------------------
  # Seurat v5: join RNA layers for SingleCellExperiment
  # conversion.
  #
  # IMPORTANT:
  # This is done ONLY on the temporary doublet-detection copy.
  # The main QC object remains layer-separated.
  # ------------------------------------------------------------
  
  cat("RNA layers before joining:\n")
  print(
    Layers(
      seurat_dbl[["RNA"]]
    )
  )
  
  cat("\nJoining RNA layers for scDblFinder...\n")
  
  seurat_dbl[["RNA"]] <- JoinLayers(
    seurat_dbl[["RNA"]]
  )
  
  cat("\nRNA layers after joining:\n")
  print(
    Layers(
      seurat_dbl[["RNA"]]
    )
  )
  
  cat("\n✓ RNA layers joined in temporary object.\n\n")
  
  # ------------------------------------------------------------
  # Convert to SingleCellExperiment
  # ------------------------------------------------------------
  
  cat("Converting to SingleCellExperiment...\n")
  
  sce_dbl <- as.SingleCellExperiment(
    seurat_dbl
  )
  
  cat("✓ Conversion complete.\n\n")
  
  # ------------------------------------------------------------
  # Run scDblFinder sample-aware
  # ------------------------------------------------------------
  
  cat("Running scDblFinder...\n")
  cat("Samples:", n_distinct(seurat_dbl$GSM), "\n\n")
  
  set.seed(12345)
  
  sce_dbl <- scDblFinder::scDblFinder(
    sce_dbl,
    samples = "GSM",
    verbose = TRUE
  )
  
  cat("\n✓ scDblFinder complete.\n\n")
  
  # ------------------------------------------------------------
  # Transfer results back to Seurat
  # ------------------------------------------------------------
  
  dbl_metadata <- as.data.frame(
    SummarizedExperiment::colData(sce_dbl)
  )
  
  common_cells <- intersect(
    colnames(seurat_dbl),
    rownames(dbl_metadata)
  )
  
  if (length(common_cells) != ncol(seurat_dbl)) {
    
    stop(
      "ERROR: Cell names did not match after scDblFinder."
    )
    
  }
  
  seurat_dbl$scDblFinder.score <- (
    dbl_metadata[
      colnames(seurat_dbl),
      "scDblFinder.score"
    ]
  )
  
  seurat_dbl$scDblFinder.class <- (
    dbl_metadata[
      colnames(seurat_dbl),
      "scDblFinder.class"
    ]
  )
  
  seurat_dbl$Potential_Doublet <- (
    seurat_dbl$scDblFinder.class == "doublet"
  )
  
  # ------------------------------------------------------------
  # Overall doublet classification
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("scDblFinder CLASSIFICATION\n")
  cat("============================================================\n\n")
  
  print(
    table(
      seurat_dbl$scDblFinder.class,
      useNA = "ifany"
    )
  )
  
  # ------------------------------------------------------------
  # Doublet rate by sample
  # ------------------------------------------------------------
  
  doublet_by_sample <- seurat_dbl@meta.data %>%
    
    rownames_to_column("Cell") %>%
    
    group_by(
      GSM,
      Donor,
      Phase
    ) %>%
    
    summarise(
      Cells = n(),
      
      Doublets = sum(
        scDblFinder.class == "doublet",
        na.rm = TRUE
      ),
      
      Doublet_Rate = round(
        100 * Doublets / Cells,
        2
      ),
      
      Median_Doublet_Score = median(
        scDblFinder.score,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    arrange(
      desc(Doublet_Rate)
    )
  
  cat("\nDoublet rate by sample:\n\n")
  
  print(
    doublet_by_sample,
    n = Inf
  )
  
  write_csv(
    doublet_by_sample,
    "results/tables/doublet_rate_by_sample.csv"
  )
  
  # ------------------------------------------------------------
  # Doublet rate by clinical state
  # ------------------------------------------------------------
  
  doublet_by_phase <- seurat_dbl@meta.data %>%
    
    group_by(
      Phase
    ) %>%
    
    summarise(
      Cells = n(),
      
      Doublets = sum(
        scDblFinder.class == "doublet",
        na.rm = TRUE
      ),
      
      Doublet_Rate = round(
        100 * Doublets / Cells,
        2
      ),
      
      .groups = "drop"
    )
  
  cat("\nDoublet rate by clinical state:\n\n")
  
  print(
    doublet_by_phase
  )
  
  write_csv(
    doublet_by_phase,
    "results/tables/doublet_rate_by_phase.csv"
  )
  
  # ------------------------------------------------------------
  # Score distribution
  # ------------------------------------------------------------
  
  p_doublet_score <- ggplot(
    seurat_dbl@meta.data,
    aes(
      x = scDblFinder.score
    )
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    theme_minimal() +
    
    labs(
      title = "scDblFinder Score Distribution",
      x = "scDblFinder score",
      y = "Cells"
    )
  
  ggsave(
    "results/figures/08_scDblFinder_score_distribution.png",
    p_doublet_score,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # Doublet rate by sample
  # ------------------------------------------------------------
  
  p_doublet_sample <- ggplot(
    doublet_by_sample,
    aes(
      x = reorder(
        GSM,
        Doublet_Rate
      ),
      y = Doublet_Rate
    )
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    theme_minimal() +
    
    labs(
      title = "Predicted Doublet Rate by Sample",
      x = "GSM",
      y = "Predicted doublets (%)"
    )
  
  ggsave(
    "results/figures/09_doublet_rate_by_sample.png",
    p_doublet_sample,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # Replace working object with annotated object
  # ------------------------------------------------------------
  
  seurat_filtered <- seurat_dbl
  
  rm(
    seurat_dbl,
    sce_dbl,
    dbl_metadata
  )
  
  gc()
  
  cat(
    "\n✓ Doublet predictions attached to Seurat metadata.\n"
  )
  
  cat(
    "IMPORTANT: Doublets have NOT been automatically removed.\n"
  )
  
} else {
  
  cat(
    "⚠ scDblFinder is not installed.\n\n"
  )
  
  seurat_filtered$Potential_Doublet <- NA
  
}

# ============================================================================
# 27. DOUBLEt SUMMARY
# ============================================================================

cat("\n=== STEP 17: FINAL DOUBLET SUMMARY ===\n\n")

if (
  "scDblFinder.class" %in%
  colnames(seurat_filtered@meta.data)
) {
  
  doublet_summary <- seurat_filtered@meta.data %>%
    
    summarise(
      
      Cells = n(),
      
      Predicted_Doublets = sum(
        scDblFinder.class == "doublet",
        na.rm = TRUE
      ),
      
      Predicted_Singlets = sum(
        scDblFinder.class == "singlet",
        na.rm = TRUE
      ),
      
      Doublet_Rate = round(
        100 *
          Predicted_Doublets /
          Cells,
        2
      )
    )
  
  print(
    doublet_summary
  )
  
  write_csv(
    doublet_summary,
    "results/tables/doublet_summary.csv"
  )
  
} else {
  
  cat(
    "scDblFinder classification unavailable.\n"
  )
  
}

# ============================================================================
# DOUBLEt VALIDATION CHECK 1
# QC METRICS BY scDblFinder CLASS
# ============================================================================

cat("\n=== DOUBLET VALIDATION 1: QC METRICS BY CLASS ===\n\n")

doublet_qc_metrics <- seurat_filtered@meta.data %>%
  
  mutate(
    Doublet_Status = ifelse(
      scDblFinder.class == "doublet",
      "Predicted doublet",
      "Predicted singlet"
    )
  ) %>%
  
  group_by(
    Doublet_Status
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Median_Genes = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_UMIs = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_MT = median(
      percent.mt,
      na.rm = TRUE
    ),
    
    Median_Doublet_Score = median(
      scDblFinder.score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(
  doublet_qc_metrics
)

write_csv(
  doublet_qc_metrics,
  "results/tables/doublet_validation_qc_metrics.csv"
)

# ============================================================================
# DOUBLEt VALIDATION CHECK 2
# GENE AND UMI DISTRIBUTIONS
# ============================================================================

doublet_plot_data <- seurat_filtered@meta.data %>%
  
  mutate(
    Doublet_Status = factor(
      ifelse(
        scDblFinder.class == "doublet",
        "Predicted doublet",
        "Predicted singlet"
      ),
      levels = c(
        "Predicted singlet",
        "Predicted doublet"
      )
    )
  )

p_doublet_genes <- ggplot(
  doublet_plot_data,
  aes(
    x = Doublet_Status,
    y = nFeature_RNA
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Detected genes by predicted doublet status",
    x = NULL,
    y = "nFeature_RNA"
  )

ggsave(
  "results/figures/10_doublet_validation_genes.png",
  p_doublet_genes,
  width = 8,
  height = 6,
  dpi = 300
)


p_doublet_umis <- ggplot(
  doublet_plot_data,
  aes(
    x = Doublet_Status,
    y = nCount_RNA
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "UMI counts by predicted doublet status",
    x = NULL,
    y = "nCount_RNA"
  )

ggsave(
  "results/figures/11_doublet_validation_umis.png",
  p_doublet_umis,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================================
# DOUBLEt VALIDATION CHECK 3
# CONTRIBUTION TO TOTAL DOUBLET BURDEN
# ============================================================================

doublet_contribution <- seurat_filtered@meta.data %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Predicted_Doublets = sum(
      scDblFinder.class == "doublet",
      na.rm = TRUE
    ),
    
    Doublet_Rate = 100 *
      Predicted_Doublets /
      Cells,
    
    .groups = "drop"
  ) %>%
  
  mutate(
    Percent_of_All_Doublets = 100 *
      Predicted_Doublets /
      sum(Predicted_Doublets)
  ) %>%
  
  arrange(
    desc(Predicted_Doublets)
  )

print(
  doublet_contribution,
  n = Inf
)

write_csv(
  doublet_contribution,
  "results/tables/doublet_contribution_by_sample.csv"
)

# ============================================================================
# DOUBLEt VALIDATION CHECK 4
# SCORE DISTRIBUTION BY CLASS
# ============================================================================

doublet_score_summary <- seurat_filtered@meta.data %>%
  
  group_by(
    scDblFinder.class
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Min_Score = min(
      scDblFinder.score,
      na.rm = TRUE
    ),
    
    Q1_Score = quantile(
      scDblFinder.score,
      0.25,
      na.rm = TRUE
    ),
    
    Median_Score = median(
      scDblFinder.score,
      na.rm = TRUE
    ),
    
    Mean_Score = mean(
      scDblFinder.score,
      na.rm = TRUE
    ),
    
    Q3_Score = quantile(
      scDblFinder.score,
      0.75,
      na.rm = TRUE
    ),
    
    Max_Score = max(
      scDblFinder.score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(
  doublet_score_summary
)

write_csv(
  doublet_score_summary,
  "results/tables/doublet_score_summary_by_class.csv"
)

p_doublet_score_by_class <- ggplot(
  doublet_plot_data,
  aes(
    x = scDblFinder.class,
    y = scDblFinder.score
  )
) +
  
  geom_violin(
    trim = FALSE
  ) +
  
  theme_minimal() +
  
  labs(
    title = "scDblFinder score by predicted class",
    x = "Predicted class",
    y = "Doublet score"
  )

ggsave(
  "results/figures/12_doublet_score_by_class.png",
  p_doublet_score_by_class,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================================
# DOUBLEt VALIDATION CHECK 5
# DOUBLEt SCORE VS RNA COMPLEXITY
# ============================================================================

p_doublet_complexity <- ggplot(
  doublet_plot_data,
  aes(
    x = scDblFinder.score,
    y = nFeature_RNA,
    color = scDblFinder.class
  )
) +
  
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Doublet score versus detected genes",
    x = "scDblFinder score",
    y = "nFeature_RNA",
    color = "Predicted class"
  )

ggsave(
  "results/figures/13_doublet_score_vs_genes.png",
  p_doublet_complexity,
  width = 9,
  height = 7,
  dpi = 300
)


p_doublet_umi_score <- ggplot(
  doublet_plot_data,
  aes(
    x = scDblFinder.score,
    y = nCount_RNA,
    color = scDblFinder.class
  )
) +
  
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Doublet score versus UMI counts",
    x = "scDblFinder score",
    y = "nCount_RNA",
    color = "Predicted class"
  )

ggsave(
  "results/figures/14_doublet_score_vs_umis.png",
  p_doublet_umi_score,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================================
# 28. REMOVE PREDICTED DOUBLETS
# ============================================================================

cat("\n=== STEP 18: REMOVING PREDICTED DOUBLETS ===\n\n")

if (
  !"scDblFinder.class" %in%
  colnames(seurat_filtered@meta.data)
) {
  
  stop(
    "scDblFinder classification unavailable. ",
    "Cannot proceed with doublet removal."
  )
  
}

cells_before_doublet_removal <- ncol(
  seurat_filtered
)

predicted_doublets <- sum(
  seurat_filtered$scDblFinder.class ==
    "doublet",
  na.rm = TRUE
)

predicted_singlets <- sum(
  seurat_filtered$scDblFinder.class ==
    "singlet",
  na.rm = TRUE
)

seurat_singlets <- subset(
  seurat_filtered,
  subset = scDblFinder.class == "singlet"
)

cells_after_doublet_removal <- ncol(
  seurat_singlets
)

cells_removed_doublets <- (
  cells_before_doublet_removal -
    cells_after_doublet_removal
)

pct_retained_after_doublet_removal <- (
  100 *
    cells_after_doublet_removal /
    cells_before_doublet_removal
)

cat(
  "Cells before doublet removal:",
  cells_before_doublet_removal,
  "\n"
)

cat(
  "Predicted doublets removed:",
  predicted_doublets,
  "\n"
)

cat(
  "Predicted singlets retained:",
  predicted_singlets,
  "\n"
)

cat(
  "Cells after doublet removal:",
  cells_after_doublet_removal,
  "\n"
)

cat(
  "Percent retained after doublet removal:",
  round(
    pct_retained_after_doublet_removal,
    2
  ),
  "%\n\n"
)


# ============================================================================
# 29. FINAL QC STATUS TABLE
# ============================================================================

cat("\n=== STEP 19: FINAL QC STATUS ===\n\n")

final_qc_status <- tibble(
  
  Metric = c(
    
    "Input cells",
    
    "Post-fixed-QC cells",
    "Cells removed by fixed QC",
    "Percent retained after fixed QC",
    
    "Predicted doublets removed",
    "Predicted singlets retained",
    "Percent retained after doublet removal",
    
    "Final atlas input cells",
    "Total cells removed",
    "Percent retained overall",
    
    "Samples",
    "Donors",
    "Clinical states",
    
    "Smallest sample before QC (cells)",
    
    "Low-gene cells removed",
    "High-gene cells removed",
    "High-MT cells removed"
  ),
  
  Value = c(
    
    cells_before,
    
    cells_after_qc,
    cells_lost_qc,
    round(
      pct_retained_qc,
      2
    ),
    
    predicted_doublets,
    predicted_singlets,
    round(
      pct_retained_after_doublet_removal,
      2
    ),
    
    cells_after_doublet_removal,
    
    cells_before -
      cells_after_doublet_removal,
    
    round(
      100 *
        cells_after_doublet_removal /
        cells_before,
      2
    ),
    
    dplyr::n_distinct(
      seurat_singlets$GSM
    ),
    
    dplyr::n_distinct(
      seurat_singlets$Donor
    ),
    
    dplyr::n_distinct(
      seurat_singlets$Phase
    ),
    
    min(
      sample_qc_before$Cells
    ),
    
    cells_low_genes,
    cells_high_genes,
    cells_high_mt
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
# 30. SAVE QC + DOUBLET-FILTERED OBJECT
# ============================================================================

cat("\n=== STEP 20: SAVING FINAL QC OBJECT ===\n\n")

saveRDS(
  seurat_singlets,
  "results/rds_objects/seurat_qc_singlets.rds"
)

cat(
  "✓ Saved:\n",
  "results/rds_objects/seurat_qc_singlets.rds\n\n",
  sep = ""
)


# ============================================================================
# 31. SAVE QC SUMMARY
# ============================================================================

qc_summary <- tibble(
  
  metric = c(
    
    "Cells before QC",
    
    "Cells after fixed QC",
    "Cells lost by fixed QC",
    "Percent retained after fixed QC",
    
    "Predicted doublets removed",
    "Predicted singlets retained",
    "Percent retained after doublet removal",
    
    "Final atlas input cells",
    "Total cells removed",
    "Percent retained overall",
    
    "Minimum genes",
    "Maximum genes",
    "Maximum mitochondrial percentage",
    
    "Samples",
    "Donors",
    "Clinical states"
  ),
  
  value = c(
    
    cells_before,
    
    cells_after_qc,
    cells_lost_qc,
    round(
      pct_retained_qc,
      2
    ),
    
    predicted_doublets,
    predicted_singlets,
    round(
      pct_retained_after_doublet_removal,
      2
    ),
    
    cells_after_doublet_removal,
    
    cells_before -
      cells_after_doublet_removal,
    
    round(
      100 *
        cells_after_doublet_removal /
        cells_before,
      2
    ),
    
    min_features,
    max_features,
    max_mt_pct,
    
    dplyr::n_distinct(
      seurat_singlets$GSM
    ),
    
    dplyr::n_distinct(
      seurat_singlets$Donor
    ),
    
    dplyr::n_distinct(
      seurat_singlets$Phase
    )
  )
)

write_csv(
  qc_summary,
  "results/tables/qc_summary.csv"
)


# ============================================================================
# 32. FINAL GATE
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PHASE 2 — SCRIPT 02 COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Input cells:",
  cells_before,
  "\n"
)

cat(
  "Post-fixed-QC cells:",
  cells_after_qc,
  "\n"
)

cat(
  "Predicted doublets removed:",
  predicted_doublets,
  "\n"
)

cat(
  "Final atlas input cells:",
  cells_after_doublet_removal,
  "\n"
)

cat(
  "Overall retention:",
  round(
    100 *
      cells_after_doublet_removal /
      cells_before,
    2
  ),
  "%\n\n"
)

cat(
  "Clinical-state retention after fixed QC:\n"
)

print(
  phase_loss
)

cat("\n")

cat(
  "Doublet decision: REMOVE PREDICTED DOUBLETS\n"
)

cat(
  "Doublet rate:",
  round(
    100 *
      predicted_doublets /
      cells_before_doublet_removal,
    2
  ),
  "%\n\n"
)

cat(
  "QC figures:",
  "results/figures/\n"
)

cat(
  "QC tables:",
  "results/tables/\n"
)

cat(
  "Final atlas object:",
  "results/rds_objects/seurat_qc_singlets.rds\n\n"
)

cat(
  "NEXT GATE:\n"
)

cat(
  "Proceed directly to Phase 3 atlas construction using:\n"
)

cat(
  "seurat_qc_singlets.rds\n"
)

cat("\n")
cat("============================================================\n")
cat("END OF SCRIPT 02\n")
cat("============================================================\n")
