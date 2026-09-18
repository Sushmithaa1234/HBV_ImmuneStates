# ============================================================
# HBV IMMUNE STATES scRNA-seq PROJECT
# PHASE 0 — PRE-ANALYSIS SAMPLE & DATA AUDIT
# ============================================================
#
# Purpose:
#   Establish exactly what samples and expression data are
#   available before preprocessing or biological analysis.
#
# SOURCE OF SAMPLE METADATA:
#   GEO GSE182159
#
# IMPORTANT:
#   Clinical-state labels are preserved EXACTLY as supplied
#   by GEO. No manual relabeling or biological reinterpretation
#   is performed in Phase 0.
#
# GEO clinical-state labels observed in this dataset:
#
#   AC = 3 samples
#   AR = 3 samples
#   IA = 5 samples
#   IT = 6 samples
#   NL = 6 samples
#
# This phase does NOT:
#   - filter cells
#   - normalize data
#   - batch-correct data
#   - cluster cells
#   - annotate cells
#   - perform biological analysis
#   - make final QC decisions
#   - create a Seurat object
#
# IMPORTANT DATA NOTE:
#   GEO states that the supplementary expression matrices
#   contain log-counts-per-10,000 values generated after
#   normalization and logarithmic transformation.
#
#   Raw count matrices are NOT provided in GSE182159.
#
# ============================================================


# ------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------

library(tidyverse)
library(GEOquery)

setwd("D:/HBV scRNA")

raw_dir <- "Data/Raw"

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


cat("============================================\n")
cat("PHASE 0 — PRE-ANALYSIS SAMPLE & DATA AUDIT\n")
cat("============================================\n\n")


# ------------------------------------------------------------
# 2. IDENTIFY DOWNLOADED LIVER FILES
# ------------------------------------------------------------

liver_files <- list.files(
  raw_dir,
  pattern = "_Liver_.*\\.txt(\\.gz)?$",
  full.names = TRUE
)

liver_gsm <- str_extract(
  basename(liver_files),
  "GSM[0-9]+"
)

cat(
  "Liver expression files found:",
  length(liver_files),
  "\n"
)

if (length(liver_files) != 23) {
  
  stop(
    "ERROR: Expected 23 liver expression files, but found ",
    length(liver_files),
    ". Check Data/Raw."
  )
  
}

if (anyDuplicated(liver_gsm) > 0) {
  
  stop(
    "ERROR: Duplicate GSM IDs detected among downloaded files."
  )
  
}

cat("✓ Exactly 23 unique liver files detected\n")
cat("✓ No duplicate GSM IDs detected\n\n")


# ------------------------------------------------------------
# 3. LOAD AUTHORITATIVE GEO METADATA
# ------------------------------------------------------------

cat("Retrieving GEO metadata for GSE182159...\n\n")

gse <- getGEO(
  "GSE182159",
  GSEMatrix = TRUE
)

if (length(gse) != 1) {
  
  stop(
    "ERROR: Expected exactly one GEO Series Matrix object."
  )
  
}

pdat <- pData(gse[[1]])


# ------------------------------------------------------------
# 4. EXTRACT DOWNLOADED LIVER SAMPLES FROM GEO
# ------------------------------------------------------------
#
# IMPORTANT:
#   GEO metadata are retained exactly as provided.
#   No clinical-state labels are manually created or changed.
# ------------------------------------------------------------

geo_liver_metadata <- pdat %>%
  
  filter(
    `tissue:ch1` == "Liver",
    geo_accession %in% liver_gsm
  ) %>%
  
  select(
    GSM = geo_accession,
    donor = `donor:ch1`,
    phase = `Stage:ch1`,
    tissue = `tissue:ch1`,
    platform = platform_id,
    instrument = instrument_model,
    library_selection,
    library_source,
    library_strategy
  )


# ------------------------------------------------------------
# 5. VERIFY GEO SAMPLE COUNT
# ------------------------------------------------------------

cat(
  "Downloaded liver samples:",
  length(liver_gsm),
  "\n"
)

cat(
  "Matching liver samples in GEO:",
  nrow(geo_liver_metadata),
  "\n\n"
)

if (nrow(geo_liver_metadata) != 23) {
  
  stop(
    "ERROR: GEO metadata does not contain exactly 23 ",
    "matching liver samples."
  )
  
}


# ------------------------------------------------------------
# 6. BIDIRECTIONAL SAMPLE MATCHING
# ------------------------------------------------------------

missing_from_geo <- setdiff(
  liver_gsm,
  geo_liver_metadata$GSM
)

missing_from_files <- setdiff(
  geo_liver_metadata$GSM,
  liver_gsm
)


if (length(missing_from_geo) > 0) {
  
  cat(
    "ERROR: Downloaded files missing from GEO metadata:\n"
  )
  
  print(missing_from_geo)
  
  stop(
    "Sample matching failed."
  )
  
}


if (length(missing_from_files) > 0) {
  
  cat(
    "ERROR: GEO liver samples missing from downloaded files:\n"
  )
  
  print(missing_from_files)
  
  stop(
    "Sample matching failed."
  )
  
}


cat(
  "✓ Every downloaded liver file has GEO metadata\n"
)

cat(
  "✓ Every GEO liver sample has a downloaded file\n\n"
)


# ------------------------------------------------------------
# 7. VERIFY DUPLICATE GSM IDs
# ------------------------------------------------------------

if (
  anyDuplicated(geo_liver_metadata$GSM) > 0
) {
  
  stop(
    "ERROR: Duplicate GSM IDs detected in GEO metadata."
  )
  
}

cat(
  "✓ GEO GSM IDs are unique\n"
)


# ------------------------------------------------------------
# 8. VERIFY DONOR UNIQUENESS
# ------------------------------------------------------------

if (
  anyDuplicated(geo_liver_metadata$donor) > 0
) {
  
  cat(
    "WARNING: Multiple liver samples are associated ",
    "with the same donor.\n"
  )
  
} else {
  
  cat(
    "✓ Each liver sample corresponds to a unique donor\n"
  )
  
}


# ------------------------------------------------------------
# 9. VERIFY TISSUE
# ------------------------------------------------------------

if (
  !all(
    geo_liver_metadata$tissue == "Liver"
  )
) {
  
  stop(
    "ERROR: Non-liver sample entered the liver cohort."
  )
  
}

cat(
  "✓ All 23 samples are annotated as liver tissue\n\n"
)


# ------------------------------------------------------------
# 10. REPORT GEO CLINICAL-STATE LABELS
# ------------------------------------------------------------
#
# These labels are preserved exactly as supplied by GEO.
# No relabeling is performed.
# ------------------------------------------------------------

cat("============================================\n")
cat("GEO CLINICAL-STATE COMPOSITION\n")
cat("============================================\n\n")

print(
  table(
    geo_liver_metadata$phase
  )
)

cat("\n")


# ------------------------------------------------------------
# 11. CHECK FOR MISSING CLINICAL METADATA
# ------------------------------------------------------------

if (
  anyNA(geo_liver_metadata$phase)
) {
  
  stop(
    "ERROR: Missing clinical-state metadata detected."
  )
  
}

if (
  anyNA(geo_liver_metadata$donor)
) {
  
  stop(
    "ERROR: Missing donor metadata detected."
  )
  
}

cat(
  "✓ No missing clinical-state labels\n"
)

cat(
  "✓ No missing donor IDs\n\n"
)


# ------------------------------------------------------------
# 12. READ AND AUDIT EACH EXPRESSION MATRIX
# ------------------------------------------------------------

cat("============================================\n")
cat("AUDITING EXPRESSION MATRICES\n")
cat("============================================\n\n")


audit_list <- list()


for (i in seq_along(liver_files)) {
  
  file <- liver_files[i]
  
  gsm <- liver_gsm[i]
  
  cat(
    "Reading:",
    gsm,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Read the actual downloaded file format
  #
  # The downloaded files are space-delimited text matrices.
  # Rows = genes
  # Columns = cells
  # ----------------------------------------------------------
  
  expression_matrix <- read.table(
    file,
    sep = " ",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  
  
  # Convert to matrix
  
  expression_matrix <- as.matrix(
    expression_matrix
  )
  
  
  # ----------------------------------------------------------
  # Basic structural checks
  # ----------------------------------------------------------
  
  n_genes <- nrow(
    expression_matrix
  )
  
  n_cells <- ncol(
    expression_matrix
  )
  
  
  if (
    !is.numeric(expression_matrix)
  ) {
    
    stop(
      "ERROR: Expression matrix for ",
      gsm,
      " is not numeric."
    )
    
  }
  
  
  if (
    anyNA(expression_matrix)
  ) {
    
    stop(
      "ERROR: NA values detected in expression matrix for ",
      gsm
    )
    
  }
  
  
  if (
    any(
      !is.finite(expression_matrix)
    )
  ) {
    
    stop(
      "ERROR: Non-finite values detected in expression matrix for ",
      gsm
    )
    
  }
  
  
  if (
    any(
      expression_matrix < 0
    )
  ) {
    
    stop(
      "ERROR: Negative expression values detected in ",
      gsm
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Gene identifier checks
  # ----------------------------------------------------------
  
  gene_ids <- rownames(
    expression_matrix
  )
  
  
  if (
    anyNA(gene_ids) ||
    any(gene_ids == "")
  ) {
    
    stop(
      "ERROR: Missing gene identifiers detected in ",
      gsm
    )
    
  }
  
  
  duplicated_genes <- sum(
    duplicated(gene_ids)
  )
  
  
  # ----------------------------------------------------------
  # Cell barcode checks
  # ----------------------------------------------------------
  
  cell_ids <- colnames(
    expression_matrix
  )
  
  
  if (
    anyNA(cell_ids) ||
    any(cell_ids == "")
  ) {
    
    stop(
      "ERROR: Missing cell identifiers detected in ",
      gsm
    )
    
  }
  
  
  duplicated_cells <- sum(
    duplicated(cell_ids)
  )
  
  
  if (
    duplicated_cells > 0
  ) {
    
    stop(
      "ERROR: Duplicate cell identifiers detected in ",
      gsm
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Genes detected per cell
  #
  # This is a descriptive metric.
  #
  # It is NOT equivalent to nFeature_RNA calculated from
  # raw counts in a newly generated Seurat object.
  # ----------------------------------------------------------
  
  genes_detected_per_cell <- colSums(
    expression_matrix > 0
  )
  
  
  # ----------------------------------------------------------
  # Expression-scale diagnostics
  # ----------------------------------------------------------
  
  min_expression <- min(
    expression_matrix
  )
  
  median_expression <- median(
    expression_matrix
  )
  
  max_expression <- max(
    expression_matrix
  )
  
  
  # ----------------------------------------------------------
  # Create audit record
  # ----------------------------------------------------------
  
  audit_list[[gsm]] <- tibble(
    
    GSM = gsm,
    
    genes_in_matrix =
      n_genes,
    
    cells_in_matrix =
      n_cells,
    
    duplicated_gene_ids =
      duplicated_genes,
    
    duplicated_cell_ids =
      duplicated_cells,
    
    min_expression_value =
      min_expression,
    
    median_expression_value =
      median_expression,
    
    max_expression_value =
      max_expression,
    
    median_genes_detected_per_cell =
      median(
        genes_detected_per_cell
      ),
    
    mean_genes_detected_per_cell =
      mean(
        genes_detected_per_cell
      ),
    
    min_genes_detected_per_cell =
      min(
        genes_detected_per_cell
      ),
    
    max_genes_detected_per_cell =
      max(
        genes_detected_per_cell
      )
    
  )
  
  
  cat(
    "  Genes:",
    n_genes,
    "| Cells:",
    n_cells,
    "| Median genes detected/cell:",
    round(
      median(
        genes_detected_per_cell
      ),
      1
    ),
    "\n"
  )
  
  
  cat(
    "  Expression range:",
    round(
      min_expression,
      3
    ),
    "to",
    round(
      max_expression,
      3
    ),
    "\n\n"
  )
  
}


# ------------------------------------------------------------
# 13. COMBINE MATRIX AUDIT RESULTS
# ------------------------------------------------------------

sample_audit <- bind_rows(
  audit_list
)


# ------------------------------------------------------------
# 14. MERGE GEO METADATA WITH MATRIX AUDIT
# ------------------------------------------------------------

final_metadata <- geo_liver_metadata %>%
  
  left_join(
    sample_audit,
    by = "GSM"
  )


# ------------------------------------------------------------
# 15. VERIFY MATRIX AUDIT COMPLETENESS
# ------------------------------------------------------------

required_audit_columns <- c(
  "genes_in_matrix",
  "cells_in_matrix",
  "median_genes_detected_per_cell"
)


for (column in required_audit_columns) {
  
  if (
    any(
      is.na(
        final_metadata[[column]]
      )
    )
  ) {
    
    stop(
      "ERROR: Missing audit information in column: ",
      column
    )
    
  }
  
}


# ------------------------------------------------------------
# 16. REPORT DUPLICATED GENE IDs
# ------------------------------------------------------------

if (
  any(
    final_metadata$duplicated_gene_ids > 0
  )
) {
  
  cat(
    "WARNING: Duplicated gene IDs detected.\n\n"
  )
  
  print(
    final_metadata %>%
      filter(
        duplicated_gene_ids > 0
      ) %>%
      select(
        GSM,
        duplicated_gene_ids
      )
  )
  
} else {
  
  cat(
    "✓ No duplicated gene IDs detected\n"
  )
  
}


# ------------------------------------------------------------
# 17. DOCUMENT DATA REPRESENTATION
# ------------------------------------------------------------

final_metadata <- final_metadata %>%
  
  mutate(
    
    data_type =
      "processed_log_normalized_expression",
    
    normalization =
      "log-counts-per-10,000",
    
    raw_counts_available =
      FALSE
    
  )


# ------------------------------------------------------------
# 18. FINAL SAMPLE AUDIT TABLE
# ------------------------------------------------------------

final_metadata <- final_metadata %>%
  
  arrange(
    phase,
    donor
  )


cat("\n============================================\n")
cat("FINAL SAMPLE AUDIT\n")
cat("============================================\n\n")

print(
  final_metadata
)


# ------------------------------------------------------------
# 19. SUMMARY BY GEO CLINICAL STATE
# ------------------------------------------------------------

phase_summary <- final_metadata %>%
  
  group_by(
    phase
  ) %>%
  
  summarise(
    
    samples =
      n(),
    
    total_cells =
      sum(
        cells_in_matrix
      ),
    
    median_cells_per_sample =
      median(
        cells_in_matrix
      ),
    
    min_cells_per_sample =
      min(
        cells_in_matrix
      ),
    
    max_cells_per_sample =
      max(
        cells_in_matrix
      ),
    
    median_genes_detected_per_cell =
      median(
        median_genes_detected_per_cell
      ),
    
    min_sample_median_genes =
      min(
        median_genes_detected_per_cell
      ),
    
    max_sample_median_genes =
      max(
        median_genes_detected_per_cell
      ),
    
    .groups = "drop"
    
  )


# ------------------------------------------------------------
# 20. SAVE AUDIT TABLES
# ------------------------------------------------------------

write_csv(
  final_metadata,
  "results/tables/phase0_final_liver_sample_metadata.csv"
)

write_csv(
  sample_audit,
  "results/tables/phase0_expression_matrix_audit.csv"
)

write_csv(
  phase_summary,
  "results/tables/phase0_geo_phase_summary.csv"
)


# ------------------------------------------------------------
# 21. PLOT — CELLS PER SAMPLE
# ------------------------------------------------------------

plot_cells <- ggplot(
  final_metadata,
  aes(
    x = reorder(
      GSM,
      cells_in_matrix
    ),
    y = cells_in_matrix,
    fill = phase
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    title =
      "Number of Cells per Liver Sample",
    x =
      "Sample",
    y =
      "Number of Cells"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom"
  )


ggsave(
  "results/figures/phase0_cells_per_sample.png",
  plot_cells,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 22. PLOT — MEDIAN GENES DETECTED PER CELL
# ------------------------------------------------------------

plot_genes <- ggplot(
  final_metadata,
  aes(
    x = reorder(
      GSM,
      median_genes_detected_per_cell
    ),
    y = median_genes_detected_per_cell,
    fill = phase
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    title =
      "Median Genes Detected per Cell",
    x =
      "Sample",
    y =
      "Median Genes Detected per Cell"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom"
  )


ggsave(
  "results/figures/phase0_median_genes_detected.png",
  plot_genes,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 23. FINAL REPORT
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PHASE 0 AUDIT COMPLETE\n")
cat("============================================\n\n")


cat(
  "Liver samples:",
  nrow(final_metadata),
  "\n"
)


cat(
  "Unique donors:",
  n_distinct(
    final_metadata$donor
  ),
  "\n"
)


cat(
  "\nGEO clinical-state composition:\n"
)

print(
  table(
    final_metadata$phase
  )
)


cat(
  "\nData type:",
  unique(
    final_metadata$data_type
  ),
  "\n"
)


cat(
  "Normalization:",
  unique(
    final_metadata$normalization
  ),
  "\n"
)


cat(
  "Raw counts available:",
  unique(
    final_metadata$raw_counts_available
  ),
  "\n"
)


cat(
  "\nPlatform:",
  paste(
    unique(
      final_metadata$platform
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Instrument:",
  paste(
    unique(
      final_metadata$instrument
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Library strategy:",
  paste(
    unique(
      final_metadata$library_strategy
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Library selection:",
  paste(
    unique(
      final_metadata$library_selection
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Library source:",
  paste(
    unique(
      final_metadata$library_source
    ),
    collapse = ", "
  ),
  "\n"
)


cat("\n============================================\n")
cat("NO CELLS HAVE BEEN FILTERED.\n")
cat("NO QC THRESHOLDS HAVE BEEN APPLIED.\n")
cat("NO NORMALIZATION HAS BEEN PERFORMED.\n")
cat("NO BATCH CORRECTION HAS BEEN PERFORMED.\n")
cat("NO CLINICAL-STATE LABELS HAVE BEEN MANUALLY ALTERED.\n")
cat("NO BIOLOGICAL ANALYSIS HAS BEEN PERFORMED.\n")
cat("============================================\n\n")


cat(
  "Phase 0 is ready for review before Phase 1.\n"
)
