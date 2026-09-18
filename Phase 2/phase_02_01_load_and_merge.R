# ============================================================================
# HBV Immune States scRNA-seq Analysis
# Phase 2 — Script 01: Load and Merge All 23 Liver Samples
# ============================================================================
#
# Purpose:
#   Load the 23 GEO-verified liver expression matrices into Seurat,
#   preserve sample/donor/clinical-state metadata, verify data integrity,
#   and merge the samples into a single Seurat object.
#
# IMPORTANT:
#   The GEO matrices contain processed log-counts-per-10,000 expression
#   values, NOT raw UMI count matrices.
#
#   Therefore:
#     - No normalization is performed here.
#     - No cells or genes are filtered here.
#     - Expression values are stored in the Seurat v5 "data" layer.
#     - Raw UMI/count-based quantities are NOT reconstructed.
#
# ============================================================================


# ============================================================================
# 1. LOAD PACKAGES
# ============================================================================

library(Seurat)
library(tidyverse)
library(here)


# ============================================================================
# 2. SET WORKING DIRECTORY
# ============================================================================

setwd(here())

cat("============================================================\n")
cat("PHASE 2: SCRIPT 01 — LOAD AND MERGE\n")
cat("============================================================\n\n")


# ============================================================================
# 3. LOAD PHASE 0 METADATA
# ============================================================================

cat("Loading Phase 0 metadata...\n")

metadata <- read.csv(
  "results/tables/phase0_final_liver_sample_metadata.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "✓ Metadata loaded:",
  nrow(metadata),
  "samples\n"
)

cat(
  "  Clinical states:",
  paste(unique(metadata$phase), collapse = ", "),
  "\n"
)

cat(
  "  Unique donors:",
  n_distinct(metadata$donor),
  "\n\n"
)


# ----------------------------------------------------------------------------
# Verify metadata structure
# ----------------------------------------------------------------------------

if (nrow(metadata) != 23) {
  stop(
    "ERROR: Expected 23 metadata records, found ",
    nrow(metadata),
    "."
  )
}

if (n_distinct(metadata$GSM) != 23) {
  stop(
    "ERROR: Expected 23 unique GSM IDs in Phase 0 metadata."
  )
}

if (n_distinct(metadata$donor) != 23) {
  stop(
    "ERROR: Expected 23 unique donors in Phase 0 metadata."
  )
}

required_metadata_columns <- c(
  "GSM",
  "donor",
  "phase",
  "tissue"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "ERROR: Required metadata columns missing: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

if (anyNA(metadata$GSM) ||
    anyNA(metadata$donor) ||
    anyNA(metadata$phase) ||
    anyNA(metadata$tissue)) {
  stop(
    "ERROR: Missing required metadata detected."
  )
}

cat("✓ Phase 0 metadata structure verified\n\n")


# ============================================================================
# 4. IDENTIFY LIVER EXPRESSION MATRICES
# ============================================================================

cat("Identifying liver expression matrices...\n\n")

liver_files <- list.files(
  "Data/Raw",
  pattern = "_Liver_.*\\.txt$",
  full.names = TRUE
)

if (length(liver_files) != 23) {
  stop(
    "ERROR: Expected 23 liver expression files, found ",
    length(liver_files),
    "."
  )
}

liver_gsm <- str_extract(
  basename(liver_files),
  "GSM[0-9]+"
)

if (any(is.na(liver_gsm))) {
  stop(
    "ERROR: Could not extract GSM ID from one or more filenames."
  )
}

if (anyDuplicated(liver_gsm) > 0) {
  stop(
    "ERROR: Duplicate GSM IDs detected among expression files."
  )
}

cat(
  "✓ Found",
  length(liver_files),
  "liver expression matrices\n"
)

cat("✓ GSM identifiers are unique\n\n")


# ============================================================================
# 5. VERIFY FILE ↔ PHASE 0 METADATA MATCH
# ============================================================================

cat("Verifying expression-file / metadata correspondence...\n")

missing_metadata_for_files <- setdiff(
  liver_gsm,
  metadata$GSM
)

missing_files_for_metadata <- setdiff(
  metadata$GSM,
  liver_gsm
)

if (length(missing_metadata_for_files) > 0) {
  stop(
    "ERROR: Expression files without Phase 0 metadata: ",
    paste(missing_metadata_for_files, collapse = ", ")
  )
}

if (length(missing_files_for_metadata) > 0) {
  stop(
    "ERROR: Phase 0 samples without expression files: ",
    paste(missing_files_for_metadata, collapse = ", ")
  )
}

cat("✓ All 23 expression files match Phase 0 metadata\n\n")


# ============================================================================
# 6. INITIALIZE OBJECT LIST AND LOAD SUMMARY
# ============================================================================

seurat_list <- list()

load_summary <- data.frame(
  GSM = character(),
  Donor = character(),
  Phase = character(),
  Genes = integer(),
  Cells = integer(),
  stringsAsFactors = FALSE
)


# ============================================================================
# 7. LOAD EACH PROCESSED EXPRESSION MATRIX
# ============================================================================

cat("============================================================\n")
cat("LOADING EXPRESSION MATRICES\n")
cat("============================================================\n\n")


for (file in liver_files) {
  
  # --------------------------------------------------------------------------
  # Extract GSM ID
  # --------------------------------------------------------------------------
  
  gsm_id <- basename(file) %>%
    str_extract("GSM[0-9]+")
  
  
  # --------------------------------------------------------------------------
  # Retrieve Phase 0 metadata
  # --------------------------------------------------------------------------
  
  sample_meta <- metadata %>%
    filter(GSM == gsm_id)
  
  if (nrow(sample_meta) != 1) {
    stop(
      "ERROR: Expected exactly one metadata record for ",
      gsm_id,
      ", found ",
      nrow(sample_meta),
      "."
    )
  }
  
  donor_id <- sample_meta$donor[[1]]
  phase <- sample_meta$phase[[1]]
  
  
  cat(
    "Loading ",
    gsm_id,
    " (",
    donor_id,
    ", ",
    phase,
    ")...\n",
    sep = ""
  )
  
  
  # --------------------------------------------------------------------------
  # Read processed expression matrix
  #
  # The local GEO supplementary files are space-delimited.
  # --------------------------------------------------------------------------
  
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
  
  expression_matrix <- as.matrix(expression_matrix)
  
  
  # --------------------------------------------------------------------------
  # Validate expression matrix
  # --------------------------------------------------------------------------
  
  if (!is.numeric(expression_matrix)) {
    stop(
      "ERROR: Expression matrix for ",
      gsm_id,
      " is not numeric."
    )
  }
  
  if (nrow(expression_matrix) == 0 ||
      ncol(expression_matrix) == 0) {
    stop(
      "ERROR: Empty expression matrix for ",
      gsm_id,
      "."
    )
  }
  
  if (anyNA(expression_matrix)) {
    stop(
      "ERROR: NA values detected in ",
      gsm_id,
      "."
    )
  }
  
  if (any(!is.finite(expression_matrix))) {
    stop(
      "ERROR: Non-finite values detected in ",
      gsm_id,
      "."
    )
  }
  
  if (any(expression_matrix < 0)) {
    stop(
      "ERROR: Negative expression values detected in ",
      gsm_id,
      "."
    )
  }
  
  if (anyDuplicated(rownames(expression_matrix)) > 0) {
    stop(
      "ERROR: Duplicate gene IDs detected in ",
      gsm_id,
      "."
    )
  }
  
  if (anyDuplicated(colnames(expression_matrix)) > 0) {
    stop(
      "ERROR: Duplicate cell IDs detected in ",
      gsm_id,
      "."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Verify expected gene count
  # --------------------------------------------------------------------------
  
  if (nrow(expression_matrix) != 24452) {
    stop(
      "ERROR: Expected 24,452 genes in ",
      gsm_id,
      ", found ",
      nrow(expression_matrix),
      "."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Create Seurat v5 assay
  #
  # IMPORTANT:
  # The matrix contains processed log-counts-per-10,000 values.
  #
  # These values are stored in the DATA layer.
  # They are NOT treated as raw counts.
  # --------------------------------------------------------------------------
  
  expression_assay <- CreateAssay5Object(
    data = expression_matrix,
    min.cells = 0,
    min.features = 0
  )
  
  
  # --------------------------------------------------------------------------
  # Create Seurat object
  # --------------------------------------------------------------------------
  
  seurat_obj <- CreateSeuratObject(
    counts = expression_assay,
    project = gsm_id
  )
  
  
  # --------------------------------------------------------------------------
  # Add sample-level metadata
  # --------------------------------------------------------------------------
  
  seurat_obj$GSM <- gsm_id
  seurat_obj$Donor <- donor_id
  seurat_obj$Phase <- phase
  
  
  # --------------------------------------------------------------------------
  # Store object
  # --------------------------------------------------------------------------
  
  seurat_list[[gsm_id]] <- seurat_obj
  
  
  # --------------------------------------------------------------------------
  # Record loading summary
  # --------------------------------------------------------------------------
  
  load_summary <- rbind(
    load_summary,
    data.frame(
      GSM = gsm_id,
      Donor = donor_id,
      Phase = phase,
      Genes = nrow(seurat_obj),
      Cells = ncol(seurat_obj),
      stringsAsFactors = FALSE
    )
  )
  
  
  cat(
    "  ✓ ",
    gsm_id,
    " — ",
    nrow(seurat_obj),
    " genes, ",
    ncol(seurat_obj),
    " cells\n",
    sep = ""
  )
  
}


cat(
  "\n✓ Successfully loaded ",
  length(seurat_list),
  " samples\n\n",
  sep = ""
)


# ============================================================================
# 8. VERIFY LOADED CELL TOTAL
# ============================================================================

cat("Verifying loaded cell total...\n")

expected_cells <- sum(
  metadata$cells_in_matrix
)

loaded_cells <- sum(
  load_summary$Cells
)

cat(
  "  Phase 0 expected:",
  expected_cells,
  "\n"
)

cat(
  "  Loaded:",
  loaded_cells,
  "\n"
)

if (loaded_cells != expected_cells) {
  stop(
    "ERROR: Loaded cell count does not match Phase 0 audit."
  )
}

cat("✓ Loaded cell total matches Phase 0: 106,592\n\n")


# ============================================================================
# 9. MERGE ALL 23 SAMPLES
# ============================================================================

cat("============================================================\n")
cat("MERGING SAMPLES\n")
cat("============================================================\n\n")

if (length(seurat_list) != 23) {
  stop(
    "ERROR: Expected 23 Seurat objects, found ",
    length(seurat_list),
    "."
  )
}

seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "HBV_Liver"
)


cat("✓ Samples merged successfully\n\n")

cat("Merged object dimensions:\n")
cat("  Genes:", nrow(seurat_merged), "\n")
cat("  Cells:", ncol(seurat_merged), "\n\n")


# ============================================================================
# 10. VERIFY MERGED OBJECT
# ============================================================================

cat("============================================================\n")
cat("VERIFYING MERGED OBJECT\n")
cat("============================================================\n\n")


# ----------------------------------------------------------------------------
# Dimensions
# ----------------------------------------------------------------------------

if (nrow(seurat_merged) != 24452) {
  stop(
    "ERROR: Expected 24,452 genes after merging."
  )
}

if (ncol(seurat_merged) != expected_cells) {
  stop(
    "ERROR: Merged cell count does not match Phase 0."
  )
}


# ----------------------------------------------------------------------------
# Cell ID uniqueness
# ----------------------------------------------------------------------------

if (anyDuplicated(colnames(seurat_merged)) > 0) {
  stop(
    "ERROR: Duplicate cell IDs detected after merging."
  )
}


# ----------------------------------------------------------------------------
# Sample metadata
# ----------------------------------------------------------------------------

if (anyNA(seurat_merged$GSM) ||
    anyNA(seurat_merged$Donor) ||
    anyNA(seurat_merged$Phase)) {
  stop(
    "ERROR: Missing sample metadata detected after merging."
  )
}


# ----------------------------------------------------------------------------
# Verify expected number of samples and donors
# ----------------------------------------------------------------------------

if (n_distinct(seurat_merged$GSM) != 23) {
  stop(
    "ERROR: Expected 23 GSM identifiers after merging."
  )
}

if (n_distinct(seurat_merged$Donor) != 23) {
  stop(
    "ERROR: Expected 23 donors after merging."
  )
}


cat("✓ Genes:", nrow(seurat_merged), "\n")
cat("✓ Cells:", ncol(seurat_merged), "\n")
cat("✓ Unique GSMs:", n_distinct(seurat_merged$GSM), "\n")
cat("✓ Unique donors:", n_distinct(seurat_merged$Donor), "\n")
cat("✓ Cell IDs are unique\n\n")


# ============================================================================
# 11. VERIFY CLINICAL-STATE METADATA
# ============================================================================

cat("=== VERIFYING CLINICAL-STATE METADATA ===\n\n")

# Cell distribution in the merged Seurat object
cell_phase_distribution <- table(seurat_merged$Phase)

cat("Clinical-state cell distribution:\n")
print(cell_phase_distribution)
cat("\n")

# Expected number of samples per clinical state from Phase 0
expected_phase_counts <- c(
  AC = 3,
  AR = 3,
  IA = 5,
  IT = 6,
  NL = 6
)

# Validate sample-level distribution using the original metadata
sample_phase_distribution <- table(metadata$phase)

cat("Clinical-state sample distribution:\n")
print(sample_phase_distribution)
cat("\n")

if (!identical(
  as.integer(sample_phase_distribution[names(expected_phase_counts)]),
  as.integer(expected_phase_counts)
)) {
  stop(
    "ERROR: Clinical-state SAMPLE distribution does not match Phase 0."
  )
}

cat("✓ Clinical-state sample distribution matches Phase 0.\n\n")

sample_cell_counts <- table(seurat_merged$GSM)

cat("Cells per GSM:\n")
print(sample_cell_counts)
cat("\n")

if (length(sample_cell_counts) != nrow(metadata)) {
  stop("ERROR: Number of GSMs in merged object does not match Phase 0.")
}

if (!all(names(sample_cell_counts) %in% metadata$GSM)) {
  stop("ERROR: Unexpected GSM detected in merged object.")
}

cat("✓ All 23 GSMs represented in merged object.\n\n")


# ============================================================================
# 12. CELLS PER CLINICAL STATE
# ============================================================================

cat("Cells per clinical state:\n")

phase_summary <- seurat_merged@meta.data %>%
  group_by(Phase) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(Donor),
    .groups = "drop"
  )

print(phase_summary)

cat("\n")


# ============================================================================
# 13. CELLS PER DONOR
# ============================================================================

cat("Cells per donor:\n")

donor_summary <- seurat_merged@meta.data %>%
  group_by(Donor) %>%
  summarise(
    GSM = first(GSM),
    Phase = first(Phase),
    n_cells = n(),
    .groups = "drop"
  )

print(donor_summary, n = 23)

cat("\n✓ Donor/sample metadata verified\n\n")


# ============================================================================
# 14. VERIFY PROCESSED-DATA REPRESENTATION
# ============================================================================

cat("============================================================\n")
cat("DATA REPRESENTATION\n")
cat("============================================================\n\n")

cat("Expression data:\n")
cat("  Type: processed expression\n")
cat("  Normalization: log-counts-per-10,000\n")
cat("  Raw counts available: FALSE\n\n")

cat(
  "No normalization, filtering, or batch correction was performed.\n\n"
)


# ============================================================================
# 15. SAVE CHECKPOINT
# ============================================================================

cat("Saving merged Seurat object...\n")

dir.create(
  "results/rds_objects",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  seurat_merged,
  "results/rds_objects/seurat_merged_processed_expression.rds"
)

cat(
  "✓ Saved to:\n",
  "  results/rds_objects/seurat_merged_processed_expression.rds\n\n",
  sep = ""
)


# ============================================================================
# 16. SAVE LOAD SUMMARY
# ============================================================================

write.csv(
  load_summary,
  "results/tables/phase2_script01_load_summary.csv",
  row.names = FALSE
)

write.csv(
  phase_summary,
  "results/tables/phase2_script01_phase_summary.csv",
  row.names = FALSE
)

cat("✓ Loading summaries saved\n\n")


# ============================================================================
# 17. FINAL REPORT
# ============================================================================

cat("============================================================\n")
cat("PHASE 2 — SCRIPT 01 COMPLETE\n")
cat("============================================================\n\n")

cat("Samples loaded:        ", length(seurat_list), "\n", sep = "")
cat("Unique donors:         ", n_distinct(seurat_merged$Donor), "\n", sep = "")
cat("Genes:                 ", nrow(seurat_merged), "\n", sep = "")
cat("Cells:                 ", ncol(seurat_merged), "\n", sep = "")

cat("\nClinical states:\n")
print(table(seurat_merged$Phase))

cat("\nProcessing performed:\n")
cat("  ✓ Expression matrices loaded\n")
cat("  ✓ Sample metadata attached\n")
cat("  ✓ Samples merged\n")
cat("  ✓ Cell/sample totals verified\n")
cat("  ✓ Clinical-state metadata verified\n")
cat("  ✓ Data representation preserved\n")

cat("\nProcessing NOT performed:\n")
cat("  ✓ No cell filtering\n")
cat("  ✓ No gene filtering\n")
cat("  ✓ No normalization\n")
cat("  ✓ No batch correction\n")
cat("  ✓ No clustering\n")
cat("  ✓ No biological analysis\n")

cat("\nNext: Phase 2 — Script 02 (Formal QC)\n")
cat("============================================================\n")
