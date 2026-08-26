# ============================================================================
# HBV Immune States scRNA-seq Analysis
# Script 01: Load and merge all 23 liver samples
# ============================================================================

library(Seurat)
library(tidyverse)
library(here)

# Set working directory
setwd(here())

cat("=== PHASE 2: Script 01 — Load and Merge ===\n\n")

# ============================================================================
# Step 1: Load metadata
# ============================================================================

cat("Loading sample metadata...\n")
metadata <- read.csv("results/tables/final_liver_sample_metadata.csv", stringsAsFactors = FALSE)

cat("✓ Metadata loaded:", nrow(metadata), "samples\n")
cat("  Clinical states:", paste(unique(metadata$phase), collapse = ", "), "\n")
cat("  Donors:", n_distinct(metadata$donor), "unique\n\n")

# ============================================================================
# Step 2: Load all liver samples
# ============================================================================

cat("Loading expression matrices...\n\n")

liver_files <- list.files("Data/Raw", pattern = "_Liver_.*\\.txt$", full.names = TRUE)

if (length(liver_files) != 23) {
  stop("ERROR: Expected 23 liver files, found", length(liver_files))
}

cat("Found", length(liver_files), "liver samples\n\n")

seurat_list <- list()
load_summary <- data.frame(
  GSM = character(),
  Donor = character(),
  Phase = character(),
  Genes = integer(),
  Cells = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

# ============================================================================
# Step 3: Create Seurat object for each sample
# ============================================================================

for (file in liver_files) {
  # Extract GSM ID from filename
  gsm_id <- basename(file) %>% str_extract("GSM[0-9]+")
  
  # Look up metadata for this sample
  sample_meta <- metadata %>% filter(GSM == gsm_id)
  
  if (nrow(sample_meta) == 0) {
    cat("⚠️  WARNING: No metadata found for", gsm_id, "\n")
    next
  }
  
  donor_id <- sample_meta$donor
  phase <- sample_meta$phase
  qc_flag <- sample_meta$sample_qc_flag
  
  cat("Loading", gsm_id, "(", donor_id, ",", phase, ")...\n")
  
  # Read count matrix (space-delimited)
  counts <- read.table(
    file,
    sep = " ",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  counts <- as.matrix(counts)
  
  # Check for valid data
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    cat("  ⚠️  WARNING: Empty matrix for", gsm_id, "— skipping\n")
    next
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = gsm_id,
    min.cells = 3,
    min.features = 0
  )
  
  # Add metadata
  seurat_obj@meta.data$GSM <- gsm_id
  seurat_obj@meta.data$Donor <- donor_id
  seurat_obj@meta.data$Phase <- phase
  seurat_obj@meta.data$QC_Flag <- qc_flag
  
  # Store in list
  seurat_list[[gsm_id]] <- seurat_obj
  
  # Log summary
  load_summary <- rbind(load_summary, data.frame(
    GSM = gsm_id,
    Donor = donor_id,
    Phase = phase,
    Genes = nrow(seurat_obj),
    Cells = ncol(seurat_obj),
    Status = ifelse(qc_flag == "PASS", "✓", "⚠️ "),
    stringsAsFactors = FALSE
  ))
  
  cat("  ✓", gsm_id, "—", nrow(seurat_obj), "genes,", ncol(seurat_obj), "cells [", qc_flag, "]\n")
}

cat("\n✓ Loaded", length(seurat_list), "samples successfully\n\n")

# ============================================================================
# Step 4: Merge all samples
# ============================================================================

cat("Merging all samples...\n")

if (length(seurat_list) < 2) {
  stop("ERROR: Need at least 2 samples to merge")
}

seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  project = "HBV_Liver"
)

cat("✓ Merged object dimensions:\n")
cat("  Cells:", ncol(seurat_merged), "\n")
cat("  Genes:", nrow(seurat_merged), "\n\n")

# ============================================================================
# Step 5: Verify metadata
# ============================================================================

cat("Verifying metadata...\n")
cat("\nPhase distribution:\n")
print(table(seurat_merged$Phase))

cat("\nQC flag distribution:\n")
print(table(seurat_merged$QC_Flag))

cat("\nCells per phase:\n")
phase_summary <- seurat_merged@meta.data %>%
  group_by(Phase) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(Donor),
    .groups = 'drop'
  )
print(phase_summary)

cat("\n✓ Metadata verified\n\n")

# ============================================================================
# Step 6: Save checkpoint
# ============================================================================

cat("Saving checkpoint...\n")

dir.create("results/rds_objects", recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_merged, "results/rds_objects/seurat_merged_raw.rds")

cat("✓ Saved to results/rds_objects/seurat_merged_raw.rds\n\n")

# ============================================================================
# Step 7: Summary report
# ============================================================================

cat("=== LOAD SUMMARY ===\n")
print(load_summary)

cat("\n=== METADATA SNAPSHOT ===\n")
print(head(seurat_merged@meta.data))

cat("\n=== PHASE 2 SCRIPT 01 COMPLETE ===\n")
cat("Next: Run Script 02 (QC filtering)\n")
