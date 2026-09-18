# ============================================================================
# PHASE 3 — GLOBAL CELL ATLAS CONSTRUCTION
# ============================================================================
#
# Script: phase_03_01_atlas_construction.R
#
# Purpose:
# Construct a representative global reference atlas from the complete
# Phase 2 QC-retained dataset.
#
# Strategy:
# 1. Load the Phase 2 processed-expression object.
# 2. Verify processed log-CP10K expression and metadata.
# 3. Do NOT re-normalize the supplied expression.
# 4. Identify 2,000 consensus highly variable genes across the 23 samples.
# 5. Construct a 20,000-cell global reference sketch:
#    - all 23 samples/donors represented;
#    - minimum 20-cell safeguard per sample;
#    - remaining cells allocated approximately proportionally;
#    - within-sample leverage-score-weighted sampling without replacement.
# 6. Validate sample, donor, and clinical-state representation.
# 7. Scale the 20,000-cell reference using the 2,000 HVGs.
# 8. Perform PCA, neighbour graph construction, clustering, and UMAP.
# 9. Generate global atlas diagnostics by cluster, clinical state, and donor.
#
# Biological framing:
# Clinical states are treated as cohort-level comparisons rather than a
# presumed linear disease trajectory.
#
# Important:
# The complete Phase 2 dataset is retained unchanged.
# The 20,000-cell atlas is a reference sketch used for exploratory
# dimensionality reduction, clustering, and global annotation.
#
# Input:
# results/rds_objects/seurat_qc_processed_expression.rds
#
# Output:
# results/rds_objects/phase3_atlas_sketch_umap_clustered.rds
#
# ============================================================================


# ============================================================================
# SECTION 1 — ENVIRONMENT SETUP
# ============================================================================

cat("============================================================\n")
cat("PHASE 3 — GLOBAL CELL ATLAS CONSTRUCTION\n")
cat("============================================================\n\n")

cat("=== SECTION 1: ENVIRONMENT SETUP ===\n\n")

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
library(Matrix)

setwd(here())

rm(list = ls())
gc()

set.seed(12345)

dir.create(
  "results/rds_objects",
  recursive = TRUE,
  showWarnings = FALSE
)

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

cat("✓ Environment initialized\n")
cat("✓ Random seed set to 12345\n")
cat("✓ Output directories verified\n\n")


# ============================================================================
# SECTION 2 — LOAD PHASE 2 OBJECT + DATASET VALIDATION
# ============================================================================

cat("=== SECTION 2: LOADING PHASE 2 QC OBJECT ===\n\n")

seurat_obj <- readRDS(
  "results/rds_objects/seurat_qc_processed_expression.rds"
)

DefaultAssay(seurat_obj) <- "RNA"

cat(
  "Cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Genes:",
  nrow(seurat_obj),
  "\n"
)

cat(
  "Assays:",
  paste(
    Assays(seurat_obj),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Samples:",
  dplyr::n_distinct(seurat_obj$GSM),
  "\n"
)

cat(
  "Donors:",
  dplyr::n_distinct(seurat_obj$Donor),
  "\n"
)

cat(
  "Clinical states:",
  dplyr::n_distinct(seurat_obj$Phase),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 2.1 — Validate expected dataset dimensions
# ----------------------------------------------------------------------------

if (ncol(seurat_obj) != 106592) {
  
  stop(
    "ERROR: Expected 106,592 cells but found ",
    ncol(seurat_obj),
    "."
  )
}

if (nrow(seurat_obj) != 24452) {
  
  stop(
    "ERROR: Expected 24,452 genes but found ",
    nrow(seurat_obj),
    "."
  )
}

if (dplyr::n_distinct(seurat_obj$GSM) != 23) {
  
  stop(
    "ERROR: Expected 23 samples."
  )
}

if (dplyr::n_distinct(seurat_obj$Donor) != 23) {
  
  stop(
    "ERROR: Expected 23 donors."
  )
}

if (dplyr::n_distinct(seurat_obj$Phase) != 5) {
  
  stop(
    "ERROR: Expected 5 clinical states."
  )
}

cat("✓ Dataset dimensions match Phase 2 checkpoint\n\n")


# ----------------------------------------------------------------------------
# 2.2 — Required metadata
# ----------------------------------------------------------------------------

required_metadata <- c(
  "GSM",
  "Donor",
  "Phase"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(seurat_obj@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "ERROR: Required metadata columns are missing: ",
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


# ----------------------------------------------------------------------------
# 2.3 — Confirm sample-specific processed data layers
# ----------------------------------------------------------------------------

rna_layers <- Layers(
  seurat_obj[["RNA"]]
)

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

cat(
  "Sample-specific processed data layers:",
  length(data_layers),
  "\n"
)

cat(
  "Raw counts layers:",
  length(counts_layers),
  "\n\n"
)

if (length(data_layers) != 23) {
  
  stop(
    "ERROR: Expected 23 sample-specific processed data layers, found ",
    length(data_layers),
    "."
  )
}

if (length(counts_layers) > 0) {
  
  stop(
    "ERROR: Unexpected raw counts layers detected. ",
    "The Phase 2 object is expected to contain processed expression only."
  )
}

cat("✓ Processed-expression layer structure confirmed\n")
cat("✓ No raw UMI counts layers present\n\n")


# ----------------------------------------------------------------------------
# 2.4 — Confirm layer-to-sample correspondence
# ----------------------------------------------------------------------------

expected_data_layers <- paste0(
  "data.",
  sort(
    unique(
      seurat_obj$GSM
    )
  )
)

observed_data_layers <- sort(
  data_layers
)

if (!identical(
  observed_data_layers,
  expected_data_layers
)) {
  
  stop(
    "ERROR: Processed data layers do not correspond exactly to the 23 GSMs."
  )
}

cat("✓ All 23 GSM-specific data layers correspond to metadata samples\n\n")


# ============================================================================
# SECTION 3 — CONFIRM FULL DATASET IS RETAINED
# ============================================================================

cat("=== SECTION 3: FULL DATASET RETENTION VALIDATION ===\n\n")

sample_snapshot <- seurat_obj@meta.data %>%
  count(
    GSM,
    Donor,
    Phase,
    name = "Cells"
  ) %>%
  arrange(
    Phase,
    GSM
  )

cat("Sample-level cell counts:\n\n")

print(
  sample_snapshot
)

cat("\n")

write_csv(
  sample_snapshot,
  "results/tables/phase3_dataset_snapshot.csv"
)

cat("✓ Dataset snapshot saved\n\n")


# ----------------------------------------------------------------------------
# 3.1 — Clinical-state counts
# ----------------------------------------------------------------------------

phase_snapshot <- seurat_obj@meta.data %>%
  count(
    Phase,
    name = "Cells"
  ) %>%
  arrange(
    Phase
  )

cat("Clinical-state cell counts:\n\n")

print(
  phase_snapshot
)

cat("\n")

write_csv(
  phase_snapshot,
  "results/tables/phase3_clinical_state_snapshot.csv"
)

cat("✓ Clinical-state snapshot saved\n\n")


# ============================================================================
# SECTION 4 — USE SUPPLIED PROCESSED EXPRESSION
# ============================================================================
#
# The GEO matrices supplied to Phase 2 are already processed log-CP10K
# expression values.
#
# Therefore:
#
# - no additional LogNormalize step is performed;
# - no raw UMI counts are reconstructed;
# - the supplied expression representation is retained.
#
# HVG selection below is performed directly from these supplied expression
# values, independently within each sample layer.
#
# ============================================================================

cat("=== SECTION 4: PROCESSED EXPRESSION VALIDATION ===\n\n")

cat(
  "Expression representation:",
  "processed log-CP10K\n"
)

cat(
  "Additional normalization:",
  "NOT PERFORMED\n\n"
)

cat("✓ Supplied processed expression will be used directly\n")
cat("✓ No double-normalization performed\n\n")


# ============================================================================
# SECTION 5 — CONSENSUS HIGHLY VARIABLE GENES
# ============================================================================

cat("=== SECTION 5: HIGHLY VARIABLE GENES ===\n\n")

rna <- seurat_obj[["RNA"]]

data_layers <- grep(
  "^data(\\.|$)",
  Layers(rna),
  value = TRUE
)

if (length(data_layers) != 23) {
  
  stop(
    "ERROR: Expected 23 data layers, found ",
    length(data_layers),
    "."
  )
}

cat(
  "Calculating HVGs independently across ",
  length(data_layers),
  " sample layers...\n\n"
)


# ----------------------------------------------------------------------------
# 5.1 — Calculate per-sample HVGs
# ----------------------------------------------------------------------------
#
# The source matrices contain processed log-CP10K values and no raw counts.
#
# A temporary legacy Assay is therefore constructed using ONLY the existing
# processed expression in its data slot.
#
# This does NOT:
# - create fake counts;
# - perform normalization;
# - alter the Phase 2 object.
#
# It simply provides FindVariableFeatures with a single data matrix rather
# than a multi-layer Assay5 whose default layer is ambiguous.
# ----------------------------------------------------------------------------

layer_hvgs <- vector(
  mode = "list",
  length = length(data_layers)
)

names(layer_hvgs) <- data_layers

for (layer_name in data_layers) {
  
  cat(
    "Processing:",
    layer_name,
    "\n"
  )
  
  layer_matrix <- LayerData(
    object = rna,
    layer = layer_name
  )
  
  if (ncol(layer_matrix) == 0) {
    
    stop(
      "ERROR: Layer ",
      layer_name,
      " contains zero cells."
    )
  }
  
  if (nrow(layer_matrix) != nrow(seurat_obj)) {
    
    stop(
      "ERROR: Layer ",
      layer_name,
      " does not contain the expected 24,452 genes."
    )
  }
  
  # Create a temporary legacy Assay with the existing processed
  # expression stored directly as the data layer.
  #
  # IMPORTANT:
  # This is NOT a counts matrix.
  temporary_assay <- CreateAssayObject(
    data = layer_matrix
  )
  
  temporary_assay <- FindVariableFeatures(
    object = temporary_assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  current_hvgs <- VariableFeatures(
    temporary_assay
  )
  
  if (length(current_hvgs) != 2000) {
    
    stop(
      "ERROR: ",
      layer_name,
      " returned ",
      length(current_hvgs),
      " HVGs instead of 2,000."
    )
  }
  
  layer_hvgs[[layer_name]] <- current_hvgs
  
  rm(
    layer_matrix,
    temporary_assay,
    current_hvgs
  )
  
  gc()
}

cat("\n✓ Per-sample HVG calculation complete\n\n")


# ----------------------------------------------------------------------------
# 5.2 — Build per-sample HVG ranking table
# ----------------------------------------------------------------------------

hvg_rank_table <- purrr::map2_dfr(
  
  layer_hvgs,
  
  names(layer_hvgs),
  
  ~ tibble(
    
    Gene = .x,
    
    Sample = .y,
    
    Rank = seq_along(.x)
    
  )
)

if (nrow(hvg_rank_table) != 23 * 2000) {
  
  stop(
    "ERROR: Expected 46,000 per-sample HVG records."
  )
}

cat(
  "Per-sample HVG records:",
  nrow(hvg_rank_table),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 5.3 — Frequency of HVG selection across samples
# ----------------------------------------------------------------------------

hvg_frequency <- hvg_rank_table %>%
  count(
    Gene,
    name = "Frequency"
  )


# ----------------------------------------------------------------------------
# 5.4 — Median HVG rank across samples
# ----------------------------------------------------------------------------

hvg_median_rank <- hvg_rank_table %>%
  group_by(
    Gene
  ) %>%
  summarise(
    Median_Rank = median(Rank),
    .groups = "drop"
  )


# ----------------------------------------------------------------------------
# 5.5 — Consensus ranking
# ----------------------------------------------------------------------------
#
# Ranking principle:
#
# 1. genes selected as HVGs in more samples are prioritized;
# 2. median HVG rank breaks ties;
# 3. gene name provides deterministic final tie-breaking.
#
# This is the same general frequency + median-rank logic documented by
# Seurat for consensus integration-feature selection.
# ----------------------------------------------------------------------------

hvg_consensus <- hvg_frequency %>%
  left_join(
    hvg_median_rank,
    by = "Gene"
  ) %>%
  arrange(
    desc(Frequency),
    Median_Rank,
    Gene
  )

if (nrow(hvg_consensus) < 2000) {
  
  stop(
    "ERROR: Fewer than 2,000 consensus HVG candidates were available."
  )
}


# ----------------------------------------------------------------------------
# 5.6 — Select final 2,000 HVGs
# ----------------------------------------------------------------------------

hvg <- hvg_consensus %>%
  slice_head(
    n = 2000
  ) %>%
  pull(
    Gene
  )

cat(
  "Number of consensus highly variable genes:",
  length(hvg),
  "\n"
)

if (length(hvg) != 2000) {
  
  stop(
    "ERROR: Expected exactly 2,000 consensus HVGs."
  )
}

VariableFeatures(seurat_obj) <- hvg

cat("✓ Exactly 2,000 consensus HVGs assigned\n\n")


# ----------------------------------------------------------------------------
# 5.7 — HVG frequency diagnostics
# ----------------------------------------------------------------------------

cat(
  "HVGs selected in all 23 samples:",
  sum(hvg_consensus$Frequency >= 23),
  "\n"
)

cat(
  "HVGs selected in >= 20 samples:",
  sum(hvg_consensus$Frequency >= 20),
  "\n"
)

cat(
  "HVGs selected in >= 15 samples:",
  sum(hvg_consensus$Frequency >= 15),
  "\n\n"
)

write_csv(
  hvg_consensus,
  "results/tables/phase3_hvg_consensus_ranking.csv"
)

write_csv(
  tibble(
    HVG = hvg,
    Frequency =
      hvg_consensus$Frequency[
        match(
          hvg,
          hvg_consensus$Gene
        )
      ],
    Median_Rank =
      hvg_consensus$Median_Rank[
        match(
          hvg,
          hvg_consensus$Gene
        )
      ]
  ),
  "results/tables/phase3_final_2000_hvgs.csv"
)

cat("✓ HVG diagnostics saved\n\n")


# ============================================================================
# SECTION 6 — GLOBAL 20,000-CELL LEVERAGE-SCORE SKETCH
# ============================================================================

cat("=== SECTION 6: GLOBAL 20,000-CELL SKETCH ===\n\n")

target_sketch_cells <- 20000L

minimum_cells_per_sample <- 20L

total_cells <- ncol(
  seurat_obj
)

cat(
  "Full dataset:",
  total_cells,
  "cells\n"
)

cat(
  "Target sketch:",
  target_sketch_cells,
  "cells\n"
)

cat(
  "Target fraction:",
  round(
    100 *
      target_sketch_cells /
      total_cells,
    2
  ),
  "%\n"
)

cat(
  "Minimum sample representation:",
  minimum_cells_per_sample,
  "cells\n\n"
)


# ----------------------------------------------------------------------------
# 6.1 — Sample sizes
# ----------------------------------------------------------------------------

sample_sizes <- seurat_obj@meta.data %>%
  count(
    GSM,
    Donor,
    Phase,
    name = "Original_Cells"
  ) %>%
  arrange(
    GSM
  )

if (nrow(sample_sizes) != 23) {
  
  stop(
    "ERROR: Expected 23 sample/donor records."
  )
}


# ----------------------------------------------------------------------------
# 6.2 — Validate minimum representation feasibility
# ----------------------------------------------------------------------------

minimum_required_cells <-
  nrow(sample_sizes) *
  minimum_cells_per_sample

if (minimum_required_cells >= target_sketch_cells) {
  
  stop(
    "ERROR: Minimum representation requirement exceeds sketch capacity."
  )
}


# ----------------------------------------------------------------------------
# 6.3 — Allocate minimum representation
# ----------------------------------------------------------------------------

sample_sizes <- sample_sizes %>%
  mutate(
    Minimum_Allocation =
      minimum_cells_per_sample
  )

remaining_cells <-
  target_sketch_cells -
  sum(
    sample_sizes$Minimum_Allocation
  )

cat(
  "Cells reserved for minimum representation:",
  sum(sample_sizes$Minimum_Allocation),
  "\n"
)

cat(
  "Cells remaining for proportional allocation:",
  remaining_cells,
  "\n\n"
)


# ----------------------------------------------------------------------------
# 6.4 — Allocate remaining cells approximately proportionally
# ----------------------------------------------------------------------------
#
# Largest-remainder allocation is used so the final allocation is exactly
# 20,000 cells.
# ----------------------------------------------------------------------------

sample_sizes <- sample_sizes %>%
  mutate(
    Proportional_Exact =
      Original_Cells /
      sum(Original_Cells) *
      remaining_cells,
    
    Proportional_Base =
      floor(
        Proportional_Exact
      ),
    
    Proportional_Remainder =
      Proportional_Exact -
      Proportional_Base
  )

cells_left_after_floor <-
  remaining_cells -
  sum(
    sample_sizes$Proportional_Base
  )

if (cells_left_after_floor > 0) {
  
  extra_indices <- order(
    sample_sizes$Proportional_Remainder,
    decreasing = TRUE
  )[seq_len(
    cells_left_after_floor
  )]
  
  sample_sizes$Proportional_Base[
    extra_indices
  ] <-
    sample_sizes$Proportional_Base[
      extra_indices
    ] + 1L
}


# ----------------------------------------------------------------------------
# 6.5 — Final allocation
# ----------------------------------------------------------------------------

sample_sizes <- sample_sizes %>%
  mutate(
    Sketch_Cells =
      Minimum_Allocation +
      Proportional_Base,
    
    Sketch_Fraction =
      Sketch_Cells /
      Original_Cells,
    
    Sketch_Percent =
      100 *
      Sketch_Fraction
  )

cat("Final sample allocation:\n\n")

print(
  sample_sizes
)

cat("\n")

if (
  sum(sample_sizes$Sketch_Cells) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Final sample allocations do not sum to exactly 20,000."
  )
}

if (
  any(
    sample_sizes$Sketch_Cells >
    sample_sizes$Original_Cells
  )
) {
  
  stop(
    "ERROR: A sample was allocated more sketch cells than available cells."
  )
}

cat("✓ Exactly 20,000 cells allocated across 23 samples\n")

cat(
  "✓ Every sample receives at least",
  minimum_cells_per_sample,
  "cells\n\n"
)


# ----------------------------------------------------------------------------
# 6.6 — Save allocation table
# ----------------------------------------------------------------------------

write_csv(
  sample_sizes,
  "results/tables/phase3_sketch_allocation_by_sample.csv"
)

cat("✓ Sketch allocation table saved\n\n")

# ============================================================================
# SECTION 6.7 — CALCULATE SAMPLE-SPECIFIC LEVERAGE SCORES
# ============================================================================
#
# The leverage score is calculated from a low-dimensional PCA representation
# of the 2,000 consensus HVGs within each sample.
#
# This avoids asking LeverageScore() to perform its internal decomposition
# directly on a 2,000-feature × many-cell matrix, which can become
# numerically unstable for this data representation.
#
# The resulting leverage scores are then used for within-sample
# probability-weighted sampling without replacement.
#
# ============================================================================

cat(
  "=== CALCULATING SAMPLE-SPECIFIC LEVERAGE SCORES ===\n\n"
)

leverage_scores <- vector(
  mode = "list",
  length = nrow(sample_sizes)
)

names(leverage_scores) <-
  sample_sizes$GSM


for (
  i in seq_len(
    nrow(sample_sizes)
  )
) {
  
  GSM_id <- sample_sizes$GSM[i]
  
  layer_name <- paste0(
    "data.",
    GSM_id
  )
  
  cat(
    "Calculating leverage scores:",
    GSM_id,
    "\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Extract consensus HVGs from the sample-specific processed-expression layer
  # --------------------------------------------------------------------------
  
  layer_matrix <- LayerData(
    object = seurat_obj[["RNA"]],
    layer = layer_name,
    features = hvg
  )
  
  if (
    ncol(layer_matrix) == 0
  ) {
    
    stop(
      "ERROR: Layer ",
      layer_name,
      " contains zero cells."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Transpose to cell × gene orientation
  # --------------------------------------------------------------------------
  
  cell_by_gene <- t(
    as.matrix(layer_matrix)
  )
  
  
  n_cells <- nrow(
    cell_by_gene
  )
  
  n_features <- ncol(
    cell_by_gene
  )
  
  
  # --------------------------------------------------------------------------
  # Determine a valid PCA dimensionality.
  #
  # We use up to 50 dimensions, but never more than the mathematical rank
  # permitted by the number of cells.
  #
  # The very small 66-cell sample therefore uses at most 50 dimensions,
  # while larger samples also use 50.
  # --------------------------------------------------------------------------
  
  n_pca_dims <- min(
    50L,
    n_cells - 1L,
    n_features
  )
  
  if (
    n_pca_dims < 2L
  ) {
    
    stop(
      "ERROR: ",
      GSM_id,
      " has insufficient cells/features for leverage calculation."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # PCA of the sample-specific expression matrix
  #
  # prcomp() is used here because these are already processed expression
  # values and we need a stable numerical representation for the leverage
  # calculation.
  # --------------------------------------------------------------------------
  
  pca_result <- tryCatch(
    
    {
      
      prcomp(
        x =
          cell_by_gene,
        
        center =
          TRUE,
        
        scale. =
          FALSE,
        
        rank. =
          n_pca_dims
      )
      
    },
    
    error = function(e) {
      
      stop(
        "ERROR: PCA failed for ",
        GSM_id,
        ".\n",
        "Original error: ",
        conditionMessage(e)
      )
      
    }
  )
  
  
  # --------------------------------------------------------------------------
  # Extract cell embeddings
  # --------------------------------------------------------------------------
  
  cell_embeddings <- pca_result$x
  
  if (
    ncol(cell_embeddings) < 2
  ) {
    
    stop(
      "ERROR: Fewer than 2 PCA dimensions were obtained for ",
      GSM_id,
      "."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Calculate approximate leverage scores.
  #
  # For an orthogonal PCA score matrix U, leverage is proportional to the
  # squared row norm of the cell coordinates after normalization by the
  # singular-value scale.
  #
  # We calculate leverage from the normalized PCA left-singular vectors.
  # --------------------------------------------------------------------------
  
  singular_values <- pca_result$sdev
  
  usable_dims <- min(
    ncol(cell_embeddings),
    sum(
      singular_values >
        sqrt(.Machine$double.eps)
    )
  )
  
  if (
    usable_dims < 2
  ) {
    
    stop(
      "ERROR: Insufficient numerical rank for leverage calculation in ",
      GSM_id,
      "."
    )
  }
  
  normalized_embeddings <-
    sweep(
      cell_embeddings[
        ,
        seq_len(usable_dims),
        drop = FALSE
      ],
      
      MARGIN = 2,
      
      STATS = singular_values[
        seq_len(usable_dims)
      ],
      
      FUN = "/"
    )
  
  
  # Row-wise squared norm
  score_vector <- rowSums(
    normalized_embeddings^2
  )
  
  
  # --------------------------------------------------------------------------
  # Add a small positive floor so every cell remains eligible for sampling.
  # --------------------------------------------------------------------------
  
  score_vector <- pmax(
    score_vector,
    .Machine$double.eps
  )
  
  names(score_vector) <-
    rownames(cell_by_gene)
  
  
  # --------------------------------------------------------------------------
  # Validate scores
  # --------------------------------------------------------------------------
  
  if (
    length(score_vector) !=
    n_cells
  ) {
    
    stop(
      "ERROR: Leverage-score vector length does not match cell count for ",
      GSM_id,
      "."
    )
  }
  
  if (
    any(
      !is.finite(score_vector)
    )
  ) {
    
    stop(
      "ERROR: Non-finite leverage scores detected in ",
      GSM_id,
      "."
    )
  }
  
  if (
    any(
      score_vector < 0
    )
  ) {
    
    stop(
      "ERROR: Negative leverage scores detected in ",
      GSM_id,
      "."
    )
  }
  
  if (
    sum(score_vector) <= 0
  ) {
    
    stop(
      "ERROR: Leverage scores sum to zero in ",
      GSM_id,
      "."
    )
  }
  
  
  # Store scores
  leverage_scores[[GSM_id]] <-
    score_vector
  
  
  cat(
    "  Cells:",
    n_cells,
    "\n"
  )
  
  cat(
    "  PCA dimensions:",
    usable_dims,
    "\n"
  )
  
  cat(
    "  Score range:",
    format(
      min(score_vector),
      scientific = TRUE
    ),
    "to",
    format(
      max(score_vector),
      scientific = TRUE
    ),
    "\n\n"
  )
  
  
  rm(
    layer_matrix,
    cell_by_gene,
    pca_result,
    cell_embeddings,
    singular_values,
    normalized_embeddings,
    score_vector
  )
  
  gc()
}


cat(
  "✓ Leverage scores calculated for all 23 samples\n\n"
)

# ============================================================================
# SECTION 6.8 — LEVERAGE-SCORE-WEIGHTED SAMPLING
# ============================================================================

cat(
  "=== SELECTING LEVERAGE-WEIGHTED REPRESENTATIVE CELLS ===\n\n"
)

selected_cells_by_sample <- vector(
  mode = "list",
  length = nrow(sample_sizes)
)

names(selected_cells_by_sample) <-
  sample_sizes$GSM

for (
  i in seq_len(
    nrow(sample_sizes)
  )
) {
  
  GSM_id <- sample_sizes$GSM[i]
  
  n_to_select <-
    sample_sizes$Sketch_Cells[i]
  
  scores <-
    leverage_scores[[GSM_id]]
  
  if (is.null(scores)) {
    
    stop(
      "ERROR: No leverage scores found for ",
      GSM_id,
      "."
    )
  }
  
  if (
    n_to_select >
    length(scores)
  ) {
    
    stop(
      "ERROR: Requested ",
      n_to_select,
      " cells from ",
      GSM_id,
      " but only ",
      length(scores),
      " are available."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Convert leverage scores to sampling probabilities.
  # --------------------------------------------------------------------------
  
  sampling_probabilities <-
    scores /
    sum(scores)
  
  if (
    any(
      !is.finite(
        sampling_probabilities
      )
    )
  ) {
    
    stop(
      "ERROR: Invalid sampling probabilities for ",
      GSM_id,
      "."
    )
  }
  
  set.seed(
    12345 + i
  )
  
  selected_cells <- sample(
    x = names(scores),
    size = n_to_select,
    replace = FALSE,
    prob = sampling_probabilities
  )
  
  selected_cells_by_sample[[GSM_id]] <-
    selected_cells
  
  cat(
    GSM_id,
    ": selected",
    length(selected_cells),
    "of",
    length(scores),
    "cells\n"
  )
}


# ============================================================================
# SECTION 6.9 — COMBINE AND VALIDATE SELECTED CELLS
# ============================================================================

sketch_cells <- unlist(
  selected_cells_by_sample,
  use.names = FALSE
)

if (
  anyDuplicated(sketch_cells) > 0
) {
  
  stop(
    "ERROR: Duplicate cell IDs were selected."
  )
}

cat(
  "\nTotal selected cells:",
  length(sketch_cells),
  "\n"
)

if (
  length(sketch_cells) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Expected exactly ",
    target_sketch_cells,
    " cells but selected ",
    length(sketch_cells),
    "."
  )
}

cat("✓ Exactly 20,000 unique cells selected\n\n")


# ============================================================================
# SECTION 6.10 — CREATE THE ATLAS SKETCH
# ============================================================================

cat(
  "Creating 20,000-cell atlas sketch...\n\n"
)

atlas_sketch <- subset(
  x = seurat_obj,
  cells = sketch_cells
)

if (
  ncol(atlas_sketch) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Atlas sketch does not contain exactly 20,000 cells."
  )
}

if (
  nrow(atlas_sketch) !=
  nrow(seurat_obj)
) {
  
  stop(
    "ERROR: Atlas sketch does not retain the full gene set."
  )
}

cat(
  "✓ Atlas sketch dimensions:",
  ncol(atlas_sketch),
  "cells ×",
  nrow(atlas_sketch),
  "genes\n\n"
)


# ============================================================================
# SECTION 7 — SKETCH REPRESENTATION AUDIT
# ============================================================================

cat("=== SECTION 7: SKETCH REPRESENTATION AUDIT ===\n\n")


# ----------------------------------------------------------------------------
# 7.1 — Sample representation
# ----------------------------------------------------------------------------

sketch_sample_counts <- atlas_sketch@meta.data %>%
  count(
    GSM,
    Donor,
    Phase,
    name = "Sketch_Cells"
  )

sketch_representation <- sample_sizes %>%
  select(
    GSM,
    Donor,
    Phase,
    Original_Cells
  ) %>%
  left_join(
    sketch_sample_counts,
    by = c(
      "GSM",
      "Donor",
      "Phase"
    )
  ) %>%
  mutate(
    Sketch_Cells =
      replace_na(
        Sketch_Cells,
        0L
      ),
    
    Sketch_Fraction =
      Sketch_Cells /
      Original_Cells,
    
    Sketch_Percent =
      100 *
      Sketch_Fraction
  ) %>%
  arrange(
    Phase,
    GSM
  )

cat("Sample representation:\n\n")

print(
  sketch_representation
)

cat("\n")

if (
  any(
    sketch_representation$Sketch_Cells < 1
  )
) {
  
  stop(
    "ERROR: One or more samples are absent from the atlas sketch."
  )
}

if (
  nrow(sketch_representation) != 23
) {
  
  stop(
    "ERROR: Expected 23 samples in the representation table."
  )
}

if (
  sum(sketch_representation$Sketch_Cells) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Sample representation does not sum to 20,000 cells."
  )
}

cat("✓ All 23 samples represented\n")
cat("✓ Sample representation totals exactly 20,000 cells\n\n")


# ----------------------------------------------------------------------------
# 7.2 — Clinical-state representation
# ----------------------------------------------------------------------------

original_phase_counts <- seurat_obj@meta.data %>%
  count(
    Phase,
    name = "Original_Cells"
  )

sketch_phase_counts <- atlas_sketch@meta.data %>%
  count(
    Phase,
    name = "Sketch_Cells"
  )

phase_representation <- original_phase_counts %>%
  left_join(
    sketch_phase_counts,
    by = "Phase"
  ) %>%
  mutate(
    Sketch_Fraction =
      Sketch_Cells /
      Original_Cells,
    
    Sketch_Percent =
      100 *
      Sketch_Fraction,
    
    Original_Fraction =
      Original_Cells /
      sum(Original_Cells),
    
    Sketch_Cohort_Fraction =
      Sketch_Cells /
      sum(Sketch_Cells)
  ) %>%
  arrange(
    Phase
  )

cat("Clinical-state representation:\n\n")

print(
  phase_representation
)

cat("\n")

write_csv(
  phase_representation,
  "results/tables/phase3_sketch_representation_by_phase.csv"
)

write_csv(
  sketch_representation,
  "results/tables/phase3_sketch_representation_by_sample.csv"
)

cat("✓ Clinical-state representation saved\n")
cat("✓ Sample representation saved\n\n")


# ----------------------------------------------------------------------------
# 7.3 — Donor representation
# ----------------------------------------------------------------------------

donor_representation <- atlas_sketch@meta.data %>%
  count(
    Donor,
    GSM,
    Phase,
    name = "Sketch_Cells"
  ) %>%
  arrange(
    Phase,
    GSM
  )

if (
  nrow(donor_representation) != 23
) {
  
  stop(
    "ERROR: Expected 23 donors in the atlas sketch."
  )
}

cat("✓ All 23 donors represented\n\n")

write_csv(
  donor_representation,
  "results/tables/phase3_sketch_representation_by_donor.csv"
)


# ----------------------------------------------------------------------------
# 7.4 — Representation plot
# ----------------------------------------------------------------------------

representation_plot_data <-
  sketch_representation %>%
  select(
    GSM,
    Original_Cells,
    Sketch_Cells
  ) %>%
  pivot_longer(
    cols = c(
      Original_Cells,
      Sketch_Cells
    ),
    names_to = "Dataset",
    values_to = "Cells"
  )

p_sketch_representation <-
  ggplot(
    representation_plot_data,
    aes(
      x = GSM,
      y = Cells,
      fill = Dataset
    )
  ) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title =
      "Original Dataset and 20,000-Cell Atlas Sketch",
    x =
      "Sample",
    y =
      "Number of Cells",
    fill =
      NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1
      ),
    
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_sketch_representation_by_sample.png",
  p_sketch_representation,
  width = 14,
  height = 7,
  dpi = 300
)

cat("✓ Sketch representation plot saved\n\n")


# ============================================================================
# SECTION 8 — SCALE THE 20,000-CELL ATLAS SKETCH
# ============================================================================

cat("=== SECTION 8: SCALING THE ATLAS SKETCH ===\n\n")

VariableFeatures(atlas_sketch) <-
  hvg

sketch_hvg <-
  VariableFeatures(atlas_sketch)

if (
  length(sketch_hvg) != 2000
) {
  
  stop(
    "ERROR: Atlas sketch does not contain exactly 2,000 HVGs."
  )
}

cat(
  "Scaling:",
  length(sketch_hvg),
  "HVGs\n"
)

cat(
  "Cells:",
  ncol(atlas_sketch),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 8.1 — Scale sketch
# ----------------------------------------------------------------------------

atlas_sketch <- ScaleData(
  object = atlas_sketch,
  assay = "RNA",
  features = sketch_hvg,
  verbose = TRUE
)

cat("\n✓ Scaling completed\n\n")


# ----------------------------------------------------------------------------
# 8.2 — Validate scaled data
# ----------------------------------------------------------------------------

scaled_layers <- Layers(
  atlas_sketch[["RNA"]],
  search = "^scale\\.data"
)

cat(
  "Scaled-data layers detected:",
  length(scaled_layers),
  "\n"
)

if (
  length(scaled_layers) == 0
) {
  
  stop(
    "ERROR: No scaled-data layer was created."
  )
}

cat("✓ Scaled-data representation detected\n\n")


# ----------------------------------------------------------------------------
# 8.3 — Save scaled checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_scaled.rds"
)

cat(
  "✓ Scaled atlas checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_scaled.rds\n\n"
)


# ============================================================================
# SECTION 9 — PCA
# ============================================================================

cat("=== SECTION 9: PCA ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_scaled.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

cat(
  "PCA input cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "PCA features:",
  length(
    VariableFeatures(atlas_sketch)
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 9.1 — Run PCA
# ----------------------------------------------------------------------------

atlas_sketch <- RunPCA(
  object = atlas_sketch,
  assay = "RNA",
  features = VariableFeatures(atlas_sketch),
  npcs = 50,
  verbose = TRUE
)

cat("\n✓ PCA completed\n")

cat(
  "Principal components calculated:",
  ncol(
    Embeddings(
      atlas_sketch,
      reduction = "pca"
    )
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 9.2 — Validate PCA
# ----------------------------------------------------------------------------

pca_embeddings <- Embeddings(
  atlas_sketch,
  reduction = "pca"
)

if (
  nrow(pca_embeddings) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: PCA embedding does not contain all 20,000 cells."
  )
}

if (
  ncol(pca_embeddings) != 50
) {
  
  stop(
    "ERROR: Expected 50 PCs."
  )
}

if (
  any(
    !is.finite(pca_embeddings)
  )
) {
  
  stop(
    "ERROR: Non-finite values detected in PCA embeddings."
  )
}

cat("✓ PCA contains exactly 20,000 cells and 50 PCs\n\n")


# ----------------------------------------------------------------------------
# 9.3 — Elbow plot
# ----------------------------------------------------------------------------

p_elbow <- ElbowPlot(
  object = atlas_sketch,
  reduction = "pca",
  ndims = 50
) +
  labs(
    title =
      "PCA Elbow Plot — 20,000-Cell Global Atlas",
    x =
      "Principal Component",
    y =
      "Standard Deviation"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_pca_elbow_plot.png",
  p_elbow,
  width = 9,
  height = 6,
  dpi = 300
)

cat("✓ PCA elbow plot saved\n\n")


# ----------------------------------------------------------------------------
# 9.4 — Save PCA checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_pca50.rds"
)

cat("✓ PCA checkpoint saved\n\n")


# ============================================================================
# SECTION 10 — SELECT PCS FOR GLOBAL ATLAS
# ============================================================================
#
# A fixed 1:20 PC range is used as a predefined starting point.
#
# The elbow plot is retained as a diagnostic rather than being used to
# retrospectively tune the analysis after biological annotations are known.
#
# ============================================================================

cat("=== SECTION 10: SELECTING PCS ===\n\n")

selected_pcs <- 1:20

available_pcs <- ncol(
  Embeddings(
    atlas_sketch,
    reduction = "pca"
  )
)

if (
  available_pcs < 20
) {
  
  stop(
    "ERROR: Fewer than 20 PCs are available."
  )
}

cat(
  "Selected PCs:",
  paste(
    selected_pcs,
    collapse = ", "
  ),
  "\n\n"
)

selected_pc_table <- tibble(
  PC = selected_pcs
)

write_csv(
  selected_pc_table,
  "results/tables/phase3_selected_pcs.csv"
)

cat("✓ PC selection saved\n\n")


# ============================================================================
# SECTION 11 — NEIGHBOUR GRAPH + CLUSTERING
# ============================================================================

cat("=== SECTION 11: NEIGHBOUR GRAPH + CLUSTERING ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_pca50.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

cat(
  "Cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "PCs used:",
  paste(
    selected_pcs,
    collapse = ", "
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 11.1 — Construct neighbour graph
# ----------------------------------------------------------------------------

atlas_sketch <- FindNeighbors(
  object = atlas_sketch,
  reduction = "pca",
  dims = selected_pcs,
  verbose = TRUE
)

cat("\n✓ Neighbour graph constructed\n\n")


# ----------------------------------------------------------------------------
# 11.2 — Validate graph
# ----------------------------------------------------------------------------

graph_names <- names(
  atlas_sketch@graphs
)

cat("Graphs stored:\n\n")

print(
  graph_names
)

cat("\n")

if (
  length(graph_names) == 0
) {
  
  stop(
    "ERROR: No neighbour graph was created."
  )
}

cat("✓ Neighbour graph detected\n\n")


# ----------------------------------------------------------------------------
# 11.3 — Cluster
# ----------------------------------------------------------------------------

set.seed(12345)

atlas_sketch <- FindClusters(
  object = atlas_sketch,
  resolution = 0.5,
  verbose = TRUE
)

cat("\n✓ Clustering completed\n\n")


# ----------------------------------------------------------------------------
# 11.4 — Validate clusters
# ----------------------------------------------------------------------------

if (
  !"seurat_clusters" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "ERROR: seurat_clusters was not created."
  )
}

cluster_table <- table(
  atlas_sketch$seurat_clusters
)

if (
  any(
    is.na(
      atlas_sketch$seurat_clusters
    )
  )
) {
  
  stop(
    "ERROR: Missing cluster assignments detected."
  )
}

if (
  sum(cluster_table) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: Cluster assignments do not account for all cells."
  )
}

cat(
  "Number of clusters:",
  length(cluster_table),
  "\n\n"
)

print(
  cluster_table
)

cat("\n")

cluster_table_df <- tibble(
  Cluster =
    names(cluster_table),
  
  Cells =
    as.integer(cluster_table)
)

write_csv(
  cluster_table_df,
  "results/tables/phase3_global_cluster_sizes.csv"
)

cat("✓ Cluster-size table saved\n\n")


# ----------------------------------------------------------------------------
# 11.5 — Save clustered checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_clustered.rds"
)

cat("✓ Clustered atlas checkpoint saved\n\n")


# ============================================================================
# SECTION 12 — UMAP
# ============================================================================

cat("=== SECTION 12: UMAP ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_clustered.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

set.seed(12345)

atlas_sketch <- RunUMAP(
  object = atlas_sketch,
  reduction = "pca",
  dims = selected_pcs,
  verbose = TRUE
)

cat("\n✓ UMAP completed\n\n")


# ----------------------------------------------------------------------------
# 12.1 — Validate UMAP
# ----------------------------------------------------------------------------

umap_embeddings <- Embeddings(
  atlas_sketch,
  reduction = "umap"
)

if (
  nrow(umap_embeddings) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: UMAP embedding does not contain all atlas cells."
  )
}

if (
  ncol(umap_embeddings) != 2
) {
  
  stop(
    "ERROR: UMAP should contain two dimensions."
  )
}

if (
  any(
    !is.finite(umap_embeddings)
  )
) {
  
  stop(
    "ERROR: Non-finite UMAP values detected."
  )
}

cat("✓ UMAP contains 20,000 cells and 2 dimensions\n\n")


# ============================================================================
# SECTION 13 — UMAP BY CLUSTER
# ============================================================================

cat("=== SECTION 13: UMAP BY CLUSTER ===\n\n")

p_global_cluster <- DimPlot(
  object = atlas_sketch,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Global UMAP — Atlas Sketch by Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold"
      )
  )

ggsave(
  "results/figures/phase3_global_umap_by_cluster.png",
  p_global_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Cluster UMAP saved\n\n")


# ============================================================================
# SECTION 14 — UMAP BY CLINICAL STATE
# ============================================================================

cat("=== SECTION 14: UMAP BY CLINICAL STATE ===\n\n")

clinical_states <- unique(
  atlas_sketch$Phase
)

if (
  length(clinical_states) != 5
) {
  
  stop(
    "ERROR: Expected all 5 clinical states in atlas sketch."
  )
}

p_umap_phase <- DimPlot(
  object = atlas_sketch,
  reduction = "umap",
  group.by = "Phase",
  raster = TRUE
) +
  ggtitle(
    "Global UMAP — Atlas Sketch by Clinical State"
  ) +
  theme_classic() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold"
      )
  )

ggsave(
  "results/figures/phase3_global_umap_by_clinical_state.png",
  p_umap_phase,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Clinical-state UMAP saved\n\n")


# ============================================================================
# SECTION 15 — UMAP BY DONOR
# ============================================================================

cat("=== SECTION 15: UMAP BY DONOR ===\n\n")

donors <- unique(
  atlas_sketch$Donor
)

if (
  length(donors) != 23
) {
  
  stop(
    "ERROR: Expected 23 donors in atlas sketch."
  )
}

p_umap_donor <- DimPlot(
  object = atlas_sketch,
  reduction = "umap",
  group.by = "Donor",
  raster = TRUE
) +
  ggtitle(
    "Global UMAP — Atlas Sketch by Donor"
  ) +
  theme_classic() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold"
      ),
    
    legend.position =
      "right"
  )

ggsave(
  "results/figures/phase3_global_umap_by_donor.png",
  p_umap_donor,
  width = 12,
  height = 9,
  dpi = 300
)

cat("✓ Donor UMAP saved\n\n")


# ============================================================================
# SECTION 16 — FINAL ATLAS CHECKPOINT
# ============================================================================

cat("=== SECTION 16: FINAL ATLAS CHECKPOINT ===\n\n")

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_umap_clustered.rds"
)

cat(
  "✓ Final atlas checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_umap_clustered.rds\n\n"
)


# ============================================================================
# FINAL VALIDATION
# ============================================================================

cat(
  "============================================================\n"
)

cat(
  "✓ PHASE 3 SCRIPT 01 — COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat("Full dataset:\n")

cat(
  "  Cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "  Genes:",
  nrow(seurat_obj),
  "\n"
)

cat(
  "  Samples:",
  dplyr::n_distinct(seurat_obj$GSM),
  "\n"
)

cat(
  "  Donors:",
  dplyr::n_distinct(seurat_obj$Donor),
  "\n\n"
)

cat("Atlas sketch:\n")

cat(
  "  Cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "  Genes:",
  nrow(atlas_sketch),
  "\n"
)

cat(
  "  HVGs:",
  length(
    VariableFeatures(atlas_sketch)
  ),
  "\n"
)

cat(
  "  Samples:",
  dplyr::n_distinct(atlas_sketch$GSM),
  "\n"
)

cat(
  "  Donors:",
  dplyr::n_distinct(atlas_sketch$Donor),
  "\n"
)

cat(
  "  Clinical states:",
  dplyr::n_distinct(atlas_sketch$Phase),
  "\n"
)

cat(
  "  Clusters:",
  length(
    unique(
      atlas_sketch$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "  UMAP dimensions:",
  ncol(
    Embeddings(
      atlas_sketch,
      "umap"
    )
  ),
  "\n\n"
)

cat("Sketch design:\n")

cat(
  "  ✓ 20,000-cell target\n",
  "  ✓ All 23 samples represented\n",
  "  ✓ All 23 donors represented\n",
  "  ✓ Approximately proportional sample allocation\n",
  "  ✓ Minimum 20-cell representation safeguard\n",
  "  ✓ Within-sample leverage-score-weighted sampling\n",
  "  ✓ Sampling without replacement\n",
  "  ✓ Full dataset retained unchanged\n\n"
)

cat("Expression handling:\n")

cat(
  "  ✓ Supplied processed log-CP10K expression retained\n",
  "  ✓ No additional LogNormalize step\n",
  "  ✓ No raw UMI counts reconstructed\n\n"
)

cat("Atlas construction:\n")

cat(
  "  ✓ 2,000 consensus HVGs identified\n",
  "  ✓ Sketch scaled\n",
  "  ✓ 50 PCs calculated\n",
  "  ✓ PCs 1–20 used for global graph\n",
  "  ✓ Neighbour graph constructed\n",
  "  ✓ Clustering completed\n",
  "  ✓ UMAP completed\n",
  "  ✓ Cluster/state/donor diagnostics generated\n\n"
)

cat("Final checkpoint:\n")

cat(
  "  results/rds_objects/phase3_atlas_sketch_umap_clustered.rds\n\n"
)

cat(
  "Ready for phase_03_02_cluster_markers_annotation.R\n"
)

cat(
  "============================================================\n"
)
