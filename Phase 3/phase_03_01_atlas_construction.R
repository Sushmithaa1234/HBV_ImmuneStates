# ============================================================================
# PHASE 3 — GLOBAL CELL ATLAS CONSTRUCTION
# ============================================================================
#
# Script: phase_03_01_atlas_construction.R
#
# Sections:
#   1  — Environment setup
#   2  — Load QC-filtered singlet object + metadata checks
#   3  — Confirm only singlets are retained
#   4  — Log-normalization
#   5  — Identify top 2,000 highly variable genes using VST
#   6  — Global 20,000-cell leverage-score sketch
#   7  — Sketch representation by sample
#   8  — Scale the 20,000-cell sketch
#   9  — Separate the 20,000-cell atlas sketch
#   10 — PCA on the atlas sketch
#   11 — Select PCs 1–20
#   12 — Neighbour graph + clustering
#   13 — UMAP
#   14 — Global UMAP by cluster
#   15 — UMAP by clinical state
#   16 — UMAP by donor
#
# Biological framing:
#   Clinical states are treated as cohort-level comparisons rather than a
#   presumed linear disease trajectory.
#
# Important:
#   The full QC-filtered dataset is retained unchanged.
#   The 20,000-cell atlas is a representative sketch used for exploratory
#   dimensionality reduction, clustering, and global annotation.
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
# SECTION 2 — LOAD QC-FILTERED SINGLET OBJECT + METADATA CHECKS
# ============================================================================

cat("=== SECTION 2: LOADING QC-FILTERED SINGLET OBJECT ===\n\n")

seurat_obj <- readRDS(
  "results/rds_objects/seurat_qc_singlets.rds"
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
# 2.1 — Required metadata
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

cat("✓ Required metadata columns present:\n")
cat(
  paste(
    required_metadata,
    collapse = ", "
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 2.2 — Metadata overview
# ----------------------------------------------------------------------------

metadata_overview <- seurat_obj@meta.data %>%
  summarise(
    Cells = n(),
    Samples = n_distinct(GSM),
    Donors = n_distinct(Donor),
    Clinical_States = n_distinct(Phase)
  )

print(metadata_overview)

cat("\n")


# ----------------------------------------------------------------------------
# 2.3 — Sample-level cell counts
# ----------------------------------------------------------------------------

dataset_snapshot <- seurat_obj@meta.data %>%
  dplyr::count(
    GSM,
    Donor,
    Phase,
    name = "Cells"
  ) %>%
  dplyr::arrange(
    Phase,
    GSM
  )

cat("Cell counts by sample:\n\n")

print(dataset_snapshot)

write_csv(
  dataset_snapshot,
  "results/tables/phase3_dataset_snapshot.csv"
)

cat("\n✓ Dataset snapshot saved\n\n")


# ============================================================================
# SECTION 3 — CONFIRM ONLY SINGLET CELLS ARE RETAINED
# ============================================================================

cat("=== SECTION 3: SINGLET VALIDATION ===\n\n")

if (
  !"scDblFinder.class" %in%
  colnames(seurat_obj@meta.data)
) {
  
  stop(
    "ERROR: 'scDblFinder.class' metadata column is not present. ",
    "Cannot verify singlet retention."
  )
  
}

singlet_table <- table(
  seurat_obj$scDblFinder.class,
  useNA = "ifany"
)

cat("scDblFinder classification:\n\n")

print(singlet_table)

cat("\n")


# ----------------------------------------------------------------------------
# 3.1 — Explicit validation
# ----------------------------------------------------------------------------

non_singlet_cells <- sum(
  seurat_obj$scDblFinder.class != "singlet",
  na.rm = TRUE
)

missing_classification <- sum(
  is.na(
    seurat_obj$scDblFinder.class
  )
)

if (non_singlet_cells > 0) {
  
  stop(
    "ERROR: ",
    non_singlet_cells,
    " non-singlet cells remain in the input object."
  )
  
}

if (missing_classification > 0) {
  
  stop(
    "ERROR: ",
    missing_classification,
    " cells have missing scDblFinder classifications."
  )
  
}

cat(
  "✓ All retained cells are classified as singlets\n"
)

cat(
  "✓ No doublets retained\n"
)

cat(
  "✓ No cells have missing doublet classifications\n\n"
)


# ----------------------------------------------------------------------------
# 3.2 — Save singlet validation summary
# ----------------------------------------------------------------------------

singlet_validation <- tibble(
  
  Total_Cells = ncol(seurat_obj),
  
  Singlet_Cells = sum(
    seurat_obj$scDblFinder.class == "singlet"
  ),
  
  Doublet_Cells = sum(
    seurat_obj$scDblFinder.class == "doublet"
  ),
  
  Missing_Classification = sum(
    is.na(
      seurat_obj$scDblFinder.class
    )
  )
  
)

write_csv(
  singlet_validation,
  "results/tables/phase3_singlet_validation.csv"
)

cat(
  "✓ Singlet validation summary saved\n\n"
)


# ============================================================================
# SECTION 4 — LOG-NORMALIZATION
# ============================================================================

cat("=== SECTION 4: LOG-NORMALIZATION ===\n\n")

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- NormalizeData(
  
  object = seurat_obj,
  
  assay = "RNA",
  
  normalization.method = "LogNormalize",
  
  scale.factor = 10000,
  
  verbose = TRUE
  
)

cat(
  "\n✓ Log-normalization complete\n\n"
)


# ============================================================================
# SECTION 5 — HIGHLY VARIABLE FEATURES
# ============================================================================

cat("=== SECTION 5: HIGHLY VARIABLE GENES ===\n\n")

seurat_obj <- FindVariableFeatures(
  
  object = seurat_obj,
  
  assay = "RNA",
  
  selection.method = "vst",
  
  nfeatures = 2000,
  
  verbose = TRUE
  
)

hvg <- VariableFeatures(
  seurat_obj
)

cat(
  "\nNumber of highly variable genes:",
  length(hvg),
  "\n"
)

if (length(hvg) != 2000) {
  
  stop(
    "ERROR: Expected exactly 2,000 highly variable genes, but ",
    length(hvg),
    " were identified."
  )
  
}

cat(
  "✓ Exactly 2,000 HVGs identified using VST\n\n"
)


# ----------------------------------------------------------------------------
# 5.1 — Save HVG list
# ----------------------------------------------------------------------------

hvg_table <- tibble(
  Rank = seq_along(hvg),
  Gene = hvg
)

write_csv(
  hvg_table,
  "results/tables/phase3_highly_variable_genes.csv"
)

cat(
  "✓ HVG list saved\n\n"
)


# ============================================================================
# SECTION 6 — GLOBAL 20,000-CELL LEVERAGE-SCORE SKETCH
# ============================================================================
#
# The QC-filtered dataset contains approximately 105,000 singlet cells
# across 23 sample/donor layers.
#
# A representative global sketch of exactly 20,000 cells is constructed.
#
# IMPORTANT:
#   Seurat's SketchData() treats a single ncells value as a per-layer
#   request. Therefore, it is NOT used here to construct the global
#   20,000-cell sketch.
#
# Instead:
#   1. Sample-specific RNA data layers are identified.
#   2. Leverage scores are calculated independently within each sample.
#   3. Exactly proportional numbers of cells are allocated to each sample.
#   4. Cells are sampled without replacement with probability proportional
#      to leverage score.
#   5. The selected cells are extracted from the original object to create
#      the actual 20,000-cell sketch object.
#
# The full object remains unchanged.
#
# ============================================================================

cat("=== SECTION 6: GLOBAL LEVERAGE-SCORE SKETCH ===\n\n")

target_sketch_cells <- 20000L

total_cells <- ncol(seurat_obj)

cat(
  "Full dataset:",
  total_cells,
  "cells\n"
)

cat(
  "Target global sketch:",
  target_sketch_cells,
  "cells\n"
)

cat(
  "Fraction represented:",
  round(
    100 *
      target_sketch_cells /
      total_cells,
    2
  ),
  "%\n\n"
)

if (
  target_sketch_cells >= total_cells
) {
  
  stop(
    "ERROR: Target sketch size must be smaller than the full dataset."
  )
}


# ----------------------------------------------------------------------------
# 6.1 — Identify sample-specific data layers
# ----------------------------------------------------------------------------

data_layers <- Layers(
  seurat_obj[["RNA"]],
  search = "^data\\."
)

if (length(data_layers) == 0) {
  
  stop(
    "ERROR: No sample-specific RNA data layers were found."
  )
}

cat(
  "Number of sample-specific data layers:",
  length(data_layers),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 6.2 — Verify that every layer corresponds to a GSM
# ----------------------------------------------------------------------------

layer_GSM <- sub(
  "^data\\.",
  "",
  data_layers
)

if (
  !all(
    layer_GSM %in%
    unique(seurat_obj$GSM)
  )
) {
  
  stop(
    "ERROR: One or more RNA data layers could not be matched to a GSM."
  )
}


# ----------------------------------------------------------------------------
# 6.3 — Calculate sample sizes
# ----------------------------------------------------------------------------

sample_sizes <- seurat_obj@meta.data %>%
  dplyr::count(
    GSM,
    name = "Original_Cells"
  ) %>%
  dplyr::arrange(GSM)

cat(
  "Sample sizes:\n\n"
)

print(sample_sizes)

cat("\n")


# ----------------------------------------------------------------------------
# 6.4 — Allocate exactly 20,000 cells across samples
# ----------------------------------------------------------------------------
#
# Proportional allocation with largest-remainder rounding.
#
# ----------------------------------------------------------------------------

sample_sizes <- sample_sizes %>%
  mutate(
    
    Exact_Allocation =
      Original_Cells /
      sum(Original_Cells) *
      target_sketch_cells,
    
    Base_Allocation =
      floor(
        Exact_Allocation
      ),
    
    Fractional_Remainder =
      Exact_Allocation -
      Base_Allocation
    
  )

remaining_cells <- target_sketch_cells -
  sum(
    sample_sizes$Base_Allocation
  )

sample_sizes$Sketch_Cells <-
  sample_sizes$Base_Allocation

if (remaining_cells > 0) {
  
  extra_indices <- order(
    sample_sizes$Fractional_Remainder,
    decreasing = TRUE
  )[seq_len(remaining_cells)]
  
  sample_sizes$Sketch_Cells[
    extra_indices
  ] <-
    sample_sizes$Sketch_Cells[
      extra_indices
    ] + 1L
}

sample_sizes <- sample_sizes %>%
  mutate(
    
    Sketch_Fraction =
      Sketch_Cells /
      Original_Cells,
    
    Sketch_Percent =
      100 *
      Sketch_Fraction
    
  )

cat(
  "Proportional sketch allocation:\n\n"
)

print(sample_sizes)

cat(
  "\nTotal allocated cells:",
  sum(sample_sizes$Sketch_Cells),
  "\n\n"
)

if (
  sum(sample_sizes$Sketch_Cells) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Sample allocations do not sum to exactly 20,000 cells."
  )
}


# ----------------------------------------------------------------------------
# 6.5 — Calculate leverage scores within each sample
# ----------------------------------------------------------------------------

cat(
  "=== CALCULATING LEVERAGE SCORES ===\n\n"
)

leverage_scores <- vector(
  mode = "list",
  length = length(data_layers)
)

names(leverage_scores) <- data_layers


for (
  i in seq_along(data_layers)
) {
  
  layer_name <- data_layers[i]
  
  GSM_id <- sub(
    "^data\\.",
    "",
    layer_name
  )
  
  cat(
    "Calculating leverage scores:",
    GSM_id,
    "\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Extract normalized expression for the 2,000 HVGs
  # --------------------------------------------------------------------------
  
  layer_matrix <- LayerData(
    object = seurat_obj[["RNA"]],
    layer = layer_name,
    features = VariableFeatures(seurat_obj)
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
  # Calculate leverage scores
  # --------------------------------------------------------------------------
  
  score_vector <- LeverageScore(
    
    object = layer_matrix,
    
    nsketch = min(
      5000L,
      ncol(layer_matrix)
    ),
    
    seed = 12345,
    
    verbose = FALSE
    
  )
  
  
  score_vector <- as.numeric(
    score_vector
  )
  
  names(score_vector) <- colnames(
    layer_matrix
  )
  
  
  if (
    length(score_vector) !=
    ncol(layer_matrix)
  ) {
    
    stop(
      "ERROR: Leverage-score vector length does not match cell count ",
      "for ",
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
  
  
  leverage_scores[[layer_name]] <-
    score_vector
  
  
  cat(
    "  Cells:",
    length(score_vector),
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
}


cat(
  "✓ Leverage scores calculated for all sample layers\n\n"
)


# ----------------------------------------------------------------------------
# 6.6 — Select cells using leverage-score weighted sampling
# ----------------------------------------------------------------------------

cat(
  "=== SELECTING REPRESENTATIVE CELLS ===\n\n"
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
  
  layer_name <- paste0(
    "data.",
    GSM_id
  )
  
  scores <-
    leverage_scores[[layer_name]]
  
  
  if (
    is.null(scores)
  ) {
    
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
      " cells are available."
    )
  }
  
  
  set.seed(
    12345 + i
  )
  
  selected_cells <- sample(
    
    x = names(scores),
    
    size = n_to_select,
    
    replace = FALSE,
    
    prob = scores
    
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


# ----------------------------------------------------------------------------
# 6.7 — Combine selected cells
# ----------------------------------------------------------------------------

sketch_cells <- unlist(
  selected_cells_by_sample,
  use.names = FALSE
)

sketch_cells <- unique(
  sketch_cells
)

cat(
  "\nTotal unique selected cells:",
  length(sketch_cells),
  "\n"
)

if (
  length(sketch_cells) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Exactly ",
    target_sketch_cells,
    " unique cells were expected, but ",
    length(sketch_cells),
    " were selected."
  )
}

cat(
  "✓ Exactly 20,000 unique cells selected\n\n"
)


# ----------------------------------------------------------------------------
# 6.8 — Construct actual 20,000-cell sketch object
# ----------------------------------------------------------------------------

cat(
  "Creating the 20,000-cell sketch object...\n\n"
)

sketch_obj <- subset(
  x = seurat_obj,
  cells = sketch_cells
)


# ----------------------------------------------------------------------------
# 6.9 — Validate sketch dimensions
# ----------------------------------------------------------------------------

if (
  ncol(sketch_obj) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Sketch object contains ",
    ncol(sketch_obj),
    " cells instead of exactly ",
    target_sketch_cells,
    "."
  )
}

cat(
  "✓ sketch_obj contains exactly",
  ncol(sketch_obj),
  "cells\n"
)

cat(
  "✓ Feature set:",
  nrow(sketch_obj),
  "genes\n\n"
)


# ============================================================================
# SECTION 7 — SKETCH REPRESENTATION BY SAMPLE
# ============================================================================

cat("=== SECTION 7: SKETCH REPRESENTATION BY SAMPLE ===\n\n")


# ----------------------------------------------------------------------------
# 7.1 — Original sample counts
# ----------------------------------------------------------------------------

original_sample_counts <- table(
  seurat_obj$GSM
)

original_sample_counts <- as.numeric(
  original_sample_counts
)

names(original_sample_counts) <-
  names(
    table(seurat_obj$GSM)
  )


# ----------------------------------------------------------------------------
# 7.2 — Sketch sample counts
# ----------------------------------------------------------------------------

sketch_sample_counts <- table(
  sketch_obj$GSM
)

sketch_sample_counts <- as.numeric(
  sketch_sample_counts
)

names(sketch_sample_counts) <-
  names(
    table(sketch_obj$GSM)
  )


# ----------------------------------------------------------------------------
# 7.3 — Construct representation table
# ----------------------------------------------------------------------------

sketch_representation <- tibble(
  
  GSM =
    names(original_sample_counts),
  
  Original_Cells =
    unname(original_sample_counts),
  
  Sketch_Cells =
    unname(
      sketch_sample_counts[
        names(original_sample_counts)
      ]
    )
  
) %>%
  mutate(
    
    Sketch_Fraction =
      Sketch_Cells /
      Original_Cells,
    
    Sketch_Percent =
      100 *
      Sketch_Fraction
    
  )


# ----------------------------------------------------------------------------
# 7.4 — Validate representation
# ----------------------------------------------------------------------------

if (
  any(
    is.na(
      sketch_representation$Sketch_Cells
    )
  )
) {
  
  stop(
    "ERROR: One or more original samples are absent from the sketch."
  )
}

if (
  any(
    sketch_representation$Sketch_Cells < 1
  )
) {
  
  stop(
    "ERROR: One or more samples have zero cells in the sketch."
  )
}

if (
  sum(
    sketch_representation$Sketch_Cells
  ) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Sketch sample counts do not sum to 20,000."
  )
}

if (
  nrow(sketch_representation) !=
  nrow(dataset_snapshot)
) {
  
  stop(
    "ERROR: Number of samples represented in sketch does not match ",
    "the original dataset."
  )
}


# ----------------------------------------------------------------------------
# 7.5 — Display representation
# ----------------------------------------------------------------------------

cat(
  "Sample-level sketch representation:\n\n"
)

print(
  sketch_representation,
  n = Inf
)

cat("\n")

cat(
  "Original dataset cells:",
  sum(sketch_representation$Original_Cells),
  "\n"
)

cat(
  "Sketch cells:",
  sum(sketch_representation$Sketch_Cells),
  "\n"
)

cat(
  "Number of samples:",
  nrow(sketch_representation),
  "\n"
)

cat(
  "Minimum sample representation:",
  round(
    min(sketch_representation$Sketch_Percent),
    2
  ),
  "%\n"
)

cat(
  "Maximum sample representation:",
  round(
    max(sketch_representation$Sketch_Percent),
    2
  ),
  "%\n\n"
)


# ----------------------------------------------------------------------------
# 7.6 — Save representation table
# ----------------------------------------------------------------------------

write_csv(
  sketch_representation,
  "results/tables/phase3_sketch_representation_by_sample.csv"
)

cat(
  "✓ Sample-level sketch representation table saved.\n\n"
)


# ----------------------------------------------------------------------------
# 7.7 — Plot original versus sketch representation
# ----------------------------------------------------------------------------

representation_plot_data <-
  sketch_representation %>%
  dplyr::select(
    GSM,
    Original_Cells,
    Sketch_Cells
  ) %>%
  tidyr::pivot_longer(
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
      "Original Dataset and 20,000-Cell Sketch Representation",
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

cat(
  "✓ Sketch representation plot saved.\n\n"
)

cat(
  "=== SECTION 7 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 8 — SCALE THE 20,000-CELL SKETCH
# ============================================================================

cat("=== SECTION 8: SCALING THE 20,000-CELL SKETCH ===\n\n")

sketch_cells <- colnames(
  sketch_obj
)

n_sketch_cells <- length(
  sketch_cells
)

cat(
  "Sketch cells identified:",
  n_sketch_cells,
  "\n"
)

cat(
  "Target sketch size:",
  target_sketch_cells,
  "\n\n"
)

if (
  n_sketch_cells !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: Expected ",
    target_sketch_cells,
    " sketch cells, but found ",
    n_sketch_cells,
    "."
  )
}

DefaultAssay(sketch_obj) <- "RNA"

hvg <- VariableFeatures(
  seurat_obj
)

if (
  length(hvg) != 2000
) {
  
  stop(
    "ERROR: Expected 2,000 HVGs from Section 5, but found ",
    length(hvg),
    "."
  )
}

sketch_hvg <- intersect(
  hvg,
  rownames(sketch_obj[["RNA"]])
)

if (
  length(sketch_hvg) != 2000
) {
  
  stop(
    "ERROR: Expected all 2,000 HVGs to be present in the sketch, ",
    "but only ",
    length(sketch_hvg),
    " were found."
  )
}

VariableFeatures(sketch_obj) <-
  sketch_hvg

cat(
  "Scaling:",
  length(sketch_hvg),
  "HVGs\n"
)

cat(
  "Cells:",
  n_sketch_cells,
  "\n\n"
)


# ----------------------------------------------------------------------------
# 8.1 — Scale sketch
# ----------------------------------------------------------------------------

sketch_obj <- ScaleData(
  
  object = sketch_obj,
  
  assay = "RNA",
  
  features = sketch_hvg,
  
  verbose = TRUE
  
)

cat(
  "\n✓ Scaling completed successfully\n\n"
)


# ----------------------------------------------------------------------------
# 8.2 — Validate scaled data
# ----------------------------------------------------------------------------

scaled_layers <- Layers(
  sketch_obj[["RNA"]],
  search = "^scale.data"
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

cat(
  "✓ Scaled-data layer successfully created\n\n"
)


# ----------------------------------------------------------------------------
# 8.3 — Save checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  sketch_obj,
  "results/rds_objects/phase3_sketch_scaled.rds"
)

cat(
  "✓ Scaled sketch checkpoint saved:\n",
  "  results/rds_objects/phase3_sketch_scaled.rds\n\n"
)

cat(
  "=== SECTION 8 COMPLETE ===\n\n"
)

# ============================================================================
# SECTION 9 — SEPARATE THE 20,000-CELL ATLAS SKETCH
# ============================================================================

cat("=== SECTION 9: SEPARATING THE 20,000-CELL ATLAS SKETCH ===\n\n")


# ----------------------------------------------------------------------------
# 9.1 — Load scaled sketch checkpoint
# ----------------------------------------------------------------------------

sketch_obj <- readRDS(
  "results/rds_objects/phase3_sketch_scaled.rds"
)

DefaultAssay(sketch_obj) <- "RNA"


# ----------------------------------------------------------------------------
# 9.2 — Identify sketch cells
# ----------------------------------------------------------------------------

sketch_cells <- colnames(
  sketch_obj
)

n_sketch_cells <- length(
  sketch_cells
)

cat(
  "Cells represented by sketch:",
  n_sketch_cells,
  "\n"
)

cat(
  "Expected sketch size:",
  target_sketch_cells,
  "\n\n"
)

if (
  n_sketch_cells !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: The sketch contains ",
    n_sketch_cells,
    " cells rather than the expected ",
    target_sketch_cells,
    "."
  )
}


# ----------------------------------------------------------------------------
# 9.3 — Create genuinely separate atlas sketch
# ----------------------------------------------------------------------------

atlas_sketch <- subset(
  x = sketch_obj,
  cells = sketch_cells
)


# ----------------------------------------------------------------------------
# 9.4 — Validate separated object
# ----------------------------------------------------------------------------

if (
  ncol(atlas_sketch) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: atlas_sketch contains ",
    ncol(atlas_sketch),
    " cells rather than exactly ",
    target_sketch_cells,
    "."
  )
}

if (
  !setequal(
    colnames(atlas_sketch),
    sketch_cells
  )
) {
  
  stop(
    "ERROR: atlas_sketch does not contain exactly the selected cells."
  )
}


# ----------------------------------------------------------------------------
# 9.5 — Confirm sample representation
# ----------------------------------------------------------------------------

atlas_sketch_sample_counts <- table(
  atlas_sketch$GSM
)

cat(
  "Atlas sketch cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "Genes:",
  nrow(atlas_sketch),
  "\n"
)

cat(
  "Samples represented:",
  length(atlas_sketch_sample_counts),
  "\n\n"
)

print(
  atlas_sketch_sample_counts
)

cat("\n")

if (
  length(atlas_sketch_sample_counts) !=
  dplyr::n_distinct(seurat_obj$GSM)
) {
  
  stop(
    "ERROR: Not all original samples are represented in the atlas sketch."
  )
}


# ----------------------------------------------------------------------------
# 9.6 — Confirm HVGs
# ----------------------------------------------------------------------------

atlas_sketch_hvg <- intersect(
  
  VariableFeatures(seurat_obj),
  
  rownames(atlas_sketch)
  
)

if (
  length(atlas_sketch_hvg) != 2000
) {
  
  stop(
    "ERROR: Expected all 2,000 HVGs to be present in atlas_sketch, ",
    "but found ",
    length(atlas_sketch_hvg),
    "."
  )
}

VariableFeatures(atlas_sketch) <-
  atlas_sketch_hvg

cat(
  "HVGs available:",
  length(atlas_sketch_hvg),
  "\n"
)

cat(
  "✓ All 2,000 HVGs retained\n\n"
)


# ----------------------------------------------------------------------------
# 9.7 — Save atlas sketch
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch.rds"
)

cat(
  "✓ Genuine 20,000-cell atlas sketch created\n"
)

cat(
  "✓ Saved to:\n",
  "  results/rds_objects/phase3_atlas_sketch.rds\n\n"
)

cat(
  "=== SECTION 9 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 10 — PCA ON THE 20,000-CELL ATLAS SKETCH
# ============================================================================

cat("=== SECTION 10: PCA — 50 PRINCIPAL COMPONENTS ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"


# ----------------------------------------------------------------------------
# 10.1 — Confirm PCA input
# ----------------------------------------------------------------------------

if (
  ncol(atlas_sketch) !=
  target_sketch_cells
) {
  
  stop(
    "ERROR: PCA input does not contain exactly 20,000 cells."
  )
}

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
# 10.2 — Run PCA
# ----------------------------------------------------------------------------

atlas_sketch <- RunPCA(
  
  object = atlas_sketch,
  
  assay = "RNA",
  
  features = VariableFeatures(atlas_sketch),
  
  npcs = 50,
  
  verbose = TRUE
  
)

cat(
  "\n✓ PCA completed\n"
)

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
# 10.3 — Validate PCA
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
    "ERROR: PCA embedding contains ",
    nrow(pca_embeddings),
    " cells rather than ",
    target_sketch_cells,
    "."
  )
}

if (
  ncol(pca_embeddings) != 50
) {
  
  stop(
    "ERROR: Expected 50 PCs, but ",
    ncol(pca_embeddings),
    " were calculated."
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

cat(
  "✓ PCA contains exactly 20,000 cells and 50 PCs\n\n"
)


# ----------------------------------------------------------------------------
# 10.4 — Generate elbow plot
# ----------------------------------------------------------------------------

p_elbow <- ElbowPlot(
  
  object = atlas_sketch,
  
  reduction = "pca",
  
  ndims = 50
  
) +
  labs(
    title =
      "PCA Elbow Plot — 20,000-Cell Atlas Sketch",
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


# ----------------------------------------------------------------------------
# 10.5 — Save elbow plot
# ----------------------------------------------------------------------------

ggsave(
  "results/figures/phase3_pca_elbow_plot.png",
  p_elbow,
  width = 9,
  height = 6,
  dpi = 300
)

cat(
  "✓ PCA elbow plot saved:\n",
  "  results/figures/phase3_pca_elbow_plot.png\n\n"
)


# ----------------------------------------------------------------------------
# 10.6 — Save PCA checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_pca50.rds"
)

cat(
  "✓ PCA checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_pca50.rds\n\n"
)

cat(
  "=== SECTION 10 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 11 — SELECT THE FIRST 20 PRINCIPAL COMPONENTS
# ============================================================================

cat("=== SECTION 11: SELECTING THE FIRST 20 PCs ===\n\n")

selected_pcs <- 1:20

cat(
  "Selected principal components:",
  paste(
    selected_pcs,
    collapse = ", "
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 11.1 — Validate selected dimensions
# ----------------------------------------------------------------------------

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

if (
  max(selected_pcs) > available_pcs
) {
  
  stop(
    "ERROR: Selected PC exceeds the number of calculated PCs."
  )
}

cat(
  "✓ First 20 PCs selected for downstream analysis\n"
)

cat(
  "✓ PC range:",
  min(selected_pcs),
  "to",
  max(selected_pcs),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 11.2 — Save selected PC information
# ----------------------------------------------------------------------------

selected_pc_table <- tibble(
  PC = selected_pcs
)

write_csv(
  selected_pc_table,
  "results/tables/phase3_selected_pcs.csv"
)

cat(
  "✓ Selected PC list saved:\n",
  "  results/tables/phase3_selected_pcs.csv\n\n"
)


# ----------------------------------------------------------------------------
# 11.3 — Save checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_pca50.rds"
)

cat(
  "=== SECTION 11 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 12 — NEIGHBOUR GRAPH + CLUSTERING
# ============================================================================

cat("=== SECTION 12: NEIGHBOUR GRAPH + CLUSTERING ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_pca50.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

selected_pcs <- 1:20

cat(
  "Cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "Available PCs:",
  ncol(
    Embeddings(
      atlas_sketch,
      reduction = "pca"
    )
  ),
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
# 12.1 — Validate PCA dimensions
# ----------------------------------------------------------------------------

available_pcs <- ncol(
  Embeddings(
    atlas_sketch,
    reduction = "pca"
  )
)

if (
  available_pcs < max(selected_pcs)
) {
  
  stop(
    "ERROR: Requested PCs exceed the number of available PCs."
  )
}


# ----------------------------------------------------------------------------
# 12.2 — Construct neighbour graph
# ----------------------------------------------------------------------------

cat(
  "Constructing neighbour graph...\n\n"
)

atlas_sketch <- FindNeighbors(
  
  object = atlas_sketch,
  
  reduction = "pca",
  
  dims = selected_pcs,
  
  verbose = TRUE
  
)

cat(
  "\n✓ Neighbour graph constructed\n\n"
)


# ----------------------------------------------------------------------------
# 12.3 — Validate neighbour graph
# ----------------------------------------------------------------------------

graph_names <- names(
  atlas_sketch@graphs
)

cat(
  "Graphs currently stored in object:\n\n"
)

print(graph_names)

cat("\n")

if (
  length(graph_names) == 0
) {
  
  stop(
    "ERROR: No neighbour graph was created."
  )
}

cat(
  "✓ Neighbour graph successfully detected\n\n"
)


# ----------------------------------------------------------------------------
# 12.4 — Cluster cells
# ----------------------------------------------------------------------------

cat(
  "Clustering cells at resolution 0.5...\n\n"
)

set.seed(12345)

atlas_sketch <- FindClusters(
  
  object = atlas_sketch,
  
  resolution = 0.5,
  
  verbose = TRUE
  
)

cat(
  "\n✓ Clustering completed\n\n"
)


# ----------------------------------------------------------------------------
# 12.5 — Validate cluster assignment
# ----------------------------------------------------------------------------

if (
  !"seurat_clusters" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "ERROR: 'seurat_clusters' was not created."
  )
}

cluster_table <- table(
  atlas_sketch$seurat_clusters
)

cat(
  "Number of clusters:",
  length(cluster_table),
  "\n\n"
)

cat(
  "Cells per cluster:\n\n"
)

print(cluster_table)

cat("\n")

if (
  sum(cluster_table) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: Cluster assignments do not account for all cells."
  )
}

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

cat(
  "✓ Every atlas cell has a cluster assignment\n\n"
)


# ----------------------------------------------------------------------------
# 12.6 — Save cluster-size table
# ----------------------------------------------------------------------------

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

cat(
  "✓ Cluster-size table saved\n\n"
)


# ----------------------------------------------------------------------------
# 12.7 — Save clustering checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_clustered.rds"
)

cat(
  "✓ Clustered atlas sketch checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_clustered.rds\n\n"
)

cat(
  "=== SECTION 12 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 13 — UMAP
# ============================================================================

cat("=== SECTION 13: UMAP ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_clustered.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

selected_pcs <- 1:20

cat(
  "UMAP input cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "UMAP input dimensions:",
  paste(
    selected_pcs,
    collapse = ", "
  ),
  "\n\n"
)


# ----------------------------------------------------------------------------
# 13.1 — Run UMAP
# ----------------------------------------------------------------------------

cat(
  "Running UMAP...\n\n"
)

set.seed(12345)

atlas_sketch <- RunUMAP(
  
  object = atlas_sketch,
  
  reduction = "pca",
  
  dims = selected_pcs,
  
  verbose = TRUE
  
)

cat(
  "\n✓ UMAP completed\n\n"
)


# ----------------------------------------------------------------------------
# 13.2 — Validate UMAP
# ----------------------------------------------------------------------------

umap_embeddings <- Embeddings(
  atlas_sketch,
  reduction = "umap"
)

cat(
  "UMAP dimensions:",
  ncol(umap_embeddings),
  "\n"
)

cat(
  "UMAP cells:",
  nrow(umap_embeddings),
  "\n\n"
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
    "ERROR: UMAP should contain exactly two dimensions."
  )
}

if (
  any(
    !is.finite(umap_embeddings)
  )
) {
  
  stop(
    "ERROR: Non-finite values detected in UMAP embeddings."
  )
}

cat(
  "✓ UMAP contains exactly 20,000 cells and 2 dimensions\n\n"
)


# ----------------------------------------------------------------------------
# 13.3 — Save UMAP checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_umap.rds"
)

cat(
  "✓ UMAP checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_umap.rds\n\n"
)

cat(
  "=== SECTION 13 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 14 — GLOBAL UMAP BY CLUSTER
# ============================================================================

cat("=== SECTION 14: GLOBAL UMAP BY CLUSTER ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_umap.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"


# ----------------------------------------------------------------------------
# 14.1 — Generate cluster UMAP
# ----------------------------------------------------------------------------

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

print(
  p_global_cluster
)


# ----------------------------------------------------------------------------
# 14.2 — Save plot
# ----------------------------------------------------------------------------

ggsave(
  "results/figures/phase3_global_umap_by_cluster.png",
  p_global_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

cat(
  "\n✓ Global cluster UMAP saved:\n",
  "  results/figures/phase3_global_umap_by_cluster.png\n\n"
)

cat(
  "=== SECTION 14 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 15 — UMAP BY CLINICAL STATE
# ============================================================================

cat("=== SECTION 15: UMAP BY CLINICAL STATE ===\n\n")


# ----------------------------------------------------------------------------
# 15.1 — Validate clinical-state metadata
# ----------------------------------------------------------------------------

if (
  !"Phase" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "ERROR: 'Phase' metadata column is missing."
  )
}

clinical_states <- unique(
  atlas_sketch$Phase
)

cat(
  "Clinical states represented:\n\n"
)

print(
  sort(clinical_states)
)

cat("\n")


# ----------------------------------------------------------------------------
# 15.2 — Generate clinical-state UMAP
# ----------------------------------------------------------------------------

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

print(
  p_umap_phase
)


# ----------------------------------------------------------------------------
# 15.3 — Save plot
# ----------------------------------------------------------------------------

ggsave(
  "results/figures/phase3_global_umap_by_clinical_state.png",
  p_umap_phase,
  width = 10,
  height = 8,
  dpi = 300
)

cat(
  "\n✓ Clinical-state UMAP saved:\n",
  "  results/figures/phase3_global_umap_by_clinical_state.png\n\n"
)

cat(
  "=== SECTION 15 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 16 — UMAP BY DONOR
# ============================================================================

cat("=== SECTION 16: UMAP BY DONOR ===\n\n")


# ----------------------------------------------------------------------------
# 16.1 — Validate donor metadata
# ----------------------------------------------------------------------------

if (
  !"Donor" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "ERROR: 'Donor' metadata column is missing."
  )
}

donors <- unique(
  atlas_sketch$Donor
)

cat(
  "Number of donors represented:",
  length(donors),
  "\n\n"
)

if (
  length(donors) != 23
) {
  
  stop(
    "ERROR: Expected 23 donors, but found ",
    length(donors),
    "."
  )
}


# ----------------------------------------------------------------------------
# 16.2 — Generate donor UMAP
# ----------------------------------------------------------------------------

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

print(
  p_umap_donor
)


# ----------------------------------------------------------------------------
# 16.3 — Save donor UMAP
# ----------------------------------------------------------------------------

ggsave(
  "results/figures/phase3_global_umap_by_donor.png",
  p_umap_donor,
  width = 12,
  height = 9,
  dpi = 300
)

cat(
  "\n✓ Donor UMAP saved:\n",
  "  results/figures/phase3_global_umap_by_donor.png\n\n"
)


# ----------------------------------------------------------------------------
# 16.4 — Save final Sections 1–16 checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_umap_clustered.rds"
)

cat(
  "✓ Final Sections 1–16 atlas checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_umap_clustered.rds\n\n"
)


# ============================================================================
# FINAL CHECKPOINT — SECTIONS 1–16
# ============================================================================

cat(
  "============================================================\n"
)

cat(
  "✓ PHASE 3 SCRIPT 01 — SECTIONS 1–16 COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Full dataset cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Full dataset genes:",
  nrow(seurat_obj),
  "\n"
)

cat(
  "Atlas sketch cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "HVGs:",
  length(
    VariableFeatures(atlas_sketch)
  ),
  "\n"
)

cat(
  "Clusters:",
  length(
    unique(
      atlas_sketch$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "UMAP dimensions:",
  ncol(
    Embeddings(
      atlas_sketch,
      "umap"
    )
  ),
  "\n"
)

cat(
  "Clinical states:",
  length(
    unique(
      atlas_sketch$Phase
    )
  ),
  "\n"
)

cat(
  "Donors:",
  length(
    unique(
      atlas_sketch$Donor
    )
  ),
  "\n\n"
)

cat(
  "Completed:\n",
  "  ✓ Environment setup\n",
  "  ✓ QC-filtered singlet object loaded\n",
  "  ✓ Metadata validated\n",
  "  ✓ Singlet retention validated\n",
  "  ✓ Log-normalization completed\n",
  "  ✓ 2,000 HVGs identified using VST\n",
  "  ✓ Exactly 20,000-cell leverage-score sketch constructed\n",
  "  ✓ All sample layers represented\n",
  "  ✓ Sketch representation validated\n",
  "  ✓ 20,000-cell sketch scaled\n",
  "  ✓ Separate atlas sketch created\n",
  "  ✓ PCA calculated to 50 components\n",
  "  ✓ PCs 1–20 selected\n",
  "  ✓ Neighbour graph constructed\n",
  "  ✓ Global clustering completed\n",
  "  ✓ UMAP calculated using PCs 1–20\n",
  "  ✓ Global UMAP by cluster generated\n",
  "  ✓ UMAP by clinical state generated\n",
  "  ✓ UMAP by donor generated\n",
  "  ✓ Final Sections 1–16 checkpoint saved\n\n"
)

cat(
  "Final checkpoint:\n",
  "  results/rds_objects/phase3_atlas_sketch_umap_clustered.rds\n\n"
)

cat(
  "Ready for phase_03_02_cluster_markers_annotation.R\n"
)

cat(
  "============================================================\n"
)

