###############################################################################
# PHASE 3 — FULL-DATASET PCA PROJECTION + LABEL TRANSFER + FINAL ATLAS
#
# Purpose:
#   1. Load the complete Phase 2 processed-expression dataset.
#   2. Project all 106,592 cells into the validated PCA space learned from the
#      20,000-cell Phase 3 atlas.
#   3. Transfer atlas cell-type labels using PCA-space k-nearest neighbours.
#   4. Quantify neighbour agreement and transfer confidence.
#   5. Perform targeted biological QC of transferred CD8_T and
#      Monocyte_Myeloid populations.
#   6. Generate a NEW full-dataset UMAP from projected PCA coordinates.
#   7. Save the final full-dataset annotated Seurat object.
#
# IMPORTANT REPRESENTATION NOTE:
#   The GEO matrices are processed log-CP10K expression values.
#   They are NOT raw UMI counts.
#
#   Therefore this script:
#     - does NOT reconstruct raw counts
#     - does NOT run NormalizeData()
#     - does NOT calculate conventional nCount_RNA
#     - does NOT perform another normalization step
#     - does NOT perform batch correction
#     - does NOT filter cells
#
# IMPORTANT PROJECTION NOTE:
#   The Phase 2 RNA assay contains 23 sample-specific data layers:
#
#       data.GSM5519469
#       data.GSM5519471
#       ...
#
#   Script 03 does NOT JoinLayers().
#
#   Each sample-specific processed-expression layer is projected
#   independently into the SAME reference PCA space.
#
#   The projection uses:
#     - the 1,999 genes actually represented in the atlas PCA loadings
#     - reference-derived feature means
#     - reference-derived feature SDs
#     - the existing processed log-CP10K values
#     - the atlas PCA loadings
#
#   The resulting 50-dimensional projections are combined only after
#   projection. The large expression matrices are never joined.
#
# Input:
#   results/rds_objects/seurat_qc_processed_expression.rds
#   results/rds_objects/phase3_atlas_sketch_final_annotation.rds
#
# Output:
#   results/rds_objects/phase3_final_full_dataset.rds
#
# Phase 3 atlas:
#   20,000 reference cells
#   20 clusters
#   PCA = 50 dimensions
#   Atlas graph = PCs 1:20
#
# Full dataset:
#   106,592 cells
#   24,452 genes
#   23 donors
#   5 clinical states
#
# Software:
#   Seurat 5.5.1
#   SeuratObject 5.4.0
###############################################################################


###############################################################################
# SECTION 25 — SETUP + LOAD FULL DATASET + ATLAS
###############################################################################

rm(list = ls())

set.seed(20260918)


###############################################################################
# Package loading
###############################################################################

suppressPackageStartupMessages({
  
  library(Seurat)
  library(SeuratObject)
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  library(Matrix)
  
  library(RANN)
  
  library(ggplot2)
  library(pheatmap)
})


###############################################################################
# Explicit package validation
###############################################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "tibble",
  "Matrix",
  "RANN",
  "ggplot2",
  "pheatmap"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0L) {
  
  stop(
    "Required packages are missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    call. = FALSE
  )
}


###############################################################################
# Version report
###############################################################################

cat("\n============================================================\n")
cat("PHASE 3 — SOFTWARE ENVIRONMENT\n")
cat("============================================================\n")

cat(
  "Seurat: ",
  as.character(packageVersion("Seurat")),
  "\n",
  sep = ""
)

cat(
  "SeuratObject: ",
  as.character(packageVersion("SeuratObject")),
  "\n",
  sep = ""
)


###############################################################################
# Version lock
###############################################################################

if (
  as.character(packageVersion("Seurat")) !=
  "5.5.1"
) {
  
  stop(
    "This script was validated for Seurat 5.5.1. Found ",
    as.character(packageVersion("Seurat")),
    ".",
    call. = FALSE
  )
}


if (
  as.character(packageVersion("SeuratObject")) !=
  "5.4.0"
) {
  
  stop(
    "This script was validated for SeuratObject 5.4.0. Found ",
    as.character(packageVersion("SeuratObject")),
    ".",
    call. = FALSE
  )
}


###############################################################################
# File paths
###############################################################################

full_data_path <- file.path(
  "results",
  "rds_objects",
  "seurat_qc_processed_expression.rds"
)


atlas_path <- file.path(
  "results",
  "rds_objects",
  "phase3_atlas_sketch_final_annotation.rds"
)


output_dir <- file.path(
  "results",
  "rds_objects"
)


results_dir <- "results"


if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


if (!dir.exists(results_dir)) {
  
  dir.create(
    results_dir,
    recursive = TRUE
  )
}


###############################################################################
# Locked constants
###############################################################################

n_full_cells_expected <- 106592L

n_reference_cells <- 20000L

n_genes_expected <- 24452L

n_pcs <- 50L

# The atlas graph was constructed using PCs 1:20.
# Label transfer therefore uses the same 20-dimensional space.
selected_pcs <- 1:20

knn_k <- 30L

knn_batch_size <- 5000L

low_agreement_threshold <- 0.80

low_margin_threshold <- 0.20

very_high_agreement_threshold <- 0.95


###############################################################################
# Load objects
###############################################################################

if (!file.exists(full_data_path)) {
  
  stop(
    "Full Phase 2 object not found:\n",
    full_data_path,
    call. = FALSE
  )
}


if (!file.exists(atlas_path)) {
  
  stop(
    "Phase 3 atlas object not found:\n",
    atlas_path,
    call. = FALSE
  )
}


cat("\nLoading Phase 2 full dataset...\n")

seurat_obj <- readRDS(
  full_data_path
)


cat("Loading Phase 3 atlas...\n")

atlas_sketch <- readRDS(
  atlas_path
)


###############################################################################
# Basic object validation
###############################################################################

if (!inherits(seurat_obj, "Seurat")) {
  
  stop(
    "Full dataset is not a Seurat object.",
    call. = FALSE
  )
}


if (!inherits(atlas_sketch, "Seurat")) {
  
  stop(
    "Atlas is not a Seurat object.",
    call. = FALSE
  )
}


###############################################################################
# Validate RNA assays
###############################################################################

if (!"RNA" %in% names(seurat_obj@assays)) {
  
  stop(
    "RNA assay not found in full Phase 2 object.",
    call. = FALSE
  )
}


if (!"RNA" %in% names(atlas_sketch@assays)) {
  
  stop(
    "RNA assay not found in atlas object.",
    call. = FALSE
  )
}


DefaultAssay(seurat_obj) <- "RNA"

DefaultAssay(atlas_sketch) <- "RNA"


###############################################################################
# Full dataset dimensions
###############################################################################

full_cell_count <- ncol(seurat_obj)

full_gene_count <- nrow(
  seurat_obj[["RNA"]]
)


cat("\n============================================================\n")
cat("PHASE 3 — FULL DATASET PROJECTION\n")
cat("============================================================\n")

cat(
  "Full cells: ",
  full_cell_count,
  "\n",
  sep = ""
)

cat(
  "Full genes: ",
  full_gene_count,
  "\n",
  sep = ""
)

cat(
  "Reference cells: ",
  ncol(atlas_sketch),
  "\n",
  sep = ""
)

cat(
  "Reference genes: ",
  nrow(atlas_sketch[["RNA"]]),
  "\n",
  sep = ""
)


###############################################################################
# Dimension checks
###############################################################################

if (
  full_cell_count !=
  n_full_cells_expected
) {
  
  stop(
    "Unexpected full-dataset cell count. Expected ",
    n_full_cells_expected,
    " but found ",
    full_cell_count,
    ".",
    call. = FALSE
  )
}


if (
  full_gene_count !=
  n_genes_expected
) {
  
  stop(
    "Unexpected full-dataset gene count. Expected ",
    n_genes_expected,
    " but found ",
    full_gene_count,
    ".",
    call. = FALSE
  )
}


if (
  ncol(atlas_sketch) !=
  n_reference_cells
) {
  
  stop(
    "Unexpected atlas size. Expected ",
    n_reference_cells,
    " but found ",
    ncol(atlas_sketch),
    ".",
    call. = FALSE
  )
}


if (
  nrow(atlas_sketch[["RNA"]]) !=
  n_genes_expected
) {
  
  stop(
    "Unexpected atlas gene count.",
    call. = FALSE
  )
}


###############################################################################
# Gene identity
###############################################################################

full_genes <- rownames(
  seurat_obj[["RNA"]]
)


atlas_genes <- rownames(
  atlas_sketch[["RNA"]]
)


if (!identical(full_genes, atlas_genes)) {
  
  common_genes <- intersect(
    full_genes,
    atlas_genes
  )
  
  if (length(common_genes) != n_genes_expected) {
    
    stop(
      "Full dataset and atlas do not contain the same gene set.",
      call. = FALSE
    )
  }
  
  stop(
    "Full dataset and atlas contain the same genes but in different order.",
    call. = FALSE
  )
}


cat(
  "Gene identity/order: PASS\n"
)


###############################################################################
# Metadata validation
###############################################################################

required_metadata <- c(
  "Phase",
  "Donor",
  "GSM"
)


missing_metadata <- setdiff(
  required_metadata,
  colnames(seurat_obj@meta.data)
)


if (length(missing_metadata) > 0L) {
  
  stop(
    "Required metadata missing: ",
    paste(
      missing_metadata,
      collapse = ", "
    ),
    call. = FALSE
  )
}


cat(
  "Clinical states: ",
  paste(
    sort(
      unique(
        as.character(
          seurat_obj$Phase
        )
      )
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "Donors: ",
  length(
    unique(
      as.character(
        seurat_obj$Donor
      )
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Samples: ",
  length(
    unique(
      as.character(
        seurat_obj$GSM
      )
    )
  ),
  "\n",
  sep = ""
)


###############################################################################
# SECTION 25A — VALIDATE PROCESSED-EXPRESSION LAYERS
###############################################################################

cat("\n============================================================\n")
cat("SECTION 25A — PROCESSED-EXPRESSION LAYERS\n")
cat("============================================================\n")


full_rna_layers <- Layers(
  seurat_obj[["RNA"]]
)


full_data_layers <- grep(
  "^data\\.",
  full_rna_layers,
  value = TRUE
)


full_counts_layers <- grep(
  "^counts\\.",
  full_rna_layers,
  value = TRUE
)


cat(
  "Detected data layers: ",
  length(full_data_layers),
  "\n",
  sep = ""
)


cat(
  "Detected counts layers: ",
  length(full_counts_layers),
  "\n",
  sep = ""
)


if (length(full_data_layers) == 0L) {
  
  stop(
    "No sample-specific RNA data layers found.",
    call. = FALSE
  )
}


if (length(full_data_layers) != 23L) {
  
  stop(
    "Expected 23 sample-specific data layers but found ",
    length(full_data_layers),
    ".",
    call. = FALSE
  )
}


if (length(full_counts_layers) > 0L) {
  
  warning(
    "Counts layers are present but will NOT be used."
  )
}


cat("\nData layers:\n")

print(
  full_data_layers
)


###############################################################################
# Validate each data layer
###############################################################################

layer_cell_counts <- integer(
  length(full_data_layers)
)

names(layer_cell_counts) <-
  full_data_layers


layer_cell_ids <- vector(
  mode = "list",
  length = length(full_data_layers)
)


for (i in seq_along(full_data_layers)) {
  
  layer_name <- full_data_layers[i]
  
  cat(
    "Validating ",
    layer_name,
    "...\n",
    sep = ""
  )
  
  layer_matrix <- LayerData(
    object = seurat_obj[["RNA"]],
    layer = layer_name
  )
  
  layer_cell_counts[i] <- ncol(
    layer_matrix
  )
  
  layer_cell_ids[[i]] <- colnames(
    layer_matrix
  )
  
  if (
    nrow(layer_matrix) !=
    n_genes_expected
  ) {
    
    stop(
      "Layer ",
      layer_name,
      " contains ",
      nrow(layer_matrix),
      " genes; expected ",
      n_genes_expected,
      ".",
      call. = FALSE
    )
  }
  
  if (
    !identical(
      rownames(layer_matrix),
      full_genes
    )
  ) {
    
    stop(
      "Gene ordering differs in layer ",
      layer_name,
      ".",
      call. = FALSE
    )
  }
  
  if (
    anyNA(layer_matrix)
  ) {
    
    stop(
      "NA values detected in layer ",
      layer_name,
      ".",
      call. = FALSE
    )
  }
  
  if (
    any(!is.finite(layer_matrix))
  ) {
    
    stop(
      "Non-finite values detected in layer ",
      layer_name,
      ".",
      call. = FALSE
    )
  }
  
  rm(
    layer_matrix
  )
}


###############################################################################
# Validate total layer cell count
###############################################################################

if (
  sum(layer_cell_counts) !=
  n_full_cells_expected
) {
  
  stop(
    "Sum of sample-specific layer cell counts is ",
    sum(layer_cell_counts),
    " but expected ",
    n_full_cells_expected,
    ".",
    call. = FALSE
  )
}


###############################################################################
# Validate unique cell IDs
###############################################################################

all_layer_cell_ids <- unlist(
  layer_cell_ids,
  use.names = FALSE
)


if (
  length(
    unique(
      all_layer_cell_ids
    )
  ) !=
  n_full_cells_expected
) {
  
  stop(
    "Cell IDs are not unique across the 23 data layers.",
    call. = FALSE
  )
}


if (
  !identical(
    sort(all_layer_cell_ids),
    sort(colnames(seurat_obj))
  )
) {
  
  stop(
    "The cells represented by the 23 data layers do not exactly match ",
    "the cells in the Phase 2 Seurat object.",
    call. = FALSE
  )
}


cat(
  "\nAll 23 processed-expression layers validated.\n"
)


cat(
  "Total cells represented: ",
  sum(layer_cell_counts),
  "\n",
  sep = ""
)


###############################################################################
# Representation lock
###############################################################################

cat("\nRepresentation locked:\n")

cat(
  "Input values = existing GEO processed log-CP10K expression.\n"
)

cat(
  "NormalizeData(): NOT RUN\n"
)

cat(
  "Raw-count reconstruction: NOT PERFORMED\n"
)

cat(
  "JoinLayers(): NOT USED\n"
)

cat(
  "Batch correction: NOT PERFORMED\n"
)

cat(
  "Cell filtering: NOT PERFORMED\n"
)


###############################################################################
# SECTION 25B — VALIDATE ATLAS PCA
###############################################################################

cat("\n============================================================\n")
cat("SECTION 25B — ATLAS PCA VALIDATION\n")
cat("============================================================\n")


if (
  !"pca" %in%
  Reductions(atlas_sketch)
) {
  
  stop(
    "Atlas does not contain a PCA reduction.",
    call. = FALSE
  )
}


atlas_pca <- Embeddings(
  atlas_sketch,
  reduction = "pca"
)


atlas_loadings <- Loadings(
  atlas_sketch,
  reduction = "pca"
)


if (
  ncol(atlas_pca) <
  n_pcs
) {
  
  stop(
    "Atlas PCA contains fewer than ",
    n_pcs,
    " dimensions.",
    call. = FALSE
  )
}


if (
  ncol(atlas_loadings) <
  n_pcs
) {
  
  stop(
    "Atlas PCA loadings contain fewer than ",
    n_pcs,
    " dimensions.",
    call. = FALSE
  )
}


###############################################################################
# Use the genes actually represented in the PCA loadings
###############################################################################

pca_features <- rownames(
  atlas_loadings
)


if (
  length(pca_features) == 0L
) {
  
  stop(
    "No PCA loading features found.",
    call. = FALSE
  )
}


missing_pca_features <- setdiff(
  pca_features,
  full_genes
)


if (
  length(missing_pca_features) > 0L
) {
  
  stop(
    "Atlas PCA contains features absent from full dataset: ",
    paste(
      head(
        missing_pca_features,
        20
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}


cat(
  "Atlas PCA dimensions: ",
  ncol(atlas_pca),
  "\n",
  sep = ""
)


cat(
  "Atlas PCA loading features: ",
  length(pca_features),
  "\n",
  sep = ""
)


cat(
  "Atlas variable features: ",
  length(
    VariableFeatures(atlas_sketch)
  ),
  "\n",
  sep = ""
)


cat(
  "PCA loading features retained: ",
  length(pca_features),
  "\n"
)


if (
  length(pca_features) !=
  nrow(atlas_loadings)
) {
  
  stop(
    "PCA feature count does not match PCA loading rows.",
    call. = FALSE
  )
}


cat(
  "PCA validation: PASS\n"
)


###############################################################################
# SECTION 25C — CALCULATE REFERENCE SCALING PARAMETERS
###############################################################################
# The reference PCA was learned from the atlas' processed RNA data after
# reference-based centering/scaling. We reproduce Seurat's projection logic
# directly here because ProjectCellEmbeddings() dispatch on a matrix/Assay5
# object is version-sensitive in the installed Seurat release.
#
# The projection is:
#   1. obtain reference means and SDs for the PCA features
#   2. scale each query layer with those reference statistics
#   3. multiply the scaled query matrix by the reference PCA loadings
#
# No normalization is performed.
###############################################################################

cat("\n============================================================\n")
cat("SECTION 25C — REFERENCE PCA SCALING PARAMETERS\n")
cat("============================================================\n")

atlas_data_layers <- Layers(
  atlas_sketch[["RNA"]],
  search = "^data$"
)

if (length(atlas_data_layers) != 1L) {
  stop(
    "Expected exactly one atlas RNA data layer named 'data'. Found: ",
    paste(atlas_data_layers, collapse = ", "),
    call. = FALSE
  )
}

atlas_data <- LayerData(
  object = atlas_sketch[["RNA"]],
  layer = atlas_data_layers[1]
)
atlas_data <- atlas_data[pca_features, , drop = FALSE]

if (!inherits(atlas_data, "dgCMatrix")) {
  atlas_data <- as(atlas_data, "dgCMatrix")
}

reference_feature_mean <- Seurat:::RowMeanSparse(atlas_data)
reference_feature_sd <- sqrt(Seurat:::RowVarSparse(atlas_data))

names(reference_feature_mean) <- pca_features
names(reference_feature_sd) <- pca_features

reference_feature_mean[!is.finite(reference_feature_mean)] <- 1
reference_feature_sd[!is.finite(reference_feature_sd)] <- 1
reference_feature_sd[reference_feature_sd == 0] <- 1

if (length(reference_feature_mean) != length(pca_features) ||
    length(reference_feature_sd) != length(pca_features)) {
  stop("Reference scaling parameters have incorrect lengths.", call. = FALSE)
}

if (any(!is.finite(reference_feature_mean)) ||
    any(!is.finite(reference_feature_sd))) {
  stop("Non-finite reference scaling parameters remain.", call. = FALSE)
}

cat("Reference scaling features: ", length(pca_features), "\n", sep = "")
cat("Reference scaling parameters: PASS\n")

gc()

###############################################################################
# SECTION 25D — VALIDATE THE PROJECTION MATHEMATICS ON THE ATLAS ITSELF
###############################################################################
# This is the critical safeguard. Before projecting 106,592 query cells, we
# project the 20,000 reference cells through the same algebra and compare the
# reconstructed coordinates with the stored atlas PCA embeddings.
###############################################################################

cat("\n============================================================\n")
cat("SECTION 25D — ATLAS SELF-PROJECTION VALIDATION\n")
cat("============================================================\n")

reference_pca_loadings <- atlas_loadings[
  pca_features,
  seq_len(n_pcs),
  drop = FALSE
]

if (!identical(rownames(reference_pca_loadings), pca_features)) {
  stop("Reference PCA loading feature order could not be established.", call. = FALSE)
}

scale_reference_matrix <- function(mat) {
  mat <- mat[pca_features, , drop = FALSE]
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  Seurat:::FastSparseRowScaleWithKnownStats(
    mat = mat,
    mu = reference_feature_mean,
    sigma = reference_feature_sd,
    display_progress = FALSE
  )
}

atlas_scaled <- scale_reference_matrix(atlas_data)

atlas_projected_check <- t(reference_pca_loadings) %*% atlas_scaled
atlas_projected_check <- as.matrix(atlas_projected_check)
atlas_projected_check <- t(atlas_projected_check)

rownames(atlas_projected_check) <- rownames(atlas_pca)
colnames(atlas_projected_check) <- colnames(atlas_pca)[seq_len(n_pcs)]

stored_atlas_pca <- atlas_pca[, seq_len(n_pcs), drop = FALSE]

# PCA signs are mathematically arbitrary. Compare both the direct correlation
# and the sign-aligned correlation for every PC.
pc_correlations <- vapply(
  seq_len(n_pcs),
  function(j) {
    cor(
      atlas_projected_check[, j],
      stored_atlas_pca[, j],
      use = "pairwise.complete.obs"
    )
  },
  numeric(1)
)

pc_signs <- ifelse(pc_correlations < 0, -1, 1)
atlas_projected_aligned <- sweep(
  atlas_projected_check,
  2,
  pc_signs,
  "*"
)

pc_correlations_aligned <- vapply(
  seq_len(n_pcs),
  function(j) {
    cor(
      atlas_projected_aligned[, j],
      stored_atlas_pca[, j],
      use = "pairwise.complete.obs"
    )
  },
  numeric(1)
)

pc_rmse <- vapply(
  seq_len(n_pcs),
  function(j) {
    sqrt(mean((
      atlas_projected_aligned[, j] - stored_atlas_pca[, j]
    )^2))
  },
  numeric(1)
)

atlas_projection_validation <- data.frame(
  PC = seq_len(n_pcs),
  Correlation = pc_correlations,
  Sign_Aligned_Correlation = pc_correlations_aligned,
  RMSE = pc_rmse,
  Sign_Used = pc_signs
)

write.csv(
  atlas_projection_validation,
  file = file.path(results_dir, "phase3_atlas_self_projection_validation.csv"),
  row.names = FALSE
)

cat("\nAtlas self-projection validation:\n")
print(atlas_projection_validation)

# Require essentially exact recovery. If this fails, STOP rather than silently
# transferring labels in a coordinate system that has not been validated.
if (any(!is.finite(pc_correlations_aligned)) ||
    any(pc_correlations_aligned < 0.999) ||
    any(!is.finite(pc_rmse))) {
  stop(
    "Atlas self-projection validation failed. The stored PCA coordinates " ,
    "cannot be reproduced reliably from the stored RNA data/scaling statistics. " ,
    "No full-dataset label transfer will be performed.",
    call. = FALSE
  )
}

cat("\nAtlas self-projection: PASS\n")
cat("Minimum sign-aligned PC correlation: ",
    sprintf("%.6f", min(pc_correlations_aligned)), "\n", sep = "")

rm(atlas_scaled, atlas_projected_check, atlas_projected_aligned)
gc()

###############################################################################
# SECTION 26 — VALIDATE FINAL ATLAS ANNOTATION
###############################################################################

###############################################################################

cat("\n============================================================\n")
cat("SECTION 26 — ATLAS ANNOTATION VALIDATION\n")
cat("============================================================\n")


if (
  !"Final_Cell_Type" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "Final_Cell_Type metadata is missing from atlas.",
    call. = FALSE
  )
}


if (
  !"Atlas_Cluster" %in%
  colnames(atlas_sketch@meta.data)
) {
  
  stop(
    "Atlas_Cluster metadata is missing from atlas.",
    call. = FALSE
  )
}


atlas_labels <- as.character(
  atlas_sketch$Final_Cell_Type
)


atlas_cluster <- as.character(
  atlas_sketch$Atlas_Cluster
)


if (
  anyNA(atlas_labels)
) {
  
  stop(
    "Atlas contains NA Final_Cell_Type labels.",
    call. = FALSE
  )
}


if (
  any(atlas_labels == "")
) {
  
  stop(
    "Atlas contains blank Final_Cell_Type labels.",
    call. = FALSE
  )
}


observed_clusters <- sort(
  unique(
    atlas_cluster
  )
)


n_atlas_clusters <- length(
  observed_clusters
)


cat(
  "Atlas clusters: ",
  n_atlas_clusters,
  "\n",
  sep = ""
)


cat(
  "Clusters: ",
  paste(
    observed_clusters,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


if (
  n_atlas_clusters !=
  20L
) {
  
  stop(
    "Expected 20 atlas clusters but found ",
    n_atlas_clusters,
    ".",
    call. = FALSE
  )
}


###############################################################################
# One final label per cluster
###############################################################################

cluster_label_check <- atlas_sketch@meta.data %>%
  
  dplyr::mutate(
    Atlas_Cluster = as.character(
      Atlas_Cluster
    ),
    Final_Cell_Type = as.character(
      Final_Cell_Type
    )
  ) %>%
  
  dplyr::group_by(
    Atlas_Cluster
  ) %>%
  
  dplyr::summarise(
    n_labels = dplyr::n_distinct(
      Final_Cell_Type
    ),
    labels = paste(
      sort(
        unique(
          Final_Cell_Type
        )
      ),
      collapse = "; "
    ),
    .groups = "drop"
  )


if (
  any(
    cluster_label_check$n_labels !=
    1L
  )
) {
  
  print(
    cluster_label_check
  )
  
  stop(
    "At least one atlas cluster has multiple final labels.",
    call. = FALSE
  )
}


###############################################################################
# Cluster → final label lookup
###############################################################################

cluster_label_lookup <- atlas_sketch@meta.data %>%
  
  dplyr::mutate(
    Atlas_Cluster = as.character(
      Atlas_Cluster
    ),
    Final_Cell_Type = as.character(
      Final_Cell_Type
    )
  ) %>%
  
  dplyr::select(
    Atlas_Cluster,
    Final_Cell_Type
  ) %>%
  
  dplyr::distinct()


cat("\nFinal atlas labels:\n")

print(
  cluster_label_lookup
)


###############################################################################
# Expected label vocabulary
###############################################################################

expected_labels <- c(
  "CD8_T",
  "NK",
  "CD4_T",
  "MAIT",
  "B_Cell",
  "Unresolved_Lymphoid",
  "Treg",
  "Monocyte_Myeloid",
  "Dendritic_Cell",
  "Plasma_Cell",
  "pDC"
)


unexpected_labels <- setdiff(
  unique(
    atlas_labels
  ),
  expected_labels
)


if (
  length(unexpected_labels) > 0L
) {
  
  stop(
    "Unexpected Final_Cell_Type labels detected: ",
    paste(
      unexpected_labels,
      collapse = ", "
    ),
    call. = FALSE
  )
}


###############################################################################
# SECTION 27 — SAMPLE-WISE PCA PROJECTION
###############################################################################

cat("\n============================================================\n")
cat("SECTION 27 — SAMPLE-WISE PCA PROJECTION\n")
cat("============================================================\n")

cat("Projecting all 23 sample-specific data layers into the validated atlas PCA.\n")

reference_pca_loadings <- atlas_loadings[
  pca_features,
  seq_len(n_pcs),
  drop = FALSE
]

project_single_layer <- function(layer_matrix) {
  query_matrix <- layer_matrix[pca_features, , drop = FALSE]
  if (!inherits(query_matrix, "dgCMatrix")) {
    query_matrix <- as(query_matrix, "dgCMatrix")
  }
  
  query_scaled <- Seurat:::FastSparseRowScaleWithKnownStats(
    mat = query_matrix,
    mu = reference_feature_mean,
    sigma = reference_feature_sd,
    display_progress = FALSE
  )
  
  projected <- t(reference_pca_loadings) %*% query_scaled
  projected <- as.matrix(t(projected))
  
  rownames(projected) <- colnames(layer_matrix)
  colnames(projected) <- colnames(atlas_pca)[seq_len(n_pcs)]
  
  if (nrow(projected) != ncol(layer_matrix)) {
    stop(
      "Projected matrix has ", nrow(projected),
      " rows but query contains ", ncol(layer_matrix), " cells.",
      call. = FALSE
    )
  }
  
  if (ncol(projected) != n_pcs) {
    stop("Projected matrix has the wrong number of PCA dimensions.", call. = FALSE)
  }
  
  if (any(!is.finite(projected))) {
    stop("Non-finite projected PCA coordinates detected.", call. = FALSE)
  }
  
  projected
}

projected_embeddings_list <- vector("list", length(full_data_layers))
names(projected_embeddings_list) <- full_data_layers

for (i in seq_along(full_data_layers)) {
  layer_name <- full_data_layers[i]
  
  cat("\n------------------------------------------------------------\n")
  cat("Projecting layer ", i, " / ", length(full_data_layers), ": ",
      layer_name, "\n", sep = "")
  
  layer_matrix <- LayerData(
    object = seurat_obj[["RNA"]],
    layer = layer_name
  )
  
  projected_embeddings_list[[i]] <- project_single_layer(layer_matrix)
  
  cat("Cells projected: ", nrow(projected_embeddings_list[[i]]), "\n", sep = "")
  cat("Dimensions: ", ncol(projected_embeddings_list[[i]]), "\n", sep = "")
  cat("Projection status: PASS\n")
  
  rm(layer_matrix)
  gc()
}

for (i in seq_along(projected_embeddings_list)) {
  projection <- projected_embeddings_list[[i]]
  if (is.null(projection)) {
    stop("Projection for layer ", full_data_layers[i], " is NULL.", call. = FALSE)
  }
  if (nrow(projection) != layer_cell_counts[i]) {
    stop("Projection cell count mismatch for ", full_data_layers[i], ".", call. = FALSE)
  }
  if (ncol(projection) != n_pcs) {
    stop("Projection dimension mismatch for ", full_data_layers[i], ".", call. = FALSE)
  }
}

cat("\nCombining projected PCA coordinates...\n")
query_pca_full <- do.call(rbind, projected_embeddings_list)

if (!is.matrix(query_pca_full) ||
    nrow(query_pca_full) != n_full_cells_expected ||
    ncol(query_pca_full) != n_pcs) {
  stop("Combined projected PCA has incorrect dimensions.", call. = FALSE)
}

if (!identical(sort(rownames(query_pca_full)), sort(colnames(seurat_obj)))) {
  stop("Projected PCA cell identities do not exactly match the Phase 2 object.", call. = FALSE)
}

query_pca_full <- query_pca_full[colnames(seurat_obj), , drop = FALSE]

if (!identical(rownames(query_pca_full), colnames(seurat_obj))) {
  stop("Projected PCA could not be aligned to Seurat cell order.", call. = FALSE)
}

if (any(!is.finite(query_pca_full))) {
  stop("Combined projected PCA contains non-finite values.", call. = FALSE)
}

cat("\nFull-dataset PCA projection: PASS\n")
cat("Cells: ", nrow(query_pca_full), "\n", sep = "")
cat("Dimensions: ", ncol(query_pca_full), "\n", sep = "")

if ("pca.full" %in% Reductions(seurat_obj)) {
  seurat_obj[["pca.full"]] <- NULL
}

seurat_obj[["pca.full"]] <- CreateDimReducObject(
  embeddings = query_pca_full,
  loadings = reference_pca_loadings,
  key = "PCFull_",
  assay = "RNA"
)

stored_pca_full <- Embeddings(seurat_obj, reduction = "pca.full")
if (!identical(rownames(stored_pca_full), colnames(seurat_obj)) ||
    ncol(stored_pca_full) != n_pcs) {
  stop("Stored pca.full failed validation.", call. = FALSE)
}

cat("pca.full stored in Seurat object: PASS\n")

gc()

###############################################################################
# SECTION 28 — KNN LABEL TRANSFER + CONFIDENCE
###############################################################################

###############################################################################

cat("\n============================================================\n")
cat("SECTION 28 — KNN LABEL TRANSFER\n")
cat("============================================================\n")


###############################################################################
# Reference PCA
#
# PCs 1:20 match the dimensions used to construct the atlas graph.
###############################################################################

reference_pca_for_transfer <- atlas_pca[
  ,
  selected_pcs,
  drop = FALSE
]


query_pca_for_transfer <- query_pca_full[
  ,
  selected_pcs,
  drop = FALSE
]


###############################################################################
# Reference labels
###############################################################################

reference_clusters <- as.character(
  atlas_sketch$Atlas_Cluster
)


reference_labels <- cluster_label_lookup$Final_Cell_Type[
  match(
    reference_clusters,
    cluster_label_lookup$Atlas_Cluster
  )
]


if (
  anyNA(reference_labels)
) {
  
  stop(
    "Some atlas cells could not be assigned Final_Cell_Type.",
    call. = FALSE
  )
}


###############################################################################
# Result vectors
###############################################################################

n_query <- nrow(
  query_pca_for_transfer
)


transferred_cluster <- rep(
  NA_character_,
  n_query
)


transferred_label <- rep(
  NA_character_,
  n_query
)


neighbour_agreement <- rep(
  NA_real_,
  n_query
)


second_best_agreement <- rep(
  NA_real_,
  n_query
)


agreement_margin <- rep(
  NA_real_,
  n_query
)


###############################################################################
# Batch-wise KNN
###############################################################################

batch_starts <- seq(
  1L,
  n_query,
  by = knn_batch_size
)


cat(
  "KNN batches: ",
  length(batch_starts),
  "\n",
  sep = ""
)


for (
  batch_start in batch_starts
) {
  
  batch_end <- min(
    batch_start +
      knn_batch_size -
      1L,
    n_query
  )
  
  
  batch_indices <- batch_start:batch_end
  
  
  query_batch <- query_pca_for_transfer[
    batch_indices,
    ,
    drop = FALSE
  ]
  
  
  nn_result <- RANN::nn2(
    data = reference_pca_for_transfer,
    query = query_batch,
    k = knn_k
  )
  
  
  neighbour_indices <- nn_result$nn.idx
  
  
  for (
    i in seq_along(batch_indices)
  ) {
    
    cell_index <- batch_indices[i]
    
    
    neighbour_reference_indices <-
      neighbour_indices[i, ]
    
    
    neighbour_clusters <-
      reference_clusters[
        neighbour_reference_indices
      ]
    
    
    neighbour_labels <-
      reference_labels[
        neighbour_reference_indices
      ]
    
    
    ###########################################################################
    # Cluster vote
    ###########################################################################
    
    cluster_counts <- table(
      neighbour_clusters
    )
    
    
    cluster_counts <- sort(
      cluster_counts,
      decreasing = TRUE
    )
    
    
    transferred_cluster[cell_index] <-
      names(
        cluster_counts
      )[1]
    
    
    ###########################################################################
    # Cell-type vote
    ###########################################################################
    
    label_counts <- table(
      neighbour_labels
    )
    
    
    label_counts <- sort(
      label_counts,
      decreasing = TRUE
    )
    
    
    transferred_label[cell_index] <-
      names(
        label_counts
      )[1]
    
    
    ###########################################################################
    # Agreement
    ###########################################################################
    
    top_count <- as.numeric(
      label_counts[1]
    )
    
    
    neighbour_agreement[cell_index] <-
      top_count /
      knn_k
    
    
    if (
      length(label_counts) >=
      2L
    ) {
      
      second_best_agreement[cell_index] <-
        as.numeric(
          label_counts[2]
        ) /
        knn_k
      
    } else {
      
      second_best_agreement[cell_index] <-
        0
    }
    
    
    agreement_margin[cell_index] <-
      neighbour_agreement[cell_index] -
      second_best_agreement[cell_index]
  }
  
  
  cat(
    "Completed cells ",
    batch_start,
    "–",
    batch_end,
    " of ",
    n_query,
    "\n",
    sep = ""
  )
}


###############################################################################
# Validate transfer
###############################################################################

if (
  anyNA(
    transferred_label
  )
) {
  
  stop(
    "Some cells did not receive a transferred label.",
    call. = FALSE
  )
}


if (
  anyNA(
    neighbour_agreement
  )
) {
  
  stop(
    "Some cells did not receive an agreement score.",
    call. = FALSE
  )
}


if (
  any(
    !is.finite(
      neighbour_agreement
    )
  )
) {
  
  stop(
    "Non-finite neighbour agreement detected.",
    call. = FALSE
  )
}


###############################################################################
# Confidence classification
#
# Low:
#   agreement < 0.80
#
# Intermediate:
#   agreement >= 0.80 but margin < 0.20
#
# High:
#   agreement >= 0.80 and margin >= 0.20
###############################################################################

transfer_confidence <- dplyr::case_when(
  
  neighbour_agreement <
    low_agreement_threshold ~
    "Low",
  
  neighbour_agreement >=
    low_agreement_threshold &
    agreement_margin <
    low_margin_threshold ~
    "Intermediate",
  
  TRUE ~
    "High"
)


###############################################################################
# Store transfer metadata
###############################################################################

seurat_obj$Transferred_Atlas_Cluster <-
  transferred_cluster


seurat_obj$Final_Cell_Type <-
  transferred_label


seurat_obj$Transfer_Agreement <-
  neighbour_agreement


seurat_obj$Transfer_SecondBest_Agreement <-
  second_best_agreement


seurat_obj$Transfer_Agreement_Margin <-
  agreement_margin


seurat_obj$Transfer_Confidence <-
  transfer_confidence


# QC flag only.
# It is NOT a filtering criterion.
seurat_obj$Low_Confidence <-
  seurat_obj$Transfer_Confidence !=
  "High"


###############################################################################
# Transfer summary
###############################################################################

transfer_summary <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Final_Cell_Type,
    Transfer_Confidence,
    name = "Cells"
  ) %>%
  
  dplyr::arrange(
    Final_Cell_Type,
    Transfer_Confidence
  )


cat("\nTransfer summary:\n")

print(
  transfer_summary
)


###############################################################################
# Overall confidence
###############################################################################

confidence_summary <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Transfer_Confidence,
    name = "Cells"
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  )


cat("\nOverall transfer confidence:\n")

print(
  confidence_summary
)


###############################################################################
# Very high agreement
###############################################################################

very_high_agreement_cells <- sum(
  seurat_obj$Transfer_Agreement >=
    very_high_agreement_threshold
)


cat(
  "\nCells with >= ",
  very_high_agreement_threshold * 100,
  "% neighbour agreement: ",
  very_high_agreement_cells,
  "\n",
  sep = ""
)


###############################################################################
# SECTION 29 — TARGETED BIOLOGICAL QC
###############################################################################

cat("\n============================================================\n")
cat("SECTION 29 — TARGETED BIOLOGICAL QC\n")
cat("============================================================\n")


###############################################################################
# Dynamic target clusters
###############################################################################

cd8_clusters <- cluster_label_lookup$Atlas_Cluster[
  cluster_label_lookup$Final_Cell_Type ==
    "CD8_T"
]


monocyte_myeloid_clusters <-
  cluster_label_lookup$Atlas_Cluster[
    cluster_label_lookup$Final_Cell_Type ==
      "Monocyte_Myeloid"
  ]


cat(
  "CD8_T atlas clusters: ",
  paste(
    cd8_clusters,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "Monocyte_Myeloid atlas clusters: ",
  paste(
    monocyte_myeloid_clusters,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


if (
  length(cd8_clusters) == 0L
) {
  
  warning(
    "No CD8_T atlas clusters identified."
  )
}


if (
  length(monocyte_myeloid_clusters) == 0L
) {
  
  warning(
    "No Monocyte_Myeloid atlas clusters identified."
  )
}


###############################################################################
# Diagnostic marker sets
###############################################################################

cd8_markers <- c(
  "CD3D",
  "CD3E",
  "CD3G",
  "TRBC1",
  "TRBC2",
  "CD8A",
  "CD8B",
  "CCL5",
  "CST7",
  "GZMK",
  "GZMA",
  "NKG7",
  "CCL4",
  "CCL4L2"
)


monocyte_myeloid_markers <- c(
  "LYZ",
  "S100A8",
  "S100A9",
  "FCN1",
  "CD14",
  "VCAN",
  "TREM1",
  "FCAR",
  "FPR2",
  "CTSS",
  "CTSD",
  "LGALS3",
  "OLR1",
  "CST3"
)


macrophage_associated_markers <- c(
  "C1QA",
  "C1QB",
  "C1QC",
  "APOE",
  "TREM2",
  "CD68",
  "MSR1",
  "MARCO",
  "LPL"
)


###############################################################################
# Available markers
###############################################################################

available_cd8_markers <- intersect(
  cd8_markers,
  full_genes
)


available_monocyte_markers <- intersect(
  monocyte_myeloid_markers,
  full_genes
)


available_macrophage_markers <- intersect(
  macrophage_associated_markers,
  full_genes
)


diagnostic_markers <- unique(
  c(
    available_cd8_markers,
    available_monocyte_markers,
    available_macrophage_markers
  )
)


###############################################################################
# Diagnostic expression summary
#
# IMPORTANT:
#   This is computed layer-by-layer so no JoinLayers() is required.
#   The expression values remain the original processed log-CP10K values.
###############################################################################

diagnostic_marker_rows <- list()
row_counter <- 1L

if (length(diagnostic_markers) > 0L) {
  for (layer_name in full_data_layers) {
    layer_matrix <- LayerData(
      object = seurat_obj[["RNA"]],
      layer = layer_name
    )
    present <- intersect(diagnostic_markers, rownames(layer_matrix))
    if (length(present) == 0L) next
    
    cells <- colnames(layer_matrix)
    meta_layer <- seurat_obj@meta.data[cells, , drop = FALSE]
    expr <- layer_matrix[present, , drop = FALSE]
    
    for (label in unique(as.character(meta_layer$Final_Cell_Type))) {
      idx <- which(as.character(meta_layer$Final_Cell_Type) == label)
      if (length(idx) == 0L) next
      vals <- Matrix::rowMeans(expr[, idx, drop = FALSE])
      row <- as.data.frame(as.list(vals), check.names = FALSE)
      row$Final_Cell_Type <- label
      row$Layer <- layer_name
      diagnostic_marker_rows[[row_counter]] <- row
      row_counter <- row_counter + 1L
    }
    rm(layer_matrix, expr, meta_layer)
    gc()
  }
  
  diagnostic_layer_means <- dplyr::bind_rows(diagnostic_marker_rows)
  
  diagnostic_means <- diagnostic_layer_means %>%
    dplyr::group_by(Final_Cell_Type) %>%
    dplyr::summarise(
      dplyr::across(all_of(intersect(diagnostic_markers, colnames(diagnostic_layer_means))),
                    ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  write.csv(
    diagnostic_means,
    file = file.path(results_dir, "phase3_full_dataset_diagnostic_marker_means.csv"),
    row.names = FALSE
  )
  
  cat("\nDiagnostic marker means saved.\n")
}

###############################################################################
# CD8 transfer QC
###############################################################################

cd8_transfer_qc <- seurat_obj@meta.data %>%
  
  dplyr::filter(
    Final_Cell_Type ==
      "CD8_T"
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Mean_Agreement =
      mean(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Median_Agreement =
      median(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Fraction_High_Confidence =
      mean(
        Transfer_Confidence ==
          "High"
      ),
    
    Fraction_Low_Agreement =
      mean(
        Transfer_Agreement <
          low_agreement_threshold
      ),
    
    Fraction_High_Agreement =
      mean(
        Transfer_Agreement >=
          very_high_agreement_threshold
      )
  )


cat(
  "\nCD8_T transfer QC:\n"
)

print(
  cd8_transfer_qc
)


###############################################################################
# Monocyte/Myeloid transfer QC
###############################################################################

monocyte_transfer_qc <- seurat_obj@meta.data %>%
  
  dplyr::filter(
    Final_Cell_Type ==
      "Monocyte_Myeloid"
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Mean_Agreement =
      mean(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Median_Agreement =
      median(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Fraction_High_Confidence =
      mean(
        Transfer_Confidence ==
          "High"
      ),
    
    Fraction_Low_Agreement =
      mean(
        Transfer_Agreement <
          low_agreement_threshold
      ),
    
    Fraction_High_Agreement =
      mean(
        Transfer_Agreement >=
          very_high_agreement_threshold
      )
  )


cat(
  "\nMonocyte_Myeloid transfer QC:\n"
)

print(
  monocyte_transfer_qc
)


###############################################################################
# Biological annotation note
###############################################################################

cat(
  "\nBiological annotation note:\n"
)

cat(
  "Monocyte_Myeloid is intentionally broad at the global atlas stage.\n"
)

cat(
  "Macrophage-associated markers are diagnostic only here.\n"
)

cat(
  "Macrophage-state resolution is reserved for Phase 5.\n"
)


###############################################################################
# SECTION 30 — FULL-DATASET CELL-TYPE COMPOSITION
###############################################################################

cat("\n============================================================\n")
cat("SECTION 30 — FULL-DATASET COMPOSITION\n")
cat("============================================================\n")


###############################################################################
# Overall composition
###############################################################################

full_Cell_Type_counts <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Cells
    )
  )


cat(
  "\nFull-dataset cell-type composition:\n"
)

print(
  full_Cell_Type_counts
)


write.csv(
  full_Cell_Type_counts,
  file = file.path(
    results_dir,
    "phase3_full_dataset_Cell_Type_composition.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Composition by clinical phase
###############################################################################

composition_by_phase <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Phase,
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  ) %>%
  
  dplyr::ungroup()


write.csv(
  composition_by_phase,
  file = file.path(
    results_dir,
    "phase3_full_dataset_composition_by_phase.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Composition by donor
###############################################################################

composition_by_donor <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    Phase,
    Final_Cell_Type,
    name = "Cells"
  )


write.csv(
  composition_by_donor,
  file = file.path(
    results_dir,
    "phase3_full_dataset_composition_by_donor.csv"
  ),
  row.names = FALSE
)


###############################################################################
# High-confidence sensitivity analysis
#
# IMPORTANT:
# This does NOT modify or filter the final object.
###############################################################################

high_confidence_composition <-
  seurat_obj@meta.data %>%
  
  dplyr::filter(
    Transfer_Confidence ==
      "High"
  ) %>%
  
  dplyr::count(
    Phase,
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  ) %>%
  
  dplyr::ungroup()


write.csv(
  high_confidence_composition,
  file = file.path(
    results_dir,
    "phase3_high_confidence_composition_sensitivity.csv"
  ),
  row.names = FALSE
)


###############################################################################
# SECTION 31 — FULL-DATASET UMAP
###############################################################################

cat("\n============================================================\n")
cat("SECTION 31 — FULL-DATASET UMAP\n")
cat("============================================================\n")


###############################################################################
# IMPORTANT:
#
# This is a NEW UMAP fitted to the full-dataset projected PCA coordinates.
#
# It is NOT a projection of the atlas UMAP.
###############################################################################

if (
  "umap.full" %in%
  Reductions(seurat_obj)
) {
  
  seurat_obj[["umap.full"]] <- NULL
}


cat(
  "Fitting full-dataset UMAP using projected PCs ",
  paste(
    selected_pcs,
    collapse = ":"
  ),
  "...\n",
  sep = ""
)


seurat_obj <- RunUMAP(
  object = seurat_obj,
  reduction = "pca.full",
  dims = selected_pcs,
  reduction.name = "umap.full",
  reduction.key = "UMAPFull_",
  seed.use = 20260918,
  verbose = TRUE
)


###############################################################################
# Validate UMAP
###############################################################################

if (
  !"umap.full" %in%
  Reductions(seurat_obj)
) {
  
  stop(
    "umap.full reduction was not created.",
    call. = FALSE
  )
}


full_umap <- Embeddings(
  seurat_obj,
  reduction = "umap.full"
)


if (
  nrow(full_umap) !=
  n_full_cells_expected
) {
  
  stop(
    "umap.full contains ",
    nrow(full_umap),
    " cells; expected ",
    n_full_cells_expected,
    ".",
    call. = FALSE
  )
}


if (
  ncol(full_umap) !=
  2L
) {
  
  stop(
    "umap.full is not two-dimensional.",
    call. = FALSE
  )
}


###############################################################################
# UMAP by final cell type
###############################################################################

umap_plot_Cell_Type <- DimPlot(
  object = seurat_obj,
  reduction = "umap.full",
  group.by = "Final_Cell_Type",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  
  ggtitle(
    "HBV Liver scRNA-seq — Full Dataset Cell-Type Annotation"
  ) +
  
  theme_classic()


ggsave(
  filename = file.path(
    results_dir,
    "phase3_full_dataset_umap_Cell_Type.pdf"
  ),
  plot = umap_plot_Cell_Type,
  width = 11,
  height = 8
)


###############################################################################
# UMAP by clinical phase
###############################################################################

umap_plot_phase <- DimPlot(
  object = seurat_obj,
  reduction = "umap.full",
  group.by = "Phase",
  raster = TRUE
) +
  
  ggtitle(
    "HBV Liver scRNA-seq — Full Dataset Clinical States"
  ) +
  
  theme_classic()


ggsave(
  filename = file.path(
    results_dir,
    "phase3_full_dataset_umap_phase.pdf"
  ),
  plot = umap_plot_phase,
  width = 10,
  height = 8
)


###############################################################################
# UMAP by transfer confidence
###############################################################################

umap_plot_confidence <- DimPlot(
  object = seurat_obj,
  reduction = "umap.full",
  group.by = "Transfer_Confidence",
  raster = TRUE
) +
  
  ggtitle(
    "HBV Liver scRNA-seq — Label Transfer Confidence"
  ) +
  
  theme_classic()


ggsave(
  filename = file.path(
    results_dir,
    "phase3_full_dataset_umap_transfer_confidence.pdf"
  ),
  plot = umap_plot_confidence,
  width = 10,
  height = 8
)


###############################################################################
# SECTION 32 — TRANSFER AGREEMENT FIGURES
###############################################################################

cat("\n============================================================\n")
cat("SECTION 32 — TRANSFER AGREEMENT\n")
cat("============================================================\n")


###############################################################################
# Agreement distribution
###############################################################################

agreement_df <- seurat_obj@meta.data %>%
  
  dplyr::mutate(
    Cell_ID =
      rownames(
        seurat_obj@meta.data
      )
  )


agreement_plot <- ggplot(
  agreement_df,
  aes(
    x =
      Transfer_Agreement
  )
) +
  
  geom_histogram(
    bins = 50
  ) +
  
  geom_vline(
    xintercept =
      low_agreement_threshold,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept =
      very_high_agreement_threshold,
    linetype = "dashed"
  ) +
  
  labs(
    title =
      "Full-Dataset KNN Label-Transfer Agreement",
    x =
      "Neighbour agreement",
    y =
      "Cells"
  ) +
  
  theme_classic()


ggsave(
  filename = file.path(
    results_dir,
    "phase3_transfer_agreement_distribution.pdf"
  ),
  plot = agreement_plot,
  width = 9,
  height = 6
)


###############################################################################
# Agreement by cell type
###############################################################################

agreement_by_Cell_Type <- agreement_df %>%
  
  dplyr::group_by(
    Final_Cell_Type
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Mean_Agreement =
      mean(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Median_Agreement =
      median(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Q1_Agreement =
      as.numeric(
        quantile(
          Transfer_Agreement,
          0.25,
          na.rm = TRUE
        )
      ),
    
    Q3_Agreement =
      as.numeric(
        quantile(
          Transfer_Agreement,
          0.75,
          na.rm = TRUE
        )
      ),
    
    Fraction_High =
      mean(
        Transfer_Confidence ==
          "High"
      ),
    
    Fraction_Intermediate =
      mean(
        Transfer_Confidence ==
          "Intermediate"
      ),
    
    Fraction_Low =
      mean(
        Transfer_Confidence ==
          "Low"
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Mean_Agreement
    )
  )


print(
  agreement_by_Cell_Type
)


write.csv(
  agreement_by_Cell_Type,
  file = file.path(
    results_dir,
    "phase3_transfer_agreement_by_Cell_Type.csv"
  ),
  row.names = FALSE
)


###############################################################################
# SECTION 33 — AGREEMENT BY CLINICAL STATE
###############################################################################

cat("\n============================================================\n")
cat("SECTION 33 — AGREEMENT BY CLINICAL STATE\n")
cat("============================================================\n")


###############################################################################
# Agreement by Phase
###############################################################################

agreement_by_phase <- seurat_obj@meta.data %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Mean_Agreement =
      mean(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Median_Agreement =
      median(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Fraction_High =
      mean(
        Transfer_Confidence ==
          "High"
      ),
    
    Fraction_Intermediate =
      mean(
        Transfer_Confidence ==
          "Intermediate"
      ),
    
    Fraction_Low =
      mean(
        Transfer_Confidence ==
          "Low"
      ),
    
    Fraction_Low_Agreement =
      mean(
        Transfer_Agreement <
          low_agreement_threshold
      ),
    
    .groups = "drop"
  )


print(
  agreement_by_phase
)


write.csv(
  agreement_by_phase,
  file = file.path(
    results_dir,
    "phase3_transfer_agreement_by_phase.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Low-agreement counts by Phase
###############################################################################

low_agreement_by_phase <- seurat_obj@meta.data %>%
  
  dplyr::mutate(
    Low_Agreement =
      Transfer_Agreement <
      low_agreement_threshold
  ) %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Low_Agreement_Cells =
      sum(
        Low_Agreement
      ),
    
    Low_Agreement_Fraction =
      mean(
        Low_Agreement
      ),
    
    .groups = "drop"
  )


print(
  low_agreement_by_phase
)


write.csv(
  low_agreement_by_phase,
  file = file.path(
    results_dir,
    "phase3_low_agreement_by_phase.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Agreement by donor
###############################################################################

agreement_by_donor <- seurat_obj@meta.data %>%
  
  dplyr::group_by(
    Donor,
    Phase
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Mean_Agreement =
      mean(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Median_Agreement =
      median(
        Transfer_Agreement,
        na.rm = TRUE
      ),
    
    Fraction_High =
      mean(
        Transfer_Confidence ==
          "High"
      ),
    
    Fraction_Low =
      mean(
        Transfer_Confidence ==
          "Low"
      ),
    
    .groups = "drop"
  )


write.csv(
  agreement_by_donor,
  file = file.path(
    results_dir,
    "phase3_transfer_agreement_by_donor.csv"
  ),
  row.names = FALSE
)


###############################################################################
# SECTION 34 — FINAL OBJECT + INTEGRITY CHECKS + SAVE
###############################################################################

cat("\n============================================================\n")
cat("SECTION 34 — FINAL PHASE 3 OBJECT\n")
cat("============================================================\n")


###############################################################################
# Final dimensions
###############################################################################

final_cells <- ncol(
  seurat_obj
)


final_genes <- nrow(
  seurat_obj[["RNA"]]
)


cat(
  "Final cells: ",
  final_cells,
  "\n",
  sep = ""
)


cat(
  "Final genes: ",
  final_genes,
  "\n",
  sep = ""
)


if (
  final_cells !=
  n_full_cells_expected
) {
  
  stop(
    "FINAL CELL COUNT FAILURE: expected ",
    n_full_cells_expected,
    " but found ",
    final_cells,
    ".",
    call. = FALSE
  )
}


if (
  final_genes !=
  n_genes_expected
) {
  
  stop(
    "FINAL GENE COUNT FAILURE: expected ",
    n_genes_expected,
    " but found ",
    final_genes,
    ".",
    call. = FALSE
  )
}


###############################################################################
# Final metadata validation
###############################################################################

required_final_metadata <- c(
  "Phase",
  "Donor",
  "GSM",
  "Transferred_Atlas_Cluster",
  "Final_Cell_Type",
  "Transfer_Agreement",
  "Transfer_SecondBest_Agreement",
  "Transfer_Agreement_Margin",
  "Transfer_Confidence",
  "Low_Confidence"
)


missing_final_metadata <- setdiff(
  required_final_metadata,
  colnames(
    seurat_obj@meta.data
  )
)


if (
  length(missing_final_metadata) > 0L
) {
  
  stop(
    "Missing final metadata columns: ",
    paste(
      missing_final_metadata,
      collapse = ", "
    ),
    call. = FALSE
  )
}


###############################################################################
# Final metadata completeness
###############################################################################

if (
  anyNA(
    seurat_obj$Final_Cell_Type
  )
) {
  
  stop(
    "Final_Cell_Type contains NA values.",
    call. = FALSE
  )
}


if (
  anyNA(
    seurat_obj$Transfer_Confidence
  )
) {
  
  stop(
    "Transfer_Confidence contains NA values.",
    call. = FALSE
  )
}


if (
  any(
    !is.finite(
      seurat_obj$Transfer_Agreement
    )
  )
) {
  
  stop(
    "Transfer_Agreement contains non-finite values.",
    call. = FALSE
  )
}


if (
  any(
    seurat_obj$Transfer_Agreement < 0 |
    seurat_obj$Transfer_Agreement > 1
  )
) {
  
  stop(
    "Transfer_Agreement contains values outside [0,1].",
    call. = FALSE
  )
}


###############################################################################
# Clinical metadata validation
###############################################################################

n_final_phases <- length(
  unique(
    as.character(
      seurat_obj$Phase
    )
  )
)


n_final_donors <- length(
  unique(
    as.character(
      seurat_obj$Donor
    )
  )
)


n_final_samples <- length(
  unique(
    as.character(
      seurat_obj$GSM
    )
  )
)


cat(
  "Clinical states: ",
  n_final_phases,
  "\n",
  sep = ""
)


cat(
  "Donors: ",
  n_final_donors,
  "\n",
  sep = ""
)


cat(
  "Samples: ",
  n_final_samples,
  "\n",
  sep = ""
)


if (
  n_final_phases !=
  5L
) {
  
  stop(
    "Expected 5 clinical states but found ",
    n_final_phases,
    ".",
    call. = FALSE
  )
}


if (
  n_final_donors !=
  23L
) {
  
  stop(
    "Expected 23 donors but found ",
    n_final_donors,
    ".",
    call. = FALSE
  )
}


if (
  n_final_samples !=
  23L
) {
  
  stop(
    "Expected 23 samples but found ",
    n_final_samples,
    ".",
    call. = FALSE
  )
}


###############################################################################
# Verify no cells were removed
###############################################################################

if (
  !identical(
    rownames(
      seurat_obj@meta.data
    ),
    colnames(
      seurat_obj
    )
  )
) {
  
  stop(
    "Metadata cell identities do not match Seurat cell identities.",
    call. = FALSE
  )
}


if (
  nrow(
    seurat_obj@meta.data
  ) !=
  n_full_cells_expected
) {
  
  stop(
    "Final metadata does not contain all 106,592 cells.",
    call. = FALSE
  )
}


###############################################################################
# Final cell-type counts
###############################################################################

final_Cell_Type_counts <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Cells
    )
  )


cat(
  "\nFINAL CELL-TYPE COUNTS:\n"
)

print(
  final_Cell_Type_counts
)


###############################################################################
# Final confidence counts
###############################################################################

final_confidence_counts <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Transfer_Confidence,
    name = "Cells"
  ) %>%
  
  dplyr::mutate(
    Fraction =
      Cells /
      sum(Cells),
    
    Percent =
      100 *
      Fraction
  )


cat(
  "\nFINAL TRANSFER-CONFIDENCE COUNTS:\n"
)

print(
  final_confidence_counts
)


###############################################################################
# Required reductions
###############################################################################

required_reductions <- c(
  "pca.full",
  "umap.full"
)


missing_reductions <- setdiff(
  required_reductions,
  Reductions(
    seurat_obj
  )
)


if (
  length(missing_reductions) > 0L
) {
  
  stop(
    "Required reductions missing: ",
    paste(
      missing_reductions,
      collapse = ", "
    ),
    call. = FALSE
  )
}


###############################################################################
# Validate pca.full
###############################################################################

final_pca <- Embeddings(
  seurat_obj,
  reduction = "pca.full"
)


if (
  nrow(final_pca) !=
  n_full_cells_expected
) {
  
  stop(
    "pca.full does not contain all full-dataset cells.",
    call. = FALSE
  )
}


if (
  ncol(final_pca) !=
  n_pcs
) {
  
  stop(
    "pca.full does not contain exactly 50 dimensions.",
    call. = FALSE
  )
}


if (
  !identical(
    rownames(final_pca),
    colnames(seurat_obj)
  )
) {
  
  stop(
    "pca.full cell ordering is inconsistent.",
    call. = FALSE
  )
}


###############################################################################
# Validate umap.full
###############################################################################

final_umap <- Embeddings(
  seurat_obj,
  reduction = "umap.full"
)


if (
  nrow(final_umap) !=
  n_full_cells_expected
) {
  
  stop(
    "umap.full does not contain all full-dataset cells.",
    call. = FALSE
  )
}


if (
  ncol(final_umap) !=
  2L
) {
  
  stop(
    "umap.full is not two-dimensional.",
    call. = FALSE
  )
}


###############################################################################
# Validate original 23 data layers remain intact
###############################################################################

final_rna_layers <- Layers(
  seurat_obj[["RNA"]]
)


final_data_layers <- grep(
  "^data\\.",
  final_rna_layers,
  value = TRUE
)


if (
  length(final_data_layers) !=
  23L
) {
  
  stop(
    "Final object does not retain all 23 sample-specific data layers.",
    call. = FALSE
  )
}


###############################################################################
# Final integrity report
###############################################################################

final_integrity <- list(
  
  Phase =
    3,
  
  Seurat_version =
    as.character(
      packageVersion(
        "Seurat"
      )
    ),
  
  SeuratObject_version =
    as.character(
      packageVersion(
        "SeuratObject"
      )
    ),
  
  Representation =
    "GEO processed log-CP10K expression",
  
  Raw_counts_available =
    FALSE,
  
  Normalization_performed_in_Phase3 =
    FALSE,
  
  Batch_correction_performed_in_Phase3 =
    FALSE,
  
  Full_cells =
    final_cells,
  
  Full_genes =
    final_genes,
  
  Reference_cells =
    ncol(
      atlas_sketch
    ),
  
  Atlas_clusters =
    n_atlas_clusters,
  
  Atlas_PCA_dimensions =
    n_pcs,
  
  Atlas_PCA_features =
    length(
      pca_features
    ),
  
  Label_transfer_dimensions =
    selected_pcs,
  
  KNN_k =
    knn_k,
  
  Low_agreement_threshold =
    low_agreement_threshold,
  
  Low_margin_threshold =
    low_margin_threshold,
  
  Very_high_agreement_threshold =
    very_high_agreement_threshold,
  
  Clinical_states =
    n_final_phases,
  
  Donors =
    n_final_donors,
  
  Samples =
    n_final_samples,
  
  Cells_removed_in_Phase3 =
    0L,
  
  Original_data_layers_retained =
    length(
      final_data_layers
    ),
  
  Required_reductions =
    required_reductions,
  
  Final_Cell_Types =
    sort(
      unique(
        as.character(
          seurat_obj$Final_Cell_Type
        )
      )
    ),
  
  Transfer_confidence_levels =
    sort(
      unique(
        as.character(
          seurat_obj$Transfer_Confidence
        )
      )
    )
)


###############################################################################
# Save integrity report
###############################################################################

saveRDS(
  final_integrity,
  file = file.path(
    output_dir,
    "phase3_final_full_dataset_integrity.rds"
  )
)


###############################################################################
# Save compact integrity CSV
###############################################################################

integrity_table <- data.frame(
  
  Full_cells =
    final_integrity$Full_cells,
  
  Full_genes =
    final_integrity$Full_genes,
  
  Reference_cells =
    final_integrity$Reference_cells,
  
  Atlas_clusters =
    final_integrity$Atlas_clusters,
  
  Atlas_PCA_dimensions =
    final_integrity$Atlas_PCA_dimensions,
  
  Atlas_PCA_features =
    final_integrity$Atlas_PCA_features,
  
  Clinical_states =
    final_integrity$Clinical_states,
  
  Donors =
    final_integrity$Donors,
  
  Samples =
    final_integrity$Samples,
  
  PCA_dimensions =
    ncol(final_pca),
  
  KNN_k =
    final_integrity$KNN_k,
  
  Cells_removed_in_Phase3 =
    final_integrity$Cells_removed_in_Phase3,
  
  Original_data_layers_retained =
    final_integrity$Original_data_layers_retained
)


write.csv(
  integrity_table,
  file = file.path(
    results_dir,
    "phase3_final_full_dataset_integrity.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Save final metadata
###############################################################################

final_metadata <- seurat_obj@meta.data %>%
  
  dplyr::mutate(
    Cell_ID =
      rownames(
        seurat_obj@meta.data
      )
  ) %>%
  
  dplyr::select(
    Cell_ID,
    everything()
  )


write.csv(
  final_metadata,
  file = file.path(
    results_dir,
    "phase3_final_full_dataset_metadata.csv"
  ),
  row.names = FALSE
)


###############################################################################
# Save final Seurat object
###############################################################################

final_output_path <- file.path(
  output_dir,
  "phase3_final_full_dataset.rds"
)


saveRDS(
  seurat_obj,
  file = final_output_path
)


###############################################################################
# FINAL VALIDATION MESSAGE
###############################################################################

cat("\n")
cat("============================================================\n")
cat("PHASE 3 COMPLETE\n")
cat("============================================================\n")


cat(
  "Final object:\n",
  final_output_path,
  "\n\n"
)


cat(
  "Cells: ",
  final_cells,
  "\n",
  sep = ""
)


cat(
  "Genes: ",
  final_genes,
  "\n",
  sep = ""
)


cat(
  "Reference cells: ",
  ncol(atlas_sketch),
  "\n",
  sep = ""
)


cat(
  "Atlas clusters: ",
  n_atlas_clusters,
  "\n",
  sep = ""
)


cat(
  "Atlas PCA features: ",
  length(pca_features),
  "\n",
  sep = ""
)


cat(
  "Clinical states: ",
  n_final_phases,
  "\n",
  sep = ""
)


cat(
  "Donors: ",
  n_final_donors,
  "\n",
  sep = ""
)


cat(
  "Samples: ",
  n_final_samples,
  "\n",
  sep = ""
)


cat(
  "Projected PCA dimensions: ",
  ncol(final_pca),
  "\n",
  sep = ""
)


cat(
  "Label-transfer dimensions: ",
  paste(
    selected_pcs,
    collapse = ":"
  ),
  "\n",
  sep = ""
)


cat(
  "KNN k: ",
  knn_k,
  "\n",
  sep = ""
)


cat(
  "Cells removed in Phase 3: 0\n"
)


cat(
  "\nNo additional normalization was performed.\n"
)


cat(
  "No batch correction was performed.\n"
)


cat(
  "No cells were filtered during label transfer.\n"
)


cat(
  "No JoinLayers() operation was used.\n"
)


cat(
  "\nThe 23 original processed-expression layers remain intact.\n"
)


cat(
  "Projected PCA, transferred labels, confidence metrics, and a NEW",
  " full-dataset UMAP were added to the final object.\n"
)


cat("============================================================\n")

