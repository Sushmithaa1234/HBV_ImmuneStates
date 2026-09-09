# ============================================================================
# PHASE 3 — FULL-DATASET PCA PROJECTION + LABEL TRANSFER + FINAL ATLAS
# ============================================================================
#
# Script: phase_03_03_full_dataset_projection_transfer.R
#
# Reworked directly from Sections 25–34 of the original full Phase 3 script.
# The uploaded pre-existing Script 03 was NOT used as source material.
#
# Inputs:
#   results/rds_objects/seurat_qc_singlets.rds
#   results/rds_objects/phase3_atlas_sketch_final_annotation.rds
#
# Primary final output:
#   results/rds_objects/phase3_final_full_dataset.rds
#
# Compatibility locks with Scripts 01–02:
#   - Atlas assay is RNA.
#   - Section 24 provides cell-level Final_CellType labels.
#   - Full dataset remains 105,220 QC-filtered singlet cells.
#   - Atlas reference remains 20,000 cells and 50 locked PCs.
#   - No new PCA is fitted to the full dataset.
#   - KNN label transfer is deterministic and batched.
#   - Transfer confidence is QC metadata only; cells are not removed.
#   - CD8_T and Inflammatory_Myeloid target clusters are identified
#     dynamically from the Section 24 annotation.
#
# ============================================================================

cat("============================================================\n")
cat("PHASE 3 — SCRIPT 03: FULL DATASET PROJECTION + LABEL TRANSFER\n")
cat("============================================================\n\n")

# ============================================================================
# ENVIRONMENT
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(here)
  library(Matrix)
  library(RANN)
})

setwd(here::here())
set.seed(12345)

dir.create("results/rds_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Locked structural parameters.
n_full_cells_expected <- 105220L
n_reference_cells <- 20000L
n_pcs <- 50L
selected_pcs <- 1:20
knn_k <- 30L
projection_batch_size <- 5000L
knn_batch_size <- 5000L
low_agreement_threshold <- 0.80
low_margin_threshold <- 0.20
very_high_agreement_threshold <- 0.95

cat("Environment initialized.\n")
cat("Full dataset target: ", n_full_cells_expected, " cells\n", sep = "")
cat("Atlas reference target: ", n_reference_cells, " cells\n", sep = "")
cat("Locked PCA dimensions: ", n_pcs, "\n", sep = "")
cat("KNN k: ", knn_k, "\n", sep = "")
cat("Low-agreement threshold: ", low_agreement_threshold, "\n", sep = "")
cat("Low-margin threshold: ", low_margin_threshold, "\n\n", sep = "")

# SECTION 25 — PROJECT FULL DATASET INTO LOCKED ATLAS PCA SPACE
# ============================================================

message("============================================================")
message("SECTION 25 — PROJECT FULL DATASET INTO LOCKED ATLAS PCA SPACE")
message("============================================================")

# ------------------------------------------------------------
# 25.1 — Load Section 01 and Section 24 checkpoints
# ------------------------------------------------------------

seurat_obj <- readRDS(
  here::here(
    "results/rds_objects/seurat_qc_singlets.rds"
  )
)

atlas_sketch <- readRDS(
  here::here(
    "results/rds_objects/phase3_atlas_sketch_final_annotation.rds"
  )
)

message("Full dataset cells: ", ncol(seurat_obj))
message("Final atlas sketch cells: ", ncol(atlas_sketch))

if (ncol(seurat_obj) != 105220) {
  stop(
    "Expected 105,220 cells in the full QC-filtered dataset; found ",
    ncol(seurat_obj), "."
  )
}

if (ncol(atlas_sketch) != 20000) {
  stop(
    "Final atlas sketch does not contain exactly 20,000 cells."
  )
}

# ------------------------------------------------------------
# 25.2 — Confirm required assays
# ------------------------------------------------------------

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop(
    "The full dataset does not contain the RNA assay."
  )
}

message("Full dataset assay: RNA")
message("Atlas reference assay: RNA")

# ------------------------------------------------------------
# 25.3 — Confirm locked atlas PCA
# ------------------------------------------------------------

if (!"pca" %in% Reductions(atlas_sketch)) {
  stop(
    "Locked PCA reduction was not found in atlas_sketch."
  )
}

sketch_pca_embeddings <- Embeddings(
  atlas_sketch,
  reduction = "pca"
)

locked_loadings <- Loadings(
  atlas_sketch,
  reduction = "pca"
)

if (nrow(sketch_pca_embeddings) != 20000) {
  stop(
    "The locked atlas PCA does not contain exactly 20,000 cells."
  )
}

if (ncol(sketch_pca_embeddings) != 50) {
  stop(
    "The locked atlas PCA does not contain exactly 50 PCs."
  )
}

if (ncol(locked_loadings) != 50) {
  stop(
    "The locked PCA loading matrix does not contain exactly 50 PCs."
  )
}

message(
  "Locked atlas PCA confirmed: ",
  nrow(sketch_pca_embeddings),
  " cells × ",
  ncol(sketch_pca_embeddings),
  " PCs."
)

# ------------------------------------------------------------
# 25.4 — Extract exact locked PCA feature set
# ------------------------------------------------------------

pca_features <- rownames(locked_loadings)

if (length(pca_features) == 0) {
  stop(
    "No features were found in the locked PCA loading matrix."
  )
}

message(
  "Number of locked PCA features: ",
  length(pca_features)
)

# ------------------------------------------------------------
# 25.5 — Verify PCA features in full RNA assay
# ------------------------------------------------------------

full_rna_features <- rownames(
  seurat_obj[["RNA"]]
)

missing_pca_features <- setdiff(
  pca_features,
  full_rna_features
)

if (length(missing_pca_features) > 0) {
  stop(
    "The full RNA assay is missing ",
    length(missing_pca_features),
    " features required by the locked PCA."
  )
}

message(
  "All locked PCA features are present in the full RNA assay."
)

# ------------------------------------------------------------
# 25.6 — Join full RNA layers
# ------------------------------------------------------------
#
# The full dataset may arrive from the QC checkpoint with only
# counts layers. Script 01 normalized the working object in memory,
# but the original QC checkpoint can still contain counts only.
#
# If normalized data layers are absent, apply the SAME LogNormalize
# procedure used in Script 01 before joining the layers. This does
# not fit a new PCA or alter the locked atlas PCA; it only recreates
# the normalized expression layer required for projection.
# ------------------------------------------------------------

rna_layers_before <- Layers(seurat_obj[["RNA"]])

message("RNA layers before normalization/join:")
print(rna_layers_before)

has_normalized_layer <- any(
  grepl("^data", rna_layers_before)
)

if (!has_normalized_layer) {
  
  message(
    "No normalized RNA data layers detected. ",
    "Applying Script 01 LogNormalize to the full dataset..."
  )
  
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
  
  rna_layers_after_normalization <- Layers(
    seurat_obj[["RNA"]]
  )
  
  message("RNA layers after normalization:")
  print(rna_layers_after_normalization)
  
  if (!any(grepl("^data", rna_layers_after_normalization))) {
    stop(
      "LogNormalize completed but no normalized RNA data layer ",
      "was created."
    )
  }
  
} else {
  
  message(
    "Normalized RNA data layer(s) already present; ",
    "no additional normalization performed."
  )
}

message("Joining full-dataset RNA layers...")

seurat_obj <- JoinLayers(
  object = seurat_obj,
  assay = "RNA"
)

message("RNA layers after joining:")
print(Layers(seurat_obj[["RNA"]]))

if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
  stop(
    "The joined RNA assay does not contain a normalized 'data' layer."
  )
}

# ------------------------------------------------------------
# 25.7 — Extract only the features required by locked PCA
# ------------------------------------------------------------
#
# IMPORTANT:
#
# We deliberately extract ONLY the locked PCA features.
# We do not load the complete 18,925-gene matrix into the
# projection calculation.
# ------------------------------------------------------------

message(
  "Extracting normalized expression for locked PCA features..."
)

full_pca_data <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer = "data"
)[pca_features, , drop = FALSE]

message(
  "Projection matrix: ",
  nrow(full_pca_data),
  " features × ",
  ncol(full_pca_data),
  " cells."
)

if (nrow(full_pca_data) != length(pca_features)) {
  stop(
    "Extracted PCA input does not contain the expected features."
  )
}

if (ncol(full_pca_data) != ncol(seurat_obj)) {
  stop(
    "Extracted PCA input does not contain all full-dataset cells."
  )
}

# ------------------------------------------------------------
# 25.8 — Verify feature order
# ------------------------------------------------------------

if (!identical(
  rownames(full_pca_data),
  pca_features
)) {
  stop(
    "Full-dataset PCA features are not in the exact locked order."
  )
}

# ------------------------------------------------------------
# 25.9 — Prepare locked PCA loadings
# ------------------------------------------------------------

locked_loadings <- locked_loadings[
  pca_features,
  1:50,
  drop = FALSE
]

if (!identical(
  rownames(locked_loadings),
  pca_features
)) {
  stop(
    "Locked PCA loading features are not aligned with the ",
    "full-dataset PCA input."
  )
}

# ------------------------------------------------------------
# 25.10 — Determine PCA scaling parameters
# ------------------------------------------------------------
#
# The PCA was generated from the scaled sketch assay.
#
# We retrieve the scaling parameters from the reference assay
# rather than estimating them again from the full dataset.
# ------------------------------------------------------------

message("Retrieving locked sketch scaling parameters...")

atlas_data <- SeuratObject::LayerData(
  object = atlas_sketch,
  assay = "RNA",
  layer = "data"
)

missing_atlas_features <- setdiff(
  pca_features,
  rownames(atlas_data)
)

if (length(missing_atlas_features) > 0) {
  stop(
    "The atlas RNA data layer is missing ",
    length(missing_atlas_features),
    " locked PCA features."
  )
}

# LayerData() does not guarantee that the returned rows follow the
# requested PCA-feature order. Explicitly reorder them using the
# locked PCA feature vector before calculating scaling parameters.
atlas_data <- atlas_data[
  pca_features,
  ,
  drop = FALSE
]

if (!identical(rownames(atlas_data), pca_features)) {
  stop(
    "Locked atlas RNA data could not be aligned to the exact PCA ",
    "feature order."
  )
}

# ------------------------------------------------------------
# 25.11 — Calculate locked feature means and standard deviations
# ------------------------------------------------------------
#
# These parameters reconstruct the exact scaling convention
# used for the atlas PCA.
# ------------------------------------------------------------

feature_means <- Matrix::rowMeans(atlas_data)

feature_sq_means <- Matrix::rowMeans(
  atlas_data ^ 2
)

feature_sds <- sqrt(
  pmax(
    feature_sq_means - feature_means ^ 2,
    0
  )
)

feature_sds[feature_sds == 0 | !is.finite(feature_sds)] <- 1

rm(
  atlas_data,
  feature_sq_means
)

# ------------------------------------------------------------
# 25.12 — Initialize projected PCA matrix
# ------------------------------------------------------------

n_full_cells <- ncol(full_pca_data)
n_pcs <- 50

projected_pca <- matrix(
  NA_real_,
  nrow = n_full_cells,
  ncol = n_pcs
)

rownames(projected_pca) <- colnames(full_pca_data)

colnames(projected_pca) <- paste0(
  "PC_",
  seq_len(n_pcs)
)

# ------------------------------------------------------------
# 25.13 — Batch-wise PCA projection
# ------------------------------------------------------------
#
# A batch size of 5,000 is deliberately conservative.
#
# This prevents the large matrix multiplication from attempting
# to allocate one enormous intermediate object.
# ------------------------------------------------------------

batch_size <- 5000

batch_starts <- seq(
  1,
  n_full_cells,
  by = batch_size
)

message(
  "Beginning batch-wise projection of ",
  n_full_cells,
  " cells..."
)

for (batch_start in batch_starts) {
  
  batch_end <- min(
    batch_start + batch_size - 1,
    n_full_cells
  )
  
  batch_cells <- batch_start:batch_end
  
  message(
    "Projecting cells ",
    batch_start,
    "–",
    batch_end,
    " of ",
    n_full_cells,
    "..."
  )
  
  batch_data <- full_pca_data[
    ,
    batch_cells,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Apply the scaling parameters learned from the sketch
  # ----------------------------------------------------------
  
  batch_data <- sweep(
    batch_data,
    MARGIN = 1,
    STATS = feature_means,
    FUN = "-"
  )
  
  batch_data <- sweep(
    batch_data,
    MARGIN = 1,
    STATS = feature_sds,
    FUN = "/"
  )
  
  # ----------------------------------------------------------
  # Project onto the locked PCA loading matrix
  # ----------------------------------------------------------
  
  batch_pca <- Matrix::t(
    locked_loadings
  ) %*% batch_data
  
  batch_pca <- as.matrix(
    Matrix::t(batch_pca)
  )
  
  projected_pca[
    batch_cells,
    1:n_pcs
  ] <- batch_pca
  
  rm(
    batch_data,
    batch_pca
  )
  
  gc(verbose = FALSE)
}

message("Batch-wise PCA projection completed.")

# ------------------------------------------------------------
# 25.14 — Projection integrity checks
# ------------------------------------------------------------

if (nrow(projected_pca) != ncol(seurat_obj)) {
  stop(
    "Projected PCA does not contain all full-dataset cells."
  )
}

if (ncol(projected_pca) != 50) {
  stop(
    "Projected PCA does not contain exactly 50 PCs."
  )
}

if (!identical(
  rownames(projected_pca),
  colnames(seurat_obj)
)) {
  stop(
    "Projected PCA cell names do not match the full dataset."
  )
}

if (any(!is.finite(projected_pca))) {
  stop(
    "Non-finite values detected in projected PCA coordinates."
  )
}

message(
  "Projected PCA passed all integrity checks."
)

# ------------------------------------------------------------
# 25.15 — Store projected PCA reduction
# ------------------------------------------------------------

seurat_obj[["pca.full"]] <- CreateDimReducObject(
  embeddings = projected_pca,
  loadings = locked_loadings,
  assay = "RNA",
  key = "PCFull_"
)

message(
  "Projected PCA reduction added as: pca.full"
)

# ------------------------------------------------------------
# 25.16 — Final projection summary
# ------------------------------------------------------------

projection_summary <- data.frame(
  Metric = c(
    "Full dataset cells",
    "Atlas sketch cells",
    "Locked PCA features",
    "Locked PCA dimensions",
    "Projected cells",
    "Projected dimensions",
    "Batch size",
    "Reference reduction",
    "Projected reduction"
  ),
  Value = c(
    ncol(seurat_obj),
    ncol(atlas_sketch),
    length(pca_features),
    50,
    nrow(projected_pca),
    ncol(projected_pca),
    batch_size,
    "atlas_sketch::pca",
    "seurat_obj::pca.full"
  )
)

print(projection_summary)

write.csv(
  projection_summary,
  file = "results/tables/phase3_full_dataset_projection_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 25.17 — Release temporary projection objects
# ------------------------------------------------------------

rm(
  full_pca_data,
  locked_loadings,
  feature_means,
  feature_sds,
  projected_pca
)

gc(verbose = FALSE)

# ------------------------------------------------------------
# 25.18 — Save checkpoint
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "results/rds_objects/phase3_full_dataset_projected.rds"
)

message("============================================================")
message("SECTION 25 COMPLETE")
message("")
message("Full dataset projected into the locked atlas PCA space.")
message("")
message("Reference reduction:")
message("  atlas_sketch::pca")
message("")
message("Projected reduction:")
message("  seurat_obj::pca.full")
message("")
message("Cells projected: ", n_full_cells)
message("PCs projected: ", n_pcs)
message("Batch size: ", batch_size)
message("")
message("No new PCA was fitted.")
message("No new UMAP was fitted.")
message("")
message("Checkpoint:")
message("results/rds_objects/phase3_full_dataset_projected.rds")
message("============================================================")

# ============================================================

# SECTION 26 — PREPARE LABEL TRANSFER / AGREEMENT ANALYSIS
# ============================================================

message("============================================================")
message("SECTION 26 — PREPARE LABEL TRANSFER / AGREEMENT ANALYSIS")
message("============================================================")

# ------------------------------------------------------------
# 26.1 — Basic object checks
# ------------------------------------------------------------

stopifnot(exists("seurat_obj"))
stopifnot(exists("atlas_sketch"))

message("Full dataset cells: ", ncol(seurat_obj))
message("Atlas sketch cells: ", ncol(atlas_sketch))

if (ncol(seurat_obj) != 105220) {
  warning(
    "Expected 105,220 full-dataset cells; found ",
    ncol(seurat_obj),
    ". Verify the input object."
  )
}

if (ncol(atlas_sketch) != 20000) {
  stop(
    "atlas_sketch does not contain exactly 20,000 cells."
  )
}

# ------------------------------------------------------------
# 26.2 — Confirm locked reference PCA
# ------------------------------------------------------------

if (!"pca" %in% Reductions(atlas_sketch)) {
  stop(
    "Locked reference PCA 'pca' not found in atlas_sketch."
  )
}

reference_pca <- Embeddings(
  atlas_sketch,
  reduction = "pca"
)

if (nrow(reference_pca) != 20000) {
  stop(
    "Reference PCA does not contain exactly 20,000 cells."
  )
}

if (ncol(reference_pca) < 50) {
  stop(
    "Reference PCA contains fewer than 50 dimensions."
  )
}

reference_pca <- reference_pca[, 1:50, drop = FALSE]

message(
  "Reference PCA prepared: ",
  nrow(reference_pca),
  " cells × ",
  ncol(reference_pca),
  " PCs."
)

# ------------------------------------------------------------
# 26.3 — Confirm full-dataset projected PCA
# ------------------------------------------------------------

if (!"pca.full" %in% Reductions(seurat_obj)) {
  stop(
    "Projected full-dataset PCA 'pca.full' not found."
  )
}

query_pca <- Embeddings(
  seurat_obj,
  reduction = "pca.full"
)

if (nrow(query_pca) != ncol(seurat_obj)) {
  stop(
    "Projected PCA does not contain all full-dataset cells."
  )
}

if (ncol(query_pca) < 50) {
  stop(
    "Projected full-dataset PCA contains fewer than 50 dimensions."
  )
}

query_pca <- query_pca[, 1:50, drop = FALSE]

message(
  "Query PCA prepared: ",
  nrow(query_pca),
  " cells × ",
  ncol(query_pca),
  " PCs."
)

# ------------------------------------------------------------
# 26.4 — Verify PCA coordinate names and ordering
# ------------------------------------------------------------

if (is.null(rownames(reference_pca))) {
  stop(
    "Reference PCA does not contain cell names."
  )
}

if (is.null(rownames(query_pca))) {
  stop(
    "Projected query PCA does not contain cell names."
  )
}

if (anyDuplicated(rownames(reference_pca))) {
  stop(
    "Duplicate cell names detected in reference PCA."
  )
}

if (anyDuplicated(rownames(query_pca))) {
  stop(
    "Duplicate cell names detected in query PCA."
  )
}

if (!identical(
  rownames(query_pca),
  colnames(seurat_obj)
)) {
  stop(
    "Query PCA cell names do not exactly match seurat_obj."
  )
}

message("PCA cell identities verified.")

# ------------------------------------------------------------
# 26.5 — Load locked final sketch annotation
# ------------------------------------------------------------
#
# Section 24 established the final manual cluster labels.
# These are cluster-level annotations, so they must be mapped
# to individual sketch cells using seurat_clusters.
# ------------------------------------------------------------

final_annotation_file <- here::here(
  "results/rds_objects/phase3_atlas_sketch_final_annotation.rds"
)

if (!file.exists(final_annotation_file)) {
  stop(
    "Section 24 final annotation checkpoint was not found: ",
    final_annotation_file
  )
}

atlas_sketch_final <- readRDS(
  final_annotation_file
)

message(
  "Loaded Section 24 final annotation checkpoint."
)

# ------------------------------------------------------------
# 26.6 — Confirm final annotation metadata
# ------------------------------------------------------------

required_annotation_metadata <- c(
  "seurat_clusters",
  "Final_CellType",
  "Annotation_Confidence"
)

missing_annotation_metadata <- setdiff(
  required_annotation_metadata,
  colnames(atlas_sketch_final@meta.data)
)

if (length(missing_annotation_metadata) > 0) {
  stop(
    "Missing Section 24 annotation metadata: ",
    paste(missing_annotation_metadata, collapse = ", ")
  )
}

if (anyNA(atlas_sketch_final$Final_CellType)) {
  stop("Section 24 contains missing Final_CellType values.")
}

if (anyNA(atlas_sketch_final$Annotation_Confidence)) {
  stop("Section 24 contains missing Annotation_Confidence values.")
}

message(
  "Section 24 cell-level final annotations detected: ",
  length(unique(as.character(atlas_sketch_final$Final_CellType))),
  " cell types."
)

# ------------------------------------------------------------
# 26.7 — Validate one final label per atlas cluster
# ------------------------------------------------------------

reference_cluster_map <- data.frame(
  Cluster = as.character(atlas_sketch_final$seurat_clusters),
  Final_CellType = as.character(atlas_sketch_final$Final_CellType),
  stringsAsFactors = FALSE
)

cluster_label_counts <- aggregate(
  Final_CellType ~ Cluster,
  data = reference_cluster_map,
  FUN = function(x) length(unique(x))
)

if (any(cluster_label_counts$Final_CellType != 1L)) {
  stop(
    "At least one atlas cluster maps to multiple final cell types."
  )
}

observed_clusters <- sort(
  unique(reference_cluster_map$Cluster),
  method = "radix"
)

if (length(observed_clusters) != 17L) {
  stop(
    "Expected 17 atlas clusters from Scripts 01–02; found ",
    length(observed_clusters), "."
  )
}

message(
  "Validated one final cell-type label for each of the ",
  length(observed_clusters),
  " atlas clusters."
)

# ------------------------------------------------------------
# 26.8 — Map final labels onto individual reference cells
# ------------------------------------------------------------

reference_labels <- as.character(
  atlas_sketch_final$Final_CellType
)

names(reference_labels) <- colnames(
  atlas_sketch_final
)

if (any(is.na(reference_labels))) {
  stop(
    "NA reference labels detected after cluster-to-cell mapping."
  )
}

message("Reference annotation distribution:")
print(table(reference_labels))

# ------------------------------------------------------------
# 26.9 — Verify reference PCA / labels alignment
# ------------------------------------------------------------

if (!identical(
  rownames(reference_pca),
  names(reference_labels)
)) {
  stop(
    "Reference PCA cell order does not match reference ",
    "annotation label order."
  )
}

message(
  "Reference PCA and final manual labels are perfectly aligned."
)

# ------------------------------------------------------------
# 26.10 — Lock KNN parameters
# ------------------------------------------------------------

knn_k <- 30

if (knn_k >= nrow(reference_pca)) {
  stop(
    "KNN parameter k must be smaller than the number ",
    "of reference cells."
  )
}

message("Locked KNN parameter: k = ", knn_k)

# ------------------------------------------------------------
# 26.11 — Create explicit reference/query objects
# ------------------------------------------------------------

label_transfer_reference <- list(
  cells = rownames(reference_pca),
  pca = reference_pca,
  labels = reference_labels
)

label_transfer_query <- list(
  cells = rownames(query_pca),
  pca = query_pca
)

# ------------------------------------------------------------
# 26.12 — Verify reference/query dimensions
# ------------------------------------------------------------

if (
  nrow(label_transfer_reference$pca) !=
  length(label_transfer_reference$labels)
) {
  stop(
    "Reference PCA cell count does not match reference ",
    "annotation label count."
  )
}

if (
  nrow(label_transfer_query$pca) !=
  length(label_transfer_query$cells)
) {
  stop(
    "Query PCA cell count does not match query cell count."
  )
}

message(
  "Reference object: ",
  nrow(label_transfer_reference$pca),
  " cells × ",
  ncol(label_transfer_reference$pca),
  " PCs."
)

message(
  "Query object: ",
  nrow(label_transfer_query$pca),
  " cells × ",
  ncol(label_transfer_query$pca),
  " PCs."
)

# ------------------------------------------------------------
# 26.13 — Save preparation objects
# ------------------------------------------------------------

saveRDS(
  label_transfer_reference,
  file = here::here(
    "results/rds_objects/phase3_label_transfer_reference.rds"
  )
)

saveRDS(
  label_transfer_query,
  file = here::here(
    "results/rds_objects/phase3_label_transfer_query.rds"
  )
)

# ------------------------------------------------------------
# 26.14 — Save parameter summary
# ------------------------------------------------------------

label_transfer_parameters <- data.frame(
  Parameter = c(
    "Reference cells",
    "Query cells",
    "Reference dimensions",
    "Query dimensions",
    "PCA dimensions used",
    "KNN k",
    "Reference reduction",
    "Query reduction",
    "Reference annotation source"
  ),
  Value = c(
    nrow(reference_pca),
    nrow(query_pca),
    ncol(reference_pca),
    ncol(query_pca),
    "PC1:PC50",
    knn_k,
    "atlas_sketch::pca",
    "seurat_obj::pca.full",
    "Section 24 final manual cluster annotations"
  )
)

print(label_transfer_parameters)

write.csv(
  label_transfer_parameters,
  file = here::here(
    "results/tables/phase3_label_transfer_parameters.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 26.15 — Final validation
# ------------------------------------------------------------

stopifnot(
  nrow(label_transfer_reference$pca) == 20000
)

stopifnot(
  nrow(label_transfer_query$pca) == 105220
)

stopifnot(
  ncol(label_transfer_reference$pca) == 50
)

stopifnot(
  ncol(label_transfer_query$pca) == 50
)

stopifnot(
  length(label_transfer_reference$labels) == 20000
)

stopifnot(
  all(is.finite(label_transfer_reference$pca))
)

stopifnot(
  all(is.finite(label_transfer_query$pca))
)

stopifnot(
  !anyNA(label_transfer_reference$labels)
)

message("============================================================")
message("SECTION 26 COMPLETE")
message("")
message("Reference: 20,000 annotated sketch cells")
message("Query: 105,220 full-dataset cells")
message("Shared space: locked PCA, PC1:PC50")
message("KNN parameter: k = ", knn_k)
message("")
message("Reference labels: Section 24 final manual annotations")
message("")
message("Ready for Section 27:")
message("KNN neighbour search + label agreement calculation")
message("============================================================")

# ============================================================

# SECTION 27 — KNN NEIGHBOUR AGREEMENT / LABEL TRANSFER
# ============================================================

message("============================================================")
message("SECTION 27 — KNN NEIGHBOUR AGREEMENT / LABEL TRANSFER")
message("============================================================")

# ------------------------------------------------------------
# 27.1 — Required package
# ------------------------------------------------------------

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop(
    "Package 'RANN' is required for memory-efficient nearest ",
    "neighbour search. Install it before continuing."
  )
}

# ------------------------------------------------------------
# 27.2 — Load Section 26 preparation objects if necessary
# ------------------------------------------------------------

if (!exists("label_transfer_reference")) {
  
  reference_file <- here::here(
    "results/rds_objects/phase3_label_transfer_reference.rds"
  )
  
  if (!file.exists(reference_file)) {
    stop(
      "Section 26 reference checkpoint not found: ",
      reference_file
    )
  }
  
  label_transfer_reference <- readRDS(
    reference_file
  )
}

if (!exists("label_transfer_query")) {
  
  query_file <- here::here(
    "results/rds_objects/phase3_label_transfer_query.rds"
  )
  
  if (!file.exists(query_file)) {
    stop(
      "Section 26 query checkpoint not found: ",
      query_file
    )
  }
  
  label_transfer_query <- readRDS(
    query_file
  )
}

# ------------------------------------------------------------
# 27.3 — Extract reference and query PCA matrices
# ------------------------------------------------------------

reference_pca <- label_transfer_reference$pca

query_pca <- label_transfer_query$pca

reference_labels <- label_transfer_reference$labels

# ------------------------------------------------------------
# 27.4 — Basic integrity checks
# ------------------------------------------------------------

if (nrow(reference_pca) != 20000) {
  stop(
    "Reference PCA does not contain exactly 20,000 cells."
  )
}

if (nrow(query_pca) != 105220) {
  stop(
    "Query PCA does not contain exactly 105,220 cells."
  )
}

if (ncol(reference_pca) != 50) {
  stop(
    "Reference PCA does not contain exactly 50 dimensions."
  )
}

if (ncol(query_pca) != 50) {
  stop(
    "Query PCA does not contain exactly 50 dimensions."
  )
}

if (length(reference_labels) != 20000) {
  stop(
    "Reference annotation vector does not contain exactly ",
    "20,000 labels."
  )
}

if (!identical(
  rownames(reference_pca),
  names(reference_labels)
)) {
  stop(
    "Reference PCA cell names and reference labels are not aligned."
  )
}

if (any(!is.finite(reference_pca))) {
  stop(
    "Non-finite values detected in reference PCA."
  )
}

if (any(!is.finite(query_pca))) {
  stop(
    "Non-finite values detected in query PCA."
  )
}

message(
  "Reference: ",
  nrow(reference_pca),
  " cells × ",
  ncol(reference_pca),
  " PCs"
)

message(
  "Query: ",
  nrow(query_pca),
  " cells × ",
  ncol(query_pca),
  " PCs"
)

# ------------------------------------------------------------
# 27.5 — Lock KNN parameter
# ------------------------------------------------------------

knn_k <- 30

if (knn_k >= nrow(reference_pca)) {
  stop(
    "KNN k must be smaller than the number of reference cells."
  )
}

message("KNN parameter: k = ", knn_k)

# ------------------------------------------------------------
# 27.6 — Prepare reference annotation codes
# ------------------------------------------------------------

reference_label_levels <- sort(
  unique(reference_labels)
)

reference_label_codes <- match(
  reference_labels,
  reference_label_levels
)

message(
  "Reference annotation classes: ",
  length(reference_label_levels)
)

print(reference_label_levels)

# ------------------------------------------------------------
# 27.7 — Initialize output vectors
# ------------------------------------------------------------

n_query <- nrow(query_pca)

transferred_label <- character(n_query)

top_label_count <- integer(n_query)

second_label_count <- integer(n_query)

agreement <- numeric(n_query)

label_margin <- numeric(n_query)

# ------------------------------------------------------------
# 27.8 — Configure memory-efficient query batches
# ------------------------------------------------------------

batch_size <- 5000

batch_starts <- seq(
  1,
  n_query,
  by = batch_size
)

message(
  "Query batch size: ",
  batch_size
)

message(
  "Number of query batches: ",
  length(batch_starts)
)

# ------------------------------------------------------------
# 27.9 — KNN search and label voting
# ------------------------------------------------------------
#
# RANN searches the 20,000 reference cells for each query
# batch and returns ONLY the k nearest neighbours.
#
# No full query × reference distance matrix is created.
# ------------------------------------------------------------

message("Beginning KNN search and label voting...")

for (batch_start in batch_starts) {
  
  batch_end <- min(
    batch_start + batch_size - 1,
    n_query
  )
  
  batch_indices <- batch_start:batch_end
  
  message(
    "Processing query cells ",
    batch_start,
    "–",
    batch_end,
    " of ",
    n_query,
    "..."
  )
  
  # ----------------------------------------------------------
  # Extract query batch
  # ----------------------------------------------------------
  
  query_batch <- query_pca[
    batch_indices,
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # KNN search
  # ----------------------------------------------------------
  
  nn_result <- RANN::nn2(
    data = reference_pca,
    query = query_batch,
    k = knn_k,
    treetype = "kd",
    searchtype = "standard"
  )
  
  # ----------------------------------------------------------
  # Explicitly convert neighbour indices to matrix
  # ----------------------------------------------------------
  
  nn_idx <- nn_result$nn.idx
  
  if (is.null(dim(nn_idx))) {
    
    nn_idx <- matrix(
      nn_idx,
      nrow = nrow(query_batch),
      ncol = knn_k,
      byrow = TRUE
    )
  }
  
  if (
    nrow(nn_idx) != nrow(query_batch) ||
    ncol(nn_idx) != knn_k
  ) {
    stop(
      "Unexpected dimensions returned by RANN::nn2()."
    )
  }
  
  # ----------------------------------------------------------
  # Convert neighbour indices to annotation codes
  # ----------------------------------------------------------
  
  neighbour_codes <- matrix(
    reference_label_codes[nn_idx],
    nrow = nrow(nn_idx),
    ncol = ncol(nn_idx)
  )
  
  # ----------------------------------------------------------
  # Process each query cell
  # ----------------------------------------------------------
  
  for (i in seq_len(nrow(neighbour_codes))) {
    
    global_index <- batch_indices[i]
    
    counts <- tabulate(
      neighbour_codes[i, ],
      nbins = length(reference_label_levels)
    )
    
    ranked_codes <- order(
      counts,
      decreasing = TRUE
    )
    
    top_code <- ranked_codes[1]
    second_code <- ranked_codes[2]
    
    top_count <- counts[top_code]
    second_count <- counts[second_code]
    
    transferred_label[global_index] <-
      reference_label_levels[top_code]
    
    top_label_count[global_index] <-
      top_count
    
    second_label_count[global_index] <-
      second_count
    
    agreement[global_index] <-
      top_count / knn_k
    
    label_margin[global_index] <-
      (top_count - second_count) / knn_k
  }
  
  # ----------------------------------------------------------
  # Release batch memory
  # ----------------------------------------------------------
  
  rm(
    query_batch,
    nn_result,
    nn_idx,
    neighbour_codes
  )
  
  gc(verbose = FALSE)
}

message("KNN search and label voting completed.")

# ------------------------------------------------------------
# 27.10 — Construct cell-level agreement table
# ------------------------------------------------------------

knn_agreement <- data.frame(
  Cell = rownames(query_pca),
  Transferred_Label = transferred_label,
  Top_Label_Neighbours = top_label_count,
  Second_Label_Neighbours = second_label_count,
  Neighbour_Agreement = agreement,
  Label_Margin = label_margin,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 27.11 — Basic result validation
# ------------------------------------------------------------

if (nrow(knn_agreement) != n_query) {
  stop(
    "KNN agreement table does not contain all query cells."
  )
}

if (!identical(
  knn_agreement$Cell,
  rownames(query_pca)
)) {
  stop(
    "KNN agreement cell order does not match query PCA order."
  )
}

if (anyNA(knn_agreement$Transferred_Label)) {
  stop(
    "NA transferred labels detected."
  )
}

if (any(!is.finite(
  knn_agreement$Neighbour_Agreement
))) {
  stop(
    "Non-finite neighbour agreement values detected."
  )
}

if (any(
  knn_agreement$Neighbour_Agreement < 0 |
  knn_agreement$Neighbour_Agreement > 1
)) {
  stop(
    "Neighbour agreement values outside [0,1] detected."
  )
}

if (any(
  knn_agreement$Top_Label_Neighbours < 1 |
  knn_agreement$Top_Label_Neighbours > knn_k
)) {
  stop(
    "Invalid top-label neighbour counts detected."
  )
}

# ------------------------------------------------------------
# 27.12 — Confirm complete cell assignment
# ------------------------------------------------------------

assigned_cells <- sum(
  !is.na(knn_agreement$Transferred_Label)
)

message(
  "Cells with transferred labels: ",
  assigned_cells,
  " / ",
  n_query
)

if (assigned_cells != n_query) {
  stop(
    "Not all full-dataset cells received a transferred label."
  )
}

# ------------------------------------------------------------
# 27.13 — Agreement summary
# ------------------------------------------------------------

agreement_summary <- data.frame(
  Metric = c(
    "Query cells",
    "Reference cells",
    "PCA dimensions",
    "KNN k",
    "Cells assigned",
    "Mean agreement",
    "Median agreement",
    "Minimum agreement",
    "Maximum agreement"
  ),
  Value = c(
    n_query,
    nrow(reference_pca),
    ncol(reference_pca),
    knn_k,
    assigned_cells,
    mean(knn_agreement$Neighbour_Agreement),
    median(knn_agreement$Neighbour_Agreement),
    min(knn_agreement$Neighbour_Agreement),
    max(knn_agreement$Neighbour_Agreement)
  )
)

print(agreement_summary)

# ------------------------------------------------------------
# 27.14 — Transferred-label composition
# ------------------------------------------------------------

label_transfer_counts <- as.data.frame(
  table(
    knn_agreement$Transferred_Label
  )
)

colnames(label_transfer_counts) <- c(
  "Transferred_Label",
  "Cells"
)

label_transfer_counts$Fraction <-
  label_transfer_counts$Cells /
  sum(label_transfer_counts$Cells)

print(label_transfer_counts)

# ------------------------------------------------------------
# 27.15 — Agreement distribution
# ------------------------------------------------------------

agreement_distribution <- as.data.frame(
  table(
    knn_agreement$Neighbour_Agreement
  )
)

colnames(agreement_distribution) <- c(
  "Neighbour_Agreement",
  "Cells"
)

# ------------------------------------------------------------
# 27.16 — Save cell-level KNN results
# ------------------------------------------------------------

write.csv(
  knn_agreement,
  file = here::here(
    "results/tables/phase3_knn_label_transfer_agreement.csv"
  ),
  row.names = FALSE
)

write.csv(
  label_transfer_counts,
  file = here::here(
    "results/tables/phase3_knn_transferred_label_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  agreement_summary,
  file = here::here(
    "results/tables/phase3_knn_agreement_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  agreement_distribution,
  file = here::here(
    "results/tables/phase3_knn_agreement_distribution.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 27.17 — Save KNN checkpoint
# ------------------------------------------------------------

saveRDS(
  knn_agreement,
  file = here::here(
    "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
  )
)

# ------------------------------------------------------------
# 27.18 — Final message
# ------------------------------------------------------------

message("============================================================")
message("SECTION 27 COMPLETE")
message("")
message("KNN label transfer completed.")
message("")
message("Reference cells: ", nrow(reference_pca))
message("Query cells: ", nrow(query_pca))
message("Dimensions: PC1:PC50")
message("KNN k: ", knn_k)
message("")
message("Every full-dataset cell received:")
message("  - transferred label")
message("  - top-label neighbour count")
message("  - second-label neighbour count")
message("  - neighbour agreement")
message("  - label margin")
message("")
message("No full query × reference distance matrix was created.")
message("")
message("Checkpoint:")
message(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)
message("============================================================")

# ============================================================

# SECTION 28 — FLAG LOW-AGREEMENT LABEL TRANSFERS
# ============================================================

message("============================================================")
message("SECTION 28 — LOW-AGREEMENT LABEL TRANSFER FLAGGING")
message("============================================================")
message("")

# ------------------------------------------------------------
# 28.1 Load KNN label-transfer results
# ------------------------------------------------------------

if (!exists("knn_agreement")) {
  
  knn_agreement <- readRDS(
    here::here(
      "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
    )
  )
  
}

# ------------------------------------------------------------
# 28.2 Validate required columns
# ------------------------------------------------------------

required_columns <- c(
  "Cell",
  "Transferred_Label",
  "Top_Label_Neighbours",
  "Second_Label_Neighbours",
  "Neighbour_Agreement",
  "Label_Margin"
)

missing_columns <- setdiff(
  required_columns,
  colnames(knn_agreement)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
  
}

# ------------------------------------------------------------
# 28.3 Validate cell count and values
# ------------------------------------------------------------

if (nrow(knn_agreement) != 105220) {
  
  stop(
    paste(
      "Unexpected number of transferred cells:",
      nrow(knn_agreement),
      "expected 105220."
    )
  )
  
}

if (anyDuplicated(knn_agreement$Cell) > 0) {
  
  stop(
    "Duplicate cell identifiers detected in KNN transfer table."
  )
  
}

if (any(!is.finite(knn_agreement$Neighbour_Agreement))) {
  
  stop(
    "Non-finite neighbour agreement values detected."
  )
  
}

if (any(!is.finite(knn_agreement$Label_Margin))) {
  
  stop(
    "Non-finite label margin values detected."
  )
  
}

if (any(
  knn_agreement$Neighbour_Agreement < 0 |
  knn_agreement$Neighbour_Agreement > 1
)) {
  
  stop(
    "Neighbour agreement contains values outside [0, 1]."
  )
  
}

if (any(
  knn_agreement$Label_Margin < 0 |
  knn_agreement$Label_Margin > 1
)) {
  
  stop(
    "Label margin contains values outside [0, 1]."
  )
  
}

# ------------------------------------------------------------
# 28.4 Define low-agreement criteria
# ------------------------------------------------------------
#
# These are QC flags only.
#
# Low agreement:
#   Neighbour_Agreement < 0.80
#
# Weak label separation:
#   Label_Margin < 0.20
#
# A cell is flagged if either criterion is met.
#
# IMPORTANT:
# No cells are removed here.
# No transferred labels are changed here.
# ------------------------------------------------------------

agreement_threshold <- 0.80
margin_threshold <- 0.20

knn_agreement$Low_Agreement <-
  knn_agreement$Neighbour_Agreement <
  agreement_threshold

knn_agreement$Low_Margin <-
  knn_agreement$Label_Margin <
  margin_threshold

knn_agreement$Low_Confidence <-
  knn_agreement$Low_Agreement |
  knn_agreement$Low_Margin

# ------------------------------------------------------------
# 28.5 Confidence category
# ------------------------------------------------------------

knn_agreement$Transfer_Confidence <- "High"

knn_agreement$Transfer_Confidence[
  knn_agreement$Low_Confidence
] <- "Low"

knn_agreement$Transfer_Confidence[
  !knn_agreement$Low_Agreement &
    knn_agreement$Low_Margin
] <- "Intermediate"

# ------------------------------------------------------------
# 28.6 Additional agreement categories
# ------------------------------------------------------------

knn_agreement$Agreement_Category <- cut(
  knn_agreement$Neighbour_Agreement,
  breaks = c(
    -Inf,
    0.50,
    0.80,
    0.95,
    Inf
  ),
  labels = c(
    "<50%",
    "50–<80%",
    "80–<95%",
    "≥95%"
  ),
  right = FALSE
)

# ------------------------------------------------------------
# 28.7 Summary statistics
# ------------------------------------------------------------

n_total <- nrow(knn_agreement)

n_low_agreement <- sum(
  knn_agreement$Low_Agreement
)

n_low_margin <- sum(
  knn_agreement$Low_Margin
)

n_low_confidence <- sum(
  knn_agreement$Low_Confidence
)

n_high_confidence <- sum(
  !knn_agreement$Low_Confidence
)

n_intermediate <- sum(
  knn_agreement$Transfer_Confidence ==
    "Intermediate"
)

summary_table <- data.frame(
  Metric = c(
    "Total transferred cells",
    "Low agreement cells",
    "Low margin cells",
    "Low-confidence cells",
    "Intermediate-confidence cells",
    "High-confidence cells",
    "Low agreement threshold",
    "Low margin threshold"
  ),
  Value = c(
    n_total,
    n_low_agreement,
    n_low_margin,
    n_low_confidence,
    n_intermediate,
    n_high_confidence,
    agreement_threshold,
    margin_threshold
  )
)

summary_table$Fraction <- NA_real_

summary_table$Fraction[2] <-
  n_low_agreement / n_total

summary_table$Fraction[3] <-
  n_low_margin / n_total

summary_table$Fraction[4] <-
  n_low_confidence / n_total

summary_table$Fraction[5] <-
  n_intermediate / n_total

summary_table$Fraction[6] <-
  n_high_confidence / n_total

print(summary_table)

# ------------------------------------------------------------
# 28.8 Confidence category counts
# ------------------------------------------------------------

confidence_counts <- as.data.frame(
  table(
    knn_agreement$Transfer_Confidence
  )
)

colnames(confidence_counts) <- c(
  "Transfer_Confidence",
  "Cells"
)

confidence_counts$Fraction <-
  confidence_counts$Cells / n_total

print(confidence_counts)

# ------------------------------------------------------------
# 28.9 Agreement category counts
# ------------------------------------------------------------

agreement_category_counts <- as.data.frame(
  table(
    knn_agreement$Agreement_Category
  )
)

colnames(agreement_category_counts) <- c(
  "Agreement_Category",
  "Cells"
)

agreement_category_counts$Fraction <-
  agreement_category_counts$Cells / n_total

print(agreement_category_counts)

# ------------------------------------------------------------
# 28.10 Transferred label × confidence
# ------------------------------------------------------------

label_confidence_table <- as.data.frame(
  table(
    knn_agreement$Transferred_Label,
    knn_agreement$Transfer_Confidence
  )
)

colnames(label_confidence_table) <- c(
  "Transferred_Label",
  "Transfer_Confidence",
  "Cells"
)

label_confidence_table$Fraction_Within_Label <- NA_real_

for (label in unique(
  label_confidence_table$Transferred_Label
)) {
  
  total_label_cells <- sum(
    label_confidence_table$Cells[
      label_confidence_table$Transferred_Label ==
        label
    ]
  )
  
  label_confidence_table$Fraction_Within_Label[
    label_confidence_table$Transferred_Label ==
      label
  ] <-
    label_confidence_table$Cells[
      label_confidence_table$Transferred_Label ==
        label
    ] / total_label_cells
  
}

print(label_confidence_table)

# ------------------------------------------------------------
# 28.11 Identify the lowest-agreement transferred cells
# ------------------------------------------------------------

lowest_agreement_cells <- knn_agreement[
  order(
    knn_agreement$Neighbour_Agreement,
    knn_agreement$Label_Margin
  ),
  required_columns
]

lowest_agreement_cells <- head(
  lowest_agreement_cells,
  100
)

# ------------------------------------------------------------
# 28.12 Save low-confidence cell table
# ------------------------------------------------------------

low_confidence_cells <- knn_agreement[
  knn_agreement$Low_Confidence,
  ,
  drop = FALSE
]

write.csv(
  low_confidence_cells,
  file = here::here(
    "results/tables/phase3_low_confidence_transfers.csv"
  ),
  row.names = FALSE
)

write.csv(
  lowest_agreement_cells,
  file = here::here(
    "results/tables/phase3_lowest_agreement_cells_top100.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 28.13 Save summary tables
# ------------------------------------------------------------

write.csv(
  summary_table,
  file = here::here(
    "results/tables/phase3_low_agreement_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  confidence_counts,
  file = here::here(
    "results/tables/phase3_transfer_confidence_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  agreement_category_counts,
  file = here::here(
    "results/tables/phase3_agreement_category_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  label_confidence_table,
  file = here::here(
    "results/tables/phase3_label_by_transfer_confidence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 28.14 Final validation
# ------------------------------------------------------------

if (nrow(knn_agreement) != n_total) {
  
  stop(
    "Unexpected row count after confidence flagging."
  )
  
}

if (anyNA(knn_agreement$Transfer_Confidence)) {
  
  stop(
    "NA values detected in Transfer_Confidence."
  )
  
}

if (anyNA(knn_agreement$Agreement_Category)) {
  
  stop(
    "NA values detected in Agreement_Category."
  )
  
}

if (
  sum(
    knn_agreement$Low_Confidence
  ) !=
  n_low_confidence
) {
  
  stop(
    "Low-confidence cell count validation failed."
  )
  
}

# ------------------------------------------------------------
# 28.15 Save checkpoint
# ------------------------------------------------------------

saveRDS(
  knn_agreement,
  file = here::here(
    "results/rds_objects/phase3_low_agreement_flagged.rds"
  )
)

# ------------------------------------------------------------
# 28.16 Report
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 28 COMPLETE")
message("============================================================")
message("")

message(
  "Total transferred cells: ",
  n_total
)

message(
  "Low-agreement cells (< ",
  agreement_threshold * 100,
  "%): ",
  n_low_agreement,
  " (",
  round(100 * n_low_agreement / n_total, 2),
  "%)"
)

message(
  "Low-margin cells (< ",
  margin_threshold * 100,
  "%): ",
  n_low_margin,
  " (",
  round(100 * n_low_margin / n_total, 2),
  "%)"
)

message(
  "Low-confidence cells: ",
  n_low_confidence,
  " (",
  round(100 * n_low_confidence / n_total, 2),
  "%)"
)

message(
  "Intermediate-confidence cells: ",
  n_intermediate,
  " (",
  round(100 * n_intermediate / n_total, 2),
  "%)"
)

message(
  "High-confidence cells: ",
  n_high_confidence,
  " (",
  round(100 * n_high_confidence / n_total, 2),
  "%)"
)

message("")
message("No cells were removed.")
message("No transferred labels were changed.")
message("Low-confidence status is retained as a QC flag.")
message("")
message("Checkpoint:")
message(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)
message("============================================================")

# ============================================================

# SECTION 29 — TARGETED CD8 + MYELOID CLUSTER QC AUDIT
# ============================================================

message("============================================================")
message("SECTION 29 — TARGETED CD8 + MYELOID CLUSTER QC AUDIT")
message("============================================================")
message("")

# ------------------------------------------------------------
# 29.1 Load required objects
# ------------------------------------------------------------

knn_agreement <- readRDS(
  here::here(
    "results/rds_objects/phase3_low_agreement_flagged.rds"
  )
)

atlas_sketch <- readRDS(
  here::here(
    "results/rds_objects/phase3_atlas_sketch_final_annotation.rds"
  )
)

# ------------------------------------------------------------
# 29.2 Locked Section 24 annotation
# ------------------------------------------------------------

reference_cluster_map <- unique(data.frame(
  Cluster = as.character(atlas_sketch$seurat_clusters),
  Final_CellType = as.character(atlas_sketch$Final_CellType),
  stringsAsFactors = FALSE
))

cluster_label_counts <- aggregate(
  Final_CellType ~ Cluster,
  data = data.frame(
    Cluster = as.character(atlas_sketch$seurat_clusters),
    Final_CellType = as.character(atlas_sketch$Final_CellType)
  ),
  FUN = function(x) length(unique(x))
)

if (any(cluster_label_counts$Final_CellType != 1L)) {
  stop("At least one atlas cluster maps to multiple final cell types.")
}

cd8_clusters <- reference_cluster_map$Cluster[
  reference_cluster_map$Final_CellType == "CD8_T"
]

myeloid_clusters <- reference_cluster_map$Cluster[
  reference_cluster_map$Final_CellType == "Inflammatory_Myeloid"
]

if (length(cd8_clusters) == 0) {
  stop("No CD8_T reference clusters found in Section 24.")
}

if (length(myeloid_clusters) == 0) {
  stop("No Inflammatory_Myeloid reference clusters found in Section 24.")
}

# ------------------------------------------------------------
# 29.3 Validate KNN transfer checkpoint
# ------------------------------------------------------------

if (nrow(knn_agreement) != 105220) {
  
  stop(
    paste(
      "KNN transfer checkpoint contains",
      nrow(knn_agreement),
      "cells; expected 105220."
    )
  )
  
}

required_columns <- c(
  "Cell",
  "Transferred_Label",
  "Top_Label_Neighbours",
  "Second_Label_Neighbours",
  "Neighbour_Agreement",
  "Label_Margin",
  "Transfer_Confidence"
)

missing_columns <- setdiff(
  required_columns,
  colnames(knn_agreement)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing KNN columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
  
}

# Explicitly convert transferred labels to character
knn_agreement$Transferred_Label <- as.character(
  knn_agreement$Transferred_Label
)

knn_agreement$Transfer_Confidence <- as.character(
  knn_agreement$Transfer_Confidence
)

# ------------------------------------------------------------
# 29.4 Verify the known Section 27 label counts
# ------------------------------------------------------------

known_label_counts <- table(
  knn_agreement$Transferred_Label
)

print(known_label_counts)

if (unname(known_label_counts["CD8_T"]) == 0) {
  stop("No CD8_T cells were transferred in Section 27.")
}

if (unname(known_label_counts["Inflammatory_Myeloid"]) == 0) {
  stop("No Inflammatory_Myeloid cells were transferred in Section 27.")
}

# ------------------------------------------------------------
# 29.5 Identify full-dataset target cells DIRECTLY
# ------------------------------------------------------------
#
# No intermediate Target_Group variable is used.
# The Section 27 transferred label is the source of truth.
# ------------------------------------------------------------

cd8_cells <- knn_agreement[
  knn_agreement$Transferred_Label ==
    "CD8_T",
  ,
  drop = FALSE
]

myeloid_cells <- knn_agreement[
  knn_agreement$Transferred_Label ==
    "Inflammatory_Myeloid",
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 29.6 Verify target populations
# ------------------------------------------------------------

if (nrow(cd8_cells) == 0) {
  stop("CD8 target contains zero transferred cells.")
}

if (nrow(myeloid_cells) == 0) {
  stop("Inflammatory-myeloid target contains zero transferred cells.")
}

# ------------------------------------------------------------
# 29.7 CD8 agreement QC
# ------------------------------------------------------------

cd8_high <- cd8_cells$Transfer_Confidence ==
  "High"

cd8_low <- cd8_cells$Transfer_Confidence ==
  "Low"

if (anyNA(cd8_high)) {
  
  stop(
    "NA values detected in CD8 transfer confidence."
  )
  
}

cd8_qc <- data.frame(
  Metric = c(
    "Transferred CD8 cells",
    "High-agreement CD8 cells",
    "Low-agreement CD8 cells",
    "High-agreement fraction",
    "Mean agreement",
    "Median agreement",
    "Minimum agreement",
    "Maximum agreement"
  ),
  Value = c(
    nrow(cd8_cells),
    sum(cd8_high),
    sum(cd8_low),
    mean(cd8_high),
    mean(
      cd8_cells$Neighbour_Agreement
    ),
    median(
      cd8_cells$Neighbour_Agreement
    ),
    min(
      cd8_cells$Neighbour_Agreement
    ),
    max(
      cd8_cells$Neighbour_Agreement
    )
  )
)

print(cd8_qc)

# ------------------------------------------------------------
# 29.8 Myeloid agreement QC
# ------------------------------------------------------------

myeloid_high <- myeloid_cells$Transfer_Confidence ==
  "High"

myeloid_low <- myeloid_cells$Transfer_Confidence ==
  "Low"

if (anyNA(myeloid_high)) {
  
  stop(
    "NA values detected in myeloid transfer confidence."
  )
  
}

myeloid_qc <- data.frame(
  Metric = c(
    "Transferred inflammatory-myeloid cells",
    "High-agreement myeloid cells",
    "Low-agreement myeloid cells",
    "High-agreement fraction",
    "Mean agreement",
    "Median agreement",
    "Minimum agreement",
    "Maximum agreement"
  ),
  Value = c(
    nrow(myeloid_cells),
    sum(myeloid_high),
    sum(myeloid_low),
    mean(myeloid_high),
    mean(
      myeloid_cells$Neighbour_Agreement
    ),
    median(
      myeloid_cells$Neighbour_Agreement
    ),
    min(
      myeloid_cells$Neighbour_Agreement
    ),
    max(
      myeloid_cells$Neighbour_Agreement
    )
  )
)

print(myeloid_qc)

# ------------------------------------------------------------
# 29.9 Combined target-population QC
# ------------------------------------------------------------

target_population_qc <- data.frame(
  Population = c(
    "CD8_T",
    "Inflammatory_Myeloid"
  ),
  Reference_Clusters = c(
    "0, 4",
    "12"
  ),
  Transferred_Cells = c(
    nrow(cd8_cells),
    nrow(myeloid_cells)
  ),
  High_Agreement_Cells = c(
    sum(cd8_high),
    sum(myeloid_high)
  ),
  Low_Agreement_Cells = c(
    sum(cd8_low),
    sum(myeloid_low)
  ),
  High_Agreement_Fraction = c(
    mean(cd8_high),
    mean(myeloid_high)
  ),
  Mean_Agreement = c(
    mean(
      cd8_cells$Neighbour_Agreement
    ),
    mean(
      myeloid_cells$Neighbour_Agreement
    )
  ),
  Median_Agreement = c(
    median(
      cd8_cells$Neighbour_Agreement
    ),
    median(
      myeloid_cells$Neighbour_Agreement
    )
  )
)

print(target_population_qc)

# ------------------------------------------------------------
# 29.10 Inspect locked sketch target clusters
# ------------------------------------------------------------

sketch_cluster_vector <- as.character(
  atlas_sketch$seurat_clusters
)

sketch_target_table <- data.frame(
  Cluster = c(
    cd8_clusters,
    myeloid_clusters
  ),
  Locked_Label = c(
    "CD8_T",
    "CD8_T",
    "Inflammatory_Myeloid"
  ),
  Sketch_Cells = c(
    sum(
      sketch_cluster_vector == "0"
    ),
    sum(
      sketch_cluster_vector == "4"
    ),
    sum(
      sketch_cluster_vector == "12"
    )
  )
)

print(sketch_target_table)

# ------------------------------------------------------------
# 29.11 Marker QC
# ------------------------------------------------------------

cd8_markers <- c(
  "CD8A",
  "CD8B",
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",
  "CCL5",
  "GZMK",
  "GZMA"
)

myeloid_markers <- c(
  "LYZ",
  "CTSS",
  "CTSD",
  "CST3",
  "FCER1G",
  "TYROBP",
  "LGALS3",
  "AIF1",
  "CTSB",
  "FCGR3A"
)

macrophage_markers <- c(
  "C1QA",
  "C1QB",
  "C1QC",
  "APOE",
  "TREM2",
  "CD68",
  "MSR1"
)

available_genes <- rownames(
  atlas_sketch[["RNA"]]
)

marker_genes <- intersect(
  unique(
    c(
      cd8_markers,
      myeloid_markers,
      macrophage_markers
    )
  ),
  available_genes
)

sketch_marker_expression <- GetAssayData(
  atlas_sketch,
  assay = "RNA",
  layer = "data"
)[
  marker_genes,
  ,
  drop = FALSE
]

target_marker_table <- data.frame(
  Cluster = c(
    "0",
    "4",
    "12"
  ),
  Locked_Label = c(
    "CD8_T",
    "CD8_T",
    "Inflammatory_Myeloid"
  )
)

for (gene in marker_genes) {
  
  target_marker_table[[gene]] <- sapply(
    c("0", "4", "12"),
    function(cl) {
      
      cells <- colnames(
        atlas_sketch
      )[
        sketch_cluster_vector == cl
      ]
      
      mean(
        sketch_marker_expression[
          gene,
          cells,
          drop = TRUE
        ],
        na.rm = TRUE
      )
      
    }
  )
  
}

print(target_marker_table)

# ------------------------------------------------------------
# 29.12 Macrophage-associated marker summary
# ------------------------------------------------------------

myeloid_sketch_cells <- colnames(
  atlas_sketch
)[
  sketch_cluster_vector == "12"
]

macrophage_marker_means <- data.frame(
  Marker = intersect(
    macrophage_markers,
    marker_genes
  ),
  Mean_Expression = sapply(
    intersect(
      macrophage_markers,
      marker_genes
    ),
    function(gene) {
      
      mean(
        sketch_marker_expression[
          gene,
          myeloid_sketch_cells,
          drop = TRUE
        ],
        na.rm = TRUE
      )
      
    }
  )
)

print(macrophage_marker_means)

# ------------------------------------------------------------
# 29.13 Save target cell identifiers + QC
# ------------------------------------------------------------

target_cells <- rbind(
  data.frame(
    Cell = cd8_cells$Cell,
    Population = "CD8_T",
    Transferred_Label =
      cd8_cells$Transferred_Label,
    Top_Label_Neighbours =
      cd8_cells$Top_Label_Neighbours,
    Second_Label_Neighbours =
      cd8_cells$Second_Label_Neighbours,
    Neighbour_Agreement =
      cd8_cells$Neighbour_Agreement,
    Label_Margin =
      cd8_cells$Label_Margin,
    Transfer_Confidence =
      cd8_cells$Transfer_Confidence,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Cell = myeloid_cells$Cell,
    Population = "Inflammatory_Myeloid",
    Transferred_Label =
      myeloid_cells$Transferred_Label,
    Top_Label_Neighbours =
      myeloid_cells$Top_Label_Neighbours,
    Second_Label_Neighbours =
      myeloid_cells$Second_Label_Neighbours,
    Neighbour_Agreement =
      myeloid_cells$Neighbour_Agreement,
    Label_Margin =
      myeloid_cells$Label_Margin,
    Transfer_Confidence =
      myeloid_cells$Transfer_Confidence,
    stringsAsFactors = FALSE
  )
)

# ------------------------------------------------------------
# 29.14 Write outputs
# ------------------------------------------------------------

write.csv(
  sketch_target_table,
  here::here(
    "results/tables/phase3_target_reference_clusters.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_marker_table,
  here::here(
    "results/tables/phase3_target_cluster_marker_qc.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_qc,
  here::here(
    "results/tables/phase3_cd8_target_qc.csv"
  ),
  row.names = FALSE
)

write.csv(
  myeloid_qc,
  here::here(
    "results/tables/phase3_myeloid_target_qc.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_population_qc,
  here::here(
    "results/tables/phase3_cd8_myeloid_population_qc.csv"
  ),
  row.names = FALSE
)

write.csv(
  macrophage_marker_means,
  here::here(
    "results/tables/phase3_myeloid_macrophage_marker_qc.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_cells,
  here::here(
    "results/tables/phase3_cd8_myeloid_transferred_cells.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 29.15 Save checkpoint
# ------------------------------------------------------------

saveRDS(
  list(
    reference_clusters =
      sketch_target_table,
    marker_qc =
      target_marker_table,
    cd8_qc =
      cd8_qc,
    myeloid_qc =
      myeloid_qc,
    population_qc =
      target_population_qc,
    macrophage_marker_qc =
      macrophage_marker_means,
    transferred_target_cells =
      target_cells
  ),
  here::here(
    "results/rds_objects/phase3_cd8_myeloid_target_qc.rds"
  )
)

# ------------------------------------------------------------
# 29.16 FINAL VALIDATION
# ------------------------------------------------------------

if (nrow(cd8_cells) == 0) {
  stop("FINAL VALIDATION FAILED: no CD8_T cells.")
}

if (nrow(myeloid_cells) == 0) {
  stop("FINAL VALIDATION FAILED: no Inflammatory_Myeloid cells.")
}

if (
  sum(cd8_high) +
  sum(cd8_low) !=
  nrow(cd8_cells)
) {
  
  stop(
    "FINAL VALIDATION FAILED: CD8 confidence counts."
  )
  
}

if (
  sum(myeloid_high) +
  sum(myeloid_low) !=
  nrow(myeloid_cells)
) {
  
  stop(
    "FINAL VALIDATION FAILED: myeloid confidence counts."
  )
  
}

if (
  anyNA(
    target_population_qc
  )
) {
  
  stop(
    "FINAL VALIDATION FAILED: NA values in target QC."
  )
  
}

# ------------------------------------------------------------
# 29.17 Final report
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 29 COMPLETE")
message("============================================================")
message("")

message("LOCKED REFERENCE CLUSTERS:")
message("  CD8_T: clusters ", paste(cd8_clusters, collapse = ", "))
message("  Inflammatory_Myeloid: clusters ", paste(myeloid_clusters, collapse = ", "))

message("")
message("CD8:")
message(
  "  Transferred cells: ",
  nrow(cd8_cells)
)

message(
  "  High-agreement cells: ",
  sum(cd8_high)
)

message(
  "  Low-agreement cells: ",
  sum(cd8_low)
)

message(
  "  High-agreement fraction: ",
  round(
    100 * mean(cd8_high),
    2
  ),
  "%"
)

message("")
message("MYELOID:")
message(
  "  Transferred cells: ",
  nrow(myeloid_cells)
)

message(
  "  High-agreement cells: ",
  sum(myeloid_high)
)

message(
  "  Low-agreement cells: ",
  sum(myeloid_low)
)

message(
  "  High-agreement fraction: ",
  round(
    100 * mean(myeloid_high),
    2
  ),
  "%"
)

message("")
message(
  "Cluster 12 remains locked as Inflammatory_Myeloid."
)

message(
  "No macrophage-specific relabeling was performed."
)

message("")
message("Checkpoint:")
message(
  "results/rds_objects/phase3_cd8_myeloid_target_qc.rds"
)

message("============================================================")

# ============================================================

# SECTION 30 — FULL DATASET COMPOSITION
# ============================================================

cat("\n============================================================\n")
cat("SECTION 30 — FULL DATASET COMPOSITION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 30.1 Load Section 28/29 outputs
# ------------------------------------------------------------

knn_agreement <- readRDS(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)

low_agreement_flagged <- readRDS(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)

cat("Loaded KNN agreement and low-confidence flag objects.\n")


# ------------------------------------------------------------
# 30.2 Validate KNN agreement object
# ------------------------------------------------------------

stopifnot(
  nrow(knn_agreement) == 105220
)

stopifnot(
  all(!is.na(knn_agreement$Transferred_Label))
)

required_knn_cols <- c(
  "Cell",
  "Transferred_Label",
  "Neighbour_Agreement",
  "Label_Margin"
)

missing_knn_cols <- setdiff(
  required_knn_cols,
  colnames(knn_agreement)
)

if (length(missing_knn_cols) > 0) {
  stop(
    "Missing required columns from KNN agreement object: ",
    paste(missing_knn_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 30.3 Recover confidence flags from Section 28
# ------------------------------------------------------------

required_flag_cols <- c(
  "Cell",
  "Low_Confidence",
  "Transfer_Confidence"
)

missing_flag_cols <- setdiff(
  required_flag_cols,
  colnames(low_agreement_flagged)
)

if (length(missing_flag_cols) > 0) {
  stop(
    "Missing required columns from Section 28 flagged object: ",
    paste(missing_flag_cols, collapse = ", ")
  )
}

cat("Section 28 confidence flags validated.\n")


# ------------------------------------------------------------
# 30.4 Align confidence flags to KNN agreement by Cell
# ------------------------------------------------------------

confidence_lookup <- low_agreement_flagged[
  match(
    knn_agreement$Cell,
    low_agreement_flagged$Cell
  ),
  c(
    "Cell",
    "Low_Confidence",
    "Transfer_Confidence"
  )
]

stopifnot(
  identical(
    confidence_lookup$Cell,
    knn_agreement$Cell
  )
)

knn_agreement$Low_Confidence <-
  confidence_lookup$Low_Confidence

knn_agreement$Transfer_Confidence <-
  confidence_lookup$Transfer_Confidence

stopifnot(
  !any(is.na(knn_agreement$Low_Confidence))
)

stopifnot(
  !any(is.na(knn_agreement$Transfer_Confidence))
)

cat(
  "Confidence flags successfully aligned to all ",
  nrow(knn_agreement),
  " cells.\n",
  sep = ""
)


# ------------------------------------------------------------
# 30.5 Construct full-dataset composition table
# ------------------------------------------------------------

# Load the full projected Seurat object
seurat_obj <- readRDS(
  "results/rds_objects/phase3_full_dataset_projected.rds"
)

cat("Loaded full projected Seurat object.\n")

stopifnot(
  ncol(seurat_obj) == 105220
)

meta <- seurat_obj@meta.data


# ------------------------------------------------------------
# 30.6 Identify sample and clinical-state metadata columns
# ------------------------------------------------------------

sample_col <- "GSM"
state_col <- "Phase"
donor_col <- "Donor"

missing_composition_metadata <- setdiff(
  c(sample_col, state_col, donor_col),
  colnames(meta)
)

if (length(missing_composition_metadata) > 0) {
  stop(
    "Missing required composition metadata: ",
    paste(missing_composition_metadata, collapse = ", ")
  )
}

cat("Sample column: ", sample_col, "\n", sep = "")
cat("Clinical-state column: ", state_col, "\n", sep = "")


# ------------------------------------------------------------
# 30.7 Construct composition table
# ------------------------------------------------------------

composition <- data.frame(
  Cell = rownames(meta),
  Sample = as.character(meta[[sample_col]]),
  Donor = as.character(meta[[donor_col]]),
  Clinical_State = as.character(meta[[state_col]]),
  Transferred_Label = knn_agreement$Transferred_Label,
  Neighbour_Agreement = knn_agreement$Neighbour_Agreement,
  Label_Margin = knn_agreement$Label_Margin,
  Low_Confidence = knn_agreement$Low_Confidence,
  Transfer_Confidence = knn_agreement$Transfer_Confidence,
  stringsAsFactors = FALSE
)

# Confirm exact cell alignment
stopifnot(
  identical(
    composition$Cell,
    knn_agreement$Cell
  )
)

# Confirm no missing composition metadata
stopifnot(
  !any(is.na(composition$Transferred_Label))
)

stopifnot(
  !any(is.na(composition$Sample))
)

stopifnot(
  !any(is.na(composition$Clinical_State))
)

stopifnot(
  !any(is.na(composition$Donor))
)

cat("\nComposition table constructed successfully.\n")
cat("Cells: ", nrow(composition), "\n", sep = "")
cat(
  "Cell types: ",
  length(unique(composition$Transferred_Label)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 30.8 Overall full-dataset composition
# ------------------------------------------------------------

overall_counts <- as.data.frame(
  table(composition$Transferred_Label),
  stringsAsFactors = FALSE
)

colnames(overall_counts) <- c(
  "Cell_Type",
  "Cell_Count"
)

overall_counts$Percentage <- (
  overall_counts$Cell_Count /
    sum(overall_counts$Cell_Count)
) * 100

overall_counts <- overall_counts[
  order(-overall_counts$Cell_Count),
]

rownames(overall_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("OVERALL FULL-DATASET COMPOSITION\n")
cat("------------------------------------------------------------\n")

print(overall_counts)

# Integrity check
stopifnot(
  sum(overall_counts$Cell_Count) == 105220
)

# ------------------------------------------------------------
# 30.9 Composition by clinical state — counts
# ------------------------------------------------------------

state_counts <- as.data.frame(
  table(
    composition$Clinical_State,
    composition$Transferred_Label
  ),
  stringsAsFactors = FALSE
)

colnames(state_counts) <- c(
  "Clinical_State",
  "Cell_Type",
  "Cell_Count"
)

state_counts$Percentage <- NA_real_

for (st in unique(state_counts$Clinical_State)) {
  
  idx <- state_counts$Clinical_State == st
  
  state_total <- sum(
    state_counts$Cell_Count[idx]
  )
  
  state_counts$Percentage[idx] <- (
    state_counts$Cell_Count[idx] /
      state_total
  ) * 100
}

state_counts <- state_counts[
  order(
    state_counts$Clinical_State,
    -state_counts$Cell_Count
  ),
]

rownames(state_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("COMPOSITION BY CLINICAL STATE\n")
cat("------------------------------------------------------------\n")

print(state_counts)

# Integrity check
stopifnot(
  sum(state_counts$Cell_Count) == 105220
)


# ------------------------------------------------------------
# 30.10 Composition by donor/sample — counts
# ------------------------------------------------------------

sample_counts <- as.data.frame(
  table(
    composition$Sample,
    composition$Transferred_Label
  ),
  stringsAsFactors = FALSE
)

colnames(sample_counts) <- c(
  "Sample",
  "Cell_Type",
  "Cell_Count"
)

sample_counts$Percentage <- NA_real_

for (samp in unique(sample_counts$Sample)) {
  
  idx <- sample_counts$Sample == samp
  
  sample_total <- sum(
    sample_counts$Cell_Count[idx]
  )
  
  sample_counts$Percentage[idx] <- (
    sample_counts$Cell_Count[idx] /
      sample_total
  ) * 100
}

sample_counts <- sample_counts[
  order(
    sample_counts$Sample,
    -sample_counts$Cell_Count
  ),
]

rownames(sample_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("COMPOSITION BY DONOR/SAMPLE\n")
cat("------------------------------------------------------------\n")

print(
  head(
    sample_counts,
    30
  )
)

# Integrity check
stopifnot(
  sum(sample_counts$Cell_Count) == 105220
)

# ------------------------------------------------------------
# 30.11 Define high-confidence sensitivity dataset
# ------------------------------------------------------------

high_confidence <- composition[
  composition$Transfer_Confidence == "High",
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("HIGH-CONFIDENCE SENSITIVITY DATASET\n")
cat("------------------------------------------------------------\n")

cat(
  "High-confidence cells: ",
  nrow(high_confidence),
  "\n",
  sep = ""
)

cat(
  "Low-confidence cells excluded from sensitivity analysis: ",
  sum(composition$Transfer_Confidence != "High"),
  "\n",
  sep = ""
)

high_confidence_n <- nrow(high_confidence)

stopifnot(
  high_confidence_n > 0,
  high_confidence_n <= 105220
)

stopifnot(
  all(high_confidence$Transfer_Confidence == "High")
)


# ------------------------------------------------------------
# 30.12 High-confidence overall composition
# ------------------------------------------------------------

high_confidence_counts <- as.data.frame(
  table(high_confidence$Transferred_Label),
  stringsAsFactors = FALSE
)

colnames(high_confidence_counts) <- c(
  "Cell_Type",
  "High_Confidence_Cell_Count"
)

high_confidence_counts$Percentage <- (
  high_confidence_counts$High_Confidence_Cell_Count /
    sum(high_confidence_counts$High_Confidence_Cell_Count)
) * 100

high_confidence_counts <- high_confidence_counts[
  order(-high_confidence_counts$High_Confidence_Cell_Count),
]

rownames(high_confidence_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("HIGH-CONFIDENCE OVERALL COMPOSITION\n")
cat("------------------------------------------------------------\n")

print(high_confidence_counts)

stopifnot(
  sum(high_confidence_counts$High_Confidence_Cell_Count) ==
    high_confidence_n
)


# ------------------------------------------------------------
# 30.13 High-confidence composition by clinical state
# ------------------------------------------------------------

high_conf_state_counts <- as.data.frame(
  table(
    high_confidence$Clinical_State,
    high_confidence$Transferred_Label
  ),
  stringsAsFactors = FALSE
)

colnames(high_conf_state_counts) <- c(
  "Clinical_State",
  "Cell_Type",
  "High_Confidence_Cell_Count"
)

high_conf_state_counts$Percentage <- NA_real_

for (st in unique(high_conf_state_counts$Clinical_State)) {
  
  idx <- high_conf_state_counts$Clinical_State == st
  
  state_total <- sum(
    high_conf_state_counts$High_Confidence_Cell_Count[idx]
  )
  
  high_conf_state_counts$Percentage[idx] <- (
    high_conf_state_counts$High_Confidence_Cell_Count[idx] /
      state_total
  ) * 100
}

high_conf_state_counts <- high_conf_state_counts[
  order(
    high_conf_state_counts$Clinical_State,
    -high_conf_state_counts$High_Confidence_Cell_Count
  ),
]

rownames(high_conf_state_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("HIGH-CONFIDENCE COMPOSITION BY CLINICAL STATE\n")
cat("------------------------------------------------------------\n")

print(high_conf_state_counts)

stopifnot(
  sum(high_conf_state_counts$High_Confidence_Cell_Count) ==
    high_confidence_n
)


# ------------------------------------------------------------
# 30.14 Final integrity checks
# ------------------------------------------------------------

# Overall counts must sum to the complete dataset
stopifnot(
  sum(overall_counts$Cell_Count) == 105220
)

# Clinical-state counts must also sum to the complete dataset
stopifnot(
  sum(state_counts$Cell_Count) == 105220
)

# Sample counts must also sum to the complete dataset
stopifnot(
  sum(sample_counts$Cell_Count) == 105220
)

# No missing transferred labels
stopifnot(
  sum(is.na(composition$Transferred_Label)) == 0
)

# No missing clinical state
stopifnot(
  sum(is.na(composition$Clinical_State)) == 0
)

# No missing sample
stopifnot(
  sum(is.na(composition$Sample)) == 0
)

cat("\n============================================================\n")
cat("SECTION 30 COMPLETE — ALL INTEGRITY CHECKS PASSED\n")
cat("============================================================\n")

# ============================================================
# SECTION 30 — FINAL SAVE / CHECKPOINT
# ============================================================

# ------------------------------------------------------------
# 30.14 Clinical-state percentage matrix
# ------------------------------------------------------------

state_percentage_matrix <- reshape(
  state_counts[, c(
    "Clinical_State",
    "Cell_Type",
    "Percentage"
  )],
  idvar = "Clinical_State",
  timevar = "Cell_Type",
  direction = "wide"
)

colnames(state_percentage_matrix) <- sub(
  "^Percentage\\.",
  "",
  colnames(state_percentage_matrix)
)

state_percentage_matrix <- state_percentage_matrix[
  order(state_percentage_matrix$Clinical_State),
]

rownames(state_percentage_matrix) <- NULL


# ------------------------------------------------------------
# 30.15 Sample percentage matrix
# ------------------------------------------------------------

sample_percentage_matrix <- reshape(
  sample_counts[, c(
    "Sample",
    "Cell_Type",
    "Percentage"
  )],
  idvar = "Sample",
  timevar = "Cell_Type",
  direction = "wide"
)

colnames(sample_percentage_matrix) <- sub(
  "^Percentage\\.",
  "",
  colnames(sample_percentage_matrix)
)

sample_percentage_matrix <- sample_percentage_matrix[
  order(sample_percentage_matrix$Sample),
]

rownames(sample_percentage_matrix) <- NULL


# ------------------------------------------------------------
# 30.16 Save composition tables
# ------------------------------------------------------------

write.csv(
  overall_counts,
  "results/tables/phase3_full_dataset_composition_overall.csv",
  row.names = FALSE
)

write.csv(
  state_counts,
  "results/tables/phase3_full_dataset_composition_by_clinical_state.csv",
  row.names = FALSE
)

write.csv(
  sample_counts,
  "results/tables/phase3_full_dataset_composition_by_sample.csv",
  row.names = FALSE
)

write.csv(
  state_percentage_matrix,
  "results/tables/phase3_full_dataset_composition_state_percentage_matrix.csv",
  row.names = FALSE
)

write.csv(
  sample_percentage_matrix,
  "results/tables/phase3_full_dataset_composition_sample_percentage_matrix.csv",
  row.names = FALSE
)

write.csv(
  high_confidence_counts,
  "results/tables/phase3_full_dataset_composition_high_confidence.csv",
  row.names = FALSE
)

write.csv(
  high_conf_state_counts,
  "results/tables/phase3_full_dataset_composition_high_confidence_by_state.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 30.17 Save Section 30 checkpoint
# ------------------------------------------------------------

phase3_full_dataset_composition <- list(
  
  composition = composition,
  
  overall_counts = overall_counts,
  
  state_counts = state_counts,
  
  sample_counts = sample_counts,
  
  state_percentage_matrix = state_percentage_matrix,
  
  sample_percentage_matrix = sample_percentage_matrix,
  
  high_confidence_counts = high_confidence_counts,
  
  high_conf_state_counts = high_conf_state_counts,
  
  metadata = list(
    total_cells = nrow(composition),
    total_cell_types =
      length(unique(composition$Transferred_Label)),
    sample_column = sample_col,
    donor_column = donor_col,
    clinical_state_column = state_col,
    primary_composition =
      "All 105220 transferred labels retained",
    high_confidence_definition =
      "Transfer_Confidence == High",
    high_confidence_cells =
      nrow(high_confidence),
    low_confidence_cells_retained =
      TRUE
  )
)

saveRDS(
  phase3_full_dataset_composition,
  "results/rds_objects/phase3_full_dataset_composition.rds"
)


# ------------------------------------------------------------
# 30.18 Final checkpoint validation
# ------------------------------------------------------------

stopifnot(
  file.exists(
    "results/rds_objects/phase3_full_dataset_composition.rds"
  )
)

stopifnot(
  sum(overall_counts$Cell_Count) == 105220
)

stopifnot(
  sum(state_counts$Cell_Count) == 105220
)

stopifnot(
  sum(sample_counts$Cell_Count) == 105220
)

stopifnot(
  sum(high_confidence_counts$High_Confidence_Cell_Count) ==
    high_confidence_n
)

cat("\n============================================================\n")
cat("SECTION 30 — FINAL CHECKPOINT SAVED\n")
cat("============================================================\n")
cat(
  "Checkpoint: results/rds_objects/phase3_full_dataset_composition.rds\n"
)
cat(
  "Primary cells: 105220\n"
)
cat(
  "High-confidence cells: ",
  high_confidence_n,
  "\n",
  sep = ""
)
cat(
  "Low-confidence cells retained in primary composition: TRUE\n"
)
cat("All final integrity checks PASSED.\n")
cat("============================================================\n")

# ============================================================

# SECTION 31 — FINAL UMAPs FOR THE FULL DATASET
# ============================================================

cat("\n============================================================\n")
cat("SECTION 31 — FINAL FULL-DATASET UMAPs\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 31.1 Load Section 25 full-dataset PCA projection
# ------------------------------------------------------------

seurat_obj <- readRDS(
  "results/rds_objects/phase3_full_dataset_projected.rds"
)

cat(
  "Loaded full projected Seurat object.\n"
)

stopifnot(
  ncol(seurat_obj) == 105220
)


# ------------------------------------------------------------
# 31.2 Load Section 27/28 label-transfer results
# ------------------------------------------------------------

knn_agreement <- readRDS(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)

low_agreement_flagged <- readRDS(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)

stopifnot(
  nrow(knn_agreement) == 105220
)

stopifnot(
  all(!is.na(knn_agreement$Transferred_Label))
)

cat(
  "Loaded full-dataset transferred labels.\n"
)


# ------------------------------------------------------------
# 31.3 Align transferred labels and confidence metadata
# ------------------------------------------------------------

flag_idx <- match(
  knn_agreement$Cell,
  low_agreement_flagged$Cell
)

stopifnot(
  !any(is.na(flag_idx))
)

knn_agreement$Low_Confidence <-
  low_agreement_flagged$Low_Confidence[flag_idx]

knn_agreement$Transfer_Confidence <-
  low_agreement_flagged$Transfer_Confidence[flag_idx]

stopifnot(
  !any(is.na(knn_agreement$Low_Confidence))
)

stopifnot(
  !any(is.na(knn_agreement$Transfer_Confidence))
)


# ------------------------------------------------------------
# 31.4 Add locked annotation metadata to full Seurat object
# ------------------------------------------------------------

cell_idx <- match(
  colnames(seurat_obj),
  knn_agreement$Cell
)

stopifnot(
  !any(is.na(cell_idx))
)

seurat_obj$Transferred_Label <-
  knn_agreement$Transferred_Label[cell_idx]

seurat_obj$Neighbour_Agreement <-
  knn_agreement$Neighbour_Agreement[cell_idx]

seurat_obj$Label_Margin <-
  knn_agreement$Label_Margin[cell_idx]

seurat_obj$Low_Confidence <-
  knn_agreement$Low_Confidence[cell_idx]

seurat_obj$Transfer_Confidence <-
  knn_agreement$Transfer_Confidence[cell_idx]


# ------------------------------------------------------------
# 31.5 Validate annotation alignment
# ------------------------------------------------------------

stopifnot(
  length(seurat_obj$Transferred_Label) == 105220
)

stopifnot(
  !any(is.na(seurat_obj$Transferred_Label))
)

high_confidence_n <- sum(
  seurat_obj$Transfer_Confidence == "High"
)

stopifnot(
  high_confidence_n > 0,
  high_confidence_n <= 105220
)

cat(
  "Transferred annotation metadata successfully aligned.\n"
)


# ------------------------------------------------------------
# 31.6 Confirm full-dataset PCA reduction
# ------------------------------------------------------------

if (!"pca.full" %in% Reductions(seurat_obj)) {
  stop(
    "pca.full reduction not found in the full dataset object."
  )
}

full_pca <- Embeddings(
  seurat_obj,
  reduction = "pca.full"
)

cat(
  "Full PCA dimensions: ",
  nrow(full_pca),
  " cells × ",
  ncol(full_pca),
  " PCs\n",
  sep = ""
)

stopifnot(
  nrow(full_pca) == 105220
)

stopifnot(
  ncol(full_pca) >= 20
)


# ------------------------------------------------------------
# 31.7 Compute full-dataset UMAP from locked PCA space
# ------------------------------------------------------------

cat(
  "\nComputing full-dataset UMAP from pca.full...\n"
)

set.seed(1234)

seurat_obj <- RunUMAP(
  object = seurat_obj,
  reduction = "pca.full",
  dims = 1:20,
  reduction.name = "umap.full",
  reduction.key = "UMAPFull_",
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  verbose = TRUE
)

cat(
  "Full-dataset UMAP successfully computed.\n"
)


# ------------------------------------------------------------
# 31.8 Validate UMAP
# ------------------------------------------------------------

full_umap <- Embeddings(
  seurat_obj,
  reduction = "umap.full"
)

cat(
  "Full UMAP dimensions: ",
  nrow(full_umap),
  " cells × ",
  ncol(full_umap),
  " dimensions\n",
  sep = ""
)

stopifnot(
  nrow(full_umap) == 105220
)

stopifnot(
  ncol(full_umap) >= 2
)

stopifnot(
  all(is.finite(full_umap[, 1]))
)

stopifnot(
  all(is.finite(full_umap[, 2]))
)


# ------------------------------------------------------------
# 31.9 Prepare UMAP plotting dataframe
# ------------------------------------------------------------

umap_df <- data.frame(
  Cell = rownames(full_umap),
  UMAP_1 = full_umap[, 1],
  UMAP_2 = full_umap[, 2],
  Transferred_Label =
    seurat_obj$Transferred_Label[
      match(
        rownames(full_umap),
        colnames(seurat_obj)
      )
    ],
  Clinical_State =
    seurat_obj$Phase[
      match(
        rownames(full_umap),
        colnames(seurat_obj)
      )
    ],
  Sample =
    seurat_obj$GSM[
      match(
        rownames(full_umap),
        colnames(seurat_obj)
      )
    ],
  Transfer_Confidence =
    seurat_obj$Transfer_Confidence[
      match(
        rownames(full_umap),
        colnames(seurat_obj)
      )
    ],
  Neighbour_Agreement =
    seurat_obj$Neighbour_Agreement[
      match(
        rownames(full_umap),
        colnames(seurat_obj)
      )
    ],
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(umap_df) == 105220
)

stopifnot(
  !any(is.na(umap_df$Transferred_Label))
)

stopifnot(
  !any(is.na(umap_df$Clinical_State))
)

stopifnot(
  !any(is.na(umap_df$Sample))
)


# ------------------------------------------------------------
# 31.10 Final UMAP — transferred cell type
# ------------------------------------------------------------

p_umap_celltype <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Transferred_Label
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — Transferred Cell Type",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_cell_type.png",
  p_umap_celltype,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.11 Final UMAP — clinical state
# ------------------------------------------------------------

p_umap_state <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Clinical_State
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — Clinical State",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Clinical State"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_clinical_state.png",
  p_umap_state,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.12 Final UMAP — donor/sample
# ------------------------------------------------------------

p_umap_sample <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Sample
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — Donor/Sample",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Sample"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_donor.png",
  p_umap_sample,
  width = 12,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.13 Final UMAP — transfer confidence
# ------------------------------------------------------------

p_umap_confidence <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Transfer_Confidence
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — Transfer Confidence",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Transfer Confidence"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_transfer_confidence.png",
  p_umap_confidence,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.14 Final UMAP — CD8_T highlighted
# ------------------------------------------------------------

umap_df$CD8_Highlight <- ifelse(
  umap_df$Transferred_Label == "CD8_T",
  "CD8_T",
  "Other"
)

p_umap_cd8 <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = CD8_Highlight
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — CD8 T Cells",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_cd8_highlighted.png",
  p_umap_cd8,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.15 Final UMAP — inflammatory myeloid highlighted
# ------------------------------------------------------------

umap_df$Myeloid_Highlight <- ifelse(
  umap_df$Transferred_Label == "Inflammatory_Myeloid",
  "Inflammatory_Myeloid",
  "Other"
)

p_umap_myeloid <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Myeloid_Highlight
  )
) +
  geom_point(
    size = 0.15,
    alpha = 0.55,
    stroke = 0
  ) +
  labs(
    title = "Full Dataset UMAP — Inflammatory Myeloid",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  "results/figures/phase3_final_umap_inflammatory_myeloid.png",
  p_umap_myeloid,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 31.16 Save UMAP coordinates
# ------------------------------------------------------------

write.csv(
  umap_df,
  "results/tables/phase3_full_dataset_umap_coordinates.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 31.17 Save final Section 31 checkpoint
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  "results/rds_objects/phase3_full_dataset_final_umap.rds"
)

saveRDS(
  umap_df,
  "results/rds_objects/phase3_full_dataset_umap_coordinates.rds"
)


# ------------------------------------------------------------
# 31.18 Final integrity checks
# ------------------------------------------------------------

stopifnot(
  file.exists(
    "results/rds_objects/phase3_full_dataset_final_umap.rds"
  )
)

stopifnot(
  file.exists(
    "results/rds_objects/phase3_full_dataset_umap_coordinates.rds"
  )
)

stopifnot(
  nrow(umap_df) == 105220
)

stopifnot(
  length(unique(umap_df$Transferred_Label)) >= 1
)

stopifnot(
  sum(
    umap_df$Transfer_Confidence == "High"
  ) == high_confidence_n
)

cat("\n============================================================\n")
cat("SECTION 31 COMPLETE — FINAL FULL-DATASET UMAPs SAVED\n")
cat("============================================================\n")

cat(
  "Cells: 105220\n"
)

cat(
  "Transferred cell types: ",
  length(unique(umap_df$Transferred_Label)),
  "\n",
  sep = ""
)

cat(
  "High-confidence cells: ",
  sum(umap_df$Transfer_Confidence == "High"),
  "\n",
  sep = ""
)

cat(
  "UMAP reduction: umap.full\n"
)

cat(
  "PCA input: pca.full, PCs 1:20\n"
)

cat(
  "Figures saved to: results/figures/\n"
)

cat(
  "Checkpoint: results/rds_objects/phase3_full_dataset_final_umap.rds\n"
)

cat("All integrity checks PASSED.\n")

cat("============================================================\n")

# ============================================================

# SECTION 32 — NEIGHBOUR AGREEMENT FIGURE
# ============================================================

cat("\n============================================================\n")
cat("SECTION 32 — NEIGHBOUR AGREEMENT FIGURE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 32.1 Load Section 27/28 label-transfer results
# ------------------------------------------------------------

knn_agreement <- readRDS(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)

low_agreement_flagged <- readRDS(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)

stopifnot(
  nrow(knn_agreement) == 105220
)

stopifnot(
  all(!is.na(knn_agreement$Neighbour_Agreement))
)

stopifnot(
  all(!is.na(knn_agreement$Transferred_Label))
)

cat(
  "Loaded neighbour agreement results for ",
  nrow(knn_agreement),
  " cells.\n",
  sep = ""
)


# ------------------------------------------------------------
# 32.2 Align confidence flags
# ------------------------------------------------------------

flag_idx <- match(
  knn_agreement$Cell,
  low_agreement_flagged$Cell
)

stopifnot(
  !any(is.na(flag_idx))
)

knn_agreement$Low_Confidence <-
  low_agreement_flagged$Low_Confidence[flag_idx]

knn_agreement$Transfer_Confidence <-
  low_agreement_flagged$Transfer_Confidence[flag_idx]

stopifnot(
  !any(is.na(knn_agreement$Low_Confidence))
)

cat(
  "Confidence flags successfully aligned.\n"
)


# ------------------------------------------------------------
# 32.3 Overall agreement summary
# ------------------------------------------------------------

agreement_summary <- data.frame(
  Metric = c(
    "N cells",
    "Mean",
    "Median",
    "Minimum",
    "Maximum",
    "SD",
    "Low agreement (<0.80)",
    "High agreement (>=0.80)",
    "Very high agreement (>=0.95)"
  ),
  Value = c(
    nrow(knn_agreement),
    mean(knn_agreement$Neighbour_Agreement),
    median(knn_agreement$Neighbour_Agreement),
    min(knn_agreement$Neighbour_Agreement),
    max(knn_agreement$Neighbour_Agreement),
    sd(knn_agreement$Neighbour_Agreement),
    sum(
      knn_agreement$Neighbour_Agreement < 0.80
    ),
    sum(
      knn_agreement$Neighbour_Agreement >= 0.80
    ),
    sum(
      knn_agreement$Neighbour_Agreement >= 0.95
    )
  )
)

cat("\n------------------------------------------------------------\n")
cat("NEIGHBOUR AGREEMENT SUMMARY\n")
cat("------------------------------------------------------------\n")

print(agreement_summary)


# ------------------------------------------------------------
# 32.4 Prepare plotting dataframe
# ------------------------------------------------------------

agreement_df <- data.frame(
  Cell = knn_agreement$Cell,
  Neighbour_Agreement =
    knn_agreement$Neighbour_Agreement,
  Transferred_Label =
    knn_agreement$Transferred_Label,
  Transfer_Confidence =
    knn_agreement$Transfer_Confidence,
  Low_Confidence =
    knn_agreement$Low_Confidence,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 32.5 Main neighbour agreement distribution
# ------------------------------------------------------------

p_agreement_distribution <- ggplot(
  agreement_df,
  aes(
    x = Neighbour_Agreement
  )
) +
  geom_histogram(
    bins = 30,
    boundary = 0,
    closed = "left"
  ) +
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 0.95,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0.80,
    y = Inf,
    label = "Low-agreement threshold = 0.80",
    vjust = 1.5,
    hjust = -0.05,
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0.95,
    y = Inf,
    label = "≥0.95",
    vjust = 3.5,
    hjust = -0.05,
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Full-Dataset Neighbour Agreement",
    subtitle = "KNN label-transfer agreement across 105,220 cells",
    x = "Neighbour Agreement",
    y = "Number of Cells"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  "results/figures/phase3_neighbour_agreement_distribution.png",
  p_agreement_distribution,
  width = 10,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 32.6 Agreement distribution by transferred cell type
# ------------------------------------------------------------

celltype_order <- names(
  sort(
    table(agreement_df$Transferred_Label),
    decreasing = TRUE
  )
)

agreement_df$Transferred_Label <- factor(
  agreement_df$Transferred_Label,
  levels = celltype_order
)

p_agreement_celltype <- ggplot(
  agreement_df,
  aes(
    x = Transferred_Label,
    y = Neighbour_Agreement
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Neighbour Agreement by Transferred Cell Type",
    subtitle = "Agreement distribution across the full dataset",
    x = "Transferred Cell Type",
    y = "Neighbour Agreement"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  "results/figures/phase3_neighbour_agreement_by_cell_type.png",
  p_agreement_celltype,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 32.7 Agreement categories
# ------------------------------------------------------------

agreement_df$Agreement_Category <- cut(
  agreement_df$Neighbour_Agreement,
  breaks = c(
    -Inf,
    0.50,
    0.80,
    0.95,
    Inf
  ),
  labels = c(
    "<50%",
    "50–<80%",
    "80–<95%",
    "≥95%"
  ),
  right = FALSE
)

agreement_category_counts <- as.data.frame(
  table(agreement_df$Agreement_Category),
  stringsAsFactors = FALSE
)

colnames(agreement_category_counts) <- c(
  "Agreement_Category",
  "Cell_Count"
)

agreement_category_counts$Percentage <- (
  agreement_category_counts$Cell_Count /
    sum(agreement_category_counts$Cell_Count)
) * 100

cat("\n------------------------------------------------------------\n")
cat("AGREEMENT CATEGORIES\n")
cat("------------------------------------------------------------\n")

print(agreement_category_counts)


# ------------------------------------------------------------
# 32.8 Agreement category figure
# ------------------------------------------------------------

p_agreement_categories <- ggplot(
  agreement_category_counts,
  aes(
    x = Agreement_Category,
    y = Percentage
  )
) +
  geom_col() +
  labs(
    title = "Neighbour Agreement Categories",
    x = "Neighbour Agreement",
    y = "Cells (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )

ggsave(
  "results/figures/phase3_neighbour_agreement_categories.png",
  p_agreement_categories,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 32.9 Save tables
# ------------------------------------------------------------

write.csv(
  agreement_summary,
  "results/tables/phase3_neighbour_agreement_summary.csv",
  row.names = FALSE
)

write.csv(
  agreement_category_counts,
  "results/tables/phase3_neighbour_agreement_categories.csv",
  row.names = FALSE
)

write.csv(
  agreement_df,
  "results/tables/phase3_neighbour_agreement_by_cell.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 32.10 Save Section 32 checkpoint
# ------------------------------------------------------------

phase3_neighbour_agreement_figure <- list(
  
  agreement_data = agreement_df,
  
  agreement_summary = agreement_summary,
  
  agreement_categories =
    agreement_category_counts,
  
  thresholds = list(
    low_agreement = 0.80,
    very_high_agreement = 0.95
  ),
  
  metadata = list(
    total_cells = nrow(agreement_df),
    source =
      "Section 27 KNN label-transfer agreement",
    confidence_flags =
      "Section 28",
    no_cells_removed = TRUE,
    no_labels_changed = TRUE
  )
)

saveRDS(
  phase3_neighbour_agreement_figure,
  "results/rds_objects/phase3_neighbour_agreement_figure.rds"
)


# ------------------------------------------------------------
# 32.11 Final integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(agreement_df) == 105220
)

stopifnot(
  sum(agreement_category_counts$Cell_Count) == 105220
)

low_agreement_n <- sum(
  agreement_df$Neighbour_Agreement < 0.80
)

high_agreement_n <- sum(
  agreement_df$Neighbour_Agreement >= 0.80
)

very_high_agreement_n <- sum(
  agreement_df$Neighbour_Agreement >= 0.95
)

stopifnot(
  low_agreement_n + high_agreement_n == 105220
)

stopifnot(
  very_high_agreement_n <= high_agreement_n
)

cat("\n============================================================\n")
cat("SECTION 32 COMPLETE — NEIGHBOUR AGREEMENT FIGURES SAVED\n")
cat("============================================================\n")

cat(
  "Cells: 105220\n"
)

cat(
  "Low agreement (<0.80): ",
  low_agreement_n,
  "\n",
  sep = ""
)

cat(
  "High agreement (>=0.80): ",
  high_agreement_n,
  "\n",
  sep = ""
)

cat(
  "Very high agreement (>=0.95): ",
  very_high_agreement_n,
  "\n",
  sep = ""
)

cat(
  "Figures saved to: results/figures/\n"
)

cat(
  "Checkpoint: results/rds_objects/phase3_neighbour_agreement_figure.rds\n"
)

cat("All integrity checks PASSED.\n")

cat("============================================================\n")

# ============================================================

# SECTION 33 — LOW AGREEMENT BY PHASE
# ============================================================

cat("\n============================================================\n")
cat("SECTION 33 — LOW AGREEMENT BY PHASE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 33.1 Load Section 27/28 agreement results
# ------------------------------------------------------------

knn_agreement <- readRDS(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)

low_agreement_flagged <- readRDS(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)

stopifnot(
  nrow(knn_agreement) == 105220
)

stopifnot(
  all(!is.na(knn_agreement$Neighbour_Agreement))
)

stopifnot(
  all(!is.na(knn_agreement$Transferred_Label))
)

cat(
  "Loaded agreement results for ",
  nrow(knn_agreement),
  " cells.\n",
  sep = ""
)


# ------------------------------------------------------------
# 33.2 Align Section 28 confidence flags
# ------------------------------------------------------------

flag_idx <- match(
  knn_agreement$Cell,
  low_agreement_flagged$Cell
)

stopifnot(
  !any(is.na(flag_idx))
)

knn_agreement$Low_Confidence <-
  low_agreement_flagged$Low_Confidence[flag_idx]

knn_agreement$Transfer_Confidence <-
  low_agreement_flagged$Transfer_Confidence[flag_idx]

stopifnot(
  !any(is.na(knn_agreement$Low_Confidence))
)

cat(
  "Section 28 confidence flags successfully aligned.\n"
)


# ------------------------------------------------------------
# 33.3 Load full projected dataset for Phase metadata
# ------------------------------------------------------------

seurat_obj <- readRDS(
  "results/rds_objects/phase3_full_dataset_projected.rds"
)

stopifnot(
  ncol(seurat_obj) == 105220
)

meta <- seurat_obj@meta.data

stopifnot(
  "Phase" %in% colnames(meta)
)

cat(
  "Clinical phase metadata found in column: Phase\n"
)


# ------------------------------------------------------------
# 33.4 Construct phase-level agreement table
# ------------------------------------------------------------

phase_df <- data.frame(
  Cell = knn_agreement$Cell,
  Phase = as.character(
    meta[
      match(
        knn_agreement$Cell,
        rownames(meta)
      ),
      "Phase"
    ]
  ),
  Neighbour_Agreement =
    knn_agreement$Neighbour_Agreement,
  Transferred_Label =
    knn_agreement$Transferred_Label,
  Low_Agreement =
    knn_agreement$Neighbour_Agreement < 0.80,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(phase_df) == 105220
)

stopifnot(
  !any(is.na(phase_df$Phase))
)

stopifnot(
  !any(is.na(phase_df$Neighbour_Agreement))
)

stopifnot(
  !any(is.na(phase_df$Transferred_Label))
)


# ------------------------------------------------------------
# 33.5 Validate the locked low-agreement definition
# ------------------------------------------------------------

low_agreement_n <- sum(phase_df$Low_Agreement)

stopifnot(
  low_agreement_n >= 0,
  low_agreement_n <= 105220
)

cat(
  "Low-agreement definition validated: ",
  sum(phase_df$Low_Agreement),
  " cells (<0.80).\n",
  sep = ""
)


# ------------------------------------------------------------
# 33.6 Overall phase cell counts
# ------------------------------------------------------------

phase_total_counts <- as.data.frame(
  table(phase_df$Phase),
  stringsAsFactors = FALSE
)

colnames(phase_total_counts) <- c(
  "Phase",
  "Total_Cell_Count"
)

phase_total_counts <- phase_total_counts[
  order(phase_total_counts$Phase),
]

rownames(phase_total_counts) <- NULL

cat("\n------------------------------------------------------------\n")
cat("TOTAL CELLS BY PHASE\n")
cat("------------------------------------------------------------\n")

print(phase_total_counts)

stopifnot(
  sum(phase_total_counts$Total_Cell_Count) == 105220
)


# ------------------------------------------------------------
# 33.7 Low-agreement counts by phase
# ------------------------------------------------------------

low_agreement_phase_counts <- as.data.frame(
  table(
    phase_df$Phase[
      phase_df$Low_Agreement
    ]
  ),
  stringsAsFactors = FALSE
)

colnames(low_agreement_phase_counts) <- c(
  "Phase",
  "Low_Agreement_Cell_Count"
)

low_agreement_phase_counts <- low_agreement_phase_counts[
  order(low_agreement_phase_counts$Phase),
]

rownames(low_agreement_phase_counts) <- NULL


# ------------------------------------------------------------
# 33.8 Combine total and low-agreement counts
# ------------------------------------------------------------

low_agreement_by_phase <- merge(
  phase_total_counts,
  low_agreement_phase_counts,
  by = "Phase",
  all = TRUE
)

low_agreement_by_phase[
  is.na(low_agreement_by_phase)
] <- 0

low_agreement_by_phase$High_Agreement_Cell_Count <-
  low_agreement_by_phase$Total_Cell_Count -
  low_agreement_by_phase$Low_Agreement_Cell_Count

low_agreement_by_phase$Low_Agreement_Percentage <- (
  low_agreement_by_phase$Low_Agreement_Cell_Count /
    low_agreement_by_phase$Total_Cell_Count
) * 100

low_agreement_by_phase$High_Agreement_Percentage <- (
  low_agreement_by_phase$High_Agreement_Cell_Count /
    low_agreement_by_phase$Total_Cell_Count
) * 100

low_agreement_by_phase <- low_agreement_by_phase[
  order(low_agreement_by_phase$Phase),
]

rownames(low_agreement_by_phase) <- NULL


# ------------------------------------------------------------
# 33.9 Print final phase-level summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LOW AGREEMENT BY PHASE\n")
cat("------------------------------------------------------------\n")

print(
  low_agreement_by_phase
)


# ------------------------------------------------------------
# 33.10 Integrity checks
# ------------------------------------------------------------

stopifnot(
  sum(
    low_agreement_by_phase$Total_Cell_Count
  ) == 105220
)

stopifnot(
  sum(
    low_agreement_by_phase$Low_Agreement_Cell_Count
  ) == low_agreement_n
)

stopifnot(
  sum(
    low_agreement_by_phase$High_Agreement_Cell_Count
  ) ==
    105220 - low_agreement_n
)

stopifnot(
  all(
    low_agreement_by_phase$Total_Cell_Count ==
      low_agreement_by_phase$Low_Agreement_Cell_Count +
      low_agreement_by_phase$High_Agreement_Cell_Count
  )
)

cat(
  "\nPhase-level integrity checks PASSED.\n"
)


# ------------------------------------------------------------
# 33.11 Figure — low agreement percentage by phase
# ------------------------------------------------------------

p_low_agreement_phase <- ggplot(
  low_agreement_by_phase,
  aes(
    x = Phase,
    y = Low_Agreement_Percentage
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        Low_Agreement_Percentage
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Low-Agreement Assignments by Clinical Phase",
    subtitle = "Neighbour Agreement < 0.80",
    x = "Clinical Phase",
    y = "Low-Agreement Cells (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  )

ggsave(
  "results/figures/phase3_low_agreement_by_phase.png",
  p_low_agreement_phase,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 33.12 Low agreement by phase and transferred cell type
# ------------------------------------------------------------

low_agreement_celltype_phase <- as.data.frame(
  table(
    phase_df$Phase,
    phase_df$Transferred_Label,
    phase_df$Low_Agreement
  ),
  stringsAsFactors = FALSE
)

colnames(low_agreement_celltype_phase) <- c(
  "Phase",
  "Cell_Type",
  "Low_Agreement",
  "Cell_Count"
)

low_agreement_celltype_phase <- low_agreement_celltype_phase[
  low_agreement_celltype_phase$Low_Agreement == TRUE,
]

low_agreement_celltype_phase$Low_Agreement <- NULL

rownames(low_agreement_celltype_phase) <- NULL


# ------------------------------------------------------------
# 33.13 Cell-type × phase total counts
# ------------------------------------------------------------

phase_celltype_total <- as.data.frame(
  table(
    phase_df$Phase,
    phase_df$Transferred_Label
  ),
  stringsAsFactors = FALSE
)

colnames(phase_celltype_total) <- c(
  "Phase",
  "Cell_Type",
  "Total_Cell_Count"
)

low_agreement_celltype_phase <- merge(
  phase_celltype_total,
  low_agreement_celltype_phase,
  by = c(
    "Phase",
    "Cell_Type"
  ),
  all.x = TRUE
)

low_agreement_celltype_phase$Cell_Count[
  is.na(
    low_agreement_celltype_phase$Cell_Count
  )
] <- 0

colnames(
  low_agreement_celltype_phase
)[
  colnames(
    low_agreement_celltype_phase
  ) == "Cell_Count"
] <- "Low_Agreement_Cell_Count"

low_agreement_celltype_phase$Low_Agreement_Percentage <- (
  low_agreement_celltype_phase$Low_Agreement_Cell_Count /
    low_agreement_celltype_phase$Total_Cell_Count
) * 100

low_agreement_celltype_phase <- low_agreement_celltype_phase[
  order(
    low_agreement_celltype_phase$Phase,
    -low_agreement_celltype_phase$Low_Agreement_Percentage
  ),
]

rownames(low_agreement_celltype_phase) <- NULL


# ------------------------------------------------------------
# 33.14 Save tables
# ------------------------------------------------------------

write.csv(
  low_agreement_by_phase,
  "results/tables/phase3_low_agreement_by_phase.csv",
  row.names = FALSE
)

write.csv(
  low_agreement_celltype_phase,
  "results/tables/phase3_low_agreement_by_phase_and_cell_type.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 33.15 Save Section 33 checkpoint
# ------------------------------------------------------------

phase3_low_agreement_by_phase <- list(
  
  phase_summary =
    low_agreement_by_phase,
  
  phase_celltype_summary =
    low_agreement_celltype_phase,
  
  definition = list(
    low_agreement_threshold = 0.80,
    low_agreement_rule =
      "Neighbour_Agreement < 0.80"
  ),
  
  metadata = list(
    total_cells = nrow(phase_df),
    total_low_agreement = low_agreement_n,
    total_high_agreement = sum(!phase_df$Low_Agreement),
    phase_column = "Phase",
    labels_changed = FALSE,
    cells_removed = FALSE
  )
)

saveRDS(
  phase3_low_agreement_by_phase,
  "results/rds_objects/phase3_low_agreement_by_phase.rds"
)


# ------------------------------------------------------------
# 33.16 Final file/integrity validation
# ------------------------------------------------------------

stopifnot(
  file.exists(
    "results/rds_objects/phase3_low_agreement_by_phase.rds"
  )
)

stopifnot(
  file.exists(
    "results/figures/phase3_low_agreement_by_phase.png"
  )
)

stopifnot(
  sum(
    low_agreement_by_phase$Low_Agreement_Cell_Count
  ) == low_agreement_n
)

cat("\n============================================================\n")
cat("SECTION 33 COMPLETE — LOW AGREEMENT BY PHASE\n")
cat("============================================================\n")

cat(
  "Total cells: 105220\n"
)

cat(
  "Low-agreement cells (<0.80): ",
  low_agreement_n,
  "\n",
  sep = ""
)

cat(
  "High-agreement cells (>=0.80): ",
  105220 - low_agreement_n,
  "\n",
  sep = ""
)

cat(
  "Figure: results/figures/phase3_low_agreement_by_phase.png\n"
)

cat(
  "Checkpoint: results/rds_objects/phase3_low_agreement_by_phase.rds\n"
)

cat("All integrity checks PASSED.\n")

cat("============================================================\n")

# ============================================================

# SECTION 34 — FINAL METADATA, FINAL OBJECT & INTEGRITY CHECK
# ============================================================

cat("\n============================================================\n")
cat("SECTION 34 — FINAL METADATA, FINAL OBJECT & INTEGRITY CHECK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 34.1 Load the final full-dataset object
# ------------------------------------------------------------

seurat_obj <- readRDS(
  "results/rds_objects/phase3_full_dataset_final_umap.rds"
)

cat("Final full-dataset Seurat object loaded.\n")


# ------------------------------------------------------------
# 34.2 Load final label-transfer results
# ------------------------------------------------------------

knn_agreement <- readRDS(
  "results/rds_objects/phase3_knn_label_transfer_agreement.rds"
)

low_agreement_flagged <- readRDS(
  "results/rds_objects/phase3_low_agreement_flagged.rds"
)

cat("Label-transfer and confidence objects loaded.\n")


# ------------------------------------------------------------
# 34.3 Basic object integrity
# ------------------------------------------------------------

stopifnot(
  ncol(seurat_obj) == 105220
)

stopifnot(
  nrow(knn_agreement) == 105220
)

stopifnot(
  all(!is.na(knn_agreement$Transferred_Label))
)

cat(
  "Cell count: ",
  ncol(seurat_obj),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.4 Validate required dimensional reductions
# ------------------------------------------------------------

required_reductions <- c(
  "pca.full",
  "umap.full"
)

missing_reductions <- setdiff(
  required_reductions,
  Reductions(seurat_obj)
)

if (length(missing_reductions) > 0) {
  stop(
    "Missing required dimensional reductions: ",
    paste(missing_reductions, collapse = ", ")
  )
}

cat(
  "Required reductions present: ",
  paste(required_reductions, collapse = ", "),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.5 Align and validate final transferred metadata
# ------------------------------------------------------------

cell_idx <- match(
  colnames(seurat_obj),
  knn_agreement$Cell
)

stopifnot(
  !any(is.na(cell_idx))
)

seurat_obj$Transferred_Label <-
  knn_agreement$Transferred_Label[cell_idx]

seurat_obj$Neighbour_Agreement <-
  knn_agreement$Neighbour_Agreement[cell_idx]

seurat_obj$Label_Margin <-
  knn_agreement$Label_Margin[cell_idx]


# ------------------------------------------------------------
# 34.6 Align Section 28 confidence flags
# ------------------------------------------------------------

flag_idx <- match(
  colnames(seurat_obj),
  low_agreement_flagged$Cell
)

stopifnot(
  !any(is.na(flag_idx))
)

seurat_obj$Low_Confidence <-
  low_agreement_flagged$Low_Confidence[flag_idx]

seurat_obj$Transfer_Confidence <-
  low_agreement_flagged$Transfer_Confidence[flag_idx]


# ------------------------------------------------------------
# 34.7 Validate final annotation metadata
# ------------------------------------------------------------

required_final_metadata <- c(
  "Transferred_Label",
  "Neighbour_Agreement",
  "Label_Margin",
  "Low_Confidence",
  "Transfer_Confidence",
  "Phase",
  "GSM",
  "Donor"
)

missing_metadata <- setdiff(
  required_final_metadata,
  colnames(seurat_obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required final metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

stopifnot(
  !any(is.na(seurat_obj$Transferred_Label))
)

stopifnot(
  !any(is.na(seurat_obj$Neighbour_Agreement))
)

stopifnot(
  !any(is.na(seurat_obj$Label_Margin))
)

stopifnot(
  !any(is.na(seurat_obj$Low_Confidence))
)

stopifnot(
  !any(is.na(seurat_obj$Transfer_Confidence))
)

stopifnot(
  !any(is.na(seurat_obj$Phase))
)

stopifnot(
  !any(is.na(seurat_obj$GSM))
)


# ------------------------------------------------------------
# 34.8 Verify transferred labels against Section 27
# ------------------------------------------------------------

stopifnot(
  identical(
    as.character(seurat_obj$Transferred_Label),
    as.character(
      knn_agreement$Transferred_Label[cell_idx]
    )
  )
)

stopifnot(
  identical(
    as.numeric(seurat_obj$Neighbour_Agreement),
    as.numeric(
      knn_agreement$Neighbour_Agreement[cell_idx]
    )
  )
)

stopifnot(
  identical(
    as.numeric(seurat_obj$Label_Margin),
    as.numeric(
      knn_agreement$Label_Margin[cell_idx]
    )
  )
)


# ------------------------------------------------------------
# 34.9 Verify confidence flags against Section 28
# ------------------------------------------------------------

stopifnot(
  identical(
    as.logical(seurat_obj$Low_Confidence),
    as.logical(
      low_agreement_flagged$Low_Confidence[flag_idx]
    )
  )
)

stopifnot(
  identical(
    as.character(seurat_obj$Transfer_Confidence),
    as.character(
      low_agreement_flagged$Transfer_Confidence[flag_idx]
    )
  )
)


# ------------------------------------------------------------
# 34.10 Verify locked low-agreement definition
# ------------------------------------------------------------

expected_low_agreement <- (
  seurat_obj$Neighbour_Agreement < 0.80
)

stopifnot(
  identical(
    as.logical(seurat_obj$Low_Confidence),
    as.logical(expected_low_agreement)
  )
)

cat(
  "Low-agreement definition verified: Neighbour_Agreement < 0.80\n"
)


# ------------------------------------------------------------
# 34.11 Verify confidence totals
# ------------------------------------------------------------

high_confidence_n <- sum(
  seurat_obj$Transfer_Confidence == "High"
)

low_confidence_n <- sum(
  seurat_obj$Transfer_Confidence != "High"
)

low_agreement_n <- sum(
  seurat_obj$Neighbour_Agreement < 0.80
)

stopifnot(
  high_confidence_n + low_confidence_n == 105220
)

stopifnot(
  low_agreement_n >= 0,
  low_agreement_n <= 105220
)

stopifnot(
  high_confidence_n + low_confidence_n == 105220
)


# ------------------------------------------------------------
# 34.12 Verify cell-type composition
# ------------------------------------------------------------

final_label_counts <- table(
  seurat_obj$Transferred_Label
)

stopifnot(
  sum(final_label_counts) == 105220
)

stopifnot(
  length(final_label_counts) >= 1
)

cat(
  "Transferred cell types: ",
  length(final_label_counts),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.13 Verify expected locked labels
# ------------------------------------------------------------

expected_labels <- sort(
  unique(
    as.character(
      knn_agreement$Transferred_Label
    )
  )
)

stopifnot(
  setequal(
    names(final_label_counts),
    expected_labels
  )
)


# ------------------------------------------------------------
# 34.14 Verify target populations
# ------------------------------------------------------------

final_cd8_n <- sum(
  seurat_obj$Transferred_Label == "CD8_T"
)

final_myeloid_n <- sum(
  seurat_obj$Transferred_Label == "Inflammatory_Myeloid"
)

stopifnot(
  final_cd8_n > 0
)

stopifnot(
  final_myeloid_n > 0
)

cat(
  "CD8_T cells: ",
  final_cd8_n,
  "\n",
  sep = ""
)

cat(
  "Inflammatory_Myeloid cells: ",
  final_myeloid_n,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.15 Verify Phase composition
# ------------------------------------------------------------

phase_counts <- table(
  seurat_obj$Phase
)

stopifnot(
  sum(phase_counts) == 105220
)

expected_phases <- c(
  "NL",
  "IT",
  "IA",
  "AR",
  "AC"
)

stopifnot(
  setequal(
    names(phase_counts),
    expected_phases
  )
)

cat(
  "Clinical phases present: ",
  paste(
    names(sort(phase_counts)),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.16 Verify donor/sample composition
# ------------------------------------------------------------

sample_counts_final <- table(
  seurat_obj$GSM
)

stopifnot(
  sum(sample_counts_final) == 105220
)

stopifnot(
  length(sample_counts_final) == 23
)

cat(
  "Donor/sample count: ",
  length(sample_counts_final),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.17 Verify full PCA
# ------------------------------------------------------------

full_pca <- Embeddings(
  seurat_obj,
  reduction = "pca.full"
)

stopifnot(
  nrow(full_pca) == 105220
)

stopifnot(
  ncol(full_pca) >= 20
)

stopifnot(
  all(is.finite(full_pca))
)

cat(
  "Full PCA: ",
  nrow(full_pca),
  " × ",
  ncol(full_pca),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.18 Verify final UMAP
# ------------------------------------------------------------

full_umap <- Embeddings(
  seurat_obj,
  reduction = "umap.full"
)

stopifnot(
  nrow(full_umap) == 105220
)

stopifnot(
  ncol(full_umap) >= 2
)

stopifnot(
  all(is.finite(full_umap))
)

cat(
  "Full UMAP: ",
  nrow(full_umap),
  " × ",
  ncol(full_umap),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.19 Create final metadata table
# ------------------------------------------------------------

final_metadata <- seurat_obj@meta.data

final_metadata$Cell <- rownames(
  final_metadata
)

# Put Cell first
final_metadata <- final_metadata[
  ,
  c(
    "Cell",
    setdiff(
      colnames(final_metadata),
      "Cell"
    )
  ),
  drop = FALSE
]

# Add final annotation fields in a clearly grouped position
final_annotation_columns <- c(
  "Transferred_Label",
  "Neighbour_Agreement",
  "Label_Margin",
  "Low_Confidence",
  "Transfer_Confidence"
)

existing_annotation_columns <- intersect(
  final_annotation_columns,
  colnames(final_metadata)
)

other_columns <- setdiff(
  colnames(final_metadata),
  c(
    "Cell",
    existing_annotation_columns
  )
)

final_metadata <- final_metadata[
  ,
  c(
    "Cell",
    existing_annotation_columns,
    other_columns
  ),
  drop = FALSE
]

stopifnot(
  nrow(final_metadata) == 105220
)

stopifnot(
  !anyDuplicated(final_metadata$Cell)
)


# ------------------------------------------------------------
# 34.20 Save final metadata table
# ------------------------------------------------------------

write.csv(
  final_metadata,
  "results/tables/phase3_final_metadata.csv",
  row.names = FALSE
)

cat(
  "Final metadata table saved:\n",
  "results/tables/phase3_final_metadata.csv\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.21 Save final Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  "results/rds_objects/phase3_final_full_dataset.rds"
)

cat(
  "Final Seurat object saved:\n",
  "results/rds_objects/phase3_final_full_dataset.rds\n",
  sep = ""
)


# ------------------------------------------------------------
# 34.22 Build final integrity report
# ------------------------------------------------------------

final_integrity_report <- data.frame(
  Check = c(
    "Total cells",
    "Unique cells",
    "Missing transferred labels",
    "Transferred cell types",
    "High-confidence cells",
    "Low-confidence cells",
    "Low-agreement cells",
    "CD8_T cells",
    "Inflammatory_Myeloid cells",
    "Clinical phases",
    "Donor/sample count",
    "PCA rows",
    "PCA finite",
    "UMAP rows",
    "UMAP finite",
    "Metadata rows"
  ),
  Value = c(
    ncol(seurat_obj),
    length(unique(colnames(seurat_obj))),
    sum(is.na(seurat_obj$Transferred_Label)),
    length(unique(seurat_obj$Transferred_Label)),
    high_confidence_n,
    low_confidence_n,
    low_agreement_n,
    final_cd8_n,
    final_myeloid_n,
    length(unique(seurat_obj$Phase)),
    length(unique(seurat_obj$GSM)),
    nrow(full_pca),
    all(is.finite(full_pca)),
    nrow(full_umap),
    all(is.finite(full_umap)),
    nrow(final_metadata)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  final_integrity_report,
  "results/tables/phase3_final_integrity_check.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 34.23 Save final integrity checkpoint
# ------------------------------------------------------------

final_integrity <- list(
  
  status = "PASSED",
  
  total_cells = 105220,
  
  transferred_cell_types =
    sort(
      table(
        seurat_obj$Transferred_Label
      ),
      decreasing = TRUE
    ),
  
  confidence = list(
    high = high_confidence_n,
    low = low_confidence_n,
    low_agreement = low_agreement_n,
    low_agreement_threshold = 0.80
  ),
  
  target_populations = list(
    CD8_T = final_cd8_n,
    Inflammatory_Myeloid = final_myeloid_n
  ),
  
  clinical_phases =
    table(seurat_obj$Phase),
  
  donor_samples =
    table(seurat_obj$GSM),
  
  reductions = list(
    pca.full = dim(full_pca),
    umap.full = dim(full_umap)
  ),
  
  files = list(
    final_metadata =
      "results/tables/phase3_final_metadata.csv",
    final_object =
      "results/rds_objects/phase3_final_full_dataset.rds",
    integrity_report =
      "results/tables/phase3_final_integrity_check.csv"
  )
)

saveRDS(
  final_integrity,
  "results/rds_objects/phase3_final_integrity_check.rds"
)


# ------------------------------------------------------------
# 34.24 Final file validation
# ------------------------------------------------------------

required_final_files <- c(
  "results/tables/phase3_final_metadata.csv",
  "results/rds_objects/phase3_final_full_dataset.rds",
  "results/tables/phase3_final_integrity_check.csv",
  "results/rds_objects/phase3_final_integrity_check.rds"
)

stopifnot(
  all(
    file.exists(required_final_files)
  )
)


# ------------------------------------------------------------
# 34.25 FINAL INTEGRITY CHECK
# ------------------------------------------------------------

stopifnot(
  ncol(seurat_obj) == 105220
)

stopifnot(
  nrow(final_metadata) == 105220
)

stopifnot(
  length(unique(final_metadata$Cell)) == 105220
)

stopifnot(
  sum(is.na(seurat_obj$Transferred_Label)) == 0
)

stopifnot(
  length(unique(seurat_obj$Transferred_Label)) ==
    length(expected_labels)
)

stopifnot(
  high_confidence_n + low_confidence_n == 105220
)

stopifnot(
  low_agreement_n >= 0,
  low_agreement_n <= 105220
)

stopifnot(
  final_cd8_n > 0
)

stopifnot(
  final_myeloid_n > 0
)

stopifnot(
  length(unique(seurat_obj$Phase)) == 5
)

stopifnot(
  length(unique(seurat_obj$GSM)) == 23
)

stopifnot(
  nrow(full_pca) == 105220
)

stopifnot(
  all(is.finite(full_pca))
)

stopifnot(
  nrow(full_umap) == 105220
)

stopifnot(
  all(is.finite(full_umap))
)

cat("\n============================================================\n")
cat("FINAL INTEGRITY CHECK — PASSED\n")
cat("============================================================\n")

cat("\nFINAL FILES:\n")

cat(
  "1. results/tables/phase3_final_metadata.csv\n"
)

cat(
  "2. results/rds_objects/phase3_final_full_dataset.rds\n"
)

cat(
  "3. results/tables/phase3_final_integrity_check.csv\n"
)

cat(
  "4. results/rds_objects/phase3_final_integrity_check.rds\n"
)

cat("\n============================================================\n")
cat("PHASE 3 GLOBAL CELL ATLAS — FINALIZED\n")
cat("============================================================\n")

