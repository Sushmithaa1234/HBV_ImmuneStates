# ============================================================================
# HBV Immune States scRNA-seq Analysis
# Script 03: Global Cell Atlas and Cell-Type Annotation
#
# Purpose:
#   1. Load the QC-filtered singlet dataset
#   2. Normalize the data and identify variable features
#   3. Construct a representative leverage-score sketch
#   4. Perform PCA, neighbor graph construction, clustering and UMAP
#   5. Evaluate donor/sample contribution to clusters
#   6. Identify cluster marker genes
#   7. Define and inspect canonical lineage marker programs
#   8. Pause for manual biological cluster annotation
#   9. Project the full dataset into sketch-derived space
#  10. Transfer cluster/cell-type labels using KNN neighbor agreement
#  11. Flag low-agreement assignments
#  12. Summarize cell-type composition by donor and clinical state
#
# IMPORTANT:
#   - Predicted doublets were evaluated during Phase 2 and removed before
#     atlas construction.
#   - This analysis begins with the QC-filtered singlet dataset.
#   - No batch correction is applied automatically.
#   - Donor/sample effects are evaluated before deciding whether integration
#     is required for downstream analyses.
#   - Manual annotation is based on cluster markers, canonical marker programs,
#     and donor/sample representation.
#   - KNN transfer agreement is not a calibrated biological probability.
#
# WORKFLOW:
#   PART A: Run Sections 0-21B.
#           Review outputs and complete phase3_cluster_annotation.csv.
#   PART B: Resume at Section 22.
# ============================================================================

# ============================================================================
# SECTION 0 — SETUP
# ============================================================================

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
library(patchwork)
library(Matrix)
library(pheatmap)
library(FNN)

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

cat("============================================================\n")
cat("PHASE 3 — GLOBAL CELL ATLAS + CELL-TYPE ANNOTATION\n")
cat("============================================================\n\n")

# ============================================================================
# SECTION 1 — LOAD QC-FILTERED SINGLET DATASET
# ============================================================================

cat("=== SECTION 1: LOADING QC-FILTERED SINGLET DATASET ===\n\n")

seurat_obj <- readRDS(
  "results/rds_objects/seurat_qc_singlets.rds"
)

DefaultAssay(seurat_obj) <- "RNA"

cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")
cat("Samples:", n_distinct(seurat_obj$GSM), "\n")
cat("Donors:", n_distinct(seurat_obj$Donor), "\n")
cat("Clinical states:", n_distinct(seurat_obj$Phase), "\n\n")

if (ncol(seurat_obj) != 105147) {
  warning(
    "Input object contains ", ncol(seurat_obj),
    " cells rather than the expected 105147 singlets. Verify the input file."
  )
}

# --------------------------------------------------------------------------
# Required metadata check
# --------------------------------------------------------------------------

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
    "ERROR: Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
  
}

cat("✓ Required metadata columns present\n\n")

# --------------------------------------------------------------------------
# Dataset snapshot
# --------------------------------------------------------------------------

dataset_snapshot <- seurat_obj@meta.data %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    Cells = n(),
    .groups = "drop"
  )

print(
  dataset_snapshot,
  n = Inf
)

write_csv(
  dataset_snapshot,
  "results/tables/phase3_dataset_snapshot.csv"
)

cat("\n✓ Dataset snapshot saved\n\n")

# ============================================================================
# SECTION 2 — INPUT OBJECT STATUS
# ============================================================================

cat("=== SECTION 2: INPUT OBJECT STATUS ===\n\n")

cat(
  "This object contains QC-filtered cells after removal of predicted\n",
  "scDblFinder doublets during Phase 2.\n\n"
)

if ("scDblFinder.class" %in% colnames(seurat_obj@meta.data)) {
  
  cat("Retained scDblFinder metadata:\n")
  
  print(
    table(
      seurat_obj$scDblFinder.class,
      useNA = "ifany"
    )
  )
  
  cat("\n")
  
}

# ============================================================================
# SECTION 3 — NORMALIZATION
# ============================================================================

cat("=== SECTION 3: NORMALIZATION ===\n\n")

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("✓ Log-normalization complete\n\n")

# ============================================================================
# SECTION 4 — HIGHLY VARIABLE FEATURES
# ============================================================================

cat("=== SECTION 4: HIGHLY VARIABLE FEATURES ===\n\n")

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

hvg <- VariableFeatures(seurat_obj)

cat(
  "Highly variable genes:",
  length(hvg),
  "\n\n"
)

write_csv(
  data.frame(
    Gene = hvg
  ),
  "results/tables/phase3_highly_variable_genes.csv"
)

cat("✓ HVG list saved\n\n")

# ============================================================================
# SECTION 5 — LEVERAGE-SCORE SKETCH
# ============================================================================

cat("=== SECTION 5: REPRESENTATIVE SKETCH CONSTRUCTION ===\n\n")

target_sketch_cells <- 20000

cat(
  "Target sketch size:",
  target_sketch_cells,
  "cells\n"
)

sample_sizes <- seurat_obj@meta.data %>%
  
  count(
    GSM,
    Donor,
    Phase,
    name = "Original_Cells"
  )

print.data.frame(
  sample_sizes,
  row.names = FALSE
)

write_csv(
  sample_sizes,
  "results/tables/phase3_sample_sizes_before_sketch.csv"
)

sketch_obj <- SketchData(
  object = seurat_obj,
  ncells = target_sketch_cells,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(sketch_obj) <- "sketch"

sketch_cells <- colnames(
  sketch_obj[["sketch"]]
)

cat("\nSketch assay created.\n")

cat(
  "Sketch cells:",
  length(sketch_cells),
  "\n\n"
)

if (length(sketch_cells) > target_sketch_cells) {
  stop(
    "ERROR: Sketch contains more cells than the requested target."
  )
}

# ============================================================================
# SECTION 6 — SKETCH REPRESENTATION BY SAMPLE
# ============================================================================

cat("=== SECTION 6: SKETCH REPRESENTATION BY SAMPLE ===\n\n")

sketch_metadata <- sketch_obj@meta.data %>%
  
  rownames_to_column("Cell") %>%
  
  mutate(
    In_Sketch = Cell %in% sketch_cells
  )

sketch_representation <- sketch_metadata %>%
  
  group_by(
    GSM,
    Donor,
    Phase
  ) %>%
  
  summarise(
    
    Original_Cells = n(),
    
    Sketch_Cells =
      sum(In_Sketch),
    
    Percent_Represented =
      round(
        100 *
          Sketch_Cells /
          Original_Cells,
        2
      ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Sketch_Cells
  )

print(
  sketch_representation,
  n = Inf
)

write_csv(
  sketch_representation,
  "results/tables/phase3_sketch_representation_by_sample.csv"
)

zero_representation <- sketch_representation %>%
  filter(Sketch_Cells == 0)

if (nrow(zero_representation) > 0) {
  
  cat("\n⚠ Samples with ZERO sketch representation:\n\n")
  
  print(
    zero_representation,
    n = Inf
  )
  
} else {
  
  cat("\n✓ Every sample has at least one sketch cell.\n")
  
}

small_samples <- sketch_representation %>%
  
  filter(
    Original_Cells < 500
  )

if (nrow(small_samples) > 0) {
  
  cat(
    "\n⚠ Small samples detected (<500 cells):\n"
  )
  
  print(
    small_samples,
    n = Inf
  )
  
} else {
  
  cat(
    "\n✓ No samples contain fewer than 500 cells.\n"
  )
  
}

cat("\n")

# ============================================================================
# SECTION 7 — SCALE DATA
# ============================================================================

cat("=== SECTION 7: SCALING SKETCH DATA ===\n\n")

sketch_obj <- ScaleData(
  sketch_obj,
  features = VariableFeatures(sketch_obj),
  verbose = FALSE
)

cat("✓ Scaling complete\n\n")

# ============================================================================
# SECTION 8 — PCA
# ============================================================================

cat("=== SECTION 8: PCA ===\n\n")

sketch_obj <- RunPCA(
  sketch_obj,
  features = VariableFeatures(sketch_obj),
  npcs = 50,
  verbose = FALSE
)

cat("✓ PCA complete\n\n")

p_elbow <- ElbowPlot(
  sketch_obj,
  ndims = 50
) +
  ggtitle(
    "Global Atlas — PCA Elbow Plot"
  )

ggsave(
  "results/figures/phase3_pca_elbow.png",
  p_elbow,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✓ PCA elbow plot saved\n\n")

# ============================================================================
# SECTION 9 — PCA SELECTION
# ============================================================================

cat("=== SECTION 9: PCA DIMENSION SELECTION ===\n\n")

dims_use <- 1:20

cat(
  "Using PCs:",
  min(dims_use),
  "to",
  max(dims_use),
  "\n\n"
)

# ============================================================================
# SECTION 10 — NEIGHBORS + CLUSTERING
# ============================================================================

cat("=== SECTION 10: NEIGHBOR GRAPH + CLUSTERING ===\n\n")

sketch_obj <- FindNeighbors(
  sketch_obj,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

sketch_obj <- FindClusters(
  sketch_obj,
  resolution = 0.5,
  verbose = FALSE
)

n_clusters <- length(
  unique(
    sketch_obj$seurat_clusters
  )
)

cat(
  "Clusters detected:",
  n_clusters,
  "\n\n"
)

# ============================================================================
# SECTION 11 — UMAP
# ============================================================================

cat("=== SECTION 11: UMAP ===\n\n")

sketch_obj <- RunUMAP(
  sketch_obj,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap.cluster",
  reduction.key = "UMAPCLUSTER_",
  verbose = FALSE
)

cat("✓ UMAP complete\n\n")

# ============================================================================
# SECTION 12 — GLOBAL UMAP BY CLUSTER
# ============================================================================

cat("=== SECTION 12: GLOBAL UMAP BY CLUSTER ===\n\n")

p_cluster <- DimPlot(
  sketch_obj,
  reduction = "umap.cluster",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle(
    "Global Liver Atlas — Preliminary Clusters"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_umap_clusters_sketch.png",
  p_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Cluster UMAP saved\n\n")

# ============================================================================
# SECTION 13 — UMAP BY CLINICAL STATE
# ============================================================================

cat("=== SECTION 13: UMAP BY CLINICAL STATE ===\n\n")

p_phase <- DimPlot(
  sketch_obj,
  reduction = "umap.cluster",
  group.by = "Phase",
  pt.size = 0.4
) +
  ggtitle(
    "Global Liver Atlas — Clinical State"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_umap_by_phase_sketch.png",
  p_phase,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Clinical-state UMAP saved\n\n")

# ============================================================================
# SECTION 14 — UMAP BY DONOR
# ============================================================================

cat("=== SECTION 14: UMAP BY DONOR ===\n\n")

p_donor <- DimPlot(
  sketch_obj,
  reduction = "umap.cluster",
  group.by = "Donor",
  pt.size = 0.4
) +
  ggtitle(
    "Global Liver Atlas — Donor Distribution"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_umap_by_donor_sketch.png",
  p_donor,
  width = 12,
  height = 10,
  dpi = 300
)

cat("✓ Donor UMAP saved\n\n")

# ============================================================================
# SECTION 15 — CLUSTER × CLINICAL STATE
# ============================================================================

cat("=== SECTION 15: CLUSTER × CLINICAL STATE ===\n\n")

cluster_phase <- table(
  Cluster = sketch_obj$seurat_clusters,
  Phase = sketch_obj$Phase
)

cluster_phase_prop <- prop.table(
  cluster_phase,
  margin = 1
)

print(cluster_phase)

cat("\nCluster composition proportions:\n")

print(
  round(
    cluster_phase_prop,
    3
  )
)

write.csv(
  as.data.frame.matrix(cluster_phase),
  "results/tables/phase3_cluster_by_phase.csv"
)

write.csv(
  as.data.frame.matrix(cluster_phase_prop),
  "results/tables/phase3_cluster_by_phase_proportions.csv"
)

cat("\n✓ Cluster × phase tables saved\n\n")

# ============================================================================
# SECTION 16 — CLUSTER × DONOR
# ============================================================================

cat("=== SECTION 16: CLUSTER × DONOR ===\n\n")

cluster_donor <- table(
  Cluster = sketch_obj$seurat_clusters,
  Donor = sketch_obj$Donor
)

donors_per_cluster <- apply(
  cluster_donor,
  1,
  function(x) sum(x > 0)
)

largest_donor_fraction <- apply(
  cluster_donor,
  1,
  function(x) {
    
    x <- x[x > 0]
    
    max(x) / sum(x)
    
  }
)

donor_dominance <- data.frame(
  
  Cluster =
    names(largest_donor_fraction),
  
  Donors_Represented =
    donors_per_cluster[
      names(largest_donor_fraction)
    ],
  
  Largest_Donor_Fraction =
    round(
      largest_donor_fraction,
      3
    )
  
) %>%
  
  arrange(
    desc(Largest_Donor_Fraction)
  )

print(
  donor_dominance,
  row.names = FALSE
)

write_csv(
  donor_dominance,
  "results/tables/phase3_donor_dominance_by_cluster.csv"
)

cat(
  "\nInterpretation guide:\n",
  "- Clusters represented across multiple donors are more consistent\n",
  "  with shared biological populations.\n",
  "- Strong single-donor dominance requires closer inspection before\n",
  "  interpreting a cluster as a general biological state.\n\n"
)

# ============================================================================
# SECTION 17 — SAMPLE DOMINANCE
# ============================================================================

cat("=== SECTION 17: CLUSTER × SAMPLE ===\n\n")

cluster_sample <- table(
  Cluster = sketch_obj$seurat_clusters,
  GSM = sketch_obj$GSM
)

samples_per_cluster <- apply(
  cluster_sample,
  1,
  function(x) sum(x > 0)
)

largest_sample_fraction <- apply(
  cluster_sample,
  1,
  function(x) {
    
    x <- x[x > 0]
    
    max(x) / sum(x)
    
  }
)

sample_dominance <- data.frame(
  
  Cluster =
    names(largest_sample_fraction),
  
  Samples_Represented =
    samples_per_cluster[
      names(largest_sample_fraction)
    ],
  
  Largest_Sample_Fraction =
    round(
      largest_sample_fraction,
      3
    )
  
) %>%
  
  arrange(
    desc(Largest_Sample_Fraction)
  )

print(
  sample_dominance,
  row.names = FALSE
)

write_csv(
  sample_dominance,
  "results/tables/phase3_sample_dominance_by_cluster.csv"
)

cat("\n✓ Sample dominance assessment saved\n\n")

# ============================================================================
# SECTION 18 — CLUSTER MARKER DISCOVERY
# ============================================================================

cat("=== SECTION 18: CLUSTER MARKER DISCOVERY ===\n\n")

DefaultAssay(sketch_obj) <- "sketch"

cluster_markers <- FindAllMarkers(
  sketch_obj,
  assay = "sketch",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write_csv(
  cluster_markers,
  "results/tables/phase3_cluster_markers.csv"
)

top_markers <- cluster_markers %>%
  
  group_by(cluster) %>%
  
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  
  ungroup()

write_csv(
  top_markers,
  "results/tables/phase3_top10_markers_per_cluster.csv"
)

cat("\nTop markers per cluster:\n")

print(
  top_markers,
  n = Inf
)

cat("\n✓ Cluster marker analysis complete\n\n")

# ============================================================================
# SECTION 19 — CANONICAL LINEAGE MARKERS
# ============================================================================

cat("=== SECTION 19: CANONICAL LINEAGE MARKER PROGRAMS ===\n\n")

lineage_markers <- list(
  
  "CD4_T" = c(
    "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "CCR7"
  ),
  
  "CD8_T" = c(
    "CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "NKG7"
  ),
  
  "NK" = c(
    "NKG7", "KLRD1", "GNLY", "PRF1", "GZMB", "FCGR3A"
  ),
  
  "B_Cell" = c(
    "MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37"
  ),
  
  "Plasma_Cell" = c(
    "JCHAIN", "MZB1", "XBP1", "SDC1", "IGHG1"
  ),
  
  "Myeloid" = c(
    "LYZ", "FCER1G", "TYROBP", "LST1", "AIF1"
  ),
  
  "Macrophage" = c(
    "C1QA", "C1QB", "C1QC", "APOE", "CD68"
  ),
  
  "Monocyte" = c(
    "S100A8", "S100A9", "LYZ", "CTSD", "FCN1"
  ),
  
  "Dendritic_Cell" = c(
    "FCER1A", "CLEC10A", "CD1C", "HLA-DRA"
  ),
  
  "Hepatocyte" = c(
    "ALB", "APOA1", "APOC3", "TTR", "FGB"
  ),
  
  "Cholangiocyte" = c(
    "KRT19", "KRT7", "EPCAM", "KRT8", "KRT18"
  ),
  
  "Endothelial" = c(
    "KDR", "ESAM", "EMCN", "PECAM1", "VWF"
  ),
  
  "Stellate" = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"
  ),
  
  "Mast_Cell" = c(
    "TPSAB1", "TPSB2", "KIT", "MS4A2"
  ),
  
  "MAIT" = c(
    "TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "DPP4", "CCR6"
  ),
  
  "GammaDelta_T" = c(
    "TRDC",
    "TRGC1",
    "TRGC2",
    "TRDV1",
    "TRDV2"
  ),
  
  "pDC" = c(
    "CLEC4C",
    "LILRA4",
    "GZMB",
    "PTCRA",
    "TCF4"
  )
  
)

marker_presence <- tibble(
  
  Cell_Type =
    names(lineage_markers),
  
  Markers =
    sapply(
      lineage_markers,
      function(x) paste(x, collapse = "; ")
    ),
  
  Markers_Present =
    sapply(
      lineage_markers,
      function(x) {
        sum(
          x %in% rownames(sketch_obj)
        )
      }
    )
  
)

print(marker_presence)

write_csv(
  marker_presence,
  "results/tables/phase3_marker_programs.csv"
)

cat("\n✓ Canonical marker programs defined\n\n")

# ============================================================================
# SECTION 20 — MARKER HEATMAP
# ============================================================================

cat("=== SECTION 20: CLUSTER MARKER HEATMAP ===\n\n")

markers_for_heatmap <- unique(
  unlist(
    lineage_markers
  )
)

markers_for_heatmap <- intersect(
  markers_for_heatmap,
  rownames(sketch_obj)
)

if (length(markers_for_heatmap) > 1) {
  
  p_marker_heatmap <- DoHeatmap(
    sketch_obj,
    features = markers_for_heatmap,
    group.by = "seurat_clusters",
    raster = TRUE
  ) +
    NoLegend()
  
  ggsave(
    "results/figures/phase3_canonical_marker_heatmap.png",
    p_marker_heatmap,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  cat("✓ Marker heatmap saved\n\n")
  
} else {
  
  cat(
    "⚠ Insufficient canonical markers found for heatmap.\n\n"
  )
  
}

# ============================================================================
# SECTION 21 — MULTI-METHOD CELL-TYPE ANNOTATION
# ============================================================================

cat("=== SECTION 21: MULTI-METHOD CELL-TYPE ANNOTATION ===\n\n")

DefaultAssay(sketch_obj) <- "sketch"

if (!requireNamespace("SingleR", quietly = TRUE)) {
  stop(
    "ERROR: Package 'SingleR' is required. Install it with BiocManager::install('SingleR')."
  )
}

if (!requireNamespace("celldex", quietly = TRUE)) {
  stop(
    "ERROR: Package 'celldex' is required. Install it with BiocManager::install('celldex')."
  )
}

cluster_ids <- sort(
  unique(
    as.character(sketch_obj$seurat_clusters)
  )
)

marker_programs <- lineage_markers

marker_programs <- lapply(
  marker_programs,
  function(x) {
    intersect(
      x,
      rownames(sketch_obj)
    )
  }
)

marker_programs <- marker_programs[
  lengths(marker_programs) > 0
]

if (length(marker_programs) == 0) {
  stop(
    "ERROR: No canonical marker programs contain genes present in the sketch object."
  )
}

program_names <- names(marker_programs)

for (program_name in program_names) {
  
  sketch_obj <- AddModuleScore(
    object = sketch_obj,
    features = list(
      marker_programs[[program_name]]
    ),
    assay = "sketch",
    name = paste0(
      "Module_",
      program_name,
      "_"
    ),
    ctrl = 100,
    seed = 12345,
    search = FALSE,
    verbose = FALSE
  )
  
  old_score_name <- paste0(
    "Module_",
    program_name,
    "_1"
  )
  
  new_score_name <- paste0(
    "Module_",
    program_name
  )
  
  sketch_obj@meta.data[[new_score_name]] <-
    sketch_obj@meta.data[[old_score_name]]
  
  sketch_obj@meta.data[[old_score_name]] <- NULL
  
}

module_score_columns <- paste0(
  "Module_",
  program_names
)

module_scores_by_cluster <- sketch_obj@meta.data %>%
  
  mutate(
    Cluster = as.character(
      seurat_clusters
    )
  ) %>%
  
  select(
    Cluster,
    all_of(module_score_columns)
  ) %>%
  
  group_by(
    Cluster
  ) %>%
  
  summarise(
    across(
      everything(),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

module_scores_by_cluster <- module_scores_by_cluster %>%
  filter(
    !is.na(Cluster)
  ) %>%
  mutate(
    Cluster = as.character(Cluster)
  )

module_score_matrix <- module_scores_by_cluster %>%
  column_to_rownames(
    "Cluster"
  ) %>%
  as.matrix()

storage.mode(module_score_matrix) <- "numeric"

colnames(module_score_matrix) <- gsub(
  "^Module_",
  "",
  colnames(module_score_matrix)
)

module_top_labels <- apply(
  module_score_matrix,
  1,
  function(x) {
    names(x)[which.max(x)]
  }
)

module_top_scores <- apply(
  module_score_matrix,
  1,
  max
)

module_second_scores <- apply(
  module_score_matrix,
  1,
  function(x) {
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    sort(
      x,
      decreasing = TRUE
    )[2]
  }
)

module_score_gap <- module_top_scores -
  module_second_scores

module_annotation <- tibble(
  
  Cluster = names(
    module_top_labels
  ),
  
  Module_Label =
    unname(
      module_top_labels
    ),
  
  Module_Top_Score =
    unname(
      module_top_scores
    ),
  
  Module_Second_Score =
    unname(
      module_second_scores
    ),
  
  Module_Score_Gap =
    unname(
      module_score_gap
    )
  
) %>%
  
  arrange(
    as.numeric(Cluster)
  )

write_csv(
  module_annotation,
  "results/tables/phase3_module_score_annotation.csv"
)

write.csv(
  module_score_matrix,
  "results/tables/phase3_module_scores_by_cluster.csv",
  row.names = TRUE
)

p_module_heatmap <- pheatmap::pheatmap(
  module_score_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "column",
  main = "Canonical Marker Program Scores by Cluster",
  filename = "results/figures/phase3_module_score_heatmap.png",
  width = 12,
  height = 8
)

top_marker_genes <- cluster_markers %>%
  
  mutate(
    cluster = as.character(cluster)
  ) %>%
  
  group_by(
    cluster
  ) %>%
  
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  
  slice_head(
    n = 50
  ) %>%
  
  summarise(
    Top_Markers = list(
      unique(gene)
    ),
    .groups = "drop"
  )

marker_overlap_matrix <- matrix(
  0,
  nrow = length(cluster_ids),
  ncol = length(marker_programs),
  dimnames = list(
    cluster_ids,
    names(marker_programs)
  )
)

for (cluster_id in cluster_ids) {
  
  cluster_marker_row <- top_marker_genes %>%
    
    filter(
      cluster == cluster_id
    )
  
  if (nrow(cluster_marker_row) == 0) {
    next
  }
  
  cluster_genes <- cluster_marker_row$Top_Markers[[1]]
  
  for (program_name in names(marker_programs)) {
    
    marker_overlap_matrix[
      cluster_id,
      program_name
    ] <-
      sum(
        cluster_genes %in%
          marker_programs[[program_name]]
      )
    
  }
  
}

marker_overlap_top_labels <- apply(
  marker_overlap_matrix,
  1,
  function(x) {
    
    if (max(x) == 0) {
      return("Unresolved")
    }
    
    names(x)[which.max(x)]
    
  }
)

marker_overlap_top_scores <- apply(
  marker_overlap_matrix,
  1,
  max
)

marker_overlap_annotation <- tibble(
  
  Cluster = names(
    marker_overlap_top_labels
  ),
  
  Marker_Overlap_Label =
    unname(
      marker_overlap_top_labels
    ),
  
  Marker_Overlap_Count =
    unname(
      marker_overlap_top_scores
    )
  
) %>%
  
  arrange(
    as.numeric(Cluster)
  )

write_csv(
  marker_overlap_annotation,
  "results/tables/phase3_marker_overlap_annotation.csv"
)

write.csv(
  marker_overlap_matrix,
  "results/tables/phase3_marker_overlap_scores.csv",
  row.names = TRUE
)

singleR_reference <- celldex::HumanPrimaryCellAtlasData()

singleR_data <- GetAssayData(
  sketch_obj,
  assay = "sketch",
  layer = "data"
)

singleR_clusters <- sketch_obj$seurat_clusters[
  colnames(singleR_data)
]

if (
  length(singleR_clusters) !=
  ncol(singleR_data)
) {
  
  stop(
    "ERROR: Number of cluster labels does not match the number of sketch cells."
  )
  
}

if (
  any(
    is.na(singleR_clusters)
  )
) {
  
  stop(
    "ERROR: Some sketch cells are missing cluster assignments."
  )
  
}

singleR_result <- SingleR::SingleR(
  test = singleR_data,
  ref = singleR_reference,
  labels = singleR_reference$label.main,
  clusters = singleR_clusters
)

singleR_annotation <- tibble(
  
  Cluster =
    rownames(
      singleR_result
    ),
  
  SingleR_Label =
    as.character(
      singleR_result$labels
    ),
  
  SingleR_Delta =
    as.numeric(
      singleR_result$delta.next
    )
  
) %>%
  
  mutate(
    Cluster = as.character(
      Cluster
    )
  ) %>%
  
  arrange(
    as.numeric(Cluster)
  )

write_csv(
  singleR_annotation,
  "results/tables/phase3_singleR_annotation.csv"
)

singleR_scores <- as.data.frame(
  singleR_result$scores
)

singleR_scores <- singleR_scores %>%
  
  rownames_to_column(
    "Cluster"
  )

write_csv(
  singleR_scores,
  "results/tables/phase3_singleR_scores_by_cluster.csv"
)

annotation_evidence <- tibble(
  
  Cluster = cluster_ids
  
) %>%
  
  left_join(
    module_annotation,
    by = "Cluster"
  ) %>%
  
  left_join(
    marker_overlap_annotation,
    by = "Cluster"
  ) %>%
  
  left_join(
    singleR_annotation,
    by = "Cluster"
  ) %>%
  
  arrange(
    as.numeric(Cluster)
  )

write_csv(
  annotation_evidence,
  "results/tables/phase3_annotation_evidence_all_methods.csv"
)

print(
  annotation_evidence,
  n = Inf
)

cluster_annotation_template <- annotation_evidence %>%
  
  mutate(
    
    Final_CellType = NA_character_,
    
    Confidence = NA_character_,
    
    Rationale = NA_character_
    
  )

write_csv(
  cluster_annotation_template,
  "results/tables/phase3_cluster_annotation_TEMPLATE.csv"
)

cat(
  "\n✓ Method 1 complete: canonical marker program scoring\n"
)

cat(
  "✓ Method 2 complete: cluster marker-program overlap\n"
)

cat(
  "✓ Method 3 complete: SingleR reference annotation\n"
)

cat(
  "✓ Consolidated annotation evidence saved\n"
)

cat(
  "✓ Final annotation template saved\n\n"
)

cat(
  "Review:\n",
  "  1. results/tables/phase3_annotation_evidence_all_methods.csv\n",
  "  2. results/tables/phase3_module_score_annotation.csv\n",
  "  3. results/tables/phase3_marker_overlap_annotation.csv\n",
  "  4. results/tables/phase3_singleR_annotation.csv\n",
  "  5. results/tables/phase3_top10_markers_per_cluster.csv\n",
  "  6. results/figures/phase3_module_score_heatmap.png\n",
  "  7. results/figures/phase3_canonical_marker_heatmap.png\n\n"
)

cat(
  "Complete:\n",
  "results/tables/phase3_cluster_annotation.csv\n\n"
)

# ============================================================================
# SECTION 21B — SAVING PRE-FINAL-ANNOTATION OBJECT
# ============================================================================

cat("=== SECTION 21B: SAVING PRE-FINAL-ANNOTATION OBJECT ===\n\n")

saveRDS(
  sketch_obj,
  "results/rds_objects/seurat_phase3_preannotation_sketch.rds"
)

cat(
  "✓ Saved annotation-evidence sketch object:\n",
  "results/rds_objects/seurat_phase3_preannotation_sketch.rds\n\n",
  sep = ""
)

cat(
  "============================================================\n"
)

cat(
  "FINAL ANNOTATION GATE REACHED\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Review the three annotation methods and complete:\n"
)

cat(
  "results/tables/phase3_cluster_annotation.csv\n\n"
)

cat(
  "Required columns:\n"
)

cat(
  "Cluster\n"
)

cat(
  "Final_CellType\n"
)

cat(
  "Confidence\n"
)

cat(
  "Rationale\n\n"
)

cluster_annotation <- read_csv(
  "results/tables/phase3_cluster_annotation_TEMPLATE.csv",
  show_col_types = FALSE
)

write_csv(
  cluster_annotation,
  "results/tables/phase3_cluster_annotation.csv"
)

cat(
  "✓ Final cluster annotation file saved:\n",
  "results/tables/phase3_cluster_annotation.csv\n"
)

cluster_annotation <- cluster_annotation %>%
  
  mutate(
    
    Final_CellType = case_when(
      
      Cluster %in% c("0", "1") ~ "CD8_T",
      
      Cluster == "2" ~ "MAIT",
      
      Cluster %in% c("3", "6") ~ "CD4_T",
      
      Cluster %in% c("4", "14") ~ "NK",
      
      Cluster == "5" ~ "NK_like",
      
      Cluster %in% c("7", "10") ~ "B_Cell",
      
      Cluster == "8" ~ "T_Cell",
      
      Cluster == "9" ~ "Monocyte",
      
      Cluster == "11" ~ "Plasma_Cell",
      
      Cluster == "12" ~ "Macrophage",
      
      Cluster == "13" ~ "Myeloid",
      
      Cluster == "15" ~ "Myeloid",
      
      Cluster == "16" ~ "pDC",
      
      Cluster == "17" ~ "GammaDelta_T",
      
      TRUE ~ "Unresolved"
      
    ),
    
    Confidence = case_when(
      
      Cluster %in% c(
        "0", "1", "3", "4", "6",
        "7", "10", "11", "12",
        "14", "16"
      ) ~ "High",
      
      Cluster %in% c(
        "2", "5", "9", "13", "17"
      ) ~ "Moderate",
      
      Cluster %in% c(
        "8", "15"
      ) ~ "Low",
      
      TRUE ~ "Low"
      
    ),
    
    Rationale = case_when(
      
      Cluster %in% c("0", "1") ~
        "CD8_T supported by module scoring and marker overlap.",
      
      Cluster == "2" ~
        "MAIT module and marker-overlap evidence support a MAIT-like T-cell identity.",
      
      Cluster %in% c("3", "6") ~
        "CD4_T supported by module scoring and marker overlap.",
      
      Cluster %in% c("4", "14") ~
        "NK identity strongly supported by module scoring and marker overlap.",
      
      Cluster == "5" ~
        "NK module signal is present but marker overlap favors myeloid lineage; retained as NK-like pending downstream refinement.",
      
      Cluster %in% c("7", "10") ~
        "B-cell identity strongly supported by module scoring and marker overlap.",
      
      Cluster == "8" ~
        "Conflicting CD8_T and CD4_T program evidence; assigned broadly as T_Cell.",
      
      Cluster == "9" ~
        "Strong myeloid program with monocyte marker overlap; assigned Monocyte.",
      
      Cluster == "11" ~
        "Plasma-cell identity strongly supported by the plasma-cell module and marker evidence.",
      
      Cluster == "12" ~
        "Myeloid and macrophage programs are both enriched; macrophage marker overlap supports Macrophage.",
      
      Cluster == "13" ~
        "Strong myeloid signal without sufficiently specific marker-overlap evidence for a narrower subtype.",
      
      Cluster == "15" ~
        "Weak and mixed NK/myeloid signals; assigned broadly as Myeloid with low confidence.",
      
      Cluster == "16" ~
        "pDC program strongly enriched and supported by canonical pDC markers.",
      
      Cluster == "17" ~
        "NK-like program is present, but gamma-delta T-cell marker overlap supports GammaDelta_T.",
      
      TRUE ~
        "Insufficient evidence for confident lineage assignment."
      
    )
    
  )

print(
  cluster_annotation %>%
    select(
      Cluster,
      Module_Label,
      Marker_Overlap_Label,
      SingleR_Label,
      Final_CellType,
      Confidence,
      Rationale
    ),
  n = Inf
)

write_csv(
  cluster_annotation,
  "results/tables/phase3_cluster_annotation.csv"
)

# ============================================================================
# SECTION 22 - FINAL CELL-TYPE ANNOTATION
# ============================================================================

cat("=== SECTION 22: FINAL CELL-TYPE ANNOTATION ===\n\n")

annotation_file <- "results/tables/phase3_cluster_annotation.csv"

if (!file.exists(annotation_file)) {
  
  stop(
    "ANNOTATION GATE: Missing file '",
    annotation_file,
    "'. Review the multi-method evidence and complete the final annotation table."
  )
  
}

cluster_annotation <- read_csv(
  annotation_file,
  show_col_types = FALSE
)

required_annotation_columns <- c(
  "Cluster",
  "Final_CellType",
  "Confidence",
  "Rationale"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(cluster_annotation)
)

if (length(missing_annotation_columns) > 0) {
  
  stop(
    "ERROR: Annotation file is missing required columns: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
  
}

cluster_annotation <- cluster_annotation %>%
  
  mutate(
    
    Cluster = as.character(
      Cluster
    ),
    
    Final_CellType = as.character(
      Final_CellType
    ),
    
    Confidence = as.character(
      Confidence
    ),
    
    Rationale = as.character(
      Rationale
    )
    
  )

if (
  any(
    is.na(
      cluster_annotation$Final_CellType
    )
  ) ||
  any(
    trimws(
      cluster_annotation$Final_CellType
    ) == ""
  )
) {
  
  stop(
    "ERROR: Every cluster must have a Final_CellType."
  )
  
}

if (
  any(
    is.na(
      cluster_annotation$Confidence
    )
  ) ||
  any(
    trimws(
      cluster_annotation$Confidence
    ) == ""
  )
) {
  
  stop(
    "ERROR: Every cluster must have a Confidence value."
  )
  
}

if (
  any(
    is.na(
      cluster_annotation$Rationale
    )
  ) ||
  any(
    trimws(
      cluster_annotation$Rationale
    ) == ""
  )
) {
  
  stop(
    "ERROR: Every cluster must have a biological Rationale."
  )
  
}

if (
  anyDuplicated(
    cluster_annotation$Cluster
  ) > 0
) {
  
  stop(
    "ERROR: Duplicate cluster IDs found in phase3_cluster_annotation.csv."
  )
  
}

expected_clusters <- sort(
  unique(
    as.character(
      sketch_obj$seurat_clusters
    )
  )
)

annotation_clusters <- sort(
  unique(
    cluster_annotation$Cluster
  )
)

if (
  !identical(
    expected_clusters,
    annotation_clusters
  )
) {
  
  missing_annotations <- setdiff(
    expected_clusters,
    annotation_clusters
  )
  
  extra_annotations <- setdiff(
    annotation_clusters,
    expected_clusters
  )
  
  stop(
    "ERROR: Cluster IDs do not exactly match the sketch clusters.\n",
    "Missing annotations: ",
    paste(
      missing_annotations,
      collapse = ", "
    ),
    "\nExtra annotations: ",
    paste(
      extra_annotations,
      collapse = ", "
    )
  )
  
}

annotation_lookup <- setNames(
  cluster_annotation$Final_CellType,
  cluster_annotation$Cluster
)

cell_annotations <- annotation_lookup[
  as.character(
    sketch_obj$seurat_clusters
  )
]

names(cell_annotations) <- colnames(
  sketch_obj
)

sketch_obj <- AddMetaData(
  object = sketch_obj,
  metadata = cell_annotations,
  col.name = "CellType.sketch"
)

confidence_lookup <- setNames(
  cluster_annotation$Confidence,
  cluster_annotation$Cluster
)

cell_confidence <- confidence_lookup[
  as.character(
    sketch_obj$seurat_clusters
  )
]

names(cell_confidence) <- colnames(
  sketch_obj
)

sketch_obj <- AddMetaData(
  object = sketch_obj,
  metadata = cell_confidence,
  col.name = "CellType_Confidence"
)

annotated_cells <- colnames(sketch_obj)[
  !is.na(
    sketch_obj$seurat_clusters
  )
]

if (
  any(
    is.na(
      sketch_obj$CellType.sketch[
        annotated_cells
      ]
    )
  )
) {
  
  stop(
    "ERROR: At least one sketch cell did not receive a final annotation."
  )
  
}

final_annotation_summary <- sketch_obj@meta.data %>%
  
  filter(
    !is.na(
      seurat_clusters
    )
  ) %>%
  
  count(
    Cluster = seurat_clusters,
    CellType = CellType.sketch,
    Confidence = CellType_Confidence,
    name = "Cells"
  ) %>%
  
  arrange(
    as.numeric(
      as.character(
        Cluster
      )
    )
  )

write_csv(
  final_annotation_summary,
  "results/tables/phase3_final_annotation_summary.csv"
)

write_csv(
  cluster_annotation,
  "results/tables/phase3_final_cluster_annotations.csv"
)

p_final_annotation <- DimPlot(
  sketch_obj,
  reduction = "umap.cluster",
  group.by = "CellType.sketch",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  
  ggtitle(
    "Global Liver Atlas — Final Cell-Type Annotation"
  ) +
  
  theme_minimal()

ggsave(
  "results/figures/phase3_umap_final_celltypes_sketch.png",
  p_final_annotation,
  width = 12,
  height = 9,
  dpi = 300
)

saveRDS(
  sketch_obj,
  "results/rds_objects/seurat_phase3_annotated_sketch.rds"
)

cat(
  "Annotated sketch cells:",
  sum(
    !is.na(
      sketch_obj$CellType.sketch[
        annotated_cells
      ]
    )
  ),
  "/",
  length(
    annotated_cells
  ),
  "\n\n"
)

cat(
  "✓ Final biological annotations attached\n"
)

cat(
  "✓ Final annotation summary saved\n"
)

cat(
  "✓ Final annotated sketch object saved\n"
)

cat(
  "✓ Final cell-type UMAP saved\n\n"
)

# ============================================================================
# SECTION 23 — PROJECT FULL DATASET INTO SKETCH SPACE
# ============================================================================

cat("=== SECTION 23: FULL-DATASET PROJECTION ===\n\n")

sketch_obj <- readRDS("results/rds_objects/seurat_phase3_annotated_sketch.rds")

dims_use <- 1:20

sketch_obj <- ProjectData(
  sketch_obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

sketch_obj <- RunUMAP(
  sketch_obj,
  reduction = "pca.full",
  dims = dims_use,
  reduction.name = "umap.full",
  reduction.key = "UMAPFULL_",
  verbose = FALSE
)

cat("✓ Full dataset projected\n\n")

# ============================================================================
# SECTION 24 — PREPARE LABEL TRANSFER
# ============================================================================

cat("=== SECTION 24: PREPARING LABEL TRANSFER ===\n\n")

reference_cells <- colnames(
  sketch_obj[["sketch"]]
)

query_cells <- setdiff(
  colnames(sketch_obj),
  reference_cells
)

cat(
  "Reference sketch cells:",
  length(reference_cells),
  "\n"
)

cat(
  "Cells requiring projection:",
  length(query_cells),
  "\n\n"
)

sketch_obj$ProjectedCluster <- NA_character_
sketch_obj$CellType.full <- NA_character_
sketch_obj$NeighborAgreement <- NA_real_

sketch_obj$ProjectedCluster[
  reference_cells
] <-
  as.character(
    sketch_obj$seurat_clusters[
      reference_cells
    ]
  )

sketch_obj$CellType.full[
  reference_cells
] <-
  sketch_obj$CellType.sketch[
    reference_cells
  ]

sketch_obj$NeighborAgreement[
  reference_cells
] <- 1

cat(
  "Sketch cells assigned direct labels.\n\n"
)

# ============================================================================
# SECTION 25 — KNN LABEL TRANSFER WITH NEIGHBOR AGREEMENT
# ============================================================================

cat("=== SECTION 25: KNN LABEL TRANSFER ===\n\n")

pca_full <- Embeddings(
  sketch_obj,
  reduction = "pca.full"
)

reference_pca <- pca_full[
  reference_cells,
  dims_use,
  drop = FALSE
]

query_pca <- pca_full[
  query_cells,
  dims_use,
  drop = FALSE
]

k_neighbors <- 30

knn_result <- get.knnx(
  data = reference_pca,
  query = query_pca,
  k = min(
    k_neighbors,
    nrow(reference_pca)
  )
)

reference_clusters <- as.character(
  sketch_obj$seurat_clusters[
    reference_cells
  ]
)

reference_celltypes <- as.character(
  sketch_obj$CellType.sketch[
    reference_cells
  ]
)

projected_clusters <- character(
  length(query_cells)
)

projected_celltypes <- character(
  length(query_cells)
)

neighbor_agreements <- numeric(
  length(query_cells)
)

for (i in seq_along(query_cells)) {
  
  neighbor_indices <- knn_result$nn.index[
    i,
  ]
  
  neighbor_clusters <- reference_clusters[
    neighbor_indices
  ]
  
  neighbor_celltypes <- reference_celltypes[
    neighbor_indices
  ]
  
  cluster_table <- table(
    neighbor_clusters
  )
  
  celltype_table <- table(
    neighbor_celltypes
  )
  
  projected_clusters[i] <- names(
    cluster_table
  )[
    which.max(
      cluster_table
    )
  ]
  
  projected_celltypes[i] <- names(
    celltype_table
  )[
    which.max(
      celltype_table
    )
  ]
  
  neighbor_agreements[i] <- max(
    celltype_table
  ) /
    sum(
      celltype_table
    )
  
}

names(projected_clusters) <- query_cells

names(projected_celltypes) <- query_cells

names(neighbor_agreements) <- query_cells

sketch_obj <- AddMetaData(
  sketch_obj,
  metadata = projected_clusters,
  col.name = "ProjectedCluster"
)

sketch_obj <- AddMetaData(
  sketch_obj,
  metadata = projected_celltypes,
  col.name = "CellType.full"
)

sketch_obj <- AddMetaData(
  sketch_obj,
  metadata = neighbor_agreements,
  col.name = "NeighborAgreement"
)

cat(
  "✓ KNN label transfer complete\n\n"
)

# ============================================================================
# SECTION 26 — FLAG LOW-AGREEMENT ASSIGNMENTS
# ============================================================================

cat("=== SECTION 26: NEIGHBOR AGREEMENT ===\n\n")

agreement_threshold <- 0.70

sketch_obj$LowAgreementTransfer <- ifelse(
  
  is.na(sketch_obj$NeighborAgreement),
  
  TRUE,
  
  sketch_obj$NeighborAgreement < agreement_threshold
  
)

sketch_obj$AnnotationStatus <-
  ifelse(
    
    sketch_obj$LowAgreementTransfer,
    
    "Low agreement",
    
    "Assigned"
    
  )

sketch_obj$AnnotationStatus[
  reference_cells
] <- "Direct sketch annotation"

agreement_summary <- sketch_obj@meta.data %>%
  
  summarise(
    
    Mean_Neighbor_Agreement =
      mean(
        NeighborAgreement,
        na.rm = TRUE
      ),
    
    Median_Neighbor_Agreement =
      median(
        NeighborAgreement,
        na.rm = TRUE
      ),
    
    Low_Agreement_Cells =
      sum(
        LowAgreementTransfer,
        na.rm = TRUE
      ),
    
    Percent_Low_Agreement =
      round(
        100 *
          mean(
            LowAgreementTransfer,
            na.rm = TRUE
          ),
        2
      )
    
  )

print(
  agreement_summary
)

write_csv(
  agreement_summary,
  "results/tables/phase3_neighbor_agreement_summary.csv"
)

cat("\n✓ Neighbor agreement summarized\n\n")

# ============================================================================
# SECTION 27 — FULL DATASET COMPOSITION
# ============================================================================

cat("=== SECTION 27: FULL DATASET COMPOSITION ===\n\n")

celltype_by_phase <- sketch_obj@meta.data %>%
  
  count(
    CellType.full,
    Phase,
    name = "Cells"
  )

write_csv(
  celltype_by_phase,
  "results/tables/phase3_celltype_by_phase.csv"
)

celltype_by_donor <- sketch_obj@meta.data %>%
  
  count(
    CellType.full,
    Donor,
    name = "Cells"
  )

write_csv(
  celltype_by_donor,
  "results/tables/phase3_celltype_by_donor.csv"
)

celltype_counts <- sketch_obj@meta.data %>%
  
  count(
    CellType.full,
    name = "Cells"
  ) %>%
  
  arrange(
    desc(Cells)
  )

print(
  celltype_counts
)

write_csv(
  celltype_counts,
  "results/tables/phase3_celltype_counts.csv"
)

# ============================================================================
# SECTION 28 — FINAL UMAPS
# ============================================================================

cat("=== SECTION 28: FINAL UMAPS ===\n\n")

p_full_clusters <- DimPlot(
  sketch_obj,
  reduction = "umap.full",
  group.by = "ProjectedCluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.1
) +
  ggtitle(
    "Full Dataset — Projected Clusters"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_full_umap_clusters.png",
  p_full_clusters,
  width = 10,
  height = 8,
  dpi = 300
)

p_full_celltypes <- DimPlot(
  sketch_obj,
  reduction = "umap.full",
  group.by = "CellType.full",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.1
) +
  ggtitle(
    "Full Dataset — Cell-Type Annotation"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_full_umap_celltypes.png",
  p_full_celltypes,
  width = 10,
  height = 8,
  dpi = 300
)

p_full_phase <- DimPlot(
  sketch_obj,
  reduction = "umap.full",
  group.by = "Phase",
  raster = TRUE,
  pt.size = 0.1
) +
  ggtitle(
    "Full Dataset — Clinical State"
  ) +
  theme_minimal()

ggsave(
  "results/figures/phase3_full_umap_phase.png",
  p_full_phase,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Final UMAPs saved\n\n")

# ============================================================================
# SECTION 29 — NEIGHBOR AGREEMENT FIGURE
# ============================================================================

p_agreement <- ggplot(
  sketch_obj@meta.data,
  aes(
    x = NeighborAgreement
  )
) +
  
  geom_histogram(
    bins = 50
  ) +
  
  geom_vline(
    xintercept = agreement_threshold,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Cell-Type KNN Neighbor Agreement",
    x = "Fraction of neighbors supporting assigned cell type",
    y = "Cells"
  )

ggsave(
  "results/figures/phase3_neighbor_agreement.png",
  p_agreement,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================================
# SECTION 30 — LOW-AGREEMENT CELLS BY CLINICAL STATE
# ============================================================================

low_agreement_by_phase <- sketch_obj@meta.data %>%
  
  group_by(
    Phase
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Low_Agreement =
      sum(
        LowAgreementTransfer,
        na.rm = TRUE
      ),
    
    Percent_Low_Agreement =
      round(
        100 *
          Low_Agreement /
          Cells,
        2
      ),
    
    .groups = "drop"
  )

print(
  low_agreement_by_phase
)

write_csv(
  low_agreement_by_phase,
  "results/tables/phase3_low_agreement_by_phase.csv"
)

# ============================================================================
# SECTION 31 — FINAL METADATA TABLE
# ============================================================================

cat("=== SECTION 31: SAVING FINAL METADATA ===\n\n")

final_metadata <- sketch_obj@meta.data %>%
  
  rownames_to_column(
    "Cell"
  )

write_csv(
  final_metadata,
  "results/tables/phase3_complete_cell_metadata.csv"
)

cat("✓ Complete metadata table saved\n\n")

# ============================================================================
# SECTION 32 — FINAL INTEGRITY CHECK
# ============================================================================

cat("============================================================\n")
cat("PHASE 3 — FINAL INTEGRITY CHECK\n")
cat("============================================================\n\n")

cat(
  "Full cells:",
  ncol(sketch_obj),
  "\n"
)

cat(
  "Genes:",
  nrow(sketch_obj),
  "\n"
)

cat(
  "Clinical states:",
  n_distinct(sketch_obj$Phase),
  "\n"
)

cat(
  "Donors:",
  n_distinct(sketch_obj$Donor),
  "\n"
)

cat(
  "Clusters:",
  length(
    unique(
      na.omit(
        sketch_obj$ProjectedCluster
      )
    )
  ),
  "\n"
)

cat(
  "Annotated cell types:",
  length(
    unique(
      na.omit(
        sketch_obj$CellType.full
      )
    )
  ),
  "\n"
)

cat(
  "Mean neighbor agreement:",
  round(
    mean(
      sketch_obj$NeighborAgreement,
      na.rm = TRUE
    ),
    3
  ),
  "\n"
)

cat(
  "Median neighbor agreement:",
  round(
    median(
      sketch_obj$NeighborAgreement,
      na.rm = TRUE
    ),
    3
  ),
  "\n"
)

cat(
  "Low-agreement cells:",
  sum(
    sketch_obj$LowAgreementTransfer,
    na.rm = TRUE
  ),
  "\n\n"
)

# ============================================================================
# SECTION 33 — SAVE FINAL OBJECT
# ============================================================================

cat("=== SECTION 33: SAVING FINAL PHASE 3 OBJECT ===\n\n")

saveRDS(
  sketch_obj,
  "results/rds_objects/seurat_phase3_global_atlas.rds"
)

cat(
  "✓ Saved final object:\n"
)

cat(
  "  results/rds_objects/seurat_phase3_global_atlas.rds\n\n"
)

# ============================================================================
# SECTION 34 — FINAL STATUS
# ============================================================================

cat("============================================================\n")
cat("PHASE 3 — SCRIPT 03 COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Completed:\n"
)

cat(
  "  ✓ QC-filtered singlet dataset loaded\n"
)

cat(
  "  ✓ Metadata validated\n"
)

cat(
  "  ✓ Normalization completed\n"
)

cat(
  "  ✓ Variable features identified\n"
)

cat(
  "  ✓ Representative sketch constructed\n"
)

cat(
  "  ✓ PCA completed\n"
)

cat(
  "  ✓ Clustering completed\n"
)

cat(
  "  ✓ Global UMAP generated\n"
)

cat(
  "  ✓ Donor/sample dominance assessed\n"
)

cat(
  "  ✓ Cluster markers identified\n"
)

cat(
  "  ✓ Canonical marker programs evaluated\n"
)

cat(
  "  ✓ Manual annotations attached\n"
)

cat(
  "  ✓ Full-dataset projection completed\n"
)

cat(
  "  ✓ Labels transferred with KNN neighbor agreement\n"
)

cat(
  "  ✓ Low-agreement assignments flagged\n"
)

cat(
  "  ✓ Composition tables saved\n"
)

cat(
  "  ✓ Final Seurat object saved\n\n"
)

cat(
  "FINAL GATE:\n"
)

cat(
  "Review biological annotations, donor/sample dominance,\n"
)

cat(
  "and low-agreement assignments before downstream analyses.\n\n"
)

cat("============================================================\n")
cat("END OF SCRIPT 03\n")
cat("============================================================\n")
