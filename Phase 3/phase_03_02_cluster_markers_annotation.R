# ==============================================================================
# HBV scRNA-seq PROJECT
# PHASE 3 — GLOBAL ATLAS CONSTRUCTION
# SCRIPT 02 — CLUSTER MARKERS, ANNOTATION EVIDENCE & FINAL CELL-TYPE LABELS
# ==============================================================================
#
# INPUT:
#   results/rds_objects/phase3_atlas_sketch_umap_clustered.rds
#
# OUTPUT:
#   results/rds_objects/phase3_atlas_sketch_annotation_evidence.rds
#   results/rds_objects/phase3_atlas_sketch_final_annotation.rds
#
# PURPOSE:
#   1. Characterize cluster composition across clinical states.
#   2. Assess donor/sample representation of every cluster.
#   3. Identify descriptive cluster marker genes.
#   4. Compare clusters against canonical lineage marker programs.
#   5. Generate marker-based module-score evidence.
#   6. Generate reference-based SingleR/Monaco evidence.
#   7. Save a complete evidence checkpoint.
#   8. Apply FINAL cell-type labels only after explicit manual adjudication.
#
# IMPORTANT:
#   - This script does NOT perform normalization.
#   - This script does NOT perform batch correction.
#   - This script does NOT perform filtering.
#   - Cluster-marker p-values are DESCRIPTIVE and are not donor-level
#     biological inference.
#   - The atlas currently contains 20 clusters. Cluster count is detected
#     dynamically rather than hard-coded.
#   - Final labels must be explicitly supplied for ALL clusters.
#
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
  library(SingleR)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(celldex)
})

set.seed(12345)

input_file <-
  "results/rds_objects/phase3_atlas_sketch_umap_clustered.rds"

evidence_output_file <-
  "results/rds_objects/phase3_atlas_sketch_annotation_evidence.rds"

final_output_file <-
  "results/rds_objects/phase3_atlas_sketch_final_annotation.rds"

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


# ==============================================================================
# 2. LOAD ATLAS
# ==============================================================================

message("Loading Phase 3 atlas...")

atlas_sketch <- readRDS(input_file)

if (!inherits(atlas_sketch, "Seurat")) {
  stop("Input object is not a Seurat object.")
}

if (!"RNA" %in% names(atlas_sketch@assays)) {
  stop("RNA assay not found.")
}

DefaultAssay(atlas_sketch) <- "RNA"


# ==============================================================================
# 3. BASIC VALIDATION
# ==============================================================================

message("Validating atlas structure...")

if (!"seurat_clusters" %in% colnames(atlas_sketch@meta.data)) {
  stop("seurat_clusters metadata is missing.")
}

cluster_ids <- sort(
  unique(
    as.character(atlas_sketch$seurat_clusters)
  )
)

n_atlas_clusters <- length(cluster_ids)

message(
  "Detected ",
  n_atlas_clusters,
  " atlas clusters: ",
  paste(cluster_ids, collapse = ", ")
)

if (n_atlas_clusters < 2) {
  stop("Fewer than two clusters detected.")
}

if (ncol(atlas_sketch) != 20000) {
  warning(
    "Atlas contains ",
    ncol(atlas_sketch),
    " cells rather than exactly 20,000."
  )
}

message(
  "Genes: ", nrow(atlas_sketch),
  "\nCells: ", ncol(atlas_sketch),
  "\nClusters: ", n_atlas_clusters
)


# ==============================================================================
# 4. CHECK CLINICAL-STATE AND DONOR METADATA
# ==============================================================================

required_metadata <- c(
  "Phase",
  "Donor",
  "GSM"
)

missing_metadata <- base::setdiff(
  required_metadata,
  colnames(atlas_sketch@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Required metadata missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}


# ==============================================================================
# 5. STANDARDIZE CLUSTER METADATA
# ==============================================================================

atlas_sketch$Atlas_Cluster <-
  factor(
    as.character(atlas_sketch$seurat_clusters),
    levels = cluster_ids
  )


# ==============================================================================
# 6. CHECKPOINT SUMMARY
# ==============================================================================

atlas_summary <- tibble(
  Metric = c(
    "Cells",
    "Genes",
    "Atlas clusters",
    "Clinical states",
    "Donors",
    "Samples"
  ),
  Value = c(
    ncol(atlas_sketch),
    nrow(atlas_sketch),
    n_atlas_clusters,
    n_distinct(atlas_sketch$Phase),
    n_distinct(atlas_sketch$Donor),
    n_distinct(atlas_sketch$GSM)
  )
)

print(atlas_summary)


# ==============================================================================
# 17. CLUSTER × CLINICAL STATE COMPOSITION
# ==============================================================================

message("Section 17: cluster × clinical state composition...")

cluster_phase_counts <- atlas_sketch@meta.data %>%
  dplyr::count(
    Atlas_Cluster,
    Phase,
    name = "Cells"
  )

cluster_phase_within_state <- cluster_phase_counts %>%
  group_by(Phase) %>%
  mutate(
    Proportion_Within_State = Cells / sum(Cells)
  ) %>%
  ungroup()

cluster_phase_within_cluster <- cluster_phase_counts %>%
  group_by(Atlas_Cluster) %>%
  mutate(
    Proportion_Within_Cluster = Cells / sum(Cells)
  ) %>%
  ungroup()

write.csv(
  cluster_phase_counts,
  "results/tables/cluster_by_phase_counts.csv",
  row.names = FALSE
)

write.csv(
  cluster_phase_within_state,
  "results/tables/cluster_by_phase_within_state.csv",
  row.names = FALSE
)

write.csv(
  cluster_phase_within_cluster,
  "results/tables/cluster_by_phase_within_cluster.csv",
  row.names = FALSE
)


# ==============================================================================
# 18. CLUSTER × DONOR REPRESENTATION
# ==============================================================================

message("Section 18: cluster × donor representation...")

cluster_donor_counts <- atlas_sketch@meta.data %>%
  dplyr::count(
    Atlas_Cluster,
    Donor,
    name = "Cells"
  )

cluster_donor_summary <- cluster_donor_counts %>%
  group_by(Atlas_Cluster) %>%
  summarise(
    Number_of_Donors = n_distinct(Donor),
    Total_Cells = sum(Cells),
    Dominant_Donor = Donor[which.max(Cells)],
    Dominant_Donor_Cells = max(Cells),
    Dominant_Donor_Percentage =
      100 * max(Cells) / sum(Cells),
    .groups = "drop"
  )

write.csv(
  cluster_donor_counts,
  "results/tables/cluster_by_donor_counts.csv",
  row.names = FALSE
)

write.csv(
  cluster_donor_summary,
  "results/tables/cluster_donor_summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 19. CLUSTER × SAMPLE REPRESENTATION
# ==============================================================================

message("Section 19: cluster × sample representation...")

cluster_sample_counts <- atlas_sketch@meta.data %>%
  dplyr::count(
    Atlas_Cluster,
    GSM,
    name = "Cells"
  )

cluster_sample_summary <- cluster_sample_counts %>%
  group_by(Atlas_Cluster) %>%
  summarise(
    Number_of_Samples = n_distinct(GSM),
    Total_Cells = sum(Cells),
    Dominant_Sample = GSM[which.max(Cells)],
    Dominant_Sample_Cells = max(Cells),
    Dominant_Sample_Percentage =
      100 * max(Cells) / sum(Cells),
    .groups = "drop"
  )

write.csv(
  cluster_sample_counts,
  "results/tables/cluster_by_sample_counts.csv",
  row.names = FALSE
)

write.csv(
  cluster_sample_summary,
  "results/tables/cluster_sample_summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 20. DESCRIPTIVE CLUSTER MARKER DISCOVERY
# ==============================================================================

message("Section 20: descriptive cluster marker discovery...")

atlas_sketch <- JoinLayers(
  object = atlas_sketch,
  assay = "RNA"
)

cluster_markers <- FindAllMarkers(
  object = atlas_sketch,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

if (nrow(cluster_markers) == 0) {
  stop("No cluster markers were detected.")
}

cluster_markers <- cluster_markers %>%
  mutate(
    cluster = as.character(cluster)
  )

write.csv(
  cluster_markers,
  "results/tables/cluster_markers_all.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Top 10 markers per cluster
# ------------------------------------------------------------------------------

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10_markers,
  "results/tables/cluster_top10_markers.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Validate marker coverage
# ------------------------------------------------------------------------------

marker_clusters <- sort(
  unique(
    cluster_markers$cluster
  )
)

if (!setequal(marker_clusters, cluster_ids)) {
  
  missing_marker_clusters <-
    base::setdiff(
      cluster_ids,
      marker_clusters
    )
  
  extra_marker_clusters <-
    base::setdiff(
      marker_clusters,
      cluster_ids
    )
  
  warning(
    "Marker output does not contain exactly the atlas cluster set.\n",
    "Missing: ",
    paste(
      missing_marker_clusters,
      collapse = ", "
    ),
    "\nExtra: ",
    paste(
      extra_marker_clusters,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 21. CANONICAL LINEAGE MARKER PROGRAMS
# ==============================================================================

message("Section 21: canonical lineage marker programs...")

marker_programs <- list(
  
  T_Cell = c(
    "CD3D",
    "CD3E",
    "CD3G",
    "TRBC1",
    "TRBC2"
  ),
  
  CD8_T = c(
    "CD8A",
    "CD8B",
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2"
  ),
  
  CD4_T = c(
    "CD4",
    "IL7R",
    "LTB",
    "CCR7",
    "SELL"
  ),
  
  Treg = c(
    "FOXP3",
    "IL7R",
    "CTLA4",
    "IL2RA",
    "TNFRSF18",
    "TNFRSF4",
    "CCR8"
  ),
  
  NK = c(
    "NKG7",
    "GNLY",
    "FCGR3A",
    "KLRD1",
    "XCL1",
    "XCL2"
  ),
  
  MAIT = c(
    "TRAV1-2",
    "SLC4A10",
    "KLRB1",
    "IL23R",
    "CCL20",
    "DPP4"
  ),
  
  GammaDelta_T = c(
    "TRDC",
    "TRGC1",
    "TRGC2",
    "CD3D",
    "CD3E"
  ),
  
  B_Cell = c(
    "CD19",
    "MS4A1",
    "CD79A",
    "CD74",
    "HLA-DRA",
    "CD22"
  ),
  
  Plasma_Cell = c(
    "MZB1",
    "JCHAIN",
    "SDC1",
    "TNFRSF17",
    "DERL3",
    "IGHG1",
    "IGHG3",
    "IGHG4"
  ),
  
  Monocyte_Myeloid = c(
    "LYZ",
    "S100A8",
    "S100A9",
    "CTSS",
    "FCN1",
    "CTSD",
    "LGALS3",
    "FCGR3A",
    "CD14",
    "TREM1"
  ),
  
  Macrophage = c(
    "LYZ",
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "TREM2",
    "CD68",
    "MSR1"
  ),
  
  Dendritic_Cell = c(
    "FCER1A",
    "CST3",
    "CLEC10A",
    "CD1C",
    "HLA-DRA",
    "FCER1G"
  ),
  
  pDC = c(
    "CLEC4C",
    "IRF7",
    "TCF4",
    "LILRA4",
    "PTCRA"
  )
)


# ------------------------------------------------------------------------------
# Keep only genes actually present in atlas
# ------------------------------------------------------------------------------

marker_programs_present <- lapply(
  marker_programs,
  function(x) {
    intersect(
      x,
      rownames(atlas_sketch)
    )
  }
)

marker_program_sizes <- tibble(
  Program = names(marker_programs_present),
  Requested_Genes = lengths(marker_programs),
  Genes_Present = lengths(marker_programs_present)
)

print(marker_program_sizes)

write.csv(
  marker_program_sizes,
  "results/tables/marker_program_gene_availability.csv",
  row.names = FALSE
)


# ==============================================================================
# 22. CANONICAL MARKER HEATMAP
# ==============================================================================

message("Section 22: canonical marker heatmap...")

genes_for_heatmap <- unique(
  unlist(marker_programs_present)
)

genes_for_heatmap <- genes_for_heatmap[
  genes_for_heatmap %in% rownames(atlas_sketch)
]

if (length(genes_for_heatmap) == 0) {
  stop(
    "No canonical marker genes are present in the atlas."
  )
}


# ------------------------------------------------------------------------------
# 22.1 AVERAGE EXPRESSION BY CLUSTER
# ------------------------------------------------------------------------------

cluster_average_expression <- AverageExpression(
  object = atlas_sketch,
  assays = "RNA",
  features = genes_for_heatmap,
  group.by = "Atlas_Cluster",
  layer = "data",
  verbose = FALSE
)$RNA


# ------------------------------------------------------------------------------
# 22.2 RESTORE ORIGINAL CLUSTER IDs
# ------------------------------------------------------------------------------

returned_cluster_ids <- colnames(
  cluster_average_expression
)

if (
  !all(
    grepl(
      "^g[0-9]+$",
      returned_cluster_ids
    )
  )
) {
  
  stop(
    "Unexpected AverageExpression cluster column names: ",
    paste(
      returned_cluster_ids,
      collapse = ", "
    )
  )
}

restored_cluster_ids <- sub(
  "^g",
  "",
  returned_cluster_ids
)

if (
  !setequal(
    restored_cluster_ids,
    cluster_ids
  )
) {
  
  missing_clusters <- base::setdiff(
    cluster_ids,
    restored_cluster_ids
  )
  
  extra_clusters <- base::setdiff(
    restored_cluster_ids,
    cluster_ids
  )
  
  stop(
    "AverageExpression cluster IDs do not match atlas clusters.\n",
    "Missing: ",
    paste(
      missing_clusters,
      collapse = ", "
    ),
    "\nExtra: ",
    paste(
      extra_clusters,
      collapse = ", "
    )
  )
}

cluster_average_expression <-
  cluster_average_expression[
    ,
    match(
      cluster_ids,
      restored_cluster_ids
    ),
    drop = FALSE
  ]

colnames(
  cluster_average_expression
) <- cluster_ids


# ------------------------------------------------------------------------------
# 22.3 FINAL CLUSTER-MATRIX VALIDATION
# ------------------------------------------------------------------------------

if (
  ncol(cluster_average_expression) !=
  n_atlas_clusters
) {
  
  stop(
    "Cluster-average matrix contains ",
    ncol(cluster_average_expression),
    " clusters, but atlas contains ",
    n_atlas_clusters,
    "."
  )
}

if (
  !identical(
    colnames(cluster_average_expression),
    cluster_ids
  )
) {
  
  stop(
    "Cluster-average matrix columns do not match detected atlas clusters."
  )
}


# ------------------------------------------------------------------------------
# 22.4 ROW-WISE Z-SCORE
# ------------------------------------------------------------------------------

heatmap_matrix <- t(
  scale(
    t(cluster_average_expression)
  )
)

heatmap_matrix[is.na(heatmap_matrix)] <- 0


# ------------------------------------------------------------------------------
# 22.5 HEATMAP
# ------------------------------------------------------------------------------

pdf(
  "results/tables/canonical_marker_heatmap.pdf",
  width = 14,
  height = 12
)

pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA,
  main = "Canonical lineage marker expression by atlas cluster"
)

dev.off()


# ==============================================================================
# 23. MULTI-METHOD ANNOTATION EVIDENCE
# ==============================================================================

message("Section 23: multi-method annotation evidence...")


# ==============================================================================
# 23.1 MARKER-PROGRAM GENE OVERLAP
# ==============================================================================

message("23.1: marker overlap...")

top50_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 50) %>%
  ungroup()


# ------------------------------------------------------------------------------
# Construct cluster × program overlap evidence
# ------------------------------------------------------------------------------

marker_overlap_table <- purrr::map_dfr(
  cluster_ids,
  function(cluster_id) {
    
    cluster_genes <- top50_markers %>%
      filter(
        cluster == cluster_id
      ) %>%
      pull(gene)
    
    purrr::map_dfr(
      names(marker_programs_present),
      function(program_name) {
        
        program_genes <- marker_programs_present[[program_name]]
        
        overlap_genes <- intersect(
          cluster_genes,
          program_genes
        )
        
        tibble(
          Cluster = cluster_id,
          Program = program_name,
          Program_Genes_Present = list(
            program_genes
          ),
          Top50_Genes = list(
            cluster_genes
          ),
          Overlap_Genes = list(
            overlap_genes
          ),
          Overlap_Count = length(
            overlap_genes
          ),
          Program_Genes_Count = length(
            program_genes
          ),
          Overlap_Fraction =
            ifelse(
              length(program_genes) > 0,
              length(overlap_genes) /
                length(program_genes),
              NA_real_
            )
        )
      }
    )
  }
)


# ------------------------------------------------------------------------------
# Flat CSV representation
# ------------------------------------------------------------------------------

marker_overlap_table_csv <- marker_overlap_table %>%
  ungroup() %>%
  mutate(
    Program_Genes_Present =
      vapply(
        Program_Genes_Present,
        paste,
        collapse = ", ",
        FUN.VALUE = character(1)
      ),
    Top50_Genes =
      vapply(
        Top50_Genes,
        paste,
        collapse = ", ",
        FUN.VALUE = character(1)
      ),
    Overlap_Genes =
      vapply(
        Overlap_Genes,
        paste,
        collapse = ", ",
        FUN.VALUE = character(1)
      )
  )

write.csv(
  marker_overlap_table_csv,
  "results/tables/marker_program_top50_overlap.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.2 ADD MODULE SCORES
# ==============================================================================

message("23.2: module scores...")

module_score_columns <- character(0)

for (program_name in names(marker_programs_present)) {
  
  genes <- marker_programs_present[[program_name]]
  
  if (length(genes) < 2) {
    
    warning(
      "Skipping module score for ",
      program_name,
      ": fewer than two genes are present."
    )
    
    next
  }
  
  atlas_sketch <- AddModuleScore(
    object = atlas_sketch,
    features = list(genes),
    name = paste0(
      "Module_",
      program_name
    ),
    assay = "RNA",
    slot = "data",
    search = FALSE
  )
  
  new_column <- paste0(
    "Module_",
    program_name,
    "1"
  )
  
  if (
    new_column %in%
    colnames(atlas_sketch@meta.data)
  ) {
    
    module_score_columns <- c(
      module_score_columns,
      new_column
    )
  }
}


# ==============================================================================
# 23.3 CLUSTER-LEVEL MODULE SCORES
# ==============================================================================

message("23.3: cluster-level module scores...")

if (length(module_score_columns) > 0) {
  
  module_score_table <- atlas_sketch@meta.data %>%
    group_by(Atlas_Cluster) %>%
    summarise(
      across(
        all_of(module_score_columns),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    )
  
} else {
  
  module_score_table <- tibble(
    Atlas_Cluster = cluster_ids
  )
}

write.csv(
  module_score_table,
  "results/tables/cluster_module_scores.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.4 BEST MODULE-SCORE PROGRAM PER CLUSTER
# ==============================================================================

message("23.4: strongest module-score programs...")

if (length(module_score_columns) > 0) {
  
  module_score_long <- module_score_table %>%
    pivot_longer(
      cols = all_of(module_score_columns),
      names_to = "Program",
      values_to = "Score"
    ) %>%
    mutate(
      Program = sub(
        "^Module_",
        "",
        Program
      ),
      Program = sub(
        "1$",
        "",
        Program
      )
    )
  
  strongest_module_program <- module_score_long %>%
    group_by(Atlas_Cluster) %>%
    arrange(
      desc(Score),
      .by_group = TRUE
    ) %>%
    mutate(
      Rank = row_number()
    ) %>%
    slice_head(n = 3) %>%
    ungroup()
  
  write.csv(
    strongest_module_program,
    "results/tables/strongest_module_programs.csv",
    row.names = FALSE
  )
  
} else {
  
  strongest_module_program <- tibble()
  
  warning(
    "No module-score programs were successfully generated."
  )
}


# ==============================================================================
# 23.5 SINGLE-R / MONACO IMMUNE REFERENCE
# ==============================================================================

message("23.5: loading Monaco immune reference...")

monaco_ref <- MonacoImmuneData()

monaco_labels <- colData(
  monaco_ref
)$label.fine

if (is.null(monaco_labels)) {
  
  stop(
    "Fine Monaco labels were not found in the reference."
  )
}


# ------------------------------------------------------------------------------
# Cluster-level expression matrix
# ------------------------------------------------------------------------------

test_expression <- SeuratObject::LayerData(
  object = atlas_sketch,
  assay = "RNA",
  layer = "data"
)

singleR_clusters <- as.character(
  atlas_sketch$Atlas_Cluster
)

if (!identical(
  colnames(test_expression),
  colnames(atlas_sketch)
)) {
  
  stop(
    "LayerData cell order does not match Seurat object cell order."
  )
}


# ------------------------------------------------------------------------------
# Run SingleR
#
# prune = FALSE intentionally retained.
#
# Therefore the expected result contains:
#   labels
#   delta.next
#
# and does NOT require:
#   pruned.labels
# ------------------------------------------------------------------------------

message(
  "Running SingleR at cluster level..."
)

singleR_results <- SingleR(
  test = test_expression,
  ref = monaco_ref,
  labels = monaco_labels,
  clusters = singleR_clusters,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  prune = FALSE
)

singleR_results <- as.data.frame(
  singleR_results
)

singleR_results$Cluster <-
  rownames(singleR_results)

rownames(singleR_results) <- NULL


# ------------------------------------------------------------------------------
# Validate number of cluster results
# ------------------------------------------------------------------------------

if (
  nrow(singleR_results) !=
  n_atlas_clusters
) {
  
  stop(
    "SingleR returned ",
    nrow(singleR_results),
    " cluster results, but atlas contains ",
    n_atlas_clusters,
    " clusters."
  )
}


# ------------------------------------------------------------------------------
# Standardize cluster IDs
# ------------------------------------------------------------------------------

singleR_results$Cluster <-
  as.character(
    singleR_results$Cluster
  )


# ------------------------------------------------------------------------------
# Confirm all expected clusters are present
# ------------------------------------------------------------------------------

if (
  !setequal(
    singleR_results$Cluster,
    cluster_ids
  )
) {
  
  missing_clusters <- base::setdiff(
    cluster_ids,
    singleR_results$Cluster
  )
  
  extra_clusters <- base::setdiff(
    singleR_results$Cluster,
    cluster_ids
  )
  
  stop(
    "SingleR cluster IDs do not match atlas clusters.\n",
    "Missing: ",
    paste(
      missing_clusters,
      collapse = ", "
    ),
    "\nExtra: ",
    paste(
      extra_clusters,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# Explicitly reorder SingleR results to atlas cluster order
# ------------------------------------------------------------------------------

singleR_results <- singleR_results[
  match(
    cluster_ids,
    singleR_results$Cluster
  ),
  ,
  drop = FALSE
]

write.csv(
  singleR_results,
  "results/tables/singler_monaco_results.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.6 COMPACT SINGLE-R SUMMARY
# ==============================================================================

message(
  "23.6: constructing compact SingleR summary..."
)

singleR_summary <- singleR_results[
  ,
  c(
    "Cluster",
    "labels",
    "delta.next"
  ),
  drop = FALSE
]

colnames(singleR_summary) <- c(
  "Atlas_Cluster",
  "SingleR_Label",
  "SingleR_Delta_Next"
)

singleR_summary$Atlas_Cluster <-
  as.character(
    singleR_summary$Atlas_Cluster
  )

write.csv(
  singleR_summary,
  "results/tables/singler_summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.7 COMBINED ANNOTATION EVIDENCE TABLE
# ==============================================================================

message(
  "23.7: constructing combined evidence table..."
)


# ------------------------------------------------------------------------------
# Best module score
# ------------------------------------------------------------------------------

if (
  nrow(strongest_module_program) > 0
) {
  
  best_module <- strongest_module_program %>%
    group_by(Atlas_Cluster) %>%
    arrange(
      Rank,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(
      Atlas_Cluster,
      Best_Module_Program = Program,
      Best_Module_Score = Score
    )
  
} else {
  
  best_module <- tibble(
    Atlas_Cluster = cluster_ids,
    Best_Module_Program = NA_character_,
    Best_Module_Score = NA_real_
  )
}


# ------------------------------------------------------------------------------
# Best marker-program overlap
# ------------------------------------------------------------------------------

best_overlap <- marker_overlap_table %>%
  group_by(Cluster) %>%
  arrange(
    desc(Overlap_Fraction),
    desc(Overlap_Count),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    Atlas_Cluster = as.character(Cluster),
    Best_Top50_Program = Program,
    Best_Top50_Overlap_Count = Overlap_Count,
    Best_Top50_Overlap_Fraction =
      Overlap_Fraction
  )


# ------------------------------------------------------------------------------
# Combined table
# ------------------------------------------------------------------------------

combined_annotation_evidence <- tibble(
  Atlas_Cluster = cluster_ids
) %>%
  
  left_join(
    best_module,
    by = "Atlas_Cluster"
  ) %>%
  
  left_join(
    best_overlap,
    by = "Atlas_Cluster"
  ) %>%
  
  left_join(
    singleR_summary,
    by = "Atlas_Cluster"
  ) %>%
  
  left_join(
    cluster_donor_summary %>%
      mutate(
        Atlas_Cluster =
          as.character(Atlas_Cluster)
      ) %>%
      select(
        Atlas_Cluster,
        Number_of_Donors,
        Dominant_Donor,
        Dominant_Donor_Percentage
      ),
    by = "Atlas_Cluster"
  ) %>%
  
  left_join(
    cluster_sample_summary %>%
      mutate(
        Atlas_Cluster =
          as.character(Atlas_Cluster)
      ) %>%
      select(
        Atlas_Cluster,
        Number_of_Samples,
        Dominant_Sample,
        Dominant_Sample_Percentage
      ),
    by = "Atlas_Cluster"
  )


# ------------------------------------------------------------------------------
# Final combined-table validation
# ------------------------------------------------------------------------------

if (
  nrow(combined_annotation_evidence) !=
  n_atlas_clusters
) {
  
  stop(
    "Combined annotation evidence contains ",
    nrow(combined_annotation_evidence),
    " rows, but atlas contains ",
    n_atlas_clusters,
    " clusters."
  )
}

if (
  !setequal(
    combined_annotation_evidence$Atlas_Cluster,
    cluster_ids
  )
) {
  
  stop(
    "Combined annotation evidence does not contain exactly the atlas clusters."
  )
}

write.csv(
  combined_annotation_evidence,
  "results/tables/combined_annotation_evidence.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.8 TOP MARKERS PER CLUSTER FOR MANUAL INSPECTION
# ==============================================================================

message(
  "23.8: generating manual annotation inspection table..."
)


# ------------------------------------------------------------------------------
# Top 20 marker genes per cluster
# ------------------------------------------------------------------------------

top20_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  summarise(
    Top20_Markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  transmute(
    Atlas_Cluster = as.character(cluster),
    Top20_Markers = Top20_Markers
  )


# ------------------------------------------------------------------------------
# Right join to ensure every atlas cluster is represented
# ------------------------------------------------------------------------------

manual_inspection_table <- top20_markers %>%
  right_join(
    combined_annotation_evidence,
    by = "Atlas_Cluster"
  )


# ------------------------------------------------------------------------------
# Explicit robust cluster ordering
# ------------------------------------------------------------------------------
#
# The previous implementation used:
#
#   match(Atlas_Cluster, cluster_ids)
#
# which produced the duplicated.default() error in the user's environment.
#
# Here the ordering key is generated directly from the actual cluster IDs.
# For the current numeric cluster IDs this gives:
#
#   0, 1, 2, ... 19
#
# ------------------------------------------------------------------------------

manual_inspection_table <- manual_inspection_table %>%
  mutate(
    .Atlas_Cluster_Order =
      match(
        as.character(Atlas_Cluster),
        as.character(cluster_ids)
      )
  ) %>%
  arrange(
    .Atlas_Cluster_Order
  ) %>%
  select(
    - .Atlas_Cluster_Order
  )


# ------------------------------------------------------------------------------
# Validate manual inspection table
# ------------------------------------------------------------------------------

if (
  nrow(manual_inspection_table) !=
  n_atlas_clusters
) {
  
  stop(
    "Manual inspection table contains ",
    nrow(manual_inspection_table),
    " rows, but atlas contains ",
    n_atlas_clusters,
    " clusters."
  )
}

if (
  !setequal(
    manual_inspection_table$Atlas_Cluster,
    cluster_ids
  )
) {
  
  stop(
    "Manual inspection table does not contain exactly the atlas clusters."
  )
}

write.csv(
  manual_inspection_table,
  "results/tables/MANUAL_ANNOTATION_INSPECTION_TABLE.csv",
  row.names = FALSE
)


# ==============================================================================
# 23.9 SAVE ANNOTATION EVIDENCE CHECKPOINT
# ==============================================================================

message(
  "23.9: saving evidence checkpoint..."
)

annotation_evidence <- list(
  
  atlas_summary =
    atlas_summary,
  
  cluster_ids =
    cluster_ids,
  
  n_atlas_clusters =
    n_atlas_clusters,
  
  cluster_phase_counts =
    cluster_phase_counts,
  
  cluster_phase_within_state =
    cluster_phase_within_state,
  
  cluster_phase_within_cluster =
    cluster_phase_within_cluster,
  
  cluster_donor_counts =
    cluster_donor_counts,
  
  cluster_donor_summary =
    cluster_donor_summary,
  
  cluster_sample_counts =
    cluster_sample_counts,
  
  cluster_sample_summary =
    cluster_sample_summary,
  
  cluster_markers =
    cluster_markers,
  
  top10_markers =
    top10_markers,
  
  top50_markers =
    top50_markers,
  
  marker_programs =
    marker_programs_present,
  
  marker_program_sizes =
    marker_program_sizes,
  
  marker_overlap_table =
    marker_overlap_table,
  
  module_score_table =
    module_score_table,
  
  strongest_module_program =
    strongest_module_program,
  
  singleR_results =
    singleR_results,
  
  singleR_summary =
    singleR_summary,
  
  combined_annotation_evidence =
    combined_annotation_evidence,
  
  manual_inspection_table =
    manual_inspection_table
)

saveRDS(
  annotation_evidence,
  evidence_output_file
)

message(
  "\nEvidence checkpoint saved to:\n",
  evidence_output_file
)


# ==============================================================================
# 23.10 STOP BEFORE FINAL ANNOTATION
# ==============================================================================

message(
  "\n============================================================\n",
  "ANNOTATION EVIDENCE CHECKPOINT COMPLETE\n",
  "============================================================\n",
  "\nDetected atlas clusters: ",
  n_atlas_clusters,
  "\nClusters: ",
  paste(
    cluster_ids,
    collapse = ", "
  ),
  "\n\nReview:\n",
  "  results/tables/MANUAL_ANNOTATION_INSPECTION_TABLE.csv\n",
  "  results/tables/combined_annotation_evidence.csv\n",
  "  results/tables/cluster_top10_markers.csv\n",
  "  results/tables/canonical_marker_heatmap.pdf\n",
  "  results/tables/singler_summary.csv\n",
  "\nFINAL LABELS HAVE NOT BEEN ASSIGNED.\n",
  "============================================================\n"
)


# ==============================================================================
# 24. FINAL CELL-TYPE ANNOTATION
# ==============================================================================

# ------------------------------------------------------------------------------
# MANUAL ADJUDICATION
# ------------------------------------------------------------------------------
#
# IMPORTANT:
#
# Before running this section, inspect the evidence generated above.
#
# Replace the NA values below with the FINAL adjudicated identity for EVERY
# detected cluster.
#
# DO NOT automatically assign labels from:
#   - highest module score
#   - highest marker overlap
#   - SingleR label
#   - top marker alone
#
# These are evidence streams that must be adjudicated together.
#
# ------------------------------------------------------------------------------


final_cluster_labels <- tibble(
  Atlas_Cluster = as.character(0:19),
  
  Final_Cell_Type = c(
    "CD8_T",               # 0
    "NK",                  # 1
    "CD4_T",               # 2
    "MAIT",                # 3
    "CD8_T",               # 4
    "NK",                  # 5
    "MAIT",                # 6
    "CD8_T",               # 7
    "CD8_T",               # 8
    "B_Cell",              # 9
    "Unresolved_Lymphoid", # 10
    "Treg",                # 11
    "NK",                  # 12
    "NK",                  # 13
    "NK",                  # 14
    "Monocyte_Myeloid",    # 15
    "Monocyte_Myeloid",    # 16
    "Dendritic_Cell",      # 17
    "Plasma_Cell",         # 18
    "pDC"                  # 19
  ),
  
  Confidence = c(
    "High",     # 0
    "High",     # 1
    "High",     # 2
    "High",     # 3
    "High",     # 4
    "High",     # 5
    "Moderate", # 6
    "Moderate", # 7
    "High",     # 8
    "High",     # 9
    "Low",      # 10
    "High",     # 11
    "High",     # 12
    "Moderate", # 13
    "High",     # 14
    "High",     # 15
    "High",     # 16
    "High",     # 17
    "High",     # 18
    "High"      # 19
  ),
  
  Rationale = c(
    "CD8A/CD8B and CCL5/GZMK support a CD8 T-cell identity; CD8 module, marker overlap, and SingleR are concordant.",
    "GNLY/FGFBP2/PRF1/GZMB/KLRD1/KLRF1 form a strong cytotoxic NK-cell program with concordant module and SingleR evidence.",
    "CCR7/IL7R/CD4/LTB/S1PR1/CCR6 support a CD4 T-cell identity; reference annotation indicates a Th17-like state.",
    "TRAV1-2/SLC4A10/KLRB1/IL23R/CCL20/DPP4 provide highly concordant MAIT evidence.",
    "CD8B/CD3D/GZMA/CCL5 support a CD8 T-cell identity with concordant module, marker and reference evidence.",
    "KLRF1/KLRC1/NCAM1/XCL1/IL2RB/CD160 support an NK-cell identity with concordant module and SingleR evidence.",
    "TRAV1-2 and IL7R/AQP3-associated lymphocyte features support MAIT; SingleR strongly agrees, although the CD4 module is also elevated.",
    "CD8A with PDCD1/TNFRSF9/ITGA4 supports an activated CD8 T-cell state; CD8 module and SingleR agree.",
    "CD8A/CD8B/CCL5/CST7 provide strong CD8 T-cell evidence across marker, module and reference approaches.",
    "CD19/MS4A1/CD79A/CD22/CD24/FCRL1/TCL1A provide a clear B-cell program with complete marker-program overlap.",
    "Conflicting lymphoid and immunoglobulin signals, together with extreme sample dominance, prevent confident lineage assignment; retained as unresolved rather than forced.",
    "FOXP3/CTLA4/TNFRSF4/TNFRSF18/ICOS/TIGIT provide a highly concordant Treg program with reference support.",
    "GNLY/FGFBP2/GZMB/GZMH/PRF1/KLRD1 provide a strong cytotoxic NK-cell program with concordant reference evidence.",
    "NK-associated cytotoxic genes and reference/module evidence favor NK, although TRDC/TRGC expression suggests a possible NK/γδ-T boundary population.",
    "Strong NK/cytotoxic program and SingleR support NK; TRDC/TRGC expression is noted but does not outweigh the broader identity evidence.",
    "FCN1/FCAR/FPR2/TREM1/IL1B/CXCL2/CXCL3/OLR1 support an inflammatory monocyte/myeloid identity.",
    "FCN1/S100A8/S100A9/CD14/VCAN/CSF3R/HK3 provide a strong classical monocyte program.",
    "CD1C/CLEC10A/FCER1A/FLT3 and C1Q-associated genes support a dendritic-cell identity with concordant reference evidence.",
    "JCHAIN/MZB1/SDC1/TNFRSF17/DERL3 and immunoglobulin genes provide a clear plasma-cell program.",
    "LILRA4/CLEC4C/TCF4/PTCRA support a plasmacytoid dendritic-cell identity with concordant module and SingleR evidence."
  )
)


# ==============================================================================
# 24.1 VALIDATE MANUAL ANNOTATION TABLE
# ==============================================================================

message(
  "Section 24.1: validating final annotation table..."
)


# ------------------------------------------------------------------------------
# Exact cluster coverage
# ------------------------------------------------------------------------------

if (
  !setequal(
    final_cluster_labels$Atlas_Cluster,
    cluster_ids
  )
) {
  
  missing_clusters <- base::setdiff(
    cluster_ids,
    final_cluster_labels$Atlas_Cluster
  )
  
  extra_clusters <- base::setdiff(
    final_cluster_labels$Atlas_Cluster,
    cluster_ids
  )
  
  stop(
    "Final annotation table does not exactly match atlas clusters.\n",
    "Missing clusters: ",
    paste(
      missing_clusters,
      collapse = ", "
    ),
    "\nExtra clusters: ",
    paste(
      extra_clusters,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# No duplicated cluster assignments
# ------------------------------------------------------------------------------

if (
  anyDuplicated(
    final_cluster_labels$Atlas_Cluster
  ) > 0
) {
  
  stop(
    "Duplicate Atlas_Cluster entries found in final annotation table."
  )
}


# ------------------------------------------------------------------------------
# No missing labels
# ------------------------------------------------------------------------------

if (
  any(
    is.na(
      final_cluster_labels$Final_Cell_Type
    ) |
    final_cluster_labels$Final_Cell_Type == ""
  )
) {
  
  stop(
    "\nFINAL ANNOTATION NOT COMPLETE.\n\n",
    "Every atlas cluster must receive an explicit Final_Cell_Type.\n",
    "This is intentional: the script will not invent labels."
  )
}


# ------------------------------------------------------------------------------
# No missing confidence
# ------------------------------------------------------------------------------

if (
  any(
    is.na(
      final_cluster_labels$Confidence
    ) |
    final_cluster_labels$Confidence == ""
  )
) {
  
  stop(
    "Every final annotation must have an explicit Confidence value."
  )
}


# ------------------------------------------------------------------------------
# No missing rationale
# ------------------------------------------------------------------------------

if (
  any(
    is.na(
      final_cluster_labels$Rationale
    ) |
    final_cluster_labels$Rationale == ""
  )
) {
  
  stop(
    "Every final annotation must have an explicit Rationale."
  )
}


# ==============================================================================
# 24.2 VALIDATE CONFIDENCE VOCABULARY
# ==============================================================================

allowed_confidence <- c(
  "High",
  "Moderate",
  "Low"
)

invalid_confidence <- base::setdiff(
  unique(
    final_cluster_labels$Confidence
  ),
  allowed_confidence
)

if (
  length(invalid_confidence) > 0
) {
  
  stop(
    "Invalid confidence value(s): ",
    paste(
      invalid_confidence,
      collapse = ", "
    ),
    "\nAllowed values: ",
    paste(
      allowed_confidence,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 24.3 VALIDATE FINAL CELL-TYPE VOCABULARY
# ==============================================================================

allowed_cell_types <- c(
  "T_Cell",
  "CD8_T",
  "CD4_T",
  "Treg",
  "NK",
  "MAIT",
  "GammaDelta_T",
  "B_Cell",
  "Plasma_Cell",
  "Monocyte_Myeloid",
  "Macrophage",
  "Dendritic_Cell",
  "pDC",
  "Unresolved_Lymphoid"
)

invalid_cell_types <- base::setdiff(
  unique(
    final_cluster_labels$Final_Cell_Type
  ),
  allowed_cell_types
)

if (
  length(invalid_cell_types) > 0
) {
  
  stop(
    "Invalid final cell-type label(s): ",
    paste(
      invalid_cell_types,
      collapse = ", "
    ),
    "\nAllowed labels: ",
    paste(
      allowed_cell_types,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 24.4 ORDER FINAL ANNOTATIONS
# ==============================================================================

final_cluster_labels <- final_cluster_labels %>%
  mutate(
    .Atlas_Cluster_Order =
      match(
        as.character(Atlas_Cluster),
        as.character(cluster_ids)
      )
  ) %>%
  arrange(
    .Atlas_Cluster_Order
  ) %>%
  select(
    - .Atlas_Cluster_Order
  )


# ==============================================================================
# 24.5 PROPAGATE CLUSTER LABELS TO ALL SKETCH CELLS
# ==============================================================================

message(
  "Propagating final cluster annotations to sketch cells..."
)

label_lookup <- setNames(
  final_cluster_labels$Final_Cell_Type,
  final_cluster_labels$Atlas_Cluster
)

confidence_lookup <- setNames(
  final_cluster_labels$Confidence,
  final_cluster_labels$Atlas_Cluster
)

rationale_lookup <- setNames(
  final_cluster_labels$Rationale,
  final_cluster_labels$Atlas_Cluster
)

atlas_sketch$Final_Cell_Type <-
  unname(
    label_lookup[
      as.character(
        atlas_sketch$Atlas_Cluster
      )
    ]
  )

atlas_sketch$Annotation_Confidence <-
  unname(
    confidence_lookup[
      as.character(
        atlas_sketch$Atlas_Cluster
      )
    ]
  )

atlas_sketch$Annotation_Rationale <-
  unname(
    rationale_lookup[
      as.character(
        atlas_sketch$Atlas_Cluster
      )
    ]
  )


# ==============================================================================
# 24.6 FINAL VALIDATION
# ==============================================================================

message(
  "Validating final annotation..."
)


# ------------------------------------------------------------------------------
# No NA labels
# ------------------------------------------------------------------------------

if (
  any(
    is.na(
      atlas_sketch$Final_Cell_Type
    )
  )
) {
  
  stop(
    "NA final cell-type labels detected after propagation."
  )
}


# ------------------------------------------------------------------------------
# All clusters represented
# ------------------------------------------------------------------------------

final_cluster_check <- atlas_sketch@meta.data %>%
  distinct(
    Atlas_Cluster,
    Final_Cell_Type
  ) %>%
  mutate(
    .Atlas_Cluster_Order =
      match(
        as.character(Atlas_Cluster),
        as.character(cluster_ids)
      )
  ) %>%
  arrange(
    .Atlas_Cluster_Order
  ) %>%
  select(
    - .Atlas_Cluster_Order
  )

if (
  nrow(final_cluster_check) !=
  n_atlas_clusters
) {
  
  stop(
    "Final annotation does not contain exactly one label per atlas cluster."
  )
}


# ------------------------------------------------------------------------------
# Exactly one cell type per cluster
# ------------------------------------------------------------------------------

cluster_label_counts <- atlas_sketch@meta.data %>%
  dplyr::distinct(
    Atlas_Cluster,
    Final_Cell_Type
  ) %>%
  dplyr::count(
    Atlas_Cluster,
    name = "Number_of_Labels"
  )

if (
  any(
    cluster_label_counts$Number_of_Labels != 1
  )
) {
  
  stop(
    "At least one cluster has multiple final cell-type labels."
  )
}


# ==============================================================================
# 24.7 FINAL ANNOTATION SUMMARY
# ==============================================================================

final_annotation_summary <- atlas_sketch@meta.data %>%
  dplyr::count(
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  dplyr::arrange(
    desc(Cells)
  )

final_cluster_annotation_summary <- atlas_sketch@meta.data %>%
  dplyr::count(
    Atlas_Cluster,
    Final_Cell_Type,
    Annotation_Confidence,
    name = "Cells"
  ) %>%
  dplyr::mutate(
    .Atlas_Cluster_Order =
      match(
        as.character(Atlas_Cluster),
        as.character(cluster_ids)
      )
  ) %>%
  dplyr::arrange(
    .Atlas_Cluster_Order
  ) %>%
  dplyr::select(
    - .Atlas_Cluster_Order
  )

write.csv(
  final_cluster_labels,
  "results/tables/final_cluster_annotation_table.csv",
  row.names = FALSE
)

write.csv(
  final_annotation_summary,
  "results/tables/final_cell_type_summary.csv",
  row.names = FALSE
)

write.csv(
  final_cluster_annotation_summary,
  "results/tables/final_cluster_annotation_summary.csv",
  row.names = FALSE
)


# ==============================================================================
# 24.8 FINAL ANNOTATION VALIDATION REPORT
# ==============================================================================

final_validation_report <- tibble(
  Metric = c(
    "Atlas cells",
    "Atlas genes",
    "Atlas clusters",
    "Final annotated clusters",
    "Clusters with NA labels",
    "Distinct final cell types"
  ),
  Value = c(
    ncol(atlas_sketch),
    nrow(atlas_sketch),
    n_atlas_clusters,
    n_distinct(atlas_sketch$Atlas_Cluster),
    sum(
      is.na(
        atlas_sketch$Final_Cell_Type
      )
    ),
    n_distinct(
      atlas_sketch$Final_Cell_Type
    )
  )
)

print(final_validation_report)

write.csv(
  final_validation_report,
  "results/tables/final_annotation_validation.csv",
  row.names = FALSE
)


# ==============================================================================
# 24.9 SAVE FINAL ANNOTATED ATLAS
# ==============================================================================

message(
  "Saving final annotated atlas..."
)

saveRDS(
  atlas_sketch,
  final_output_file
)

message(
  "\n============================================================\n",
  "PHASE 3 SCRIPT 02 COMPLETE\n",
  "============================================================\n",
  "\nFinal atlas:\n",
  final_output_file,
  "\n\nCells: ",
  ncol(atlas_sketch),
  "\nGenes: ",
  nrow(atlas_sketch),
  "\nClusters: ",
  n_atlas_clusters,
  "\nAnnotated clusters: ",
  n_distinct(
    atlas_sketch$Atlas_Cluster
  ),
  "\nNA cell-type labels: ",
  sum(
    is.na(
      atlas_sketch$Final_Cell_Type
    )
  ),
  "\n\nFinal cell-type counts:\n"
)

print(
  final_annotation_summary
)

message(
  "\n============================================================\n"
)
