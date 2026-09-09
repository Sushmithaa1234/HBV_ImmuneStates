# ============================================================================
# PHASE 3 — CLUSTER MARKERS + CELL-TYPE ANNOTATION
# ============================================================================
#
# Script: phase_03_02_cluster_markers_annotation.R
#
# Sections:
#   17 — Cluster × clinical state
#   18 — Cluster × donor dominance
#   19 — Cluster × sample dominance
#   20 — Cluster marker discovery
#   21 — Canonical lineage markers
#   22 — Canonical marker heatmap
#   23 — Multi-method annotation
#   24 — Final cell-type annotation of the sketch
#
# Input checkpoint:
#   results/rds_objects/phase3_atlas_sketch_umap_clustered.rds
#
# Final output:
#   results/rds_objects/phase3_atlas_sketch_final_annotation.rds
#
# Biological framing:
#   Clinical states are treated as cohort-level comparisons rather than a
#   presumed linear disease trajectory.
#
# IMPORTANT:
#   The atlas contains 20,000 representative sketch cells.
#   Annotation is performed at cluster level and then propagated to cells.
#   Module scores, marker overlap, and SingleR are evidence streams rather
#   than automatic biological truth.
#
# ============================================================================


# ============================================================================
# ENVIRONMENT
# ============================================================================

cat("============================================================\n")
cat("PHASE 3 — CLUSTER MARKERS + CELL-TYPE ANNOTATION\n")
cat("============================================================\n\n")

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
library(Matrix)

setwd(here())

set.seed(12345)

cat("✓ Environment initialized\n\n")


# ============================================================================
# LOAD ATLAS CHECKPOINT
# ============================================================================

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_umap_clustered.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

cat(
  "Atlas cells:",
  ncol(atlas_sketch),
  "\n"
)

cat(
  "Atlas genes:",
  nrow(atlas_sketch),
  "\n"
)

cat(
  "Clusters:",
  length(unique(atlas_sketch$seurat_clusters)),
  "\n\n"
)

if (ncol(atlas_sketch) != 20000) {
  
  stop(
    "ERROR: Expected 20,000 atlas sketch cells."
  )
}

n_atlas_clusters <- length(
  unique(
    as.character(
      atlas_sketch$seurat_clusters
    )
  )
)

if (n_atlas_clusters < 2) {
  stop(
    "ERROR: Atlas clustering produced fewer than 2 clusters."
  )
}

cat(
  "Atlas clustering produced ",
  n_atlas_clusters,
  " clusters.\n",
  sep = ""
)


# ============================================================================
# SECTION 17 — CLUSTER × CLINICAL STATE
# ============================================================================
#
# Purpose:
#   Assess the distribution of the 17 global atlas clusters across the five
#   HBV clinical states.
#
# Outputs:
#   1. Raw cell-count table
#   2. Within-clinical-state cluster proportion table
#   3. Within-cluster clinical-state composition table
#   4. Cluster × clinical-state heatmap
#
# ============================================================================

cat("=== SECTION 17: CLUSTER × CLINICAL STATE ===\n\n")


# ----------------------------------------------------------------------------
# 17.1 — Validate metadata
# ----------------------------------------------------------------------------

required_columns <- c(
  "seurat_clusters",
  "Phase"
)

missing_columns <- setdiff(
  required_columns,
  colnames(atlas_sketch@meta.data)
)

if (length(missing_columns) > 0) {
  
  stop(
    "ERROR: Required metadata columns are missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 17.2 — Raw cluster × clinical-state counts
# ----------------------------------------------------------------------------

cluster_phase_counts <-
  atlas_sketch@meta.data %>%
  dplyr::count(
    seurat_clusters,
    Phase,
    name = "Cells"
  ) %>%
  dplyr::arrange(
    as.numeric(as.character(seurat_clusters)),
    Phase
  )

cat(
  "Cluster × clinical-state cell counts:\n\n"
)

print(
  cluster_phase_counts
)


# ----------------------------------------------------------------------------
# 17.3 — Cluster proportions within each clinical state
# ----------------------------------------------------------------------------

cluster_phase_proportions <-
  cluster_phase_counts %>%
  group_by(Phase) %>%
  mutate(
    Proportion =
      Cells / sum(Cells),
    Percent =
      100 * Proportion
  ) %>%
  ungroup()


# ----------------------------------------------------------------------------
# 17.4 — Clinical-state composition within each cluster
# ----------------------------------------------------------------------------

cluster_phase_composition <-
  cluster_phase_counts %>%
  group_by(seurat_clusters) %>%
  mutate(
    Proportion =
      Cells / sum(Cells),
    Percent =
      100 * Proportion
  ) %>%
  ungroup()


# ----------------------------------------------------------------------------
# 17.5 — Save tables
# ----------------------------------------------------------------------------

write_csv(
  cluster_phase_counts,
  "results/tables/phase3_cluster_x_clinical_state_counts.csv"
)

write_csv(
  cluster_phase_proportions,
  "results/tables/phase3_cluster_x_clinical_state_proportions.csv"
)

write_csv(
  cluster_phase_composition,
  "results/tables/phase3_cluster_x_clinical_state_composition.csv"
)

cat(
  "\n✓ Cluster × clinical-state tables saved\n\n"
)


# ----------------------------------------------------------------------------
# 17.6 — Heatmap
# ----------------------------------------------------------------------------

cluster_phase_long <-
  cluster_phase_composition %>%
  mutate(
    Cluster = factor(
      seurat_clusters,
      levels = sort(
        unique(seurat_clusters)
      )
    )
  )

p_cluster_phase <-
  ggplot(
    cluster_phase_long,
    aes(
      x = Phase,
      y = Cluster,
      fill = Percent
    )
  ) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        Percent
      )
    ),
    size = 3
  ) +
  labs(
    title =
      "Global Atlas — Cluster Composition by Clinical State",
    x =
      "Clinical State",
    y =
      "Cluster",
    fill =
      "Percent"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_cluster_x_clinical_state_heatmap.png",
  p_cluster_phase,
  width = 9,
  height = 8,
  dpi = 300
)

cat(
  "✓ Cluster × clinical-state heatmap saved\n\n"
)

cat(
  "=== SECTION 17 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 18 — CLUSTER × DONOR DOMINANCE
# ============================================================================
#
# Purpose:
#   Determine whether clusters are disproportionately contributed by one
#   or a small number of donors.
#
# Outputs:
#   1. Cluster × donor counts
#   2. Cluster × donor composition
#   3. Dominant donor summary
#   4. Cluster × donor heatmap
#   5. Dominant donor plot
#
# ============================================================================

cat("=== SECTION 18: CLUSTER × DONOR DOMINANCE ===\n\n")


# ----------------------------------------------------------------------------
# 18.1 — Validate metadata
# ----------------------------------------------------------------------------

required_columns <- c(
  "seurat_clusters",
  "Donor"
)

missing_columns <- setdiff(
  required_columns,
  colnames(atlas_sketch@meta.data)
)

if (length(missing_columns) > 0) {
  
  stop(
    "ERROR: Required metadata columns are missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 18.2 — Cluster × donor counts
# ----------------------------------------------------------------------------

cluster_donor_counts <-
  atlas_sketch@meta.data %>%
  dplyr::count(
    seurat_clusters,
    Donor,
    name = "Cells"
  ) %>%
  dplyr::arrange(
    as.numeric(as.character(seurat_clusters)),
    desc(Cells)
  )


# ----------------------------------------------------------------------------
# 18.3 — Donor composition within each cluster
# ----------------------------------------------------------------------------

cluster_donor_composition <-
  cluster_donor_counts %>%
  group_by(seurat_clusters) %>%
  mutate(
    Cluster_Total =
      sum(Cells),
    Donor_Proportion =
      Cells / Cluster_Total,
    Donor_Percent =
      100 * Donor_Proportion
  ) %>%
  ungroup()


# ----------------------------------------------------------------------------
# 18.4 — Dominant donor
# ----------------------------------------------------------------------------

cluster_donor_dominance <-
  cluster_donor_composition %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::slice_max(
    order_by = Cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    seurat_clusters,
    Dominant_Donor = Donor,
    Dominant_Donor_Cells = Cells,
    Dominant_Donor_Percent = Donor_Percent
  )

# ----------------------------------------------------------------------------
# 18.5 — Cluster-level donor summary
# ----------------------------------------------------------------------------

cluster_donor_summary <-
  cluster_donor_composition %>%
  group_by(seurat_clusters) %>%
  summarise(
    Cluster_Cells =
      sum(Cells),
    Donors_Represented =
      n_distinct(Donor),
    .groups = "drop"
  ) %>%
  left_join(
    cluster_donor_dominance,
    by = "seurat_clusters"
  ) %>%
  arrange(
    as.numeric(
      as.character(seurat_clusters)
    )
  )

cat(
  "Dominant donor by cluster:\n\n"
)

print(
  cluster_donor_summary,
  n = Inf
)


# ----------------------------------------------------------------------------
# 18.6 — Save tables
# ----------------------------------------------------------------------------

write_csv(
  cluster_donor_counts,
  "results/tables/phase3_cluster_x_donor_counts.csv"
)

write_csv(
  cluster_donor_composition,
  "results/tables/phase3_cluster_x_donor_composition.csv"
)

write_csv(
  cluster_donor_summary,
  "results/tables/phase3_cluster_donor_dominance.csv"
)

cat(
  "\n✓ Cluster × donor tables saved\n\n"
)


# ----------------------------------------------------------------------------
# 18.7 — Cluster × donor heatmap
# ----------------------------------------------------------------------------

p_cluster_donor <-
  ggplot(
    cluster_donor_composition,
    aes(
      x = Donor,
      y = factor(
        seurat_clusters,
        levels = sort(
          unique(seurat_clusters)
        )
      ),
      fill = Donor_Percent
    )
  ) +
  geom_tile(
    color = "white"
  ) +
  labs(
    title =
      "Global Atlas — Cluster Composition by Donor",
    x =
      "Donor",
    y =
      "Cluster",
    fill =
      "Percent"
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        size = 7
      ),
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_cluster_x_donor_heatmap.png",
  p_cluster_donor,
  width = 14,
  height = 8,
  dpi = 300
)

cat(
  "✓ Cluster × donor heatmap saved\n\n"
)


# ----------------------------------------------------------------------------
# 18.8 — Dominant donor plot
# ----------------------------------------------------------------------------

p_donor_dominance <-
  ggplot(
    cluster_donor_summary,
    aes(
      x = factor(
        seurat_clusters,
        levels = sort(
          unique(seurat_clusters)
        )
      ),
      y = Dominant_Donor_Percent
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = Dominant_Donor
    ),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title =
      "Global Atlas — Dominant Donor Contribution by Cluster",
    x =
      "Cluster",
    y =
      "Cells from Dominant Donor (%)"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_cluster_donor_dominance.png",
  p_donor_dominance,
  width = 11,
  height = 7,
  dpi = 300
)

cat(
  "✓ Dominant-donor plot saved\n\n"
)

cat(
  "=== SECTION 18 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 19 — CLUSTER × SAMPLE DOMINANCE
# ============================================================================
#
# Purpose:
#   Evaluate whether individual GSM samples disproportionately contribute
#   particular global atlas clusters.
#
# ============================================================================

cat("=== SECTION 19: CLUSTER × SAMPLE DOMINANCE ===\n\n")


# ----------------------------------------------------------------------------
# 19.1 — Validate metadata
# ----------------------------------------------------------------------------

required_columns <- c(
  "seurat_clusters",
  "GSM"
)

missing_columns <- setdiff(
  required_columns,
  colnames(atlas_sketch@meta.data)
)

if (length(missing_columns) > 0) {
  
  stop(
    "ERROR: Required metadata columns are missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 19.2 — Cluster × sample counts
# ----------------------------------------------------------------------------

cluster_sample_counts <-
  atlas_sketch@meta.data %>%
  dplyr::count(
    seurat_clusters,
    GSM,
    name = "Cells"
  ) %>%
  dplyr::arrange(
    as.numeric(as.character(seurat_clusters)),
    desc(Cells)
  )


# ----------------------------------------------------------------------------
# 19.3 — Sample composition within each cluster
# ----------------------------------------------------------------------------

cluster_sample_composition <-
  cluster_sample_counts %>%
  group_by(seurat_clusters) %>%
  mutate(
    Cluster_Total =
      sum(Cells),
    Sample_Proportion =
      Cells / Cluster_Total,
    Sample_Percent =
      100 * Sample_Proportion
  ) %>%
  ungroup()


# ----------------------------------------------------------------------------
# 19.4 — Dominant sample
# ----------------------------------------------------------------------------

cluster_sample_dominance <-
  cluster_sample_composition %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::slice_max(
    order_by = Cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    seurat_clusters,
    Dominant_Sample = GSM,
    Dominant_Sample_Cells = Cells,
    Dominant_Sample_Percent = Sample_Percent
  )

# ----------------------------------------------------------------------------
# 19.5 — Cluster-level sample summary
# ----------------------------------------------------------------------------

cluster_sample_summary <-
  cluster_sample_composition %>%
  group_by(seurat_clusters) %>%
  summarise(
    Cluster_Cells =
      sum(Cells),
    Samples_Represented =
      n_distinct(GSM),
    .groups = "drop"
  ) %>%
  left_join(
    cluster_sample_dominance,
    by = "seurat_clusters"
  ) %>%
  arrange(
    as.numeric(
      as.character(seurat_clusters)
    )
  )

cat(
  "Dominant sample by cluster:\n\n"
)

print(
  cluster_sample_summary,
  n = Inf
)


# ----------------------------------------------------------------------------
# 19.6 — Save tables
# ----------------------------------------------------------------------------

write_csv(
  cluster_sample_counts,
  "results/tables/phase3_cluster_x_sample_counts.csv"
)

write_csv(
  cluster_sample_composition,
  "results/tables/phase3_cluster_x_sample_composition.csv"
)

write_csv(
  cluster_sample_summary,
  "results/tables/phase3_cluster_sample_dominance.csv"
)

cat(
  "\n✓ Cluster × sample tables saved\n\n"
)


# ----------------------------------------------------------------------------
# 19.7 — Cluster × sample heatmap
# ----------------------------------------------------------------------------

p_cluster_sample <-
  ggplot(
    cluster_sample_composition,
    aes(
      x = GSM,
      y = factor(
        seurat_clusters,
        levels = sort(
          unique(seurat_clusters)
        )
      ),
      fill = Sample_Percent
    )
  ) +
  geom_tile(
    color = "white"
  ) +
  labs(
    title =
      "Global Atlas — Cluster Composition by Sample",
    x =
      "Sample",
    y =
      "Cluster",
    fill =
      "Percent"
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        size = 7
      ),
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_cluster_x_sample_heatmap.png",
  p_cluster_sample,
  width = 14,
  height = 8,
  dpi = 300
)

cat(
  "✓ Cluster × sample heatmap saved\n\n"
)


# ----------------------------------------------------------------------------
# 19.8 — Dominant sample plot
# ----------------------------------------------------------------------------

p_sample_dominance <-
  ggplot(
    cluster_sample_summary,
    aes(
      x = factor(
        seurat_clusters,
        levels = sort(
          unique(seurat_clusters)
        )
      ),
      y = Dominant_Sample_Percent
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = Dominant_Sample
    ),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title =
      "Global Atlas — Dominant Sample Contribution by Cluster",
    x =
      "Cluster",
    y =
      "Cells from Dominant Sample (%)"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(
        hjust = 0.5
      )
  )

ggsave(
  "results/figures/phase3_cluster_sample_dominance.png",
  p_sample_dominance,
  width = 11,
  height = 7,
  dpi = 300
)

cat(
  "✓ Dominant-sample plot saved\n\n"
)

cat(
  "=== SECTION 19 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 20 — CLUSTER MARKER DISCOVERY
# ============================================================================

cat("=== SECTION 20: CLUSTER MARKER DISCOVERY ===\n\n")


# ----------------------------------------------------------------------------
# 20.1 — Activate RNA assay
# ----------------------------------------------------------------------------

DefaultAssay(atlas_sketch) <- "RNA"


# ----------------------------------------------------------------------------
# 20.2 — Join RNA layers
# ----------------------------------------------------------------------------

cat(
  "Joining RNA assay layers...\n"
)

atlas_sketch <- JoinLayers(
  object = atlas_sketch,
  assay = "RNA"
)

cat(
  "✓ RNA assay layers joined\n\n"
)


# ----------------------------------------------------------------------------
# 20.3 — Find cluster markers
# ----------------------------------------------------------------------------

cat(
  "Running FindAllMarkers...\n\n"
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

cat(
  "\n✓ FindAllMarkers completed\n\n"
)

# ----------------------------------------------------------------------------
# 20.4 — Validate marker table
# ----------------------------------------------------------------------------

if (nrow(cluster_markers) == 0) {
  
  stop(
    "ERROR: No cluster markers were identified."
  )
}

required_marker_columns <- c(
  "gene",
  "cluster",
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(cluster_markers)
)

if (length(missing_marker_columns) > 0) {
  
  stop(
    "ERROR: Marker table is missing required columns: ",
    paste(
      missing_marker_columns,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 20.5 — Confirm all 17 clusters have markers
# ----------------------------------------------------------------------------

expected_clusters <- sort(
  unique(
    as.character(
      atlas_sketch$seurat_clusters
    )
  )
)

cluster_markers$cluster <- as.character(
  cluster_markers$cluster
)

marker_clusters <- sort(
  unique(
    cluster_markers$cluster
  )
)

missing_clusters <- setdiff(
  expected_clusters,
  marker_clusters
)

if (length(missing_clusters) > 0) {
  
  stop(
    "ERROR: No markers identified for cluster(s): ",
    paste(
      missing_clusters,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 20.6 — Add marker rank
# ----------------------------------------------------------------------------

cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  mutate(
    Marker_Rank = row_number()
  ) %>%
  ungroup()


# ----------------------------------------------------------------------------
# 20.7 — Save complete marker table
# ----------------------------------------------------------------------------

write_csv(
  cluster_markers,
  "results/tables/phase3_cluster_markers_complete.csv"
)

cat(
  "✓ Complete marker table saved\n",
  "  Rows: ",
  nrow(cluster_markers),
  "\n",
  sep = ""
)


# ----------------------------------------------------------------------------
# 20.8 — Top 10 markers
# ----------------------------------------------------------------------------

top10_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()


top10_counts <- top10_cluster_markers %>%
  dplyr::count(
    cluster,
    name = "Marker_Count"
  )

if (
  nrow(top10_counts) != 17
) {
  
  stop(
    "ERROR: Top-10 marker table does not contain all 17 clusters."
  )
}


# ----------------------------------------------------------------------------
# 20.9 — Save top 10 table
# ----------------------------------------------------------------------------

write_csv(
  top10_cluster_markers,
  "results/tables/phase3_cluster_markers_top10.csv"
)


# ----------------------------------------------------------------------------
# 20.10 — Compact top 10 summary
# ----------------------------------------------------------------------------

top10_marker_summary <- top10_cluster_markers %>%
  group_by(cluster) %>%
  summarise(
    Top10_Markers =
      paste(
        gene,
        collapse = ", "
      ),
    .groups = "drop"
  ) %>%
  arrange(
    as.numeric(
      as.character(cluster)
    )
  )

write_csv(
  top10_marker_summary,
  "results/tables/phase3_cluster_markers_top10_summary.csv"
)

cat(
  "✓ Top-10 marker tables saved\n\n"
)

print(
  top10_marker_summary,
  n = Inf
)

cat("\n")


# ----------------------------------------------------------------------------
# 20.11 — Save marker checkpoint
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_markers.rds"
)

cat(
  "✓ Marker checkpoint saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_markers.rds\n\n"
)

cat(
  "=== SECTION 20 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 21 — CANONICAL LINEAGE MARKERS
# ============================================================================

cat("=== SECTION 21: CANONICAL LINEAGE MARKERS ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_markers.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"


# ----------------------------------------------------------------------------
# 21.1 — Define canonical marker programs
# ----------------------------------------------------------------------------

canonical_markers <- list(
  
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
    "MALAT1"
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
    "TRBC1",
    "XCL1",
    "XCL2"
  ),
  
  MAIT = c(
    "TRAV1-2",
    "SLC4A10",
    "KLRB1",
    "IL7R",
    "CCL20",
    "IL23R"
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
    "CD14"
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
    "GZMB",
    "CLEC4C",
    "IRF7",
    "TCF4",
    "LILRA4",
    "PTCRA"
  )
)


# ----------------------------------------------------------------------------
# 21.2 — Determine marker availability
# ----------------------------------------------------------------------------

canonical_markers_present <- lapply(
  canonical_markers,
  function(x) {
    unique(
      intersect(
        x,
        rownames(atlas_sketch)
      )
    )
  }
)

cat(
  "Canonical marker availability:\n\n"
)

for (cell_type in names(canonical_markers_present)) {
  
  cat(
    sprintf(
      "%-20s %2d / %2d markers present\n",
      cell_type,
      length(
        canonical_markers_present[[cell_type]]
      ),
      length(
        canonical_markers[[cell_type]]
      )
    )
  )
}

cat("\n")


# ----------------------------------------------------------------------------
# 21.3 — Validate marker programs
# ----------------------------------------------------------------------------

program_sizes <- lengths(
  canonical_markers_present
)

if (
  any(
    program_sizes < 3
  )
) {
  
  bad_programs <- names(
    program_sizes[
      program_sizes < 3
    ]
  )
  
  stop(
    "ERROR: Insufficient canonical markers in: ",
    paste(
      bad_programs,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------------
# 21.4 — Save canonical marker table
# ----------------------------------------------------------------------------

canonical_marker_table <- do.call(
  rbind,
  lapply(
    names(canonical_markers_present),
    function(cell_type) {
      
      data.frame(
        CellType =
          cell_type,
        Marker =
          canonical_markers_present[[cell_type]],
        stringsAsFactors = FALSE
      )
      
    }
  )
)

write_csv(
  canonical_marker_table,
  "results/tables/phase3_canonical_lineage_markers.csv"
)

saveRDS(
  canonical_markers_present,
  "results/rds_objects/phase3_canonical_lineage_markers.rds"
)

cat(
  "✓ Canonical marker definitions saved\n\n"
)

cat(
  "=== SECTION 21 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 22 — CANONICAL MARKER HEATMAP
# ============================================================================
#
# Cluster-level average expression is used rather than plotting all 20,000
# cells. This makes the canonical lineage programs interpretable across the
# 17 atlas clusters.
#
# ============================================================================

cat("=== SECTION 22: CANONICAL MARKER HEATMAP ===\n\n")

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_markers.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"

canonical_markers_present <- readRDS(
  "results/rds_objects/phase3_canonical_lineage_markers.rds"
)

expected_clusters <- sort(
  unique(
    as.character(
      atlas_sketch$seurat_clusters
    )
  )
)

cluster_ids <- as.character(
  Idents(atlas_sketch)
)


# ----------------------------------------------------------------------------
# 22.1 — Confirm clusters
# ----------------------------------------------------------------------------

if (
  !all(
    expected_clusters %in%
    unique(cluster_ids)
  )
) {
  
  stop(
    "ERROR: Clusters 0–17 were not all found."
  )
}


# ----------------------------------------------------------------------------
# 22.2 — Prepare marker list
# ----------------------------------------------------------------------------

canonical_marker_genes <- unique(
  unlist(
    canonical_markers_present
  )
)

canonical_marker_genes <- intersect(
  canonical_marker_genes,
  rownames(atlas_sketch)
)

if (
  length(canonical_marker_genes) == 0
) {
  
  stop(
    "ERROR: No canonical marker genes are present."
  )
}


# ----------------------------------------------------------------------------
# 22.3 — Extract normalized expression
# ----------------------------------------------------------------------------

expression_matrix <- SeuratObject::LayerData(
  object = atlas_sketch,
  assay = "RNA",
  layer = "data"
)

expression_matrix <- expression_matrix[
  canonical_marker_genes,
  ,
  drop = FALSE
]


# ----------------------------------------------------------------------------
# 22.4 — Calculate cluster-level average expression
# ----------------------------------------------------------------------------

cluster_average <- matrix(
  
  NA_real_,
  
  nrow =
    length(canonical_marker_genes),
  
  ncol =
    length(expected_clusters),
  
  dimnames = list(
    canonical_marker_genes,
    expected_clusters
  )
)

for (cluster in expected_clusters) {
  
  cells_in_cluster <- which(
    cluster_ids == cluster
  )
  
  if (
    length(cells_in_cluster) == 0
  ) {
    
    stop(
      "ERROR: Cluster ",
      cluster,
      " contains no cells."
    )
  }
  
  cluster_average[, cluster] <-
    Matrix::rowMeans(
      expression_matrix[
        ,
        cells_in_cluster,
        drop = FALSE
      ]
    )
}


# ----------------------------------------------------------------------------
# 22.5 — Order markers by lineage program
# ----------------------------------------------------------------------------

ordered_marker_genes <- unique(
  unlist(
    lapply(
      names(canonical_markers_present),
      function(cell_type) {
        
        intersect(
          canonical_markers_present[[cell_type]],
          rownames(cluster_average)
        )
        
      }
    ),
    use.names = FALSE
  )
)

cluster_average <- cluster_average[
  ordered_marker_genes,
  expected_clusters,
  drop = FALSE
]


# ----------------------------------------------------------------------------
# 22.6 — Row-wise Z-score
# ----------------------------------------------------------------------------

heatmap_matrix <- t(
  apply(
    cluster_average,
    1,
    function(x) {
      
      if (
        sd(x) == 0
      ) {
        
        return(
          rep(
            0,
            length(x)
          )
        )
      }
      
      as.numeric(
        scale(x)
      )
    }
  )
)

rownames(heatmap_matrix) <-
  rownames(cluster_average)

colnames(heatmap_matrix) <-
  colnames(cluster_average)


# ----------------------------------------------------------------------------
# 22.7 — Prepare plotting data
# ----------------------------------------------------------------------------

heatmap_df <- expand.grid(
  Marker =
    rownames(heatmap_matrix),
  Cluster =
    colnames(heatmap_matrix),
  KEEP.OUT.ATTRS =
    FALSE,
  stringsAsFactors =
    FALSE
)

heatmap_df$Expression <-
  as.vector(
    heatmap_matrix
  )

heatmap_df$Marker <- factor(
  heatmap_df$Marker,
  levels =
    rev(
      rownames(heatmap_matrix)
    )
)

heatmap_df$Cluster <- factor(
  heatmap_df$Cluster,
  levels =
    expected_clusters
)


# ----------------------------------------------------------------------------
# 22.8 — Generate heatmap
# ----------------------------------------------------------------------------

heatmap_plot <- ggplot(
  heatmap_df,
  aes(
    x = Cluster,
    y = Marker,
    fill = Expression
  )
) +
  geom_tile(
    width = 0.95,
    height = 0.95
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    title =
      "Canonical Lineage Marker Expression",
    subtitle =
      "Cluster-level average expression",
    x =
      "Cluster",
    y =
      "Canonical Marker"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid =
      element_blank(),
    axis.text.x =
      element_text(
        size = 10
      ),
    axis.text.y =
      element_text(
        size = 8
      ),
    axis.title =
      element_text(
        face = "bold"
      ),
    plot.title =
      element_text(
        face = "bold",
        size = 14
      )
  )


# ----------------------------------------------------------------------------
# 22.9 — Save heatmap
# ----------------------------------------------------------------------------

output_png <-
  "results/figures/phase3_canonical_marker_heatmap.png"

ggsave(
  filename =
    output_png,
  plot =
    heatmap_plot,
  width =
    10,
  height =
    12,
  units =
    "in",
  dpi =
    300,
  limitsize =
    FALSE
)

saveRDS(
  heatmap_plot,
  "results/rds_objects/phase3_canonical_marker_heatmap.rds"
)

cat(
  "✓ Canonical marker heatmap saved\n\n"
)


# ----------------------------------------------------------------------------
# 22.10 — Validate heatmap
# ----------------------------------------------------------------------------

if (
  !file.exists(output_png)
) {
  
  stop(
    "ERROR: Canonical marker heatmap was not created."
  )
}

if (
  ncol(heatmap_matrix) != 17
) {
  
  stop(
    "ERROR: Heatmap does not contain 17 clusters."
  )
}

if (
  any(
    !is.finite(
      heatmap_matrix
    )
  )
) {
  
  stop(
    "ERROR: Heatmap contains non-finite values."
  )
}

cat(
  "✓ Canonical marker heatmap validated\n\n"
)

cat(
  "=== SECTION 22 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 23 — MULTI-METHOD CELL-TYPE ANNOTATION
# ============================================================================
#
# Evidence streams:
#   1. Canonical marker programs
#   2. Module scores
#   3. Top-50 marker/program overlap
#   4. SingleR / Monaco reference
#
# Final biological labels are manually adjudicated in Section 24.
#
# ============================================================================

cat("=== SECTION 23: MULTI-METHOD CELL-TYPE ANNOTATION ===\n\n")


# ----------------------------------------------------------------------------
# 23.1 — Load SingleR and celldex
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleR)
  library(celldex)
})

cat(
  "✓ SingleR loaded\n"
)

cat(
  "✓ celldex loaded\n\n"
)


# ----------------------------------------------------------------------------
# 23.2 — Load Monaco immune reference
# ----------------------------------------------------------------------------

monaco_ref <- celldex::MonacoImmuneData(
  cell.ont = "all"
)

if (
  !inherits(
    monaco_ref,
    "SummarizedExperiment"
  )
) {
  
  stop(
    "ERROR: MonacoImmuneData did not return a SummarizedExperiment."
  )
}

cat(
  "✓ MonacoImmuneData reference loaded\n",
  "  Reference genes: ",
  nrow(monaco_ref),
  "\n",
  "  Reference samples: ",
  ncol(monaco_ref),
  "\n\n",
  sep = ""
)


# ----------------------------------------------------------------------------
# 23.3 — Prepare marker programs
# ----------------------------------------------------------------------------

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
    "IL2RA",
    "CTLA4",
    "TNFRSF18",
    "TNFRSF4",
    "CCR8"
  ),
  
  NK = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "FCGR3A",
    "XCL1",
    "XCL2",
    "FGFBP2"
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
    "FCN1",
    "CD14",
    "CTSS",
    "CTSD",
    "LGALS3",
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
    "GZMB",
    "CLEC4C",
    "IRF7",
    "TCF4",
    "LILRA4",
    "PTCRA"
  )
)


# ----------------------------------------------------------------------------
# 23.4 — Retain genes present in atlas
# ----------------------------------------------------------------------------

marker_programs_present <- lapply(
  marker_programs,
  function(x) {
    unique(
      intersect(
        x,
        rownames(atlas_sketch)
      )
    )
  }
)

program_sizes <- lengths(
  marker_programs_present
)

if (
  any(
    program_sizes < 3
  )
) {
  
  bad_programs <- names(
    program_sizes[
      program_sizes < 3
    ]
  )
  
  stop(
    "ERROR: Insufficient marker genes in: ",
    paste(
      bad_programs,
      collapse = ", "
    )
  )
}


# Save marker programs.

annotation_marker_table <- do.call(
  rbind,
  lapply(
    names(marker_programs_present),
    function(program_name) {
      
      data.frame(
        Program =
          program_name,
        Gene =
          marker_programs_present[[program_name]],
        stringsAsFactors =
          FALSE
      )
      
    }
  )
)

write_csv(
  annotation_marker_table,
  "results/tables/phase3_annotation_marker_programs.csv"
)

cat(
  "✓ Annotation marker programs prepared\n\n"
)


# ============================================================================
# 23.5 — MODULE SCORE CALCULATION
# ============================================================================

cat(
  "Calculating module scores...\n\n"
)

DefaultAssay(atlas_sketch) <- "RNA"

set.seed(1234)

atlas_sketch <- AddModuleScore(
  
  object =
    atlas_sketch,
  
  features =
    marker_programs_present,
  
  assay =
    "RNA",
  
  name =
    "AnnotationModule",
  
  ctrl =
    50,
  
  seed =
    1234,
  
  search =
    FALSE,
  
  slot =
    "data"
)


# ----------------------------------------------------------------------------
# Rename numbered AddModuleScore columns
# ----------------------------------------------------------------------------

module_score_columns <- paste0(
  "AnnotationModule",
  seq_along(
    marker_programs_present
  )
)

if (
  !all(
    module_score_columns %in%
    colnames(
      atlas_sketch[[]]
    )
  )
) {
  
  stop(
    "ERROR: Expected AddModuleScore columns were not created."
  )
}

for (
  i in seq_along(
    marker_programs_present
  )
) {
  
  old_name <-
    module_score_columns[i]
  
  new_name <- paste0(
    names(
      marker_programs_present
    )[i],
    "_Score"
  )
  
  atlas_sketch[[new_name]] <-
    atlas_sketch[[
      old_name
    ]]
}


# Remove temporary numbered columns.

atlas_sketch[[
  module_score_columns
]] <- NULL


module_score_columns_named <- paste0(
  names(
    marker_programs_present
  ),
  "_Score"
)

cat(
  "✓ Module scores calculated\n\n"
)


# ============================================================================
# 23.6 — VERIFY MODULE SCORES
# ============================================================================

module_score_matrix <- atlas_sketch[[]][
  ,
  module_score_columns_named,
  drop = FALSE
]

if (
  ncol(module_score_matrix) !=
  length(
    marker_programs_present
  )
) {
  
  stop(
    "ERROR: Incorrect number of module-score columns."
  )
}

if (
  nrow(module_score_matrix) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: Module-score row count does not match atlas cells."
  )
}

if (
  any(
    !is.finite(
      as.matrix(
        module_score_matrix
      )
    )
  )
) {
  
  stop(
    "ERROR: Module scores contain non-finite values."
  )
}

write_csv(
  tibble(
    Cell = rownames(
      module_score_matrix
    )
  ) %>%
    bind_cols(
      as.data.frame(
        module_score_matrix
      )
    ),
  "results/tables/phase3_annotation_module_scores_cell_level.csv"
)

cat(
  "✓ Module-score matrix verified\n",
  "  Cells: ",
  nrow(module_score_matrix),
  "\n",
  "  Programs: ",
  ncol(module_score_matrix),
  "\n\n",
  sep = ""
)


# ============================================================================
# 23.7 — MODULE SCORES BY CLUSTER
# ============================================================================

expected_clusters <- sort(
  unique(
    as.character(
      atlas_sketch$seurat_clusters
    )
  )
)

cluster_ids <- as.character(
  Idents(atlas_sketch)
)

cluster_score_matrix <- matrix(
  
  NA_real_,
  
  nrow =
    length(
      module_score_columns_named
    ),
  
  ncol =
    length(
      expected_clusters
    ),
  
  dimnames = list(
    module_score_columns_named,
    expected_clusters
  )
)

for (
  cluster in expected_clusters
) {
  
  cells_in_cluster <- which(
    cluster_ids == cluster
  )
  
  if (
    length(cells_in_cluster) == 0
  ) {
    
    stop(
      "ERROR: Cluster ",
      cluster,
      " contains no cells."
    )
  }
  
  cluster_score_matrix[
    ,
    cluster
  ] <-
    colMeans(
      module_score_matrix[
        cells_in_cluster,
        ,
        drop = FALSE
      ]
    )
}

write.csv(
  cluster_score_matrix,
  "results/tables/phase3_annotation_module_scores_by_cluster.csv",
  row.names = TRUE
)

cat(
  "✓ Cluster-level module scores calculated and saved\n\n"
)


# ============================================================================
# 23.8 — MODULE-BASED ANNOTATION
# ============================================================================

module_based_annotation <- data.frame(
  Cluster =
    expected_clusters,
  stringsAsFactors =
    FALSE
)

best_program <- character(
  length(
    expected_clusters
  )
)

best_score <- numeric(
  length(
    expected_clusters
  )
)

second_score <- numeric(
  length(
    expected_clusters
  )
)

score_gap <- numeric(
  length(
    expected_clusters
  )
)


for (
  i in seq_along(
    expected_clusters
  )
) {
  
  scores <- cluster_score_matrix[
    ,
    expected_clusters[i]
  ]
  
  scores_sorted <- sort(
    scores,
    decreasing = TRUE
  )
  
  best_program[i] <-
    sub(
      "_Score$",
      "",
      names(
        scores_sorted
      )[1]
    )
  
  best_score[i] <-
    scores_sorted[1]
  
  second_score[i] <-
    scores_sorted[2]
  
  score_gap[i] <-
    scores_sorted[1] -
    scores_sorted[2]
}


module_based_annotation$Module_Label <-
  best_program

module_based_annotation$Top_Score <-
  best_score

module_based_annotation$Second_Score <-
  second_score

module_based_annotation$Score_Gap <-
  score_gap

module_based_annotation$Confidence <-
  ifelse(
    module_based_annotation$Score_Gap >= 0.20,
    "High",
    ifelse(
      module_based_annotation$Score_Gap >= 0.10,
      "Moderate",
      "Low"
    )
  )


write_csv(
  module_based_annotation,
  "results/tables/phase3_module_based_annotation.csv"
)

cat(
  "Module-based annotation:\n\n"
)

print(
  module_based_annotation,
  row.names = FALSE
)

cat("\n")


# ============================================================================
# 23.9 — MODULE SCORE HEATMAP
# ============================================================================

module_heatmap_matrix <- t(
  apply(
    cluster_score_matrix,
    1,
    function(x) {
      
      if (
        sd(x) == 0
      ) {
        
        return(
          rep(
            0,
            length(x)
          )
        )
      }
      
      as.numeric(
        scale(x)
      )
    }
  )
)

rownames(module_heatmap_matrix) <-
  rownames(cluster_score_matrix)

colnames(module_heatmap_matrix) <-
  colnames(cluster_score_matrix)


module_heatmap_df <- expand.grid(
  Program =
    rownames(
      module_heatmap_matrix
    ),
  Cluster =
    colnames(
      module_heatmap_matrix
    ),
  KEEP.OUT.ATTRS =
    FALSE,
  stringsAsFactors =
    FALSE
)

module_heatmap_df$Score <-
  as.vector(
    module_heatmap_matrix
  )

module_heatmap_df$Program <- factor(
  module_heatmap_df$Program,
  levels =
    rev(
      rownames(
        module_heatmap_matrix
      )
    )
)

module_heatmap_df$Cluster <- factor(
  module_heatmap_df$Cluster,
  levels =
    expected_clusters
)


module_heatmap_plot <- ggplot(
  module_heatmap_df,
  aes(
    x = Cluster,
    y = Program,
    fill = Score
  )
) +
  geom_tile(
    width = 0.95,
    height = 0.95
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Module\nZ-score"
  ) +
  labs(
    title =
      "Canonical Lineage Module Scores",
    subtitle =
      "Cluster-level average module scores",
    x =
      "Cluster",
    y =
      "Cell-type Program"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid =
      element_blank(),
    axis.text.x =
      element_text(
        size = 10
      ),
    axis.text.y =
      element_text(
        size = 9
      ),
    axis.title =
      element_text(
        face = "bold"
      ),
    plot.title =
      element_text(
        face = "bold",
        size = 14
      )
  )

ggsave(
  "results/figures/phase3_module_score_heatmap.png",
  module_heatmap_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

cat(
  "✓ Module-score heatmap saved\n\n"
)


# ============================================================================
# 23.10 — TOP 50 MARKERS PER CLUSTER
# ============================================================================

cat(
  "Extracting top 50 markers per cluster...\n\n"
)

cluster_markers <- read_csv(
  "results/tables/phase3_cluster_markers_complete.csv",
  show_col_types = FALSE
)

cluster_markers$cluster <- as.character(
  cluster_markers$cluster
)

if (
  !"Marker_Rank" %in%
  colnames(cluster_markers)
) {
  
  stop(
    "ERROR: Marker_Rank is missing from complete marker table."
  )
}

top50_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    Marker_Rank,
    .by_group = TRUE
  ) %>%
  slice_head(n = 50) %>%
  ungroup()

if (
  nrow(top50_markers) != 850
) {
  
  stop(
    "ERROR: Expected exactly 850 Top-50 marker rows."
  )
}

write_csv(
  top50_markers,
  "results/tables/phase3_cluster_markers_top50.csv"
)

cat(
  "✓ Top-50 marker table saved\n",
  "  Rows: 850\n\n"
)


# ============================================================================
# 23.11 — MARKER / PROGRAM OVERLAP
# ============================================================================

cat(
  "Calculating marker/program overlap...\n\n"
)

marker_programs_present <- readRDS(
  "results/rds_objects/phase3_canonical_lineage_markers.rds"
)

program_names <- names(
  marker_programs_present
)

overlap_counts <- matrix(
  
  0L,
  
  nrow =
    length(expected_clusters),
  
  ncol =
    length(program_names),
  
  dimnames = list(
    expected_clusters,
    program_names
  )
)


for (
  cluster in expected_clusters
) {
  
  cluster_genes <- unique(
    top50_markers$gene[
      top50_markers$cluster == cluster
    ]
  )
  
  for (
    program_name in program_names
  ) {
    
    overlap_counts[
      cluster,
      program_name
    ] <-
      length(
        intersect(
          cluster_genes,
          marker_programs_present[
            program_name
          ]
        )
      )
  }
}


marker_program_overlap <- as.data.frame(
  overlap_counts,
  check.names = FALSE
)

marker_program_overlap$Cluster <-
  rownames(
    marker_program_overlap
  )

marker_program_overlap <-
  marker_program_overlap[
    ,
    c(
      "Cluster",
      program_names
    ),
    drop = FALSE
  ]

overlap_matrix <- as.matrix(
  marker_program_overlap[
    ,
    program_names,
    drop = FALSE
  ]
)

marker_program_overlap$Best_Overlap_Program <-
  apply(
    overlap_matrix,
    1,
    function(x) {
      
      if (
        max(x) == 0
      ) {
        
        return(
          NA_character_
        )
      }
      
      names(x)[
        which.max(x)
      ]
    }
  )

marker_program_overlap$Best_Overlap_Count <-
  apply(
    overlap_matrix,
    1,
    max
  )

write_csv(
  marker_program_overlap,
  "results/tables/phase3_marker_program_overlap.csv"
)

cat(
  "✓ Marker/program overlap saved\n\n"
)


# ============================================================================
# 23.12 — SingleR CLUSTER-LEVEL ANNOTATION
# ============================================================================

cat(
  "Running SingleR at cluster level...\n\n"
)


# ----------------------------------------------------------------------------
# Extract normalized expression
# ----------------------------------------------------------------------------

test_expression <- SeuratObject::LayerData(
  object =
    atlas_sketch,
  assay =
    "RNA",
  layer =
    "data"
)

if (
  ncol(test_expression) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: SingleR test matrix does not match atlas cell count."
  )
}


# ----------------------------------------------------------------------------
# Cluster identities aligned to expression columns
# ----------------------------------------------------------------------------

singleR_clusters <- as.character(
  Idents(atlas_sketch)
)

names(singleR_clusters) <-
  colnames(
    test_expression
  )

if (
  length(singleR_clusters) !=
  ncol(test_expression)
) {
  
  stop(
    "ERROR: SingleR cluster vector length mismatch."
  )
}


# ----------------------------------------------------------------------------
# Verify Monaco labels
# ----------------------------------------------------------------------------

if (
  !"label.fine" %in%
  colnames(
    SummarizedExperiment::colData(
      monaco_ref
    )
  )
) {
  
  stop(
    "ERROR: Monaco reference does not contain label.fine."
  )
}

monaco_labels <- as.character(
  monaco_ref$label.fine
)

if (
  length(monaco_labels) !=
  ncol(monaco_ref)
) {
  
  stop(
    "ERROR: Monaco labels do not match reference samples."
  )
}


# ----------------------------------------------------------------------------
# Run SingleR
# ----------------------------------------------------------------------------

set.seed(1234)

singleR_results <- SingleR::SingleR(
  
  test =
    test_expression,
  
  ref =
    monaco_ref,
  
  labels =
    monaco_labels,
  
  clusters =
    singleR_clusters,
  
  assay.type.test =
    "logcounts",
  
  assay.type.ref =
    "logcounts",
  
  prune =
    FALSE
)


if (
  nrow(singleR_results) != 17
) {
  
  stop(
    "ERROR: SingleR did not return exactly 17 cluster results."
  )
}

cat(
  "✓ SingleR cluster-level annotation completed\n\n"
)


# ----------------------------------------------------------------------------
# Extract SingleR results
# ----------------------------------------------------------------------------

singleR_cluster_labels <- as.character(
  singleR_results$labels
)

singleR_delta <- as.numeric(
  singleR_results$delta.next
)

singleR_best_score <- as.numeric(
  singleR_results$scores[
    cbind(
      seq_len(
        nrow(
          singleR_results
        )
      ),
      match(
        singleR_results$labels,
        colnames(
          singleR_results$scores
        )
      )
    )
  ]
)


singleR_annotation <- data.frame(
  
  Cluster =
    rownames(
      singleR_results
    ),
  
  SingleR_Label =
    singleR_cluster_labels,
  
  SingleR_Best_Score =
    singleR_best_score,
  
  SingleR_Delta_Next =
    singleR_delta,
  
  stringsAsFactors =
    FALSE
)

singleR_annotation$Cluster <-
  as.character(
    singleR_annotation$Cluster
  )

singleR_annotation <-
  singleR_annotation[
    match(
      expected_clusters,
      singleR_annotation$Cluster
    ),
    ,
    drop = FALSE
  ]

if (
  any(
    is.na(
      singleR_annotation$SingleR_Label
    )
  )
) {
  
  stop(
    "ERROR: SingleR produced missing labels."
  )
}

write_csv(
  singleR_annotation,
  "results/tables/phase3_singleR_cluster_annotation.csv"
)

cat(
  "✓ SingleR annotation table saved\n\n"
)

print(
  singleR_annotation,
  row.names = FALSE
)

cat("\n")


# ============================================================================
# 23.13 — CONSOLIDATE ALL ANNOTATION EVIDENCE
# ============================================================================

cat(
  "Consolidating annotation evidence...\n\n"
)

module_based_annotation <- read_csv(
  "results/tables/phase3_module_based_annotation.csv",
  show_col_types = FALSE
)

module_based_annotation$Cluster <-
  as.character(
    module_based_annotation$Cluster
  )

marker_overlap_summary <-
  marker_program_overlap[
    ,
    c(
      "Cluster",
      "Best_Overlap_Program",
      "Best_Overlap_Count"
    ),
    drop = FALSE
  ]


annotation_evidence <- data.frame(
  Cluster =
    expected_clusters,
  stringsAsFactors =
    FALSE
)

annotation_evidence <-
  annotation_evidence %>%
  left_join(
    module_based_annotation,
    by = "Cluster"
  ) %>%
  left_join(
    marker_overlap_summary,
    by = "Cluster"
  ) %>%
  left_join(
    singleR_annotation,
    by = "Cluster"
  )

annotation_evidence <-
  annotation_evidence[
    match(
      expected_clusters,
      annotation_evidence$Cluster
    ),
    ,
    drop = FALSE
  ]

if (
  nrow(annotation_evidence) != 17
) {
  
  stop(
    "ERROR: Consolidated annotation table does not contain 17 clusters."
  )
}

write_csv(
  annotation_evidence,
  "results/tables/phase3_annotation_evidence_consolidated.csv"
)

cat(
  "✓ Consolidated annotation evidence saved\n\n"
)

print(
  annotation_evidence,
  row.names = FALSE
)

cat("\n")


# ============================================================================
# 23.14 — FINAL MANUAL ANNOTATION TEMPLATE
# ============================================================================

cat(
  "Creating final manual annotation template...\n\n"
)

final_annotation_template <- data.frame(
  
  Cluster =
    expected_clusters,
  
  Module_Label =
    annotation_evidence$Module_Label,
  
  Module_Score =
    annotation_evidence$Top_Score,
  
  Module_Score_Gap =
    annotation_evidence$Score_Gap,
  
  Marker_Overlap_Label =
    annotation_evidence$Best_Overlap_Program,
  
  Marker_Overlap_Count =
    annotation_evidence$Best_Overlap_Count,
  
  SingleR_Label =
    annotation_evidence$SingleR_Label,
  
  SingleR_Best_Score =
    annotation_evidence$SingleR_Best_Score,
  
  SingleR_Delta_Next =
    annotation_evidence$SingleR_Delta_Next,
  
  Final_CellType =
    NA_character_,
  
  Confidence =
    NA_character_,
  
  Rationale =
    NA_character_,
  
  stringsAsFactors =
    FALSE
)


write_csv(
  final_annotation_template,
  "results/tables/phase3_final_manual_annotation.csv"
)

cat(
  "✓ Final manual annotation template saved\n\n"
)


# ============================================================================
# 23.15 — ANNOTATION VALIDATION
# ============================================================================

cat(
  "Running annotation validation...\n\n"
)


if (
  ncol(atlas_sketch) != 20000
) {
  
  stop(
    "ERROR: Atlas does not contain 20,000 cells."
  )
}

if (
  length(
    unique(
      cluster_ids
    )
  ) != 17
) {
  
  stop(
    "ERROR: Atlas does not contain 17 clusters."
  )
}

if (
  nrow(cluster_markers) == 0
) {
  
  stop(
    "ERROR: Marker table is empty."
  )
}

if (
  nrow(top50_markers) != 850
) {
  
  stop(
    "ERROR: Top-50 marker table does not contain 850 rows."
  )
}

if (
  ncol(module_score_matrix) != 13
) {
  
  stop(
    "ERROR: Expected 13 module-score programs."
  )
}

if (
  nrow(singleR_annotation) != 17
) {
  
  stop(
    "ERROR: SingleR annotation does not contain 17 clusters."
  )
}

if (
  nrow(annotation_evidence) != 17
) {
  
  stop(
    "ERROR: Consolidated annotation does not contain 17 clusters."
  )
}

if (
  nrow(final_annotation_template) != 17
) {
  
  stop(
    "ERROR: Final annotation template does not contain 17 clusters."
  )
}

cat(
  "✓ Annotation validation passed\n\n"
)


# ============================================================================
# SAVE PRE-FINAL ANNOTATION OBJECT
# ============================================================================

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_pre_final_annotation.rds"
)

cat(
  "✓ Pre-final annotation object saved:\n",
  "  results/rds_objects/phase3_atlas_sketch_pre_final_annotation.rds\n\n"
)


cat(
  "=== SECTION 23 COMPLETE ===\n\n"
)


# ============================================================================
# SECTION 24 — FINAL CELL-TYPE ANNOTATION OF THE SKETCH
# ============================================================================

cat("=== SECTION 24: FINAL CELL-TYPE ANNOTATION ===\n\n")


# ----------------------------------------------------------------------------
# 24.1 — Reload pre-final checkpoint
# ----------------------------------------------------------------------------

atlas_sketch <- readRDS(
  "results/rds_objects/phase3_atlas_sketch_pre_final_annotation.rds"
)

DefaultAssay(atlas_sketch) <- "RNA"


# ----------------------------------------------------------------------------
# 24.2 — Confirm cluster structure
# ----------------------------------------------------------------------------

expected_clusters <- sort(
  unique(
    as.character(
      atlas_sketch$seurat_clusters
    )
  )
)

actual_clusters <- unique(
  as.character(
    Idents(atlas_sketch)
  )
)

if (
  !setequal(
    actual_clusters,
    expected_clusters
  )
) {
  
  stop(
    "ERROR: Atlas cluster identities are not exactly 0–17."
  )
}


# ----------------------------------------------------------------------------
# 24.3 — Load consolidated evidence
# ----------------------------------------------------------------------------

annotation_evidence <- read_csv(
  "results/tables/phase3_annotation_evidence_consolidated.csv",
  show_col_types = FALSE
)

annotation_evidence$Cluster <-
  as.character(
    annotation_evidence$Cluster
  )

annotation_evidence <-
  annotation_evidence[
    match(
      expected_clusters,
      annotation_evidence$Cluster
    ),
    ,
    drop = FALSE
  ]

if (
  nrow(annotation_evidence) != 17
) {
  
  stop(
    "ERROR: Annotation evidence does not contain 17 clusters."
  )
}

# ----------------------------------------------------------------------------
# 24.4 — FINAL MANUAL CLUSTER LABELS
# ----------------------------------------------------------------------------

final_cluster_labels <- c(
  
  "0"  = "CD8_T",
  "1"  = "GammaDelta_T",
  "2"  = "CD4_T",
  "3"  = "CD8_T",
  "4"  = "NK",
  "5"  = "MAIT",
  "6"  = "CD8_T",
  "7"  = "CD8_T",
  "8"  = "CD4_T",
  "9"  = "CD8_T",
  "10" = "MAIT",
  "11" = "Inflammatory_Myeloid",
  "12" = "B_Cell",
  "13" = "NK",
  "14" = "NK",
  "15" = "Plasma_Cell",
  "16" = "pDC"
)


final_cluster_confidence <- c(
  
  "0"  = "High",
  "1"  = "Moderate",
  "2"  = "High",
  "3"  = "High",
  "4"  = "High",
  "5"  = "High",
  "6"  = "Moderate",
  "7"  = "Moderate",
  "8"  = "Low",
  "9"  = "Low",
  "10" = "High",
  "11" = "High",
  "12" = "High",
  "13" = "High",
  "14" = "High",
  "15" = "High",
  "16" = "High"
)


final_cluster_rationale <- c(
  
  "0"  = "CD8_T module and effector-memory CD8 T-cell SingleR annotation agree.",
  "1"  = "SingleR supports non-Vd2 gamma-delta T cells, while the module score favors NK; retained as GammaDelta_T with moderate confidence.",
  "2"  = "CD4_T module and Th1/Th17 SingleR annotation support CD4 T-cell identity.",
  "3"  = "CD8_T module and effector-memory CD8 T-cell SingleR annotation agree.",
  "4"  = "NK module and natural-killer-cell SingleR annotation agree.",
  "5"  = "MAIT module and MAIT-cell SingleR annotation agree.",
  "6"  = "CD8_T module and effector-memory CD8 T-cell SingleR annotation agree, with a modest module-score gap.",
  "7"  = "CD8_T module and effector-memory CD8 T-cell SingleR annotation agree, with moderate module-score separation.",
  "8"  = "CD4_T module is weakly favored, while SingleR favors MAIT; retained as CD4_T with low confidence.",
  "9"  = "CD8_T module and effector-memory CD8 T-cell SingleR annotation agree, but the SingleR separation is weak.",
  "10" = "MAIT module and MAIT-cell SingleR annotation agree.",
  "11" = "Monocyte_Myeloid module and classical-monocyte SingleR annotation agree.",
  "12" = "B-cell module and non-switched memory B-cell SingleR annotation agree.",
  "13" = "NK module and natural-killer-cell SingleR annotation agree.",
  "14" = "NK module and natural-killer-cell SingleR annotation agree.",
  "15" = "Plasma-cell module and plasmablast SingleR annotation agree.",
  "16" = "pDC module and plasmacytoid dendritic-cell SingleR annotation agree."
)


# ----------------------------------------------------------------------------
# 24.5 — Validate final cluster labels
# ----------------------------------------------------------------------------

expected_clusters <- sort(
  unique(
    as.character(
      Idents(atlas_sketch)
    )
  ),
  method = "radix"
)

if (
  length(final_cluster_labels) != length(expected_clusters)
) {
  
  stop(
    "ERROR: Number of final annotations does not match number of atlas clusters."
  )
}

if (
  !setequal(
    names(final_cluster_labels),
    expected_clusters
  )
) {
  
  stop(
    "ERROR: Final annotation clusters do not match atlas clusters."
  )
}

if (
  any(
    is.na(final_cluster_labels)
  )
) {
  
  stop(
    "ERROR: Missing final cluster labels."
  )
}


# ----------------------------------------------------------------------------
# 24.6 — Construct final annotation table
# ----------------------------------------------------------------------------

cluster_counts <- table(
  factor(
    as.character(
      Idents(atlas_sketch)
    ),
    levels = expected_clusters
  )
)

final_annotation <- data.frame(
  
  Cluster = expected_clusters,
  
  Final_CellType =
    unname(
      final_cluster_labels[
        expected_clusters
      ]
    ),
  
  Cell_Count =
    as.integer(
      cluster_counts[
        expected_clusters
      ]
    ),
  
  stringsAsFactors = FALSE
)

final_annotation$Cell_Percent <-
  round(
    100 *
      final_annotation$Cell_Count /
      ncol(atlas_sketch),
    2
  )


# ----------------------------------------------------------------------------
# 24.7 — Add supporting evidence
# ----------------------------------------------------------------------------

evidence_by_cluster <- annotation_evidence

rownames(evidence_by_cluster) <-
  as.character(
    evidence_by_cluster$Cluster
  )

final_annotation$Module_Label <-
  as.character(
    evidence_by_cluster[
      final_annotation$Cluster,
      "Module_Label"
    ][[1]]
  )

final_annotation$Module_Score_Gap <-
  as.numeric(
    evidence_by_cluster[
      final_annotation$Cluster,
      "Score_Gap"
    ][[1]]
  )

final_annotation$Marker_Overlap_Label <-
  as.character(
    evidence_by_cluster[
      final_annotation$Cluster,
      "Best_Overlap_Program"
    ][[1]]
  )

final_annotation$Marker_Overlap_Count <-
  as.numeric(
    evidence_by_cluster[
      final_annotation$Cluster,
      "Best_Overlap_Count"
    ][[1]]
  )

final_annotation$SingleR_Label <-
  as.character(
    evidence_by_cluster[
      final_annotation$Cluster,
      "SingleR_Label"
    ][[1]]
  )

final_annotation$SingleR_Delta_Next <-
  as.numeric(
    evidence_by_cluster[
      final_annotation$Cluster,
      "SingleR_Delta_Next"
    ][[1]]
  )

final_annotation$Confidence <-
  as.character(
    final_cluster_confidence[
      final_annotation$Cluster
    ]
  )

final_annotation$Rationale <-
  as.character(
    final_cluster_rationale[
      final_annotation$Cluster
    ]
  )

# ----------------------------------------------------------------------------
# 24.8 — Propagate final cluster labels to cells
# ----------------------------------------------------------------------------

cell_cluster_ids <- as.character(
  Idents(atlas_sketch)
)

atlas_sketch$Final_CellType <-
  unname(
    final_cluster_labels[
      cell_cluster_ids
    ]
  )

atlas_sketch$Annotation_Confidence <-
  unname(
    final_cluster_confidence[
      cell_cluster_ids
    ]
  )


# ----------------------------------------------------------------------------
# 24.9 — Validate cell-level annotation
# ----------------------------------------------------------------------------

if (
  length(
    atlas_sketch$Final_CellType
  ) !=
  ncol(atlas_sketch)
) {
  
  stop(
    "ERROR: Final cell-type annotation does not match cell count."
  )
}

if (
  any(
    is.na(
      atlas_sketch$Final_CellType
    )
  )
) {
  
  stop(
    "ERROR: Some cells have missing final cell-type labels."
  )
}

cluster_label_check <- tapply(
  atlas_sketch$Final_CellType,
  cell_cluster_ids,
  function(x) {
    length(
      unique(x)
    )
  }
)

if (
  any(
    cluster_label_check != 1
  )
) {
  
  stop(
    "ERROR: At least one cluster contains multiple final labels."
  )
}


# ----------------------------------------------------------------------------
# 24.10 — Display final annotation
# ----------------------------------------------------------------------------

cat(
  "FINAL CLUSTER ANNOTATION:\n\n"
)

print(
  final_annotation,
  row.names = FALSE
)

cat("\n")

cat(
  "FINAL CELL-TYPE COUNTS:\n\n"
)

print(
  sort(
    table(
      atlas_sketch$Final_CellType
    ),
    decreasing = TRUE
  )
)

cat("\n")


# ----------------------------------------------------------------------------
# 24.11 — Save final annotation table
# ----------------------------------------------------------------------------

write_csv(
  final_annotation,
  "results/tables/phase3_final_cell_type_annotation.csv"
)

cat(
  "✓ Final cell-type annotation table saved\n"
)


# ----------------------------------------------------------------------------
# 24.12 — Save cluster-to-cell-type mapping
# ----------------------------------------------------------------------------

cluster_annotation_mapping <- data.frame(
  
  Cluster = expected_clusters,
  
  Final_CellType =
    unname(
      final_cluster_labels[
        expected_clusters
      ]
    ),
  
  Confidence =
    unname(
      final_cluster_confidence[
        expected_clusters
      ]
    ),
  
  stringsAsFactors = FALSE
)

write_csv(
  cluster_annotation_mapping,
  "results/tables/phase3_cluster_to_celltype_mapping.csv"
)

cat(
  "✓ Cluster-to-cell-type mapping saved\n"
)


# ----------------------------------------------------------------------------
# 24.13 — Save final annotated atlas
# ----------------------------------------------------------------------------

saveRDS(
  atlas_sketch,
  "results/rds_objects/phase3_atlas_sketch_final_annotation.rds"
)

cat(
  "✓ Final annotated atlas saved\n\n"
)


# ============================================================================
# FINAL VALIDATION — SECTIONS 17–24
# ============================================================================

if (
  ncol(atlas_sketch) != 20000
) {
  
  stop(
    "FINAL ERROR: Atlas does not contain exactly 20,000 cells."
  )
}

if (
  length(
    unique(
      as.character(
        Idents(atlas_sketch)
      )
    )
  ) != 17
) {
  
  stop(
    "FINAL ERROR: Atlas does not contain exactly 17 clusters."
  )
}

if (
  nrow(final_annotation) != 17
) {
  
  stop(
    "FINAL ERROR: Final annotation table does not contain 17 clusters."
  )
}

if (
  any(
    is.na(
      atlas_sketch$Final_CellType
    )
  )
) {
  
  stop(
    "FINAL ERROR: Missing cell-type annotations detected."
  )
}

cat(
  "============================================================\n"
)

cat(
  "✓ PHASE 3 SCRIPT 02 — SECTIONS 17–24 COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Atlas sketch cells:",
  ncol(atlas_sketch),
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
  "Final cell types:",
  length(
    unique(
      atlas_sketch$Final_CellType
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
  "\n"
)

cat(
  "Samples:",
  length(
    unique(
      atlas_sketch$GSM
    )
  ),
  "\n\n"
)

cat(
  "✓ Final manual annotation complete\n",
  "✓ Cell-level labels propagated\n",
  "✓ Annotation confidence recorded\n",
  "✓ Final annotated atlas checkpoint saved\n\n"
)

cat(
  "Final checkpoint:\n",
  "  results/rds_objects/phase3_atlas_sketch_final_annotation.rds\n"
)

cat(
  "============================================================\n"
)

