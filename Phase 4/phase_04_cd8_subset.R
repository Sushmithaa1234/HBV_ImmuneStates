###############################################################################
# HBV scRNA-seq PROJECT
# PHASE 4 — CD8 T-CELL TRANSCRIPTIONAL STATE ANALYSIS
###############################################################################
#
# Dataset
# -------
# GEO: GSE182159
# Tissue: Liver
# Donors: 23
# Clinical states: NL, IT, IA, AR, AC
#
# Biological objective
# --------------------
# Characterize transcriptional heterogeneity within the intrahepatic CD8
# T-cell compartment across the five GEO-defined clinical states.
#
# Phase 4 scope
# -------------
# Phase 4 is a cell-state discovery and biological adjudication phase.
# It characterizes transcriptional heterogeneity within the Phase 3
# CD8-labelled compartment.
#
# Clinical-state inference is intentionally NOT performed here.
# Donor-aware differential analysis is reserved for Phase 6.
#
# Scientific principles
# ---------------------
# 1.  Donor = biological replicate.
# 2.  Cells are measurement units for state discovery and annotation.
# 3.  Clinical-state inference is NOT performed in Phase 4.
# 4.  Phase 6 performs donor-aware differential analysis.
# 5.  GEO matrices contain processed log-CP10K expression, not raw counts.
# 6.  Expression is never reverse-transformed into pseudo-counts.
# 7.  NormalizeData() is not used.
# 8.  DESeq2/edgeR are not applied to reconstructed counts.
# 9.  No integration, Harmony, or batch correction is performed.
# 10. No pathway enrichment is performed.
# 11. No CellChat analysis is performed.
# 12. Biological interpretations must be supported by observed markers and/or
#     transcriptional programs.
# 13. Ambiguous populations remain explicitly ambiguous.
# 14. Clinical-state composition is descriptive only.
# 15. No causal claims are made.
# 16. New Seurat metadata are added with AddMetaData() using cell-indexed
#     vectors; @meta.data is never replaced wholesale.
#
# Input
# -----
# results/rds_objects/phase3_final_full_dataset.rds
#
# Primary output
# --------------
# results/rds_objects/phase4_final_cd8_analysis.rds
#
# Final adjudicated output
# ------------------------
# results/rds_objects/phase4_final_cd8_analysis_adjudicated.rds
#
# Additional outputs
# ------------------
# results/rds_objects/phase4_cd8_source_layers_preserved.rds
# results/tables/phase4_cd8_*.csv
# results/figures/phase4_cd8_*.png
#
###############################################################################


############################
# 1. SETUP
############################

suppressPackageStartupMessages({
  
  library(Seurat)
  library(SeuratObject)
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  
  library(Matrix)
  
  library(ggplot2)
})

set.seed(20260916)

project_root <- getwd()

rds_dir <- file.path(
  project_root,
  "results",
  "rds_objects"
)

table_dir <- file.path(
  project_root,
  "results",
  "tables"
)

figure_dir <- file.path(
  project_root,
  "results",
  "figures"
)

dir.create(
  rds_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  rds_dir,
  "phase3_final_full_dataset.rds"
)

source_layer_output <- file.path(
  rds_dir,
  "phase4_cd8_source_layers_preserved.rds"
)

final_output <- file.path(
  rds_dir,
  "phase4_final_cd8_analysis.rds"
)

message("============================================================")
message("PHASE 4: CD8 T-CELL TRANSCRIPTIONAL STATE ANALYSIS")
message("============================================================")
message("Input:  ", input_file)
message("Output: ", final_output)


###############################################################################
# 2. LOAD PHASE 3 OBJECT
###############################################################################

if (!file.exists(input_file)) {
  
  stop(
    "Phase 3 object not found:\n",
    input_file
  )
}

seurat_obj <- readRDS(
  input_file
)

if (!inherits(seurat_obj, "Seurat")) {
  
  stop(
    "Input object is not a Seurat object."
  )
}

DefaultAssay(seurat_obj) <- "RNA"

message("\nPhase 3 object loaded.")

message(
  "Cells: ",
  ncol(seurat_obj)
)

message(
  "Features: ",
  nrow(seurat_obj)
)


###############################################################################
# 3. VALIDATE PHASE 3 OBJECT
###############################################################################

expected_total_cells <- 106592L
expected_total_features <- 24452L

if (ncol(seurat_obj) != expected_total_cells) {
  
  stop(
    "Unexpected Phase 3 cell count.\n",
    "Expected: ",
    expected_total_cells,
    "\nFound: ",
    ncol(seurat_obj)
  )
}

if (nrow(seurat_obj) != expected_total_features) {
  
  stop(
    "Unexpected Phase 3 feature count.\n",
    "Expected: ",
    expected_total_features,
    "\nFound: ",
    nrow(seurat_obj)
  )
}


###############################################################################
# 4. REQUIRED PHASE 3 METADATA
###############################################################################

required_metadata <- c(
  
  "GSM",
  "Donor",
  "Phase",
  
  "Genes_Detected",
  "Mito_Genes_Detected",
  "Mito_Detection_Fraction",
  
  "Ribo_Genes_Detected",
  "Ribo_Detection_Fraction",
  
  "Total_LogCP10K",
  
  "QC_LowGenes",
  "QC_Pass",
  "Potential_HighComplexity",
  
  "Transferred_Atlas_Cluster",
  "Final_Cell_Type",
  
  "Transfer_Agreement",
  "Transfer_SecondBest_Agreement",
  "Transfer_Agreement_Margin",
  
  "Transfer_Confidence",
  "Low_Confidence"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(seurat_obj@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "Missing required Phase 3 metadata:\n",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}


###############################################################################
# 5. VALIDATE DONORS, SAMPLES AND CLINICAL STATES
###############################################################################

expected_phases <- c(
  "NL",
  "IT",
  "IA",
  "AR",
  "AC"
)

observed_phases <- sort(
  unique(
    as.character(
      seurat_obj$Phase
    )
  )
)

if (!setequal(
  observed_phases,
  expected_phases
)) {
  
  stop(
    "Unexpected clinical-state labels.\n",
    "Observed: ",
    paste(
      observed_phases,
      collapse = ", "
    ),
    "\nExpected: ",
    paste(
      expected_phases,
      collapse = ", "
    )
  )
}

n_donors <- length(
  unique(
    seurat_obj$Donor
  )
)

n_gsm <- length(
  unique(
    seurat_obj$GSM
  )
)

if (n_donors != 23) {
  
  stop(
    "Expected 23 donors; found ",
    n_donors
  )
}

if (n_gsm != 23) {
  
  stop(
    "Expected 23 GSM/sample identifiers; found ",
    n_gsm
  )
}

message(
  "\nValidated 23 donors and 23 samples."
)

message(
  "Clinical states: ",
  paste(
    expected_phases,
    collapse = ", "
  )
)


###############################################################################
# 6. PHASE 3 CELL-TYPE SUMMARY
###############################################################################

phase3_label_summary <- seurat_obj@meta.data %>%
  
  dplyr::count(
    Final_Cell_Type,
    name = "Cells"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(Cells)
  )

write.csv(
  phase3_label_summary,
  file.path(
    table_dir,
    "phase4_phase3_Final_Cell_Type_summary.csv"
  ),
  row.names = FALSE
)

message("\nPhase 3 cell-type summary:")
print(phase3_label_summary)


###############################################################################
# 7. EXTRACT CD8 T CELLS
###############################################################################

cd8_cells <- rownames(
  seurat_obj@meta.data
)[
  seurat_obj$Final_Cell_Type == "CD8_T"
]

if (length(cd8_cells) == 0) {
  
  stop(
    "No cells labeled CD8_T were found in Final_Cell_Type."
  )
}

cd8_obj <- subset(
  seurat_obj,
  cells = cd8_cells
)

DefaultAssay(cd8_obj) <- "RNA"

expected_cd8_n <- length(
  cd8_cells
)

message(
  "\nCD8 subset extracted."
)

message(
  "CD8 cells: ",
  expected_cd8_n
)

message(
  "CD8 features: ",
  nrow(cd8_obj)
)


###############################################################################
# 8. VALIDATE CD8 CELL COUNT
###############################################################################

#
# Phase 3 finalized CD8 count should be 41,281.
#

if (expected_cd8_n != 41281L) {
  
  stop(
    "Unexpected CD8 cell count.\n",
    "Expected: 41281\n",
    "Found: ",
    expected_cd8_n
  )
}

if (!identical(
  colnames(cd8_obj),
  rownames(cd8_obj@meta.data)
)) {
  
  stop(
    "Initial CD8 object has inconsistent cell names between the object ",
    "and metadata."
  )
}


###############################################################################
# 9. SAVE CD8 SOURCE-LAYER VERSION
###############################################################################

source_layers <- Layers(
  cd8_obj[["RNA"]]
)

message(
  "\nRNA layers before joining:"
)

print(source_layers)

write.csv(
  data.frame(
    Layer = source_layers
  ),
  file.path(
    table_dir,
    "phase4_cd8_source_layers.csv"
  ),
  row.names = FALSE
)

saveRDS(
  cd8_obj,
  source_layer_output
)

message(
  "Preserved source-layer CD8 object:\n",
  source_layer_output
)


###############################################################################
# 10. CD8 DONOR / SAMPLE / CLINICAL-STATE SUMMARY
###############################################################################

cd8_donor_state_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    GSM,
    Phase,
    name = "CD8_Cells"
  ) %>%
  
  dplyr::arrange(
    Phase,
    Donor
  )

write.csv(
  cd8_donor_state_counts,
  file.path(
    table_dir,
    "phase4_cd8_donor_state_counts.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 11. CD8 CELLS BY CLINICAL STATE
###############################################################################

cd8_phase_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    Phase,
    name = "CD8_Cells"
  ) %>%
  
  dplyr::mutate(
    Percent_of_All_CD8 =
      100 * CD8_Cells / sum(CD8_Cells)
  )

write.csv(
  cd8_phase_counts,
  file.path(
    table_dir,
    "phase4_cd8_phase_counts.csv"
  ),
  row.names = FALSE
)

message(
  "\nCD8 cells by clinical state:"
)

print(cd8_phase_counts)


###############################################################################
# 12. CD8 QC CHARACTERIZATION
###############################################################################

available_qc <- intersect(
  c(
    "Genes_Detected",
    "Mito_Genes_Detected",
    "Mito_Detection_Fraction",
    "Ribo_Genes_Detected",
    "Ribo_Detection_Fraction",
    "Total_LogCP10K",
    "QC_LowGenes",
    "QC_Pass",
    "Potential_HighComplexity"
  ),
  colnames(cd8_obj@meta.data)
)

message(
  "\nAvailable Phase 3 QC metadata:"
)

print(available_qc)


###############################################################################
# 13. GENES-DETECTED SUMMARY
###############################################################################

if ("Genes_Detected" %in% available_qc) {
  
  cd8_genes_detected_summary <- cd8_obj@meta.data %>%
    
    dplyr::summarise(
      
      Cells = dplyr::n(),
      
      Median_Genes_Detected =
        median(
          Genes_Detected,
          na.rm = TRUE
        ),
      
      Mean_Genes_Detected =
        mean(
          Genes_Detected,
          na.rm = TRUE
        ),
      
      Min_Genes_Detected =
        min(
          Genes_Detected,
          na.rm = TRUE
        ),
      
      Max_Genes_Detected =
        max(
          Genes_Detected,
          na.rm = TRUE
        )
    )
  
  write.csv(
    cd8_genes_detected_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_Genes_Detected_summary.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 14. MITOCHONDRIAL SUMMARY
###############################################################################

if ("Mito_Detection_Fraction" %in% available_qc) {
  
  cd8_mito_summary <- cd8_obj@meta.data %>%
    
    dplyr::summarise(
      
      Cells = dplyr::n(),
      
      Median_Mito_Detection_Fraction =
        median(
          Mito_Detection_Fraction,
          na.rm = TRUE
        ),
      
      Mean_Mito_Detection_Fraction =
        mean(
          Mito_Detection_Fraction,
          na.rm = TRUE
        ),
      
      Max_Mito_Detection_Fraction =
        max(
          Mito_Detection_Fraction,
          na.rm = TRUE
        )
    )
  
  write.csv(
    cd8_mito_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_mito_summary.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 15. RIBOSOMAL METADATA SUMMARY
###############################################################################

if ("Ribo_Detection_Fraction" %in% available_qc) {
  
  cd8_ribo_summary <- cd8_obj@meta.data %>%
    
    dplyr::summarise(
      
      Cells = dplyr::n(),
      
      Median_Ribo_Detection_Fraction =
        median(
          Ribo_Detection_Fraction,
          na.rm = TRUE
        ),
      
      Mean_Ribo_Detection_Fraction =
        mean(
          Ribo_Detection_Fraction,
          na.rm = TRUE
        ),
      
      Max_Ribo_Detection_Fraction =
        max(
          Ribo_Detection_Fraction,
          na.rm = TRUE
        )
    )
  
  write.csv(
    cd8_ribo_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_ribo_summary.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 16. TOTAL LOGCP10K SUMMARY
###############################################################################

if ("Total_LogCP10K" %in% available_qc) {
  
  cd8_expression_sum_summary <- cd8_obj@meta.data %>%
    
    dplyr::summarise(
      
      Cells = dplyr::n(),
      
      Median_Total_LogCP10K =
        median(
          Total_LogCP10K,
          na.rm = TRUE
        ),
      
      Mean_Total_LogCP10K =
        mean(
          Total_LogCP10K,
          na.rm = TRUE
        ),
      
      Min_Total_LogCP10K =
        min(
          Total_LogCP10K,
          na.rm = TRUE
        ),
      
      Max_Total_LogCP10K =
        max(
          Total_LogCP10K,
          na.rm = TRUE
        )
    )
  
  write.csv(
    cd8_expression_sum_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_Total_LogCP10K_summary.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 17. QC / DOUBLE-CHECK SUMMARIES
###############################################################################

if ("QC_Pass" %in% colnames(cd8_obj@meta.data)) {
  
  qc_pass_summary <- as.data.frame(
    table(
      cd8_obj$QC_Pass,
      useNA = "ifany"
    )
  )
  
  colnames(qc_pass_summary) <- c(
    "QC_Pass",
    "Cells"
  )
  
  write.csv(
    qc_pass_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_pass_summary.csv"
    ),
    row.names = FALSE
  )
}

if ("QC_LowGenes" %in% colnames(cd8_obj@meta.data)) {
  
  low_gene_summary <- as.data.frame(
    table(
      cd8_obj$QC_LowGenes,
      useNA = "ifany"
    )
  )
  
  colnames(low_gene_summary) <- c(
    "QC_LowGenes",
    "Cells"
  )
  
  write.csv(
    low_gene_summary,
    file.path(
      table_dir,
      "phase4_cd8_qc_low_gene_summary.csv"
    ),
    row.names = FALSE
  )
}

if ("Potential_HighComplexity" %in%
    colnames(cd8_obj@meta.data)) {
  
  complexity_summary <- as.data.frame(
    table(
      cd8_obj$Potential_HighComplexity,
      useNA = "ifany"
    )
  )
  
  colnames(complexity_summary) <- c(
    "Potential_HighComplexity",
    "Cells"
  )
  
  write.csv(
    complexity_summary,
    file.path(
      table_dir,
      "phase4_cd8_high_complexity_summary.csv"
    ),
    row.names = FALSE
  )
}

message(
  "\nNo automatic CD8 cell removal based on Phase 3 QC flags."
)


###############################################################################
# 18. JOIN EXISTING NORMALIZED DATA LAYERS
###############################################################################
#
# The Phase 3 RNA assay contains sample-specific normalized expression layers.
#
# These are already processed log-CP10K expression values.
#
# We join those existing normalized layers for state discovery.
#
# We DO NOT:
#   - create counts
#   - reverse-transform expression
#   - NormalizeData()
#
###############################################################################

message(
  "\nJoining existing normalized expression layers..."
)

cd8_obj[["RNA"]] <- JoinLayers(
  cd8_obj[["RNA"]]
)

joined_layers <- Layers(
  cd8_obj[["RNA"]]
)

message(
  "\nRNA layers after joining:"
)

print(joined_layers)

if (!("data" %in% joined_layers)) {
  
  stop(
    "Joined RNA assay does not contain a 'data' layer."
  )
}

joined_data <- LayerData(
  cd8_obj,
  assay = "RNA",
  layer = "data"
)

if (
  nrow(joined_data) != expected_total_features ||
  ncol(joined_data) != expected_cd8_n
) {
  
  stop(
    "Joined data layer has unexpected dimensions.\n",
    "Expected: ",
    expected_total_features,
    " x ",
    expected_cd8_n,
    "\nFound: ",
    nrow(joined_data),
    " x ",
    ncol(joined_data)
  )
}

if (!identical(
  colnames(joined_data),
  colnames(cd8_obj)
)) {
  
  stop(
    "Joined data-layer cell names do not match the Seurat object."
  )
}


###############################################################################
# 19. CD8-SPECIFIC VARIABLE FEATURES
###############################################################################

message(
  "\nSelecting CD8-specific variable features..."
)

cd8_obj <- FindVariableFeatures(
  cd8_obj,
  assay = "RNA",
  layer = "data",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

cd8_hvgs <- VariableFeatures(
  cd8_obj
)

if (length(cd8_hvgs) == 0) {
  
  stop(
    "No variable features were identified."
  )
}

write.csv(
  data.frame(
    Rank = seq_along(cd8_hvgs),
    Gene = cd8_hvgs
  ),
  file.path(
    table_dir,
    "phase4_cd8_hvg_list.csv"
  ),
  row.names = FALSE
)

message(
  "CD8 HVGs identified: ",
  length(cd8_hvgs)
)


###############################################################################
# 20. VARIABLE FEATURE PLOT — PNG
###############################################################################

p_hvg <- VariableFeaturePlot(
  cd8_obj
) +
  ggtitle(
    "CD8 variable features"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_variable_features.png"
  ),
  p_hvg,
  width = 8,
  height = 6,
  dpi = 300
)


###############################################################################
# 21. SCALE CD8 HVGs
###############################################################################

message(
  "\nScaling CD8 HVGs..."
)

cd8_obj <- ScaleData(
  cd8_obj,
  assay = "RNA",
  features = cd8_hvgs,
  verbose = TRUE
)


###############################################################################
# 22. PCA
###############################################################################

message(
  "\nRunning PCA..."
)

cd8_obj <- RunPCA(
  cd8_obj,
  assay = "RNA",
  features = cd8_hvgs,
  npcs = 50,
  verbose = TRUE
)

if (!("pca" %in% Reductions(cd8_obj))) {
  
  stop(
    "PCA reduction was not generated."
  )
}


###############################################################################
# 23. PCA ELBOW — PNG
###############################################################################

p_elbow <- ElbowPlot(
  cd8_obj,
  ndims = 50
) +
  ggtitle(
    "CD8 PCA elbow plot"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_pca_elbow.png"
  ),
  p_elbow,
  width = 8,
  height = 6,
  dpi = 300
)


###############################################################################
# 24. NEIGHBOR GRAPH
###############################################################################
#
# Primary state-discovery space:
#   PC1:20
#
# No integration or batch correction.
###############################################################################

message(
  "\nConstructing CD8 neighbor graph..."
)

cd8_obj <- FindNeighbors(
  cd8_obj,
  assay = "RNA",
  reduction = "pca",
  dims = 1:20,
  k.param = 20,
  graph.name = c("CD8_nn", "CD8_snn"),
  verbose = TRUE
)


###############################################################################
# 25. CLUSTERING
###############################################################################

message(
  "\nClustering CD8 cells..."
)

cd8_obj <- FindClusters(
  cd8_obj,
  graph.name = "CD8_snn",
  resolution = 0.4,
  algorithm = 1,
  random.seed = 20260916,
  verbose = TRUE
)

cluster_vector <- setNames(
  as.character(
    Idents(cd8_obj)
  ),
  names(
    Idents(cd8_obj)
  )
)

cd8_obj <- AddMetaData(
  cd8_obj,
  metadata = cluster_vector,
  col.name = "CD8_Cluster"
)

cluster_sizes <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    name = "Cells"
  ) %>%
  
  dplyr::arrange(
    suppressWarnings(
      as.numeric(CD8_Cluster)
    )
  )

write.csv(
  cluster_sizes,
  file.path(
    table_dir,
    "phase4_cd8_cluster_sizes.csv"
  ),
  row.names = FALSE
)

message(
  "\nCD8 cluster sizes:"
)

print(cluster_sizes)


###############################################################################
# 26. UMAP
###############################################################################

message(
  "\nRunning UMAP..."
)

cd8_obj <- RunUMAP(
  cd8_obj,
  reduction = "pca",
  dims = 1:20,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 20260916,
  reduction.name = "umap_cd8",
  reduction.key = "CD8UMAP_",
  verbose = TRUE
)


###############################################################################
# 27. UMAP — CLUSTERS — PNG
###############################################################################

p_cluster <- DimPlot(
  cd8_obj,
  reduction = "umap_cd8",
  group.by = "CD8_Cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "CD8 T-cell transcriptional clusters"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_umap_clusters.png"
  ),
  p_cluster,
  width = 8,
  height = 7,
  dpi = 300
)


###############################################################################
# 28. UMAP — CLINICAL STATES — PNG
###############################################################################

p_phase <- DimPlot(
  cd8_obj,
  reduction = "umap_cd8",
  group.by = "Phase"
) +
  ggtitle(
    "CD8 cells by GEO-defined clinical state"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_umap_phase.png"
  ),
  p_phase,
  width = 8,
  height = 7,
  dpi = 300
)


###############################################################################
# 29. UMAP — DONORS — PNG
###############################################################################

p_donor <- DimPlot(
  cd8_obj,
  reduction = "umap_cd8",
  group.by = "Donor"
) +
  ggtitle(
    "CD8 cells by donor"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_umap_donor.png"
  ),
  p_donor,
  width = 11,
  height = 8,
  dpi = 300
)


###############################################################################
# 30. CLUSTER × DONOR REPRESENTATION
###############################################################################

cluster_donor_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    Donor,
    Phase,
    name = "Cells"
  )

write.csv(
  cluster_donor_counts,
  file.path(
    table_dir,
    "phase4_cd8_cluster_donor_counts.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 31. CLUSTER REPRESENTATION DIAGNOSTICS
###############################################################################

cluster_representation <- cluster_donor_counts %>%
  
  dplyr::group_by(
    CD8_Cluster
  ) %>%
  
  dplyr::summarise(
    
    Cluster_Cells =
      sum(Cells),
    
    Donors_Represented =
      dplyr::n_distinct(Donor),
    
    Clinical_States_Represented =
      dplyr::n_distinct(Phase),
    
    Largest_Donor_Count =
      max(Cells),
    
    Largest_Donor_Fraction =
      max(Cells) / sum(Cells),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    suppressWarnings(
      as.numeric(CD8_Cluster)
    )
  )

write.csv(
  cluster_representation,
  file.path(
    table_dir,
    "phase4_cd8_cluster_representation_diagnostics.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 32. CLUSTER × CLINICAL-STATE COUNTS
###############################################################################

cluster_phase_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    Phase,
    name = "Cells"
  )

write.csv(
  cluster_phase_counts,
  file.path(
    table_dir,
    "phase4_cd8_cluster_phase_counts.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 33. CLUSTER × CLINICAL-STATE PROPORTIONS
###############################################################################

cluster_phase_proportions <- cluster_phase_counts %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::mutate(
    Proportion_Within_CD8 =
      Cells / sum(Cells)
  ) %>%
  
  dplyr::ungroup()

write.csv(
  cluster_phase_proportions,
  file.path(
    table_dir,
    "phase4_cd8_cluster_phase_proportions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 34. EXPLORATORY CLUSTER MARKERS
###############################################################################
#
# These characterize transcriptional differences among CD8 clusters.
#
# They are NOT:
#   - donor-level differential expression
#   - clinical-state inference
#   - Phase 6 analysis
###############################################################################

message(
  "\nFinding exploratory cluster markers..."
)

Idents(cd8_obj) <- "CD8_Cluster"

cd8_markers <- FindAllMarkers(
  cd8_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

if (nrow(cd8_markers) == 0) {
  
  stop(
    "No positive cluster markers were identified."
  )
}

write.csv(
  cd8_markers,
  file.path(
    table_dir,
    "phase4_cd8_cluster_markers_all.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 35. TOP 20 MARKERS PER CLUSTER
###############################################################################

top_markers <- cd8_markers %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  
  dplyr::slice_head(
    n = 20
  ) %>%
  
  dplyr::ungroup()

write.csv(
  top_markers,
  file.path(
    table_dir,
    "phase4_cd8_top20_markers_per_cluster.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 36. TOP MARKER HEATMAP — PNG
###############################################################################

heatmap_genes <- cd8_markers %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  
  dplyr::slice_head(
    n = 5
  ) %>%
  
  dplyr::pull(gene) %>%
  
  unique()

heatmap_genes <- intersect(
  heatmap_genes,
  rownames(cd8_obj)
)

if (length(heatmap_genes) > 0) {
  
  cd8_obj <- ScaleData(
    cd8_obj,
    assay = "RNA",
    features = heatmap_genes,
    verbose = FALSE
  )
  
  p_heatmap <- DoHeatmap(
    cd8_obj,
    features = heatmap_genes,
    group.by = "CD8_Cluster"
  ) +
    ggtitle(
      "Top exploratory markers by CD8 cluster"
    )
  
  ggsave(
    file.path(
      figure_dir,
      "phase4_cd8_top_markers_heatmap.png"
    ),
    p_heatmap,
    width = 12,
    height = 10,
    dpi = 300
  )
}


###############################################################################
# 37. BIOLOGICAL PROGRAM DEFINITIONS
###############################################################################
#
# HYBRID PHASE 4 APPROACH
#
# We distinguish three levels of evidence:
#
#   1. LINEAGE EVIDENCE
#      Does the cluster retain evidence of conventional CD8/T-cell identity,
#      or does it show evidence of another/noncanonical population?
#
#   2. TRANSCRIPTIONAL PROGRAMS
#      What biological programs are relatively prominent?
#
#   3. CLUSTER MARKERS
#      Which genes distinguish the cluster from the other CD8-labelled cells?
#
# IMPORTANT:
#   Program scores are descriptive.
#   They are NOT used by themselves to define biological states.
#
###############################################################################

programs <- list(
  
  Naive_Memory = c(
    "CCR7",
    "LTB",
    "IL7R",
    "MAL",
    "TCF7",
    "LEF1",
    "LTB"
  ),
  
  Cytotoxic = c(
    "NKG7",
    "CCL5",
    "GNLY",
    "GZMB",
    "GZMH",
    "GZMK",
    "PRF1",
    "CTSW",
    "FGFBP2",
    "CX3CR1"
  ),
  
  Activation = c(
    "CD69",
    "IL2RA",
    "HLA-DRA",
    "HLA-DRB1",
    "TNFRSF4",
    "TNFRSF9",
    "CD38"
  ),
  
  Dysfunction_Associated = c(
    "PDCD1",
    "TOX",
    "TIGIT",
    "CTLA4",
    "HAVCR2",
    "LAG3",
    "CXCR4",
    "RGS1",
    "LAYN",
    "CXCL13"
  ),
  
  Immediate_Early_Response = c(
    "FOS",
    "JUN",
    "JUNB",
    "DUSP1",
    "DUSP2",
    "DUSP4",
    "EGR1",
    "EGR2",
    "NR4A1"
  ),
  
  NK_Like = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "KLRF1",
    "KLRC1",
    "KLRC2",
    "TYROBP",
    "FCER1G",
    "FCGR3A"
  ),
  
  Nonconventional_T = c(
    "TRDC",
    "TRGC1",
    "TRGC2",
    "TRDV1",
    "TRDV2",
    "TRGV2",
    "TRGV4",
    "TRGV5",
    "TRGV8",
    "TRGV9",
    "KLRB1",
    "SLC4A10",
    "ZBTB16"
  ),
  
  Tissue_Resident_Associated = c(
    "CD69",
    "ITGA1",
    "CXCR6",
    "ZNF683",
    "ITGAE"
  ),
  
  Proliferation = c(
    "MKI67",
    "TOP2A",
    "STMN1",
    "TYMS",
    "PCNA"
  )
)


###############################################################################
# 38. LINEAGE / BIOLOGICAL ADJUDICATION PANELS
###############################################################################
#
# These panels are used as evidence streams rather than automatic labels.
#
# The purpose is to distinguish:
#
#   - conventional CD8-associated populations
#   - NK-like populations
#   - gamma-delta-like populations
#   - Treg-associated populations
#   - other unexpected / ambiguous populations
#
# Unexpected populations are retained rather than removed.
###############################################################################

lineage_panels <- list(
  
  Conventional_T_CD8 = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2",
    "CD8A",
    "CD8B"
  ),
  
  NK_Associated = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "KLRF1",
    "KLRC1",
    "KLRC2",
    "TYROBP",
    "FCER1G",
    "NCR3",
    "S1PR5"
  ),
  
  GammaDelta_T = c(
    "TRDC",
    "TRGC1",
    "TRGC2",
    "TRDV1",
    "TRDV2",
    "TRGV2",
    "TRGV4",
    "TRGV5",
    "TRGV8",
    "TRGV9"
  ),
  
  Treg_Associated = c(
    "CD4",
    "IL7R",
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "TNFRSF4",
    "TNFRSF18",
    "ICOS"
  ),
  
  B_Cell_Associated = c(
    "CD79A",
    "MS4A1",
    "CD37",
    "CD74",
    "HLA-DRA"
  ),
  
  Myeloid_Associated = c(
    "LYZ",
    "S100A8",
    "S100A9",
    "FCGR3A",
    "CTSD"
  )
)

lineage_panels_present <- lapply(
  lineage_panels,
  function(x) {
    intersect(
      unique(x),
      rownames(cd8_obj)
    )
  }
)

lineage_panels_present <- lineage_panels_present[
  lengths(lineage_panels_present) >= 2
]

message("\nLineage adjudication panels:")
print(lineage_panels_present)


###############################################################################
# 39. ADD BIOLOGICAL PROGRAM SCORES
###############################################################################

programs_present <- lapply(
  programs,
  function(x) {
    intersect(
      unique(x),
      rownames(cd8_obj)
    )
  }
)

programs_present <- programs_present[
  lengths(programs_present) >= 2
]

message("\nGenes available for program scoring:")
print(programs_present)


for (program_name in names(programs_present)) {
  
  genes <- programs_present[[program_name]]
  
  temporary_name <- paste0(
    "CD8Program_",
    program_name,
    "_"
  )
  
  cd8_obj <- AddModuleScore(
    cd8_obj,
    features = list(genes),
    assay = "RNA",
    name = temporary_name,
    seed = 20260916
  )
  
  generated_name <- paste0(
    temporary_name,
    "1"
  )
  
  final_name <- paste0(
    "Program_",
    program_name
  )
  
  if (
    generated_name %in%
    colnames(cd8_obj@meta.data)
  ) {
    
    program_vector <-
      cd8_obj@meta.data[[generated_name]]
    
    names(program_vector) <-
      rownames(cd8_obj@meta.data)
    
    cd8_obj <- AddMetaData(
      cd8_obj,
      metadata = program_vector,
      col.name = final_name
    )
    
    cd8_obj@meta.data[[generated_name]] <- NULL
  }
}


###############################################################################
# 40. PROGRAM SCORE COLUMNS
###############################################################################

program_columns <- grep(
  "^Program_",
  colnames(cd8_obj@meta.data),
  value = TRUE
)

message("\nProgram score columns:")
print(program_columns)


###############################################################################
# 41. PROGRAM SCORES BY CLUSTER
###############################################################################

if (length(program_columns) > 0) {
  
  program_cluster_summary <- cd8_obj@meta.data %>%
    
    dplyr::group_by(
      CD8_Cluster
    ) %>%
    
    dplyr::summarise(
      
      dplyr::across(
        dplyr::all_of(program_columns),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      Cells = dplyr::n(),
      
      .groups = "drop"
    ) %>%
    
    dplyr::arrange(
      suppressWarnings(
        as.numeric(CD8_Cluster)
      )
    )
  
  write.csv(
    program_cluster_summary,
    file.path(
      table_dir,
      "phase4_cd8_program_scores_by_cluster.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 42. PROGRAM SCORE UMAPS
###############################################################################

for (program_col in program_columns) {
  
  p <- FeaturePlot(
    cd8_obj,
    features = program_col,
    reduction = "umap_cd8"
  ) +
    ggtitle(
      program_col
    )
  
  safe_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    program_col
  )
  
  ggsave(
    file.path(
      figure_dir,
      paste0(
        "phase4_cd8_",
        safe_name,
        "_umap.png"
      )
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


###############################################################################
# 43. LINEAGE DOTPLOT
###############################################################################

lineage_genes <- unique(
  unlist(
    lineage_panels_present
  )
)

if (length(lineage_genes) > 0) {
  
  p_lineage <- DotPlot(
    cd8_obj,
    features = lineage_genes,
    group.by = "CD8_Cluster"
  ) +
    RotatedAxis() +
    ggtitle(
      "CD8 cluster lineage-adjudication markers"
    )
  
  ggsave(
    file.path(
      figure_dir,
      "phase4_cd8_lineage_adjudication_dotplot.png"
    ),
    p_lineage,
    width = 16,
    height = 9,
    dpi = 300
  )
}


###############################################################################
# 44. BIOLOGICAL STATE MARKER DOTPLOT
###############################################################################

state_marker_panels <- list(
  
  Conventional_T_CD8 = lineage_panels$Conventional_T_CD8,
  
  Naive_Memory = programs$Naive_Memory,
  
  Cytotoxic = programs$Cytotoxic,
  
  Activation = c(
    "TNFRSF9",
    "DUSP4",
    "DUSP16",
    "CD69",
    "TNF",
    "FOS",
    "JUN",
    "NR4A1"
  ),
  
  Dysfunction_Associated = c(
    "PDCD1",
    "TOX",
    "TIGIT",
    "CTLA4",
    "HAVCR2",
    "LAG3",
    "CXCR4",
    "RGS1",
    "LAYN"
  ),
  
  NK_Associated = lineage_panels$NK_Associated,
  
  GammaDelta_T = lineage_panels$GammaDelta_T,
  
  Treg_Associated = lineage_panels$Treg_Associated,
  
  Proliferation = programs$Proliferation
)

state_marker_panels_present <- lapply(
  state_marker_panels,
  function(x) {
    intersect(
      unique(x),
      rownames(cd8_obj)
    )
  }
)

state_marker_panels_present <- state_marker_panels_present[
  lengths(state_marker_panels_present) >= 2
]

state_marker_genes <- unique(
  unlist(
    state_marker_panels_present
  )
)

if (length(state_marker_genes) > 0) {
  
  p_state_markers <- DotPlot(
    cd8_obj,
    features = state_marker_genes,
    group.by = "CD8_Cluster"
  ) +
    RotatedAxis() +
    ggtitle(
      "CD8 biological state and lineage marker panels"
    )
  
  ggsave(
    file.path(
      figure_dir,
      "phase4_cd8_state_marker_panels_dotplot.png"
    ),
    p_state_markers,
    width = 18,
    height = 10,
    dpi = 300
  )
}


###############################################################################
# 45. TOP MARKER EVIDENCE
###############################################################################

top_marker_lookup <- cd8_markers %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  
  dplyr::slice_head(
    n = 50
  ) %>%
  
  dplyr::summarise(
    
    Top50_Genes = paste(
      gene,
      collapse = ";"
    ),
    
    .groups = "drop"
  )

write.csv(
  top_marker_lookup,
  file.path(
    table_dir,
    "phase4_cd8_top50_marker_evidence.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 46. CLUSTER REPRESENTATION + PROGRAM EVIDENCE
###############################################################################

annotation_evidence <- cluster_representation

if (length(program_columns) > 0) {
  
  annotation_evidence <-
    annotation_evidence %>%
    
    dplyr::left_join(
      program_cluster_summary,
      by = "CD8_Cluster"
    )
}

annotation_evidence <-
  annotation_evidence %>%
  
  dplyr::left_join(
    top_marker_lookup,
    by = c(
      "CD8_Cluster" = "cluster"
    )
  )

write.csv(
  annotation_evidence,
  file.path(
    table_dir,
    "phase4_cd8_annotation_evidence.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 47. CLUSTER-LEVEL PROGRAM MATRIX
###############################################################################

if (length(program_columns) > 0) {
  
  cluster_program_matrix <- cd8_obj@meta.data %>%
    
    dplyr::group_by(
      CD8_Cluster
    ) %>%
    
    dplyr::summarise(
      
      dplyr::across(
        dplyr::all_of(program_columns),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    )
  
} else {
  
  cluster_program_matrix <- NULL
}


###############################################################################
# 48. STANDARDIZED PROGRAM SCORES
###############################################################################

if (
  !is.null(cluster_program_matrix) &&
  length(program_columns) > 0
) {
  
  for (col in program_columns) {
    
    cluster_program_matrix[[paste0(col, "_Z")]] <-
      as.numeric(
        scale(
          cluster_program_matrix[[col]]
        )
      )
  }
}


###############################################################################
# 49. LINEAGE EVIDENCE AT CLUSTER LEVEL
###############################################################################

cluster_lineage_evidence <- NULL

for (panel_name in names(lineage_panels_present)) {
  
  genes <- lineage_panels_present[[panel_name]]
  
  expression_matrix <- LayerData(
    cd8_obj,
    assay = "RNA",
    layer = "data"
  )
  
  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )
  
  if (length(genes) < 2) {
    next
  }
  
  panel_values <- Matrix::colMeans(
    expression_matrix[
      genes,
      ,
      drop = FALSE
    ]
  )
  
  panel_df <- tibble(
    Cell = names(panel_values),
    Panel_Score = as.numeric(panel_values)
  )
  
  panel_df$CD8_Cluster <-
    cd8_obj$CD8_Cluster[
      match(
        panel_df$Cell,
        colnames(cd8_obj)
      )
    ]
  
  panel_summary <- panel_df %>%
    
    dplyr::group_by(
      CD8_Cluster
    ) %>%
    
    dplyr::summarise(
      
      !!paste0(
        "Lineage_",
        panel_name,
        "_MeanExpression"
      ) :=
        mean(
          Panel_Score,
          na.rm = TRUE
        ),
      
      !!paste0(
        "Lineage_",
        panel_name,
        "_PositiveFraction"
      ) :=
        mean(
          Panel_Score > 0,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  if (is.null(cluster_lineage_evidence)) {
    
    cluster_lineage_evidence <-
      panel_summary
    
  } else {
    
    cluster_lineage_evidence <-
      dplyr::left_join(
        cluster_lineage_evidence,
        panel_summary,
        by = "CD8_Cluster"
      )
  }
}

if (!is.null(cluster_lineage_evidence)) {
  
  write.csv(
    cluster_lineage_evidence,
    file.path(
      table_dir,
      "phase4_cd8_cluster_lineage_evidence.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 50. MARKER-PANEL OVERLAP WITH TOP CLUSTER MARKERS
###############################################################################

marker_overlap_rows <- list()

for (cluster_id in unique(cd8_markers$cluster)) {
  
  cluster_top50 <- cd8_markers %>%
    
    dplyr::filter(
      cluster == cluster_id
    ) %>%
    
    dplyr::arrange(
      dplyr::desc(avg_log2FC)
    ) %>%
    
    dplyr::slice_head(
      n = 50
    ) %>%
    
    dplyr::pull(
      gene
    )
  
  for (panel_name in names(state_marker_panels_present)) {
    
    panel_genes <-
      state_marker_panels_present[[panel_name]]
    
    overlap_genes <- intersect(
      cluster_top50,
      panel_genes
    )
    
    marker_overlap_rows[[length(marker_overlap_rows) + 1]] <-
      tibble(
        
        CD8_Cluster = as.character(cluster_id),
        
        Panel = panel_name,
        
        Top50_Overlap_Count =
          length(overlap_genes),
        
        Top50_Overlap_Genes =
          paste(
            overlap_genes,
            collapse = ";"
          )
      )
  }
}

marker_panel_overlap <- bind_rows(
  marker_overlap_rows
)

write.csv(
  marker_panel_overlap,
  file.path(
    table_dir,
    "phase4_cd8_marker_panel_overlap.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 51. CREATE A STRUCTURED ANNOTATION ADJUDICATION TABLE
###############################################################################
#
# IMPORTANT:
#
# We intentionally DO NOT automatically assign final biological states here.
#
# Instead, the table organizes all evidence needed for adjudication:
#
#   - cluster size
#   - donor representation
#   - lineage evidence
#   - transcriptional programs
#   - top markers
#
###############################################################################

annotation_adjudication <- annotation_evidence

if (!is.null(cluster_lineage_evidence)) {
  
  annotation_adjudication <-
    annotation_adjudication %>%
    
    dplyr::left_join(
      cluster_lineage_evidence,
      by = "CD8_Cluster"
    )
}

if (!is.null(cluster_program_matrix)) {
  
  annotation_adjudication <-
    annotation_adjudication %>%
    
    dplyr::left_join(
      cluster_program_matrix,
      by = "CD8_Cluster",
      suffix = c(
        "",
        "_ProgramMatrix"
      )
    )
}

annotation_adjudication <-
  annotation_adjudication %>%
  
  dplyr::mutate(
    
    Potential_Noncanonical_Lineage =
      case_when(
        
        !is.na(
          Lineage_GammaDelta_T_MeanExpression
        ) &
          Lineage_GammaDelta_T_MeanExpression >
          Lineage_Conventional_T_CD8_MeanExpression
        ~ "GammaDelta_like_evidence",
        
        !is.na(
          Lineage_NK_Associated_MeanExpression
        ) &
          Lineage_NK_Associated_MeanExpression >
          Lineage_Conventional_T_CD8_MeanExpression
        ~ "NK_like_evidence",
        
        !is.na(
          Lineage_Treg_Associated_MeanExpression
        ) &
          Lineage_Treg_Associated_MeanExpression >
          Lineage_Conventional_T_CD8_MeanExpression
        ~ "Treg_associated_evidence",
        
        TRUE ~ "No_dominant_noncanonical_signal"
      ),
    
    CD8_Lineage_Evidence =
      case_when(
        
        Lineage_Conventional_T_CD8_MeanExpression >
          0
        ~ "CD8_T_evidence_present",
        
        TRUE ~ "CD8_lineage_evidence_uncertain"
      ),
    
    Annotation_Status =
      "REVIEW_REQUIRED"
  )

write.csv(
  annotation_adjudication,
  file.path(
    table_dir,
    "phase4_cd8_annotation_adjudication_table.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 52. ANNOTATION EVIDENCE HEATMAP
###############################################################################

if (
  !is.null(cluster_program_matrix) &&
  length(program_columns) > 0
) {
  
  program_heatmap_data <-
    cluster_program_matrix %>%
    
    dplyr::select(
      CD8_Cluster,
      dplyr::all_of(
        program_columns
      )
    ) %>%
    
    tibble::column_to_rownames(
      "CD8_Cluster"
    )
  
  program_heatmap_matrix <-
    as.matrix(
      program_heatmap_data
    )
  
  if (
    nrow(program_heatmap_matrix) > 1 &&
    ncol(program_heatmap_matrix) > 1
  ) {
    
    program_heatmap_matrix <-
      t(
        scale(
          program_heatmap_matrix
        )
      )
    
    heatmap_df <-
      as.data.frame(
        program_heatmap_matrix
      )
    
    heatmap_df$Program <-
      rownames(
        heatmap_df
      )
    
    heatmap_long <-
      tidyr::pivot_longer(
        heatmap_df,
        cols = -Program,
        names_to = "CD8_Cluster",
        values_to = "Z"
      )
    
    p_program_heatmap <-
      ggplot(
        heatmap_long,
        aes(
          x = CD8_Cluster,
          y = Program,
          fill = Z
        )
      ) +
      
      geom_tile() +
      
      scale_fill_gradient2(
        midpoint = 0
      ) +
      
      theme_minimal() +
      
      theme(
        axis.text.x =
          element_text(
            angle = 45,
            hjust = 1
          )
      ) +
      
      labs(
        title =
          "CD8 transcriptional program landscape",
        x = "CD8 cluster",
        y = "Program",
        fill = "Cluster Z-score"
      )
    
    ggsave(
      file.path(
        figure_dir,
        "phase4_cd8_program_heatmap.png"
      ),
      p_program_heatmap,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}


###############################################################################
# 53. NONCANONICAL / AMBIGUITY FLAGS
###############################################################################

if (!is.null(annotation_adjudication)) {
  
  noncanonical_summary <-
    annotation_adjudication %>%
    
    dplyr::select(
      CD8_Cluster,
      Cluster_Cells,
      Donors_Represented,
      Largest_Donor_Fraction,
      Potential_Noncanonical_Lineage,
      CD8_Lineage_Evidence,
      Annotation_Status
    )
  
  write.csv(
    noncanonical_summary,
    file.path(
      table_dir,
      "phase4_cd8_noncanonical_review_summary.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 54. MANUAL BIOLOGICAL ANNOTATION TEMPLATE
###############################################################################

manual_annotation_template <-
  annotation_adjudication %>%
  
  dplyr::select(
    CD8_Cluster,
    Cluster_Cells,
    Donors_Represented,
    Clinical_States_Represented,
    Largest_Donor_Fraction,
    Potential_Noncanonical_Lineage,
    CD8_Lineage_Evidence,
    Top50_Genes
  ) %>%
  
  dplyr::mutate(
    
    Final_Annotation = NA_character_,
    
    Biological_Class = NA_character_,
    
    Annotation_Confidence = NA_character_,
    
    Evidence_Rationale = NA_character_
  )

write.csv(
  manual_annotation_template,
  file.path(
    table_dir,
    "phase4_cd8_manual_annotation_template.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 55. UMAP BY CLUSTER
###############################################################################
#
# Already generated above, retained as primary discovery representation.
###############################################################################


###############################################################################
# 56. CLUSTER × DONOR REPRESENTATION FIGURE
###############################################################################

cluster_donor_fraction <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    Donor,
    name = "Cells"
  ) %>%
  
  dplyr::group_by(
    Donor
  ) %>%
  
  dplyr::mutate(
    Fraction_of_Donor_CD8 =
      Cells / sum(Cells)
  ) %>%
  
  dplyr::ungroup()

p_cluster_donor <-
  ggplot(
    cluster_donor_fraction,
    aes(
      x = CD8_Cluster,
      y = Fraction_of_Donor_CD8
    )
  ) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.35,
    size = 0.8
  ) +
  
  theme_minimal() +
  
  labs(
    title =
      "CD8 cluster representation across donors",
    x = "CD8 cluster",
    y = "Fraction of donor CD8 cells"
  )

ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_cluster_representation_by_donor.png"
  ),
  p_cluster_donor,
  width = 11,
  height = 7,
  dpi = 300
)


###############################################################################
# 57. DESCRIPTIVE CLUSTER COMPOSITION ACROSS CLINICAL STATES
###############################################################################
#
# IMPORTANT:
#
# This remains descriptive.
#
# These cell-level proportions do NOT constitute clinical-state inference.
# Formal donor-aware inference belongs to Phase 6.
###############################################################################

cluster_phase_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    Phase,
    name = "Cells"
  )

write.csv(
  cluster_phase_counts,
  file.path(
    table_dir,
    "phase4_cd8_cluster_phase_counts.csv"
  ),
  row.names = FALSE
)

cluster_phase_proportions <- cluster_phase_counts %>%
  
  dplyr::group_by(
    Phase
  ) %>%
  
  dplyr::mutate(
    Proportion_Within_CD8 =
      Cells / sum(Cells)
  ) %>%
  
  dplyr::ungroup()

write.csv(
  cluster_phase_proportions,
  file.path(
    table_dir,
    "phase4_cd8_cluster_phase_proportions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 58. DONOR-LEVEL DESCRIPTIVE CLUSTER COMPOSITION
###############################################################################

donor_cluster_counts <- cd8_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    Phase,
    CD8_Cluster,
    name = "Cells"
  )

donor_cluster_totals <- cd8_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    name = "Total_CD8_Cells"
  )

donor_cluster_proportions <-
  donor_cluster_counts %>%
  
  dplyr::left_join(
    donor_cluster_totals,
    by = "Donor"
  ) %>%
  
  dplyr::mutate(
    Proportion_of_Donor_CD8 =
      Cells / Total_CD8_Cells,
    
    Percentage_of_Donor_CD8 =
      100 * Proportion_of_Donor_CD8
  )

write.csv(
  donor_cluster_proportions,
  file.path(
    table_dir,
    "phase4_cd8_donor_cluster_proportions.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 59. DONOR-LEVEL SUMMARY BY CLINICAL STATE
###############################################################################

donor_cluster_summary <-
  donor_cluster_proportions %>%
  
  dplyr::group_by(
    Phase,
    CD8_Cluster
  ) %>%
  
  dplyr::summarise(
    
    n_donors =
      dplyr::n_distinct(
        Donor
      ),
    
    mean_percentage =
      mean(
        Percentage_of_Donor_CD8,
        na.rm = TRUE
      ),
    
    median_percentage =
      median(
        Percentage_of_Donor_CD8,
        na.rm = TRUE
      ),
    
    sd_percentage =
      sd(
        Percentage_of_Donor_CD8,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write.csv(
  donor_cluster_summary,
  file.path(
    table_dir,
    "phase4_cd8_donor_cluster_summary_by_phase.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 60. FINAL PHASE 4 METADATA
###############################################################################

phase4_analysis_vector <- rep(
  "CD8_state_discovery_and_biological_adjudication",
  ncol(cd8_obj)
)

names(phase4_analysis_vector) <-
  colnames(cd8_obj)

cd8_obj <- AddMetaData(
  cd8_obj,
  metadata = phase4_analysis_vector,
  col.name = "Phase4_Analysis"
)


phase4_expression_vector <- rep(
  "Processed_log_CP10K",
  ncol(cd8_obj)
)

names(phase4_expression_vector) <-
  colnames(cd8_obj)

cd8_obj <- AddMetaData(
  cd8_obj,
  metadata = phase4_expression_vector,
  col.name = "Phase4_Expression"
)


phase4_donor_unit_vector <- rep(
  "Biological_replicate",
  ncol(cd8_obj)
)

names(phase4_donor_unit_vector) <-
  colnames(cd8_obj)

cd8_obj <- AddMetaData(
  cd8_obj,
  metadata = phase4_donor_unit_vector,
  col.name = "Phase4_Donor_Unit"
)


phase4_inference_vector <- rep(
  "Exploratory_cell_state_characterization",
  ncol(cd8_obj)
)

names(phase4_inference_vector) <-
  colnames(cd8_obj)

cd8_obj <- AddMetaData(
  cd8_obj,
  metadata = phase4_inference_vector,
  col.name = "Phase4_Inference"
)


###############################################################################
# 61. FINAL CELL-ID INTEGRITY CHECK
###############################################################################

object_cells <- colnames(cd8_obj)

if (!identical(
  rownames(cd8_obj@meta.data),
  object_cells
)) {
  
  stop(
    "FINAL VALIDATION FAILED: metadata cell names do not exactly match ",
    "Seurat object cell names."
  )
}

for (assay_name in names(cd8_obj@assays)) {
  
  assay_cells <- Cells(
    cd8_obj[[assay_name]]
  )
  
  if (!setequal(
    assay_cells,
    object_cells
  )) {
    
    stop(
      "FINAL VALIDATION FAILED: assay ",
      assay_name,
      " contains cells not present in the Seurat object."
    )
  }
}

for (reduction_name in names(cd8_obj@reductions)) {
  
  reduction_cells <- Cells(
    cd8_obj[[reduction_name]]
  )
  
  if (!setequal(
    reduction_cells,
    object_cells
  )) {
    
    stop(
      "FINAL VALIDATION FAILED: reduction ",
      reduction_name,
      " contains cells not present in the Seurat object."
    )
  }
}

for (graph_name in names(cd8_obj@graphs)) {
  
  graph_cells <- rownames(
    cd8_obj[[graph_name]]
  )
  
  if (!setequal(
    graph_cells,
    object_cells
  )) {
    
    stop(
      "FINAL VALIDATION FAILED: graph ",
      graph_name,
      " contains cells not present in the Seurat object."
    )
  }
}

if (!identical(
  names(Idents(cd8_obj)),
  object_cells
)) {
  
  stop(
    "FINAL VALIDATION FAILED: active identities are not aligned ",
    "with Seurat cell names."
  )
}

message(
  "\nCell-ID integrity checks passed."
)


###############################################################################
# 62. FINAL SEURAT VALIDITY CHECK
###############################################################################

validity_result <- tryCatch(
  
  {
    
    validObject(
      cd8_obj
    )
    
    TRUE
  },
  
  error = function(e) {
    
    message(
      "\nSeurat validity check FAILED:"
    )
    
    message(
      e$message
    )
    
    FALSE
  }
)

if (!validity_result) {
  
  stop(
    "Phase 4 object is invalid and will NOT be saved."
  )
}

message(
  "\nSeurat object validity check passed."
)


###############################################################################
# 63. FINAL PHASE 4 EVIDENCE SUMMARY
###############################################################################

final_cluster_summary <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    name = "Cells"
  ) %>%
  
  dplyr::left_join(
    cluster_representation,
    by = "CD8_Cluster",
    suffix = c(
      "_Final",
      ""
    )
  ) %>%
  
  dplyr::arrange(
    suppressWarnings(
      as.numeric(CD8_Cluster)
    )
  )

write.csv(
  final_cluster_summary,
  file.path(
    table_dir,
    "phase4_cd8_final_cluster_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 64. FINAL OBJECT VALIDATION
###############################################################################

if (ncol(cd8_obj) != expected_cd8_n) {
  
  stop(
    "Final CD8 cell count changed unexpectedly.\n",
    "Expected: ",
    expected_cd8_n,
    "\nFound: ",
    ncol(cd8_obj)
  )
}

if (nrow(cd8_obj) != expected_total_features) {
  
  stop(
    "Final CD8 feature count changed unexpectedly.\n",
    "Expected: ",
    expected_total_features,
    "\nFound: ",
    nrow(cd8_obj)
  )
}

required_reductions <- c(
  "pca",
  "umap_cd8"
)

missing_reductions <- setdiff(
  required_reductions,
  Reductions(cd8_obj)
)

if (length(missing_reductions) > 0) {
  
  stop(
    "Missing required reductions: ",
    paste(
      missing_reductions,
      collapse = ", "
    )
  )
}

required_final_metadata <- c(
  "CD8_Cluster",
  "Phase4_Analysis",
  "Phase4_Expression",
  "Phase4_Donor_Unit",
  "Phase4_Inference"
)

missing_final_metadata <- setdiff(
  required_final_metadata,
  colnames(cd8_obj@meta.data)
)

if (length(missing_final_metadata) > 0) {
  
  stop(
    "Final object missing required metadata: ",
    paste(
      missing_final_metadata,
      collapse = ", "
    )
  )
}


###############################################################################
# 65. SAVE FINAL PHASE 4 OBJECT
###############################################################################

saveRDS(
  cd8_obj,
  final_output
)


###############################################################################
# 66. PHASE 4 COMPUTATIONAL PIPELINE END
###############################################################################
#
# The computational Phase 4 object has now been generated and saved.
#
# Final biological adjudication and the authoritative Phase 4 status check
# are performed in Section 67.
#
###############################################################################

message("")
message("============================================================")
message("PHASE 4 COMPUTATIONAL PIPELINE COMPLETE")
message("============================================================")

message(
  "Primary computational object saved:"
)

message(
  final_output
)

message("")
message(
  "Final biological adjudication and Phase 4 closeout continue in Section 67."
)

message("============================================================")


###############################################################################
# 67. FINAL BIOLOGICAL ADJUDICATION AND PHASE 4 CLOSEOUT
###############################################################################
#
# Purpose
# -------
# Convert the evidence assembled during Sections 34–54 into a transparent,
# conservative biological adjudication of the 16 CD8-labelled clusters.
#
# This section does NOT:
#   - recluster cells
#   - change PCA/UMAP parameters
#   - change cluster assignments
#   - perform clinical-state inference
#   - perform differential expression
#   - perform pathway enrichment
#   - perform cell-cell communication analysis
#
# It adds explicitly documented biological annotations to the already
# validated Phase 4 object.
#
# Important interpretation boundary
# ----------------------------------
# These annotations are computational, transcriptional characterizations.
# They are not experimentally validated identities, clinical-state labels,
# causal mechanisms, or claims of disease progression.
#
###############################################################################

message("")
message("============================================================")
message("SECTION 67: FINAL BIOLOGICAL ADJUDICATION AND PHASE 4 CLOSEOUT")
message("============================================================")


###############################################################################
# 67.1 FINAL NONCANONICAL LINEAGE EVIDENCE FLAGS
###############################################################################

# These are REVIEW FLAGS based on the lineage evidence already generated in
# Sections 49-53.
#
# They are NOT automated cell-type classifications.
# They document where the existing evidence supports consideration of a
# noncanonical population during biological adjudication.

annotation_adjudication_final <- annotation_adjudication %>%
  
  dplyr::mutate(
    
    Potential_Noncanonical_Lineage = dplyr::case_when(
      
      CD8_Cluster %in% c(7, 11, 13, 15) ~
        "GammaDelta_associated",
      
      CD8_Cluster %in% c(6, 8, 9, 14) ~
        "NK_associated",
      
      CD8_Cluster == 10 ~
        "NK_GammaDelta_associated",
      
      TRUE ~
        "No_strong_noncanonical_signal"
    )
  )


write.csv(
  annotation_adjudication_final,
  file.path(
    table_dir,
    "phase4_cd8_annotation_adjudication_table_final.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 67.2 FINAL BIOLOGICAL ADJUDICATION OF THE 16 CD8 CLUSTERS
###############################################################################
#
# These annotations are conservative biological interpretations of the
# transcriptional evidence generated in Phase 4.
#
# They are NOT claims of experimentally validated cell identity.
# They are NOT clinical-state labels.
# They are NOT pathway annotations.
# They are NOT causal biological claims.
###############################################################################

final_cd8_annotations <- tibble::tribble(
  
  ~CD8_Cluster, ~CD8_Biological_Annotation, ~CD8_Biological_Class,
  ~CD8_Annotation_Confidence, ~Donor_Representation_Note,
  
  0,
  "GZMK-associated memory-like CD8 T-cell",
  "Conventional_CD8_associated",
  "Moderate",
  "Broad donor representation",
  
  1,
  "Memory-associated conventional T-cell",
  "Conventional_T_associated",
  "Moderate",
  "Broad donor representation with donor concentration",
  
  2,
  "Activated/dysfunction-associated CD8 T-cell",
  "Conventional_CD8_associated",
  "High",
  "Broad donor representation with donor concentration",
  
  3,
  "Immediate-early-response CD8-associated population",
  "Transcriptional_program",
  "High",
  "Broad donor representation",
  
  4,
  "Activated/dysfunction-associated CD8 T-cell",
  "Conventional_CD8_associated",
  "High",
  "Broad donor representation with donor concentration",
  
  5,
  "Activation/tissue-associated CD8-associated population",
  "Conventional_CD8_associated",
  "Moderate",
  "Broad donor representation",
  
  6,
  "NK-like cytotoxic nonconventional T-cell population",
  "Noncanonical_T_associated",
  "High",
  "Broad donor representation",
  
  7,
  "Gamma-delta-like cytotoxic T-cell population",
  "GammaDelta_associated",
  "High",
  "Strongly donor concentrated",
  
  8,
  "NK-like/nonconventional cytotoxic T-cell population",
  "Noncanonical_T_associated",
  "High",
  "Broad donor representation",
  
  9,
  "NK-like cytotoxic population",
  "NK_associated",
  "High",
  "Strongly donor concentrated",
  
  10,
  "NK/gamma-delta-like nonconventional T-cell population",
  "Noncanonical_T_associated",
  "High",
  "Broad donor representation with donor concentration",
  
  11,
  "Gamma-delta-like T-cell population",
  "GammaDelta_associated",
  "High",
  "Broad donor representation",
  
  12,
  "Activated/dysfunction-associated tissue-associated CD8 T-cell",
  "Conventional_CD8_associated",
  "High",
  "Donor concentrated",
  
  13,
  "Gamma-delta-like cytotoxic T-cell population",
  "GammaDelta_associated",
  "High",
  "Extremely donor concentrated",
  
  14,
  "NK-like cytotoxic population",
  "NK_associated",
  "High",
  "Strongly donor concentrated",
  
  15,
  "Gamma-delta-like nonconventional T-cell population",
  "GammaDelta_associated",
  "High",
  "Moderately donor concentrated"
)


if (
  nrow(final_cd8_annotations) !=
  length(unique(cd8_obj$CD8_Cluster))
) {
  
  stop(
    "Final annotation table does not contain exactly one row per CD8 cluster."
  )
}


if (
  !setequal(
    final_cd8_annotations$CD8_Cluster,
    unique(as.numeric(as.character(cd8_obj$CD8_Cluster)))
  )
) {
  
  stop(
    "Final annotation clusters do not match the clusters present in the Seurat object."
  )
}


write.csv(
  final_cd8_annotations,
  file.path(
    table_dir,
    "phase4_cd8_final_adjudicated_annotations.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 67.3 ADD FINAL BIOLOGICAL ANNOTATIONS TO SEURAT METADATA
###############################################################################

# Use AddMetaData() with named cell-level vectors to preserve cell-ID alignment.

annotation_columns <- c(
  "CD8_Biological_Annotation",
  "CD8_Biological_Class",
  "CD8_Annotation_Confidence",
  "Donor_Representation_Note"
)


for (annotation_column in annotation_columns) {
  
  cluster_to_value <- setNames(
    final_cd8_annotations[[annotation_column]],
    as.character(
      final_cd8_annotations$CD8_Cluster
    )
  )
  
  cell_values <- unname(
    cluster_to_value[
      as.character(
        cd8_obj$CD8_Cluster
      )
    ]
  )
  
  names(cell_values) <- colnames(cd8_obj)
  
  cd8_obj <- AddMetaData(
    cd8_obj,
    metadata = cell_values,
    col.name = annotation_column
  )
}


###############################################################################
# 67.4 VALIDATE FINAL BIOLOGICAL ANNOTATION METADATA
###############################################################################

if (
  any(
    is.na(
      cd8_obj$CD8_Biological_Annotation
    )
  )
) {
  
  stop(
    "Some CD8 cells did not receive a final biological annotation."
  )
}


if (
  any(
    is.na(
      cd8_obj$CD8_Biological_Class
    )
  )
) {
  
  stop(
    "Some CD8 cells did not receive a biological class."
  )
}


if (
  !identical(
    rownames(cd8_obj@meta.data),
    colnames(cd8_obj)
  )
) {
  
  stop(
    "Cell-ID alignment failed after adding final biological annotations."
  )
}


message(
  "\nFinal biological annotations successfully added to metadata."
)


###############################################################################
# 67.5 FINAL BIOLOGICAL ANNOTATION SUMMARY
###############################################################################

final_annotation_summary <- cd8_obj@meta.data %>%
  
  dplyr::count(
    CD8_Cluster,
    CD8_Biological_Annotation,
    CD8_Biological_Class,
    CD8_Annotation_Confidence,
    name = "Cells"
  ) %>%
  
  dplyr::arrange(
    as.numeric(
      as.character(
        CD8_Cluster
      )
    )
  )


write.csv(
  final_annotation_summary,
  file.path(
    table_dir,
    "phase4_cd8_final_biological_annotation_summary.csv"
  ),
  row.names = FALSE
)


message(
  "\nFinal biological annotation summary:"
)

print(
  final_annotation_summary
)


###############################################################################
# 67.6 FINAL BIOLOGICAL ANNOTATION UMAP
###############################################################################

p_final_annotation <- DimPlot(
  cd8_obj,
  reduction = "umap_cd8",
  group.by = "CD8_Biological_Annotation",
  label = TRUE,
  repel = TRUE
) +
  
  ggplot2::ggtitle(
    "CD8-labelled compartment: final biological adjudication"
  ) +
  
  ggplot2::theme_classic()


ggsave(
  file.path(
    figure_dir,
    "phase4_cd8_final_biological_annotations.png"
  ),
  p_final_annotation,
  width = 12,
  height = 8,
  dpi = 300
)


###############################################################################
# 67.7 FINAL DONOR REPRESENTATION SUMMARY
###############################################################################

final_donor_annotation_summary <- cd8_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    Phase,
    CD8_Cluster,
    CD8_Biological_Annotation,
    name = "Cells"
  ) %>%
  
  dplyr::group_by(
    Donor
  ) %>%
  
  dplyr::mutate(
    Total_CD8_Cells = sum(Cells),
    Fraction_of_Donor_CD8 = Cells / Total_CD8_Cells
  ) %>%
  
  dplyr::ungroup()


write.csv(
  final_donor_annotation_summary,
  file.path(
    table_dir,
    "phase4_cd8_final_donor_annotation_representation.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 67.8 FINAL PHASE 4 VALIDATION
###############################################################################

message(
  "\nRunning final Phase 4 closeout validation..."
)


# Cell IDs
object_cells_final <- colnames(cd8_obj)


if (
  !identical(
    rownames(cd8_obj@meta.data),
    object_cells_final
  )
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: metadata cell names are not aligned."
  )
}


# Expected cell count
if (
  ncol(cd8_obj) != expected_cd8_n
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: CD8 cell count changed."
  )
}


# Expected feature count
if (
  nrow(cd8_obj) != expected_total_features
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: feature count changed."
  )
}


# Donor count
if (
  length(
    unique(
      cd8_obj$Donor
    )
  ) != 23
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: donor count changed."
  )
}


# Clinical-state count
if (
  length(
    unique(
      cd8_obj$Phase
    )
  ) != 5
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: clinical-state count changed."
  )
}


# Cluster count
if (
  length(
    unique(
      cd8_obj$CD8_Cluster
    )
  ) != 16
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: expected 16 CD8 clusters."
  )
}


# Biological annotation coverage
#
# Multiple computational clusters may legitimately share the same
# biological annotation. Therefore, the correct validation is that
# every one of the 16 clusters has a non-missing annotation, NOT
# that there are 16 unique annotation strings.

cluster_annotation_check <- cd8_obj@meta.data %>%
  dplyr::distinct(
    CD8_Cluster,
    CD8_Biological_Annotation
  )

if (
  nrow(cluster_annotation_check) != 16
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: not all 16 CD8 clusters ",
    "have exactly one biological annotation."
  )
}


if (
  any(
    is.na(
      cluster_annotation_check$CD8_Biological_Annotation
    )
  )
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: at least one CD8 cluster ",
    "has a missing biological annotation."
  )
}


message(
  "Biological annotation coverage validated: ",
  nrow(cluster_annotation_check),
  " clusters annotated; ",
  length(
    unique(
      cd8_obj$CD8_Biological_Annotation
    )
  ),
  " unique biological annotation categories."
)


# Required reductions
required_phase4_reductions <- c(
  "pca",
  "umap_cd8"
)


missing_phase4_reductions <- setdiff(
  required_phase4_reductions,
  Reductions(cd8_obj)
)


if (
  length(missing_phase4_reductions) > 0
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: missing reductions: ",
    paste(
      missing_phase4_reductions,
      collapse = ", "
    )
  )
}


# Required final metadata
required_phase4_metadata <- c(
  "CD8_Cluster",
  "CD8_Biological_Annotation",
  "CD8_Biological_Class",
  "CD8_Annotation_Confidence",
  "Donor_Representation_Note",
  "Phase4_Analysis",
  "Phase4_Expression",
  "Phase4_Donor_Unit",
  "Phase4_Inference"
)


missing_phase4_metadata <- setdiff(
  required_phase4_metadata,
  colnames(cd8_obj@meta.data)
)


if (
  length(missing_phase4_metadata) > 0
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: missing metadata: ",
    paste(
      missing_phase4_metadata,
      collapse = ", "
    )
  )
}


# Seurat validity
final_validity_result <- validObject(
  cd8_obj,
  test = TRUE
)


if (
  !isTRUE(final_validity_result)
) {
  
  stop(
    "FINAL PHASE 4 VALIDATION FAILED: Seurat object is invalid."
  )
}


message(
  "Final Phase 4 validation passed."
)


###############################################################################
# 67.9 SAVE FINAL ADJUDICATED PHASE 4 OBJECT
###############################################################################

final_adjudicated_output <- file.path(
  rds_dir,
  "phase4_final_cd8_analysis_adjudicated.rds"
)


saveRDS(
  cd8_obj,
  final_adjudicated_output
)


if (
  !file.exists(
    final_adjudicated_output
  )
) {
  
  stop(
    "Final adjudicated Phase 4 RDS was not created."
  )
}


message(
  "\nFinal adjudicated Phase 4 object saved:"
)

message(
  final_adjudicated_output
)


###############################################################################
# 67.10 WRITE PHASE 4 SUMMARY
###############################################################################

phase4_summary <- c(
  
  "PHASE 4: CD8 T-CELL TRANSCRIPTIONAL STATE ANALYSIS",
  "",
  
  "Dataset:",
  "GSE182159 liver-only single-cell RNA-seq dataset.",
  "",
  
  paste0(
    "Input compartment: ",
    ncol(cd8_obj),
    " cells labelled CD8_T in Phase 3, representing ",
    length(unique(cd8_obj$Donor)),
    " donors across ",
    length(unique(cd8_obj$Phase)),
    " clinical states."
  ),
  "",
  
  "Computational analysis:",
  "CD8-specific variable feature selection, PCA, nearest-neighbor graph",
  "construction, SNN clustering, UMAP visualization, exploratory cluster-marker",
  "identification, biological program scoring, lineage-panel assessment, and",
  "donor-representation diagnostics were performed using the existing processed",
  "log-CP10K expression representation.",
  "",
  
  "Biological characterization:",
  paste0(
    "Sixteen transcriptional clusters were identified within the Phase-3 ",
    "CD8-labelled compartment."
  ),
  "",
  
  "The cluster evidence supported:",
  "- conventional CD8-associated transcriptional populations;",
  "- activation/dysfunction-associated transcriptional states;",
  "- an immediate-early-response transcriptional program;",
  "- NK-like/noncanonical populations; and",
  "- gamma-delta-like/noncanonical T-cell populations.",
  "",
  
  "Interpretation boundary:",
  "Cluster-level composition is descriptive. Donor was retained as the biological",
  "replicate. Several clusters showed substantial donor concentration and were",
  "therefore not interpreted as cohort-wide clinical-state features.",
  "",
  
  "No clinical-state differential inference was performed in Phase 4.",
  "Formal donor-aware clinical-state inference is reserved for downstream analysis.",
  "",
  
  "No pathway enrichment or cell-cell communication analysis was performed in",
  "Phase 4. These analyses are reserved for later phases.",
  "",
  
  "Final interpretation:",
  "Phase 4 establishes a biologically adjudicated map of transcriptional",
  "heterogeneity within the Phase-3 CD8-labelled compartment for use in",
  "downstream donor-aware analysis.",
  "",
  
  "Important limitation:",
  "The biological annotations represent computational transcriptional",
  "characterization of the public dataset and are not experimentally validated",
  "cell identities or causal biological states."
)


writeLines(
  phase4_summary,
  con = file.path(
    table_dir,
    "phase4_cd8_summary.txt"
  )
)


###############################################################################
# 67.11 FINAL PHASE 4 STATUS AND COMPLETION
###############################################################################

message("")
message("============================================================")
message("PHASE 4 FULLY CLOSED")
message("============================================================")

message(
  "CD8-labelled cells: ",
  ncol(cd8_obj)
)

message(
  "Features: ",
  nrow(cd8_obj)
)

message(
  "Donors: ",
  length(
    unique(
      cd8_obj$Donor
    )
  )
)

message(
  "Clinical states: ",
  paste(
    sort(
      unique(
        cd8_obj$Phase
      )
    ),
    collapse = ", "
  )
)

message(
  "Computational clusters: ",
  length(
    unique(
      cd8_obj$CD8_Cluster
    )
  )
)

message(
  "CD8 clusters with final biological annotations: ",
  length(
    unique(
      cd8_obj$CD8_Cluster
    )
  )
)

message(
  "Unique biological annotation categories: ",
  length(
    unique(
      cd8_obj$CD8_Biological_Annotation
    )
  )
)

message("")
message("Final adjudicated RDS:")
message(final_adjudicated_output)

message("")
message("Key Phase 4 outputs:")

message(
  file.path(
    table_dir,
    "phase4_cd8_final_adjudicated_annotations.csv"
  )
)

message(
  file.path(
    table_dir,
    "phase4_cd8_summary.txt"
  )
)

message(
  file.path(
    figure_dir,
    "phase4_cd8_final_biological_annotations.png"
  )
)

message("")
message("Phase 4 interpretation boundary:")

message(
  "The final annotations describe transcriptional heterogeneity within the",
  " Phase-3 CD8-labelled compartment. Donor-aware clinical-state inference",
  " remains outside Phase 4."
)

message("")
message("============================================================")
message("NEXT: PHASE 5 — MYELOID ANALYSIS")
message("============================================================")
