###############################################################################
# HBV scRNA-seq PROJECT
# PHASE 5 — MYELOID SUBSET ANALYSIS
#
# Goal:
#   Characterize transcriptional heterogeneity within the Phase 3
#   Monocyte_Myeloid compartment using the same analytical architecture
#   established for Phase 4 CD8 analysis.
#
# INPUT:
#   results/rds_objects/phase3_final_full_dataset.rds
#
# TARGET:
#   Final_Cell_Type == "Monocyte_Myeloid"
#
# EXPECTED TARGET:
#   2,489 cells
#   24,452 genes
#   23 donors
#   5 GEO-defined clinical states:
#       AC, AR, IA, IT, NL
#
# ANALYTICAL PRINCIPLES:
#   - Donor = biological replicate.
#   - Cells = units for transcriptional state discovery.
#   - No clinical-state inference in Phase 5.
#   - Clinical-state differential analysis is reserved for Phase 6.
#   - Existing processed log-CP10K expression is used as supplied.
#   - No reverse transformation to counts.
#   - No NormalizeData().
#   - No arbitrary cell removal based on subset QC flags.
#   - No causal interpretation.
#   - No claim that computational clusters are automatically biological states.
#   - Biological annotation requires multiple evidence streams.
#
# PHASE 5 OUTPUT:
#   results/rds_objects/phase5_final_myeloid_analysis_adjudicated.rds
#
# Software validated against:
#   Seurat 5.5.1
#   SeuratObject 5.4.0
###############################################################################


###############################################################################
# SECTION 1 — SETUP + REPRODUCIBILITY
###############################################################################

rm(list = ls())

set.seed(20260916)


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


###############################################################################
# 1.1 Explicit package validation
###############################################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "Matrix",
  "ggplot2"
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
# 1.2 Software version validation
###############################################################################

message("============================================================")
message("PHASE 5 — MYELOID SUBSET ANALYSIS")
message("============================================================")

message(
  "Seurat: ",
  as.character(packageVersion("Seurat"))
)

message(
  "SeuratObject: ",
  as.character(packageVersion("SeuratObject"))
)


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
# 1.3 Project paths
###############################################################################

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


###############################################################################
# 1.4 Input/output files
###############################################################################

input_file <- file.path(
  rds_dir,
  "phase3_final_full_dataset.rds"
)


source_layer_output <- file.path(
  rds_dir,
  "phase5_myeloid_source_layers_preserved.rds"
)


final_output <- file.path(
  rds_dir,
  "phase5_final_myeloid_analysis_adjudicated.rds"
)


###############################################################################
# 1.5 Locked constants
###############################################################################

expected_total_cells <- 106592L

expected_total_features <- 24452L

expected_myeloid_cells <- 2489L

expected_donors <- 23L

expected_gsm <- 23L

expected_hvgs <- 2000L

n_pcs <- 50L

analysis_dims <- 1:20

neighbor_k <- 20L

cluster_resolution <- 0.4

analysis_seed <- 20260916

target_cell_type <- "Monocyte_Myeloid"

expected_phases <- c(
  "AC",
  "AR",
  "IA",
  "IT",
  "NL"
)


message(
  "Target population: ",
  target_cell_type
)

message(
  "Expected target cells: ",
  expected_myeloid_cells
)

message(
  "Analysis seed: ",
  analysis_seed
)


###############################################################################
# SECTION 2 — LOAD + VALIDATE PHASE 3 FINAL DATASET
###############################################################################

message("")
message("============================================================")
message("SECTION 2 — LOAD + VALIDATE PHASE 3")
message("============================================================")


if (!file.exists(input_file)) {
  
  stop(
    "Phase 3 final object not found:\n",
    input_file,
    call. = FALSE
  )
  
}


seurat_obj <- readRDS(
  input_file
)


if (!inherits(seurat_obj, "Seurat")) {
  
  stop(
    "Input object is not a Seurat object.",
    call. = FALSE
  )
  
}


DefaultAssay(seurat_obj) <- "RNA"


###############################################################################
# 2.1 Global object dimensions
###############################################################################

message(
  "Phase 3 cells: ",
  ncol(seurat_obj)
)

message(
  "Phase 3 features: ",
  nrow(seurat_obj)
)


if (
  ncol(seurat_obj) !=
  expected_total_cells
) {
  
  stop(
    "Unexpected Phase 3 cell count.\n",
    "Expected: ",
    expected_total_cells,
    "\nFound: ",
    ncol(seurat_obj),
    call. = FALSE
  )
  
}


if (
  nrow(seurat_obj) !=
  expected_total_features
) {
  
  stop(
    "Unexpected Phase 3 feature count.\n",
    "Expected: ",
    expected_total_features,
    "\nFound: ",
    nrow(seurat_obj),
    call. = FALSE
  )
  
}


###############################################################################
# 2.2 Required metadata
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


if (length(missing_metadata) > 0L) {
  
  stop(
    "Missing required Phase 3 metadata:\n",
    paste(
      missing_metadata,
      collapse = ", "
    ),
    call. = FALSE
  )
  
}


###############################################################################
# 2.3 Validate clinical-state labels
###############################################################################

observed_phases <- sort(
  unique(
    as.character(
      seurat_obj$Phase
    )
  )
)


if (
  !setequal(
    observed_phases,
    expected_phases
  )
) {
  
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
    ),
    call. = FALSE
  )
  
}


###############################################################################
# 2.4 Validate donors and samples
###############################################################################

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


if (n_donors != expected_donors) {
  
  stop(
    "Expected 23 donors; found ",
    n_donors,
    call. = FALSE
  )
  
}


if (n_gsm != expected_gsm) {
  
  stop(
    "Expected 23 GSM/sample identifiers; found ",
    n_gsm,
    call. = FALSE
  )
  
}


message(
  "Validated ",
  n_donors,
  " donors and ",
  n_gsm,
  " samples."
)


###############################################################################
# 2.5 Phase 3 cell-type composition
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
  file = file.path(
    table_dir,
    "phase5_phase3_Final_Cell_Type_summary.csv"
  ),
  row.names = FALSE
)


message("")
message("Phase 3 cell-type composition:")
print(phase3_label_summary)


###############################################################################
# SECTION 3 — EXTRACT MYELOID SUBSET
###############################################################################

message("")
message("============================================================")
message("SECTION 3 — MYELOID SUBSET EXTRACTION")
message("============================================================")


myeloid_cells <- rownames(
  seurat_obj@meta.data
)[
  seurat_obj$Final_Cell_Type ==
    target_cell_type
]


if (length(myeloid_cells) == 0L) {
  
  stop(
    "No cells labelled Monocyte_Myeloid were found.",
    call. = FALSE
  )
  
}


myeloid_obj <- subset(
  seurat_obj,
  cells = myeloid_cells
)


DefaultAssay(myeloid_obj) <- "RNA"


expected_myeloid_n <- length(
  myeloid_cells
)


message(
  "Myeloid cells extracted: ",
  expected_myeloid_n
)


message(
  "Myeloid features: ",
  nrow(myeloid_obj)
)


if (
  expected_myeloid_n !=
  expected_myeloid_cells
) {
  
  stop(
    "Unexpected Monocyte_Myeloid cell count.\n",
    "Expected: ",
    expected_myeloid_cells,
    "\nFound: ",
    expected_myeloid_n,
    call. = FALSE
  )
  
}


if (
  !identical(
    colnames(myeloid_obj),
    rownames(myeloid_obj@meta.data)
  )
) {
  
  stop(
    "Cell names and metadata row names do not match.",
    call. = FALSE
  )
  
}


###############################################################################
# 3.1 Confirm exact target identity
###############################################################################

if (
  !all(
    myeloid_obj$Final_Cell_Type ==
    target_cell_type
  )
) {
  
  stop(
    "Extracted object contains cells outside Monocyte_Myeloid.",
    call. = FALSE
  )
  
}


###############################################################################
# 3.2 Preserve source layers BEFORE joining
###############################################################################

source_layers <- Layers(
  myeloid_obj[["RNA"]]
)


message("")
message("RNA layers before joining:")
print(source_layers)


write.csv(
  data.frame(
    Layer = source_layers
  ),
  file = file.path(
    table_dir,
    "phase5_myeloid_source_layers.csv"
  ),
  row.names = FALSE
)


saveRDS(
  myeloid_obj,
  source_layer_output
)


message(
  "Preserved source-layer object:\n",
  source_layer_output
)


###############################################################################
# 3.3 Donor × clinical-state representation
###############################################################################

myeloid_donor_state_counts <- myeloid_obj@meta.data %>%
  
  dplyr::count(
    Donor,
    GSM,
    Phase,
    name = "Myeloid_Cells"
  ) %>%
  
  dplyr::arrange(
    match(
      Phase,
      expected_phases
    ),
    Donor
  )


write.csv(
  myeloid_donor_state_counts,
  file = file.path(
    table_dir,
    "phase5_myeloid_donor_state_counts.csv"
  ),
  row.names = FALSE
)


myeloid_phase_counts <- myeloid_obj@meta.data %>%
  
  dplyr::count(
    Phase,
    name = "Myeloid_Cells"
  ) %>%
  
  dplyr::mutate(
    Percent_of_All_Myeloid =
      100 *
      Myeloid_Cells /
      sum(Myeloid_Cells)
  )


write.csv(
  myeloid_phase_counts,
  file = file.path(
    table_dir,
    "phase5_myeloid_phase_counts.csv"
  ),
  row.names = FALSE
)


message("")
message("Myeloid cells by clinical state:")
print(myeloid_phase_counts)


###############################################################################
# SECTION 4 — SUBSET QC CHARACTERIZATION
###############################################################################

message("")
message("============================================================")
message("SECTION 4 — MYELOID QC CHARACTERIZATION")
message("============================================================")


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
  
  colnames(
    myeloid_obj@meta.data
  )
  
)


message("Available Phase 3 QC metadata:")
print(available_qc)


###############################################################################
# 4.1 Genes detected
###############################################################################

if (
  "Genes_Detected" %in%
  available_qc
) {
  
  myeloid_genes_summary <-
    myeloid_obj@meta.data %>%
    
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
    myeloid_genes_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_Genes_Detected_summary.csv"
    ),
    row.names = FALSE
  )
  
}


###############################################################################
# 4.2 Mitochondrial detection
###############################################################################

if (
  "Mito_Detection_Fraction" %in%
  available_qc
) {
  
  myeloid_mito_summary <-
    myeloid_obj@meta.data %>%
    
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
    myeloid_mito_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_mito_summary.csv"
    ),
    row.names = FALSE
  )
  
}


###############################################################################
# 4.3 Ribosomal detection
###############################################################################

if (
  "Ribo_Detection_Fraction" %in%
  available_qc
) {
  
  myeloid_ribo_summary <-
    myeloid_obj@meta.data %>%
    
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
    myeloid_ribo_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_ribo_summary.csv"
    ),
    row.names = FALSE
  )
  
}


###############################################################################
# 4.4 Existing processed-expression sum
###############################################################################

if (
  "Total_LogCP10K" %in%
  available_qc
) {
  
  myeloid_expression_summary <-
    myeloid_obj@meta.data %>%
    
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
    myeloid_expression_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_Total_LogCP10K_summary.csv"
    ),
    row.names = FALSE
  )
  
}


###############################################################################
# 4.5 QC flags — descriptive only
###############################################################################

if (
  "QC_Pass" %in%
  colnames(myeloid_obj@meta.data)
) {
  
  qc_pass_summary <- as.data.frame(
    table(
      myeloid_obj$QC_Pass,
      useNA = "ifany"
    )
  )
  
  
  colnames(qc_pass_summary) <- c(
    "QC_Pass",
    "Cells"
  )
  
  
  write.csv(
    qc_pass_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_pass_summary.csv"
    ),
    row.names = FALSE
  )
  
}


if (
  "QC_LowGenes" %in%
  colnames(myeloid_obj@meta.data)
) {
  
  low_gene_summary <- as.data.frame(
    table(
      myeloid_obj$QC_LowGenes,
      useNA = "ifany"
    )
  )
  
  
  colnames(low_gene_summary) <- c(
    "QC_LowGenes",
    "Cells"
  )
  
  
  write.csv(
    low_gene_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_qc_low_gene_summary.csv"
    ),
    row.names = FALSE
  )
  
}


if (
  "Potential_HighComplexity" %in%
  colnames(myeloid_obj@meta.data)
) {
  
  complexity_summary <- as.data.frame(
    table(
      myeloid_obj$Potential_HighComplexity,
      useNA = "ifany"
    )
  )
  
  
  colnames(complexity_summary) <- c(
    "Potential_HighComplexity",
    "Cells"
  )
  
  
  write.csv(
    complexity_summary,
    file = file.path(
      table_dir,
      "phase5_myeloid_high_complexity_summary.csv"
    ),
    row.names = FALSE
  )
  
}


message("")
message(
  "No automatic Myeloid cell removal based on Phase 3 QC flags."
)


###############################################################################
# SECTION 5 — JOIN EXISTING NORMALIZED EXPRESSION LAYERS
###############################################################################

message("")
message("============================================================")
message("SECTION 5 — JOIN EXISTING PROCESSED EXPRESSION")
message("============================================================")


message(
  "Joining existing normalized expression layers..."
)


myeloid_obj[["RNA"]] <- JoinLayers(
  myeloid_obj[["RNA"]]
)


joined_layers <- Layers(
  myeloid_obj[["RNA"]]
)


message("RNA layers after joining:")
print(joined_layers)


if (
  !("data" %in% joined_layers)
) {
  
  stop(
    "Joined RNA assay does not contain a 'data' layer.",
    call. = FALSE
  )
  
}


joined_data <- LayerData(
  myeloid_obj,
  assay = "RNA",
  layer = "data"
)


if (
  nrow(joined_data) !=
  expected_total_features ||
  ncol(joined_data) !=
  expected_myeloid_n
) {
  
  stop(
    "Joined data layer has unexpected dimensions.\n",
    "Expected: ",
    expected_total_features,
    " x ",
    expected_myeloid_n,
    "\nFound: ",
    nrow(joined_data),
    " x ",
    ncol(joined_data),
    call. = FALSE
  )
  
}


if (
  !identical(
    colnames(joined_data),
    colnames(myeloid_obj)
  )
) {
  
  stop(
    "Joined data-layer cell names do not match Seurat object.",
    call. = FALSE
  )
  
}


message(
  "Joined data layer validated: ",
  nrow(joined_data),
  " genes × ",
  ncol(joined_data),
  " cells."
)


rm(joined_data)


###############################################################################
# SECTION 6 — MYELOID VARIABLE FEATURES + PCA
###############################################################################

message("")
message("============================================================")
message("SECTION 6 — VARIABLE FEATURES + PCA")
message("============================================================")


###############################################################################
# 6.1 Variable features
###############################################################################

message(
  "Selecting ",
  expected_hvgs,
  " myeloid variable features..."
)


myeloid_obj <- FindVariableFeatures(
  
  myeloid_obj,
  
  assay = "RNA",
  
  layer = "data",
  
  selection.method = "vst",
  
  nfeatures = expected_hvgs,
  
  verbose = TRUE
  
)


myeloid_hvgs <- VariableFeatures(
  myeloid_obj
)


if (
  length(myeloid_hvgs) !=
  expected_hvgs
) {
  
  warning(
    "Expected ",
    expected_hvgs,
    " HVGs but found ",
    length(myeloid_hvgs),
    "."
  )
  
}


write.csv(
  
  data.frame(
    Rank = seq_along(myeloid_hvgs),
    Gene = myeloid_hvgs
  ),
  
  file = file.path(
    table_dir,
    "phase5_myeloid_hvg_list.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 6.2 Variable-feature plot
###############################################################################

p_hvg <- VariableFeaturePlot(
  myeloid_obj
) +
  
  ggtitle(
    "Myeloid variable features"
  )


ggsave(
  
  file.path(
    figure_dir,
    "phase5_myeloid_variable_features.png"
  ),
  
  p_hvg,
  
  width = 8,
  height = 6,
  dpi = 300
  
)


###############################################################################
# 6.3 Scale HVGs
###############################################################################

message(
  "Scaling myeloid HVGs..."
)


myeloid_obj <- ScaleData(
  
  myeloid_obj,
  
  assay = "RNA",
  
  features = myeloid_hvgs,
  
  verbose = TRUE
  
)


###############################################################################
# 6.4 PCA
###############################################################################

message(
  "Running ",
  n_pcs,
  "-component PCA..."
)


myeloid_obj <- RunPCA(
  
  myeloid_obj,
  
  assay = "RNA",
  
  features = myeloid_hvgs,
  
  npcs = n_pcs,
  
  verbose = TRUE
  
)


if (
  !("pca" %in% Reductions(myeloid_obj))
) {
  
  stop(
    "PCA reduction was not generated.",
    call. = FALSE
  )
  
}


###############################################################################
# 6.5 PCA elbow plot
###############################################################################

p_elbow <- ElbowPlot(
  
  myeloid_obj,
  
  ndims = n_pcs
  
) +
  
  ggtitle(
    "Myeloid PCA elbow plot"
  )


ggsave(
  
  file.path(
    figure_dir,
    "phase5_myeloid_pca_elbow.png"
  ),
  
  p_elbow,
  
  width = 8,
  height = 6,
  dpi = 300
  
)


###############################################################################
# SECTION 7 — NEIGHBOR GRAPH + CLUSTERING
###############################################################################

message("")
message("============================================================")
message("SECTION 7 — NEIGHBOR GRAPH + CLUSTERING")
message("============================================================")


###############################################################################
# 7.1 Neighbor graph
###############################################################################

message(
  "Constructing myeloid neighbor graph..."
)


myeloid_obj <- FindNeighbors(
  
  myeloid_obj,
  
  assay = "RNA",
  
  reduction = "pca",
  
  dims = analysis_dims,
  
  k.param = neighbor_k,
  
  graph.name = c(
    "Myeloid_nn",
    "Myeloid_snn"
  ),
  
  verbose = TRUE
  
)


###############################################################################
# 7.2 Validate graphs
###############################################################################

required_graphs <- c(
  "Myeloid_nn",
  "Myeloid_snn"
)


available_graphs <- names(
  myeloid_obj@graphs
)


missing_graphs <- setdiff(
  required_graphs,
  available_graphs
)


if (
  length(missing_graphs) > 0L
) {
  
  stop(
    "Required graphs missing: ",
    paste(
      missing_graphs,
      collapse = ", "
    ),
    call. = FALSE
  )
  
}


###############################################################################
# 7.3 Clustering
###############################################################################

message(
  "Clustering myeloid cells..."
)


myeloid_obj <- FindClusters(
  
  myeloid_obj,
  
  graph.name = "Myeloid_snn",
  
  resolution = cluster_resolution,
  
  algorithm = 1,
  
  random.seed = analysis_seed,
  
  verbose = TRUE
  
)


###############################################################################
# 7.4 Explicit cluster metadata
###############################################################################

cluster_vector <- setNames(
  
  as.character(
    Idents(myeloid_obj)
  ),
  
  names(
    Idents(myeloid_obj)
  )
  
)


myeloid_obj <- AddMetaData(
  
  myeloid_obj,
  
  metadata = cluster_vector,
  
  col.name = "Myeloid_Cluster"
  
)


cluster_sizes <- myeloid_obj@meta.data %>%
  
  dplyr::count(
    Myeloid_Cluster,
    name = "Cells"
  ) %>%
  
  dplyr::arrange(
    suppressWarnings(
      as.numeric(
        Myeloid_Cluster
      )
    )
  )


write.csv(
  
  cluster_sizes,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_cluster_sizes.csv"
  ),
  
  row.names = FALSE
  
)


message("")
message("Myeloid cluster sizes:")
print(cluster_sizes)


###############################################################################
# SECTION 8 — UMAP
###############################################################################

message("")
message("============================================================")
message("SECTION 8 — UMAP")
message("============================================================")


myeloid_obj <- RunUMAP(
  
  myeloid_obj,
  
  reduction = "pca",
  
  dims = analysis_dims,
  
  n.neighbors = 30,
  
  min.dist = 0.3,
  
  seed.use = analysis_seed,
  
  reduction.name = "umap_myeloid",
  
  reduction.key = "MYELOIDUMAP_",
  
  verbose = TRUE
  
)


if (
  !(
    "umap_myeloid" %in%
    Reductions(myeloid_obj)
  )
) {
  
  stop(
    "Myeloid UMAP reduction was not generated.",
    call. = FALSE
  )
  
}


umap_embeddings <- Embeddings(
  myeloid_obj,
  reduction = "umap_myeloid"
)


if (
  nrow(umap_embeddings) !=
  ncol(myeloid_obj)
) {
  
  stop(
    "UMAP cell count does not match myeloid object.",
    call. = FALSE
  )
  
}


if (
  ncol(umap_embeddings) !=
  2L
) {
  
  stop(
    "Expected two-dimensional UMAP.",
    call. = FALSE
  )
  
}


if (
  !identical(
    rownames(umap_embeddings),
    colnames(myeloid_obj)
  )
) {
  
  stop(
    "UMAP cell order does not match myeloid object.",
    call. = FALSE
  )
  
}


if (
  any(
    !is.finite(
      umap_embeddings
    )
  )
) {
  
  stop(
    "UMAP contains non-finite values.",
    call. = FALSE
  )
  
}


###############################################################################
# 8.1 UMAP — clusters
###############################################################################

p_cluster <- DimPlot(
  
  myeloid_obj,
  
  reduction = "umap_myeloid",
  
  group.by = "Myeloid_Cluster",
  
  label = TRUE,
  
  repel = TRUE,
  
  raster = TRUE
  
) +
  
  ggtitle(
    "Myeloid transcriptional clusters"
  ) +
  
  theme_classic()


ggsave(
  
  file.path(
    figure_dir,
    "phase5_myeloid_umap_clusters.png"
  ),
  
  p_cluster,
  
  width = 9,
  height = 7,
  dpi = 300
  
)


###############################################################################
# 8.2 UMAP — clinical state
###############################################################################

p_phase <- DimPlot(
  
  myeloid_obj,
  
  reduction = "umap_myeloid",
  
  group.by = "Phase",
  
  raster = TRUE
  
) +
  
  ggtitle(
    "Myeloid cells by clinical state"
  ) +
  
  theme_classic()


ggsave(
  
  file.path(
    figure_dir,
    "phase5_myeloid_umap_clinical_state.png"
  ),
  
  p_phase,
  
  width = 9,
  height = 7,
  dpi = 300
  
)


###############################################################################
# SECTION 9 — DONOR REPRESENTATION OF CLUSTERS
###############################################################################

message("")
message("============================================================")
message("SECTION 9 — DONOR REPRESENTATION")
message("============================================================")


cluster_donor_counts <- myeloid_obj@meta.data %>%
  
  dplyr::count(
    
    Myeloid_Cluster,
    
    Donor,
    
    GSM,
    
    Phase,
    
    name = "Cells"
    
  )


write.csv(
  
  cluster_donor_counts,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_cluster_donor_counts.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 9.1 Number of donors represented in each cluster
###############################################################################

cluster_donor_representation <- cluster_donor_counts %>%
  
  dplyr::group_by(
    Myeloid_Cluster
  ) %>%
  
  dplyr::summarise(
    
    Donors_Represented =
      dplyr::n_distinct(
        Donor
      ),
    
    Samples_Represented =
      dplyr::n_distinct(
        GSM
      ),
    
    Clinical_States_Represented =
      dplyr::n_distinct(
        Phase
      ),
    
    Total_Cells =
      sum(
        Cells
      ),
    
    Largest_Donor_Cells =
      max(
        Cells
      ),
    
    Largest_Donor_Fraction =
      max(
        Cells
      ) /
      sum(
        Cells
      ),
    
    .groups = "drop"
    
  )


write.csv(
  
  cluster_donor_representation,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_cluster_donor_representation.csv"
  ),
  
  row.names = FALSE
  
)


message("")
message("Cluster donor representation:")
print(cluster_donor_representation)


###############################################################################
# SECTION 10 — MARKER DISCOVERY
###############################################################################

message("")
message("============================================================")
message("SECTION 10 — MARKER DISCOVERY")
message("============================================================")


###############################################################################
# 10.1 All positive markers
###############################################################################

message(
  "Running marker discovery..."
)


myeloid_markers <- FindAllMarkers(
  
  myeloid_obj,
  
  assay = "RNA",
  
  only.pos = TRUE,
  
  min.pct = 0.10,
  
  logfc.threshold = 0.25,
  
  test.use = "wilcox",
  
  verbose = TRUE
  
)


if (
  nrow(myeloid_markers) == 0L
) {
  
  stop(
    "No positive markers were identified.",
    call. = FALSE
  )
  
}


write.csv(
  
  myeloid_markers,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_all_positive_markers.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 10.2 Significant markers
###############################################################################

marker_p_column <- if (
  "p_val_adj" %in%
  colnames(myeloid_markers)
) {
  
  "p_val_adj"
  
} else {
  
  "p_val"
  
}


myeloid_markers_sig <- myeloid_markers %>%
  
  dplyr::filter(
    
    .data[[marker_p_column]] <
      0.05
    
  )


write.csv(
  
  myeloid_markers_sig,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_significant_markers.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 10.3 Top 20 markers per cluster
###############################################################################

myeloid_top20 <- myeloid_markers %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::slice_max(
    
    order_by = avg_log2FC,
    
    n = 20,
    
    with_ties = FALSE
    
  ) %>%
  
  dplyr::ungroup()


write.csv(
  
  myeloid_top20,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_top20_markers.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 10.4 Top 50 markers per cluster
###############################################################################

myeloid_top50 <- myeloid_markers %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::slice_max(
    
    order_by = avg_log2FC,
    
    n = 50,
    
    with_ties = FALSE
    
  ) %>%
  
  dplyr::ungroup()


write.csv(
  
  myeloid_top50,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_top50_markers.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# SECTION 11 — MYELOID BIOLOGICAL PROGRAMS
###############################################################################

message("")
message("============================================================")
message("SECTION 11 — MYELOID BIOLOGICAL PROGRAMS")
message("============================================================")


###############################################################################
# 11.1 Program definitions
#
# These are descriptive transcriptional programs.
# They are NOT automatic biological-state labels.
###############################################################################

myeloid_programs <- list(
  
  Monocyte_Inflammatory = c(
    
    "S100A8",
    "S100A9",
    "FCN1",
    "VCAN",
    "CTSS",
    "CTSD",
    "FCGR3A",
    "TYMP",
    "LGALS3",
    "CTSD"
    
  ),
  
  
  Macrophage_Associated = c(
    
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "TREM2",
    "CD68",
    "MSR1",
    "MARCO",
    "LPL",
    "FOLR2"
    
  ),
  
  
  Antigen_Presentation = c(
    
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-DPA1",
    "HLA-DPB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "CD74",
    "CIITA"
    
  ),
  
  
  Interferon_Response = c(
    
    "ISG15",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFI6",
    "MX1",
    "OAS1",
    "OAS2",
    "OAS3",
    "IFITM3"
    
  ),
  
  
  Phagolysosomal = c(
    
    "CTSD",
    "CTSS",
    "CTSB",
    "LAMP1",
    "LAMP2",
    "FCER1G",
    "TYROBP",
    "AIF1",
    "LGALS3"
    
  ),
  
  
  Complement_Lipid_Associated = c(
    
    "C3",
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "APOC1",
    "LPL",
    "FABP5",
    "GPNMB"
    
  ),
  
  
  Activation_Immediate_Early = c(
    
    "FOS",
    "JUN",
    "JUNB",
    "FOSB",
    "DUSP1",
    "DUSP2",
    "EGR1",
    "EGR2",
    "NR4A1",
    "NR4A2"
    
  ),
  
  
  Proliferation = c(
    
    "MKI67",
    "TOP2A",
    "TYMS",
    "STMN1",
    "TUBA1B",
    "UBE2C",
    "HMGB2",
    "CENPF"
    
  )
  
)


###############################################################################
# 11.2 Keep only genes actually present
###############################################################################

all_features <- rownames(
  myeloid_obj
)


myeloid_programs_present <- lapply(
  
  myeloid_programs,
  
  function(x) {
    
    intersect(
      x,
      all_features
    )
    
  }
  
)


myeloid_programs_present <-
  myeloid_programs_present[
    lengths(
      myeloid_programs_present
    ) > 0
  ]


message("")
message("Myeloid program coverage:")


for (
  program_name in
  names(myeloid_programs_present)
) {
  
  message(
    
    program_name,
    ": ",
    
    length(
      myeloid_programs_present[[program_name]]
    ),
    
    " genes present"
    
  )
  
}


###############################################################################
# 11.3 Program scoring
###############################################################################

if (
  length(myeloid_programs_present) > 0L
) {
  
  for (
    program_name in
    names(myeloid_programs_present)
  ) {
    
    features <- myeloid_programs_present[[program_name]]
    
    
    if (
      length(features) >= 2L
    ) {
      
      score_name <- paste0(
        "Program_",
        program_name
      )
      
      
      myeloid_obj <- AddModuleScore(
        
        myeloid_obj,
        
        features = list(
          features
        ),
        
        name = score_name,
        
        assay = "RNA"
        
      )
      
      
      generated_score <- paste0(
        score_name,
        "1"
      )
      
      
      if (
        generated_score %in%
        colnames(
          myeloid_obj@meta.data
        )
      ) {
        
        colnames(
          myeloid_obj@meta.data
        )[
          colnames(
            myeloid_obj@meta.data
          ) ==
            generated_score
        ] <- score_name
        
      }
      
    }
    
  }
  
}


###############################################################################
# 11.4 Cluster-level program summaries
###############################################################################

program_columns <- grep(
  
  "^Program_",
  
  colnames(
    myeloid_obj@meta.data
  ),
  
  value = TRUE
  
)


if (
  length(program_columns) > 0L
) {
  
  cluster_program_summary <-
    myeloid_obj@meta.data %>%
    
    dplyr::group_by(
      Myeloid_Cluster
    ) %>%
    
    dplyr::summarise(
      
      dplyr::across(
        all_of(
          program_columns
        ),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
      
    )
  
  
  write.csv(
    
    cluster_program_summary,
    
    file = file.path(
      table_dir,
      "phase5_myeloid_cluster_program_scores.csv"
    ),
    
    row.names = FALSE
    
  )
  
}


###############################################################################
# SECTION 12 — LINEAGE / IDENTITY EVIDENCE
###############################################################################

message("")
message("============================================================")
message("SECTION 12 — LINEAGE / IDENTITY EVIDENCE")
message("============================================================")


###############################################################################
# 12.1 Marker panels
#
# These panels are evidence streams, not classifiers.
###############################################################################

marker_panels <- list(
  
  Core_Myeloid = c(
    
    "LYZ",
    "TYROBP",
    "FCER1G",
    "CTSS",
    "CTSD",
    "AIF1",
    "CST3"
    
  ),
  
  
  Monocyte_Associated = c(
    
    "S100A8",
    "S100A9",
    "FCN1",
    "CD14",
    "VCAN",
    "TREM1",
    "FCAR",
    "FPR2",
    "OLR1"
    
  ),
  
  
  Macrophage_Associated = c(
    
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "TREM2",
    "CD68",
    "MSR1",
    "MARCO",
    "LPL",
    "FOLR2"
    
  ),
  
  
  Antigen_Presenting_Myeloid = c(
    
    "HLA-DRA",
    "HLA-DRB1",
    "HLA-DPA1",
    "HLA-DPB1",
    "CD74",
    "CIITA"
    
  ),
  
  
  Conventional_DC_Associated = c(
    
    "FCER1A",
    "CD1C",
    "CLEC10A",
    "CST3",
    "HLA-DRA"
    
  ),
  
  
  pDC_Associated = c(
    
    "GZMB",
    "GZMB",
    "JCHAIN",
    "GZMB",
    "IRF7",
    "GZMB",
    "TCF4",
    "CLEC4C",
    "GZMB"
    
  ),
  
  
  Neutrophil_Associated = c(
    
    "FCGR3B",
    "CSF3R",
    "S100A8",
    "S100A9",
    "FPR1",
    "MNDA"
    
  ),
  
  
  Lymphoid_Contamination = c(
    
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2",
    "CD8A",
    "CD8B",
    "MS4A1",
    "CD79A",
    "NKG7",
    "GNLY"
    
  )
  
)


###############################################################################
# 12.2 Available panel genes
###############################################################################

marker_panels_present <- lapply(
  
  marker_panels,
  
  function(x) {
    
    intersect(
      x,
      all_features
    )
    
  }
  
)


marker_panels_present <-
  marker_panels_present[
    lengths(
      marker_panels_present
    ) > 0
  ]


message("")
message("Marker-panel coverage:")


for (
  panel_name in
  names(marker_panels_present)
) {
  
  message(
    
    panel_name,
    ": ",
    
    length(
      marker_panels_present[[panel_name]]
    ),
    
    " genes present"
    
  )
  
}


###############################################################################
# 12.3 Marker panel table
###############################################################################

marker_panel_table <- dplyr::bind_rows(
  
  lapply(
    
    names(marker_panels_present),
    
    function(panel_name) {
      
      data.frame(
        
        Panel =
          panel_name,
        
        Gene =
          marker_panels_present[[panel_name]],
        
        stringsAsFactors =
          FALSE
        
      )
      
    }
    
  )
  
)


write.csv(
  
  marker_panel_table,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_marker_validation_panels.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 12.4 Marker validation DotPlot
###############################################################################

marker_features <- unique(
  unlist(
    marker_panels_present,
    use.names = FALSE
  )
)


if (
  length(marker_features) > 0L
) {
  
  p_marker_dotplot <- DotPlot(
    
    myeloid_obj,
    
    assay = "RNA",
    
    features = marker_features,
    
    group.by = "Myeloid_Cluster",
    
    dot.scale = 6
    
  ) +
    
    ggtitle(
      "Myeloid marker validation"
    ) +
    
    theme_classic() +
    
    theme(
      
      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )
      
    )
  
  
  ggsave(
    
    file.path(
      figure_dir,
      "phase5_myeloid_marker_validation_dotplot.png"
    ),
    
    p_marker_dotplot,
    
    width = 16,
    height = 9,
    dpi = 300
    
  )
  
}


###############################################################################
# 12.5 Cluster-level panel expression
###############################################################################

if (
  length(marker_features) > 0L
) {
  
  panel_expression <- FetchData(
    
    myeloid_obj,
    
    vars = marker_features
    
  )
  
  
  panel_expression$Myeloid_Cluster <-
    myeloid_obj$Myeloid_Cluster
  
  
  panel_cluster_means <-
    panel_expression %>%
    
    dplyr::group_by(
      Myeloid_Cluster
    ) %>%
    
    dplyr::summarise(
      
      dplyr::across(
        all_of(
          marker_features
        ),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
      
    )
  
  
  write.csv(
    
    panel_cluster_means,
    
    file = file.path(
      table_dir,
      "phase5_myeloid_marker_panel_cluster_means.csv"
    ),
    
    row.names = FALSE
    
  )
  
}


###############################################################################
# SECTION 13 — STRUCTURED CLUSTER EVIDENCE
###############################################################################

message("")
message("============================================================")
message("SECTION 13 — STRUCTURED CLUSTER EVIDENCE")
message("============================================================")


###############################################################################
# 13.1 Basic cluster statistics
###############################################################################

cluster_basic_summary <- myeloid_obj@meta.data %>%
  
  dplyr::group_by(
    Myeloid_Cluster
  ) %>%
  
  dplyr::summarise(
    
    Cells =
      dplyr::n(),
    
    Donors =
      dplyr::n_distinct(
        Donor
      ),
    
    Samples =
      dplyr::n_distinct(
        GSM
      ),
    
    Clinical_States =
      dplyr::n_distinct(
        Phase
      ),
    
    Largest_Donor_Cells =
      max(
        table(
          Donor
        )
      ),
    
    Largest_Donor_Fraction =
      max(
        table(
          Donor
        )
      ) /
      dplyr::n(),
    
    .groups = "drop"
    
  )


###############################################################################
# 13.2 Cluster marker summaries
###############################################################################

top_marker_summary <- myeloid_top20 %>%
  
  dplyr::group_by(
    cluster
  ) %>%
  
  dplyr::summarise(
    
    Top_Markers =
      paste(
        gene,
        collapse = ", "
      ),
    
    .groups = "drop"
    
  ) %>%
  
  dplyr::rename(
    Myeloid_Cluster = cluster
  )


###############################################################################
# 13.3 Program summaries
###############################################################################

if (
  exists(
    "cluster_program_summary"
  )
) {
  
  structured_evidence <- cluster_basic_summary %>%
    
    dplyr::left_join(
      top_marker_summary,
      by = "Myeloid_Cluster"
    ) %>%
    
    dplyr::left_join(
      cluster_program_summary,
      by = "Myeloid_Cluster"
    )
  
} else {
  
  structured_evidence <- cluster_basic_summary %>%
    
    dplyr::left_join(
      top_marker_summary,
      by = "Myeloid_Cluster"
    )
  
}


###############################################################################
# 13.4 Initial adjudication status
###############################################################################

structured_evidence <-
  structured_evidence %>%
  
  dplyr::mutate(
    
    Biological_Annotation =
      "REVIEW_REQUIRED",
    
    Annotation_Status =
      "REVIEW_REQUIRED",
    
    Interpretation_Boundary =
      paste(
        
        "Computational cluster requiring multi-evidence biological",
        "adjudication. Do not interpret as an established cell state",
        "without marker, program, lineage, and donor-representation",
        "evidence."
        
      )
    
  )


write.csv(
  
  structured_evidence,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_structured_cluster_evidence.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# SECTION 14 — MANUAL BIOLOGICAL ADJUDICATION TEMPLATE
###############################################################################

message("")
message("============================================================")
message("SECTION 14 — BIOLOGICAL ADJUDICATION TEMPLATE")
message("============================================================")


###############################################################################
# 14.1 Create one row per computational cluster
#
# IMPORTANT:
# This is deliberately NOT an automatic classifier.
#
# Biological labels will be assigned only after inspection of:
#   1. top markers
#   2. marker-panel evidence
#   3. program scores
#   4. donor representation
#   5. potential contaminating lineage signals
#
# This prevents cluster-number-based assumptions.
###############################################################################

annotation_template <- cluster_basic_summary %>%
  
  dplyr::left_join(
    top_marker_summary,
    by = "Myeloid_Cluster"
  ) %>%
  
  dplyr::mutate(
    
    Biological_Annotation =
      "REVIEW_REQUIRED",
    
    Lineage_Category =
      "REVIEW_REQUIRED",
    
    Annotation_Confidence =
      "REVIEW_REQUIRED",
    
    Donor_Representation_Interpretation =
      dplyr::case_when(
        
        Largest_Donor_Fraction >= 0.90 ~
          "Strongly donor concentrated",
        
        Largest_Donor_Fraction >= 0.75 ~
          "Donor concentrated",
        
        TRUE ~
          "Broad donor representation"
        
      ),
    
    Interpretation_Notes =
      "",
    
    Evidence_Summary =
      ""
    
  )


write.csv(
  
  annotation_template,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_manual_annotation_template.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 14.2 Save an explicit adjudication framework
###############################################################################

adjudication_framework <- tibble(
  
  Evidence_Stream = c(
    
    "Top differential markers",
    
    "Core myeloid identity",
    
    "Monocyte-associated program",
    
    "Macrophage-associated program",
    
    "Antigen-presentation program",
    
    "Interferon-response program",
    
    "Phagolysosomal program",
    
    "Complement/lipid-associated program",
    
    "Immediate-early activation program",
    
    "Non-myeloid contamination",
    
    "Donor representation",
    
    "Clinical-state distribution"
    
  ),
  
  
  Role = c(
    
    "Identify genes enriched within the computational cluster.",
    
    "Confirm myeloid identity rather than relying on cluster position.",
    
    "Evaluate inflammatory/classical-monocyte-associated transcriptional evidence.",
    
    "Evaluate macrophage-associated transcriptional evidence.",
    
    "Evaluate antigen-presentation-associated transcriptional evidence.",
    
    "Evaluate interferon-stimulated transcriptional evidence.",
    
    "Evaluate lysosomal/phagocytic machinery.",
    
    "Evaluate complement/lipid-handling-associated programs.",
    
    "Identify immediate-early transcriptional activation.",
    
    "Check whether apparent identity may instead reflect lymphoid/DC/neutrophil-associated signals.",
    
    "Determine whether a cluster is broadly represented across donors or dominated by one donor.",
    
    "Describe distribution only; no clinical-state inference is made in Phase 5."
    
  ),
  
  
  Interpretation_Rule = c(
    
    "Descriptive evidence; not sufficient alone for a biological-state label.",
    
    "Required for confident myeloid interpretation.",
    
    "Supports monocyte-associated interpretation when coherent with other evidence.",
    
    "Supports macrophage-associated interpretation when coherent with other evidence.",
    
    "Describes an antigen-presentation program and does not by itself define lineage.",
    
    "Describes an interferon-response program and does not by itself define lineage.",
    
    "Describes cellular machinery rather than an automatically named state.",
    
    "Describes metabolic/complement-associated biology rather than a predetermined phenotype.",
    
    "Treat as a transcriptional program, not necessarily a stable lineage/state.",
    
    "Conflicting lineage evidence requires conservative annotation.",
    
    "Strong donor concentration limits cohort-level generalization.",
    
    "Clinical-state comparisons are reserved for donor-aware Phase 6 analysis."
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  adjudication_framework,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_adjudication_framework.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# SECTION 15 — DONOR-LEVEL DESCRIPTIVE CLUSTER COMPOSITION
###############################################################################

message("")
message("============================================================")
message("SECTION 15 — DONOR-LEVEL DESCRIPTIVE COMPOSITION")
message("============================================================")


###############################################################################
# Calculate donor-level cluster proportions
###############################################################################

donor_cluster_counts <- myeloid_obj@meta.data %>%
  
  dplyr::count(
    
    Donor,
    
    GSM,
    
    Phase,
    
    Myeloid_Cluster,
    
    name = "Cells"
    
  )


donor_totals <- donor_cluster_counts %>%
  
  dplyr::group_by(
    Donor,
    GSM,
    Phase
  ) %>%
  
  dplyr::summarise(
    
    Total_Myeloid_Cells =
      sum(
        Cells
      ),
    
    .groups = "drop"
    
  )


donor_cluster_proportions <-
  donor_cluster_counts %>%
  
  dplyr::left_join(
    
    donor_totals,
    
    by = c(
      "Donor",
      "GSM",
      "Phase"
    )
    
  ) %>%
  
  dplyr::mutate(
    
    Fraction_of_Donor_Myeloid =
      Cells /
      Total_Myeloid_Cells,
    
    Percent_of_Donor_Myeloid =
      100 *
      Fraction_of_Donor_Myeloid
    
  )


write.csv(
  
  donor_cluster_proportions,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_donor_cluster_proportions.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# SECTION 16 — UMAP PROGRAM VISUALIZATION
###############################################################################

message("")
message("============================================================")
message("SECTION 16 — PROGRAM VISUALIZATION")
message("============================================================")


###############################################################################
# Plot available programs
###############################################################################

if (
  length(program_columns) > 0L
) {
  
  program_features_for_plot <- program_columns
  
  
  for (
    program_name in
    program_features_for_plot
  ) {
    
    p_program <- FeaturePlot(
      
      myeloid_obj,
      
      reduction = "umap_myeloid",
      
      features = program_name,
      
      raster = TRUE
      
    ) +
      
      ggtitle(
        program_name
      ) +
      
      theme_classic()
    
    
    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      program_name
    )
    
    
    ggsave(
      
      file.path(
        figure_dir,
        paste0(
          "phase5_myeloid_",
          safe_name,
          "_umap.png"
        )
      ),
      
      p_program,
      
      width = 8,
      height = 7,
      dpi = 300
      
    )
    
  }
  
}


###############################################################################
# SECTION 17 — SAVE INTERMEDIATE CHECKPOINT
###############################################################################

message("")
message("============================================================")
message("SECTION 17 — CHECKPOINT SAVE")
message("============================================================")


checkpoint_output <- file.path(
  
  rds_dir,
  
  "phase5_myeloid_analysis_checkpoint.rds"
  
)


saveRDS(
  
  myeloid_obj,
  
  checkpoint_output
  
)


message(
  "Checkpoint saved:\n",
  checkpoint_output
)


###############################################################################
# SECTION 18 — FINAL BIOLOGICAL ANNOTATION OBJECT
#
# At this stage the computational analysis is complete.
# Biological annotations remain REVIEW_REQUIRED until the actual
# cluster-specific evidence is inspected.
###############################################################################

message("")
message("============================================================")
message("SECTION 18 — FINAL BIOLOGICAL ANNOTATION OBJECT")
message("============================================================")


###############################################################################
# 18.1 Attach structured evidence metadata
###############################################################################

annotation_lookup <- annotation_template %>%
  
  dplyr::select(
    
    Myeloid_Cluster,
    
    Biological_Annotation,
    
    Lineage_Category,
    
    Annotation_Confidence,
    
    Donor_Representation_Interpretation,
    
    Interpretation_Notes,
    
    Evidence_Summary
    
  )


annotation_lookup_vector <- setNames(
  
  annotation_lookup$Biological_Annotation,
  
  annotation_lookup$Myeloid_Cluster
  
)


myeloid_biological_annotation <-
  annotation_lookup_vector[
    as.character(
      myeloid_obj$Myeloid_Cluster
    )
  ]


names(
  myeloid_biological_annotation
) <-
  colnames(
    myeloid_obj
  )


myeloid_obj <- AddMetaData(
  
  myeloid_obj,
  
  metadata =
    myeloid_biological_annotation,
  
  col.name =
    "Myeloid_Biological_Annotation"
  
)


###############################################################################
# 18.2 Attach lineage category
###############################################################################

lineage_lookup_vector <- setNames(
  
  annotation_lookup$Lineage_Category,
  
  annotation_lookup$Myeloid_Cluster
  
)


myeloid_lineage_category <-
  lineage_lookup_vector[
    as.character(
      myeloid_obj$Myeloid_Cluster
    )
  ]


names(
  myeloid_lineage_category
) <-
  colnames(
    myeloid_obj
  )


myeloid_obj <- AddMetaData(
  
  myeloid_obj,
  
  metadata =
    myeloid_lineage_category,
  
  col.name =
    "Myeloid_Lineage_Category"
  
)


###############################################################################
# 18.3 Metadata documenting analysis boundaries
###############################################################################

myeloid_obj$Phase5_Analysis <-
  "Myeloid_state_discovery_and_biological_adjudication"


myeloid_obj$Phase5_Expression <-
  "Processed_log_CP10K"


myeloid_obj$Phase5_Donor_Unit <-
  "Biological_replicate"


myeloid_obj$Phase5_Inference <-
  "Exploratory_cell_state_characterization"


myeloid_obj$Phase5_Clinical_State_Inference <-
  "Reserved_for_Phase6"


###############################################################################
# SECTION 19 — FINAL VALIDATION
#
# This is intentionally the ONLY final status/closeout section.
###############################################################################

message("")
message("============================================================")
message("SECTION 19 — FINAL VALIDATION")
message("============================================================")


###############################################################################
# 19.1 Cell identity and count
###############################################################################

if (
  ncol(myeloid_obj) !=
  expected_myeloid_cells
) {
  
  stop(
    "Final myeloid cell count changed unexpectedly.",
    call. = FALSE
  )
  
}


if (
  nrow(myeloid_obj) !=
  expected_total_features
) {
  
  stop(
    "Final myeloid feature count changed unexpectedly.",
    call. = FALSE
  )
  
}


if (
  !all(
    myeloid_obj$Final_Cell_Type ==
    target_cell_type
  )
) {
  
  stop(
    "Final object contains cells outside Monocyte_Myeloid.",
    call. = FALSE
  )
  
}


###############################################################################
# 19.2 Cell-ID integrity
###############################################################################

if (
  !identical(
    colnames(myeloid_obj),
    rownames(myeloid_obj@meta.data)
  )
) {
  
  stop(
    "Final cell IDs and metadata row names do not match.",
    call. = FALSE
  )
  
}


###############################################################################
# 19.3 Required metadata
###############################################################################

required_final_metadata <- c(
  
  "GSM",
  "Donor",
  "Phase",
  
  "Final_Cell_Type",
  
  "Myeloid_Cluster",
  
  "Myeloid_Biological_Annotation",
  
  "Myeloid_Lineage_Category",
  
  "Phase5_Analysis",
  
  "Phase5_Expression",
  
  "Phase5_Donor_Unit",
  
  "Phase5_Inference"
  
)


missing_final_metadata <- setdiff(
  
  required_final_metadata,
  
  colnames(
    myeloid_obj@meta.data
  )
  
)


if (
  length(missing_final_metadata) > 0L
) {
  
  stop(
    
    "Missing required final metadata:\n",
    
    paste(
      missing_final_metadata,
      collapse = ", "
    ),
    
    call. = FALSE
    
  )
  
}


###############################################################################
# 19.4 No missing annotations
###############################################################################

if (
  any(
    is.na(
      myeloid_obj$Myeloid_Biological_Annotation
    )
  )
) {
  
  stop(
    "Missing Myeloid_Biological_Annotation values.",
    call. = FALSE
  )
  
}


###############################################################################
# 19.5 Cluster annotation coverage
###############################################################################

cluster_annotation_check <-
  myeloid_obj@meta.data %>%
  
  dplyr::distinct(
    
    Myeloid_Cluster,
    
    Myeloid_Biological_Annotation
    
  )


n_computational_clusters <-
  dplyr::n_distinct(
    myeloid_obj$Myeloid_Cluster
  )


n_annotated_clusters <-
  nrow(
    cluster_annotation_check
  )


if (
  n_annotated_clusters !=
  n_computational_clusters
) {
  
  stop(
    
    "Not every computational cluster has a biological annotation row.\n",
    
    "Computational clusters: ",
    n_computational_clusters,
    
    "\nAnnotation rows: ",
    n_annotated_clusters,
    
    call. = FALSE
    
  )
  
}


message(
  
  "Biological annotation coverage validated: ",
  
  n_annotated_clusters,
  
  " computational clusters represented."
  
)


###############################################################################
# 19.6 Required reductions
###############################################################################

required_reductions <- c(
  
  "pca",
  
  "umap_myeloid"
  
)


missing_reductions <- setdiff(
  
  required_reductions,
  
  Reductions(
    myeloid_obj
  )
  
)


if (
  length(missing_reductions) > 0L
) {
  
  stop(
    
    "Missing required reductions: ",
    
    paste(
      missing_reductions,
      collapse = ", "
    ),
    
    call. = FALSE
    
  )
  
}


###############################################################################
# 19.7 PCA dimensions
###############################################################################

pca_embeddings_final <- Embeddings(
  
  myeloid_obj,
  
  reduction = "pca"
  
)


if (
  ncol(pca_embeddings_final) <
  n_pcs
) {
  
  stop(
    "Final PCA contains fewer than 50 PCs.",
    call. = FALSE
  )
  
}


###############################################################################
# 19.8 UMAP dimensions
###############################################################################

umap_embeddings_final <- Embeddings(
  
  myeloid_obj,
  
  reduction = "umap_myeloid"
  
)


if (
  ncol(umap_embeddings_final) !=
  2L
) {
  
  stop(
    "Final UMAP is not two-dimensional.",
    call. = FALSE
  )
  
}


if (
  any(
    !is.finite(
      umap_embeddings_final
    )
  )
) {
  
  stop(
    "Final UMAP contains non-finite values.",
    call. = FALSE
  )
  
}


###############################################################################
# 19.9 Seurat object validity
###############################################################################

validObject(
  
  myeloid_obj,
  
  test = TRUE
  
)


###############################################################################
# 19.10 Final structured evidence save
###############################################################################

write.csv(
  
  annotation_template,
  
  file = file.path(
    table_dir,
    "phase5_myeloid_final_annotation_adjudication_table.csv"
  ),
  
  row.names = FALSE
  
)


###############################################################################
# 19.11 Save final RDS
###############################################################################

saveRDS(
  
  myeloid_obj,
  
  final_output
  
)


###############################################################################
# 19.12 Reload final object and validate persistence
###############################################################################

myeloid_final_check <- readRDS(
  
  final_output
  
)


if (
  !inherits(
    myeloid_final_check,
    "Seurat"
  )
) {
  
  stop(
    "Reloaded final object is not a Seurat object.",
    call. = FALSE
  )
  
}


if (
  ncol(myeloid_final_check) !=
  ncol(myeloid_obj)
) {
  
  stop(
    "Reloaded final object has a different cell count.",
    call. = FALSE
  )
  
}


if (
  !identical(
    
    colnames(
      myeloid_final_check
    ),
    
    colnames(
      myeloid_obj
    )
    
  )
) {
  
  stop(
    "Reloaded final object cell order differs.",
    call. = FALSE
  )
  
}


if (
  !identical(
    
    myeloid_final_check$
    Myeloid_Biological_Annotation,
    
    myeloid_obj$
    Myeloid_Biological_Annotation
    
  )
) {
  
  stop(
    "Reloaded biological annotations differ from in-memory object.",
    call. = FALSE
  )
  
}


validObject(
  
  myeloid_final_check,
  
  test = TRUE
  
)


###############################################################################
# SECTION 20 — FINAL CLOSEOUT
###############################################################################

message("")
message("============================================================")
message("PHASE 5 FINAL CHECKPOINT")
message("============================================================")


message(
  "Status: COMPLETE"
)


message(
  "Myeloid-labelled cells: ",
  ncol(myeloid_final_check)
)


message(
  "Features: ",
  nrow(myeloid_final_check)
)


message(
  "Donors: ",
  dplyr::n_distinct(
    myeloid_final_check$Donor
  )
)


message(
  "Clinical states: ",
  paste(
    sort(
      unique(
        as.character(
          myeloid_final_check$Phase
        )
      )
    ),
    collapse = ", "
  )
)


message(
  "Computational clusters: ",
  dplyr::n_distinct(
    myeloid_final_check$Myeloid_Cluster
  )
)


message(
  "Clusters with biological annotation rows: ",
  nrow(
    cluster_annotation_check
  )
)


message(
  "Unique biological annotation categories currently assigned: ",
  dplyr::n_distinct(
    myeloid_final_check$
      Myeloid_Biological_Annotation
  )
)


message(
  "HVGs: ",
  length(
    VariableFeatures(
      myeloid_final_check,
      assay = "RNA"
    )
  )
)


message(
  "PCA: PRESENT"
)


message(
  "UMAP: PRESENT"
)


message(
  "Myeloid neighbor graph: PRESENT"
)


message(
  "Marker analysis: COMPLETE"
)


message(
  "Biological program scoring: COMPLETE"
)


message(
  "Lineage/state evidence: COMPLETE"
)


message(
  "Donor representation analysis: COMPLETE"
)


message(
  "Clinical-state inference: RESERVED FOR PHASE 6"
)


message(
  "Differential-state/expression analysis: RESERVED FOR PHASE 6"
)


message(
  "Pathway enrichment: RESERVED FOR PHASE 7"
)


message(
  "CellChat: RESERVED FOR PHASE 8"
)


message("")
message("FINAL RDS:")
message(final_output)


message("")
message("KEY TABLES:")
message(
  file.path(
    table_dir,
    "phase5_myeloid_top20_markers.csv"
  )
)


message(
  file.path(
    table_dir,
    "phase5_myeloid_cluster_program_scores.csv"
  )
)


message(
  file.path(
    table_dir,
    "phase5_myeloid_structured_cluster_evidence.csv"
  )
)


message(
  file.path(
    table_dir,
    "phase5_myeloid_manual_annotation_template.csv"
  )
)


message(
  file.path(
    table_dir,
    "phase5_myeloid_cluster_donor_representation.csv"
  )
)


message("")
message("KEY FIGURES:")
message(
  file.path(
    figure_dir,
    "phase5_myeloid_umap_clusters.png"
  )
)


message(
  file.path(
    figure_dir,
    "phase5_myeloid_umap_clinical_state.png"
  )
)


message(
  file.path(
    figure_dir,
    "phase5_myeloid_marker_validation_dotplot.png"
  )
)


message("")
message("============================================================")
message("PHASE 5 COMPUTATIONAL PIPELINE IS WRAPPED.")
message("============================================================")
message("")
message(
  "IMPORTANT: Biological annotations remain REVIEW_REQUIRED ",
  "until the actual Phase 5 cluster evidence is inspected."
)
message(
  "Do not interpret computational cluster numbers as biological ",
  "states before adjudication."
)
message("============================================================")

# ============================================================
# PHASE 5 — FINAL BIOLOGICAL ANNOTATION TABLE
# ============================================================
# Purpose:
#   Create the final cluster-level biological adjudication table
#   from the evidence generated during Phase 5.
#
# Important:
#   These are biological annotations of computational clusters.
#   They are NOT clinical-state conclusions.
#
# Donor = biological replicate.
# Cells = units for state discovery / annotation.
# ============================================================

library(dplyr)
library(tibble)
library(readr)

# ------------------------------------------------------------
# 1. Define final biological adjudication
# ------------------------------------------------------------

phase5_final_annotations <- tibble(
  
  Myeloid_Cluster = as.character(0:4),
  
  Biological_Annotation = c(
    "Inflammatory monocyte-associated population",
    "Activated/antigen-presenting monocyte-associated population",
    "Antigen-presenting myeloid population",
    "FCGR3A-associated myeloid population",
    "Platelet-associated population"
  ),
  
  Lineage_Category = c(
    "Monocyte_associated",
    "Monocyte_associated",
    "Myeloid_associated",
    "Myeloid_associated",
    "Non_myeloid_associated"
  ),
  
  Biological_Annotation_Confidence = c(
    "High",
    "High",
    "Moderate-High",
    "Moderate",
    "High"
  ),
  
  Donor_Representation_Interpretation = c(
    "Strongly donor concentrated",
    "Broad donor representation",
    "Broad donor representation",
    "Strongly donor concentrated",
    "Extremely donor concentrated"
  ),
  
  Interpretation_Eligible_As_Myeloid_State = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),
  
  Evidence_Summary = c(
    paste(
      "Strong inflammatory monocyte-associated evidence:",
      "S100A8/S100A12, CCR2, FCN1, CD14 and VCAN.",
      "Coherent core myeloid and antigen-presentation programs.",
      "No substantial lymphoid signal."
    ),
    
    paste(
      "Activated monocyte-associated population with strong",
      "HLA-II/CD74 antigen-presentation signal and markers including",
      "IL1A, SERPINB2, C15orf48, GPR183 and IL1R2."
    ),
    
    paste(
      "Strong core myeloid and HLA-II/CD74 antigen-presentation",
      "programs with comparatively weak macrophage-associated",
      "program. Evidence supports a myeloid population but does",
      "not justify forcing a specific macrophage identity."
    ),
    
    paste(
      "FCGR3A/CX3CR1/MS4A4A-associated myeloid population with",
      "strong HLA-II and core myeloid signal. Mixed granulocyte-",
      "associated signal and donor concentration warrant conservative",
      "annotation."
    ),
    
    paste(
      "Dominated by platelet-associated markers including PF4, PPBP,",
      "TUBB1, ITGA2B, GP9, TREML1, MPIG6B and PF4V1.",
      "This population should not be interpreted as a myeloid state.",
      "Because the dataset contains processed log-CP10K expression",
      "rather than raw counts, a definitive doublet call is not made."
    )
  ),
  
  Interpretation_Boundary = c(
    "Do not interpret as a clinical-state-specific inflammatory population without donor-aware Phase 6 analysis.",
    
    "Activation and antigen-presentation are transcriptional descriptors; do not infer causality, pathogenicity or clinical-state specificity.",
    
    "Do not relabel as macrophage, Kupffer cell, DC or disease-associated macrophage without stronger lineage evidence.",
    
    "Do not force a neutrophil or dendritic-cell identity despite some associated markers; donor concentration limits cohort-level interpretation.",
    
    "Treat as a platelet-associated/non-myeloid population. Do not interpret as a biological myeloid state or make a definitive doublet claim."
  )
)

# ------------------------------------------------------------
# 2. Add cluster sizes from the actual object
# ------------------------------------------------------------

cluster_sizes <- myeloid_obj@meta.data %>%
  count(Myeloid_Cluster, name = "Cells") %>%
  mutate(Myeloid_Cluster = as.character(Myeloid_Cluster))

phase5_final_annotations <- phase5_final_annotations %>%
  left_join(
    cluster_sizes,
    by = "Myeloid_Cluster"
  ) %>%
  select(
    Myeloid_Cluster,
    Cells,
    Biological_Annotation,
    Lineage_Category,
    Biological_Annotation_Confidence,
    Donor_Representation_Interpretation,
    Interpretation_Eligible_As_Myeloid_State,
    Evidence_Summary,
    Interpretation_Boundary
  )

# ------------------------------------------------------------
# 3. Validate that every computational cluster has exactly one
#    biological annotation
# ------------------------------------------------------------

stopifnot(
  nrow(phase5_final_annotations) == 5
)

stopifnot(
  setequal(
    phase5_final_annotations$Myeloid_Cluster,
    as.character(sort(unique(myeloid_obj$Myeloid_Cluster)))
  )
)

stopifnot(
  !anyDuplicated(phase5_final_annotations$Myeloid_Cluster)
)

stopifnot(
  !any(is.na(phase5_final_annotations$Biological_Annotation))
)

# ------------------------------------------------------------
# 4. Print final table
# ------------------------------------------------------------

print(phase5_final_annotations, n = Inf)

# ------------------------------------------------------------
# 5. Save final table
# ------------------------------------------------------------

write_csv(
  phase5_final_annotations,
  file.path(
    table_dir,
    "phase5_myeloid_final_biological_annotations.csv"
  )
)

cat("\n============================================================\n")
cat("PHASE 5 FINAL BIOLOGICAL ANNOTATION TABLE SAVED\n")
cat("============================================================\n")
cat("Clusters:", nrow(phase5_final_annotations), "\n")
cat(
  "Myeloid-interpretable clusters:",
  sum(phase5_final_annotations$Interpretation_Eligible_As_Myeloid_State),
  "\n"
)
cat(
  "Non-myeloid-associated cluster:",
  sum(!phase5_final_annotations$Interpretation_Eligible_As_Myeloid_State),
  "\n"
)
cat(
  "Output:",
  file.path(
    table_dir,
    "phase5_myeloid_final_biological_annotations.csv"
  ),
  "\n"
)
cat("============================================================\n")
