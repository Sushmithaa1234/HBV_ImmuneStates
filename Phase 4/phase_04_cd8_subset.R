# ============================================================
# PHASE 4: CD8 T-CELL ANALYSIS
# Script 04: CD8 subset
#
# Goal:
# Isolate and characterize CD8 T-cell states across
# HBV clinical states.
#
# Section 1: Setup and reproducibility
# ============================================================


# ------------------------------------------------------------
# 1.1 Load required libraries
# ------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(here)
library(Matrix)


# ------------------------------------------------------------
# 1.2 Set project root
# ------------------------------------------------------------

setwd(here())


# ------------------------------------------------------------
# 1.3 Reproducibility
# ------------------------------------------------------------

set.seed(12345)


# ------------------------------------------------------------
# 1.4 Define input/output paths
# ------------------------------------------------------------

phase3_final_path <- here(
  "results",
  "rds_objects",
  "phase3_final_full_dataset.rds"
)

output_rds_dir <- here(
  "results",
  "rds_objects"
)

output_table_dir <- here(
  "results",
  "tables"
)

output_figure_dir <- here(
  "results",
  "figures"
)


# ------------------------------------------------------------
# 1.5 Create output directories if needed
# ------------------------------------------------------------

dir.create(
  output_rds_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1.6 Confirm project structure
# ------------------------------------------------------------

stopifnot(
  dir.exists(here("results")),
  dir.exists(output_rds_dir),
  dir.exists(output_table_dir),
  dir.exists(output_figure_dir),
  file.exists(phase3_final_path)
)


# ------------------------------------------------------------
# 1.7 Analysis constants
# ------------------------------------------------------------

target_cell_type <- "CD8_T"

analysis_seed <- 12345

message("============================================================")
message("PHASE 4: CD8 T-CELL ANALYSIS")
message("============================================================")
message("Project root: ", here())
message("Input object: ", phase3_final_path)
message("Target population: ", target_cell_type)
message("Random seed: ", analysis_seed)
message("============================================================")

# ============================================================
# 2. LOAD AND VALIDATE PHASE 3 FINAL DATASET
# ============================================================


# ------------------------------------------------------------
# 2.1 Load final Phase 3 object
# ------------------------------------------------------------

seurat_obj <- readRDS(phase3_final_path)


# ------------------------------------------------------------
# 2.2 Basic object validation
# ------------------------------------------------------------

message("Checking Phase 3 final dataset...")

stopifnot(
  inherits(seurat_obj, "Seurat"),
  ncol(seurat_obj) == 105220,
  nrow(seurat_obj) == 18925
)


# ------------------------------------------------------------
# 2.3 Required metadata validation
# ------------------------------------------------------------

# ------------------------------------------------------------
# 2.3 Required metadata validation
# ------------------------------------------------------------

required_metadata <- c(
  "GSM",
  "Donor",
  "Phase",
  "Transferred_Label",
  "Transfer_Confidence",
  "Neighbour_Agreement",
  "Label_Margin",
  "Low_Confidence"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(seurat_obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2.4 Validate clinical phases
# ------------------------------------------------------------

expected_phases <- c(
  "NL",
  "IT",
  "IA",
  "AR",
  "AC"
)

observed_phases <- sort(unique(seurat_obj$Phase))

if (!setequal(observed_phases, expected_phases)) {
  stop(
    "Unexpected clinical phases detected: ",
    paste(observed_phases, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2.5 Validate donor/sample structure
# ------------------------------------------------------------

n_donors <- dplyr::n_distinct(seurat_obj$Donor)
n_gsm <- dplyr::n_distinct(seurat_obj$GSM)

if (n_donors != 23) {
  stop("Expected 23 donors; found ", n_donors)
}

if (n_gsm != 23) {
  stop("Expected 23 GSM/sample identifiers; found ", n_gsm)
}


# ------------------------------------------------------------
# 2.6 Validate transferred cell types
# ------------------------------------------------------------

transferred_labels <- sort(
  unique(seurat_obj$Transferred_Label)
)

message(
  "Transferred cell types: ",
  paste(transferred_labels, collapse = ", ")
)

if (!(target_cell_type %in% transferred_labels)) {
  stop(
    "Target population '",
    target_cell_type,
    "' was not found in Transferred_Label."
  )
}


# ------------------------------------------------------------
# 2.7 Count CD8 T cells
# ------------------------------------------------------------

cd8_n <- sum(
  seurat_obj$Transferred_Label == target_cell_type,
  na.rm = TRUE
)

if (cd8_n == 0) {
  stop("No CD8 T cells detected.")
}


# ------------------------------------------------------------
# 2.8 Confirm no missing transferred labels
# ------------------------------------------------------------

missing_transferred_labels <- sum(
  is.na(seurat_obj$Transferred_Label)
)

if (missing_transferred_labels > 0) {
  stop(
    "Found ",
    missing_transferred_labels,
    " cells with missing Transferred_Label."
  )
}


# ------------------------------------------------------------
# 2.9 Entry checkpoint summary
# ------------------------------------------------------------

message("")
message("============================================================")
message("PHASE 4 ENTRY CHECKPOINT")
message("============================================================")
message("Total cells: ", ncol(seurat_obj))
message("Total genes: ", nrow(seurat_obj))
message("Donors: ", n_donors)
message("GSM/sample IDs: ", n_gsm)
message("Clinical phases: ", paste(expected_phases, collapse = ", "))
message("CD8 T cells: ", cd8_n)
message(
  "CD8 fraction: ",
  round(100 * cd8_n / ncol(seurat_obj), 2),
  "%"
)
message(
  "Missing transferred labels: ",
  missing_transferred_labels
)
message("============================================================")

# ============================================================
# 3. EXTRACT CD8 T CELLS
# ============================================================


# ------------------------------------------------------------
# 3.1 Identify CD8 T-cell barcodes
# ------------------------------------------------------------

cd8_cells <- rownames(seurat_obj@meta.data)[
  seurat_obj$Transferred_Label == target_cell_type
]


# ------------------------------------------------------------
# 3.2 Validate CD8 cell selection
# ------------------------------------------------------------

if (length(cd8_cells) != cd8_n) {
  stop(
    "CD8 cell count mismatch: expected ",
    cd8_n,
    " but identified ",
    length(cd8_cells),
    " cells."
  )
}

if (anyDuplicated(cd8_cells) > 0) {
  stop("Duplicated cell barcodes detected in CD8 selection.")
}


# ------------------------------------------------------------
# 3.3 Extract CD8 subset
# ------------------------------------------------------------

cd8_obj <- subset(
  x = seurat_obj,
  cells = cd8_cells
)


# ------------------------------------------------------------
# 3.4 Validate extracted object
# ------------------------------------------------------------

stopifnot(
  inherits(cd8_obj, "Seurat"),
  ncol(cd8_obj) == cd8_n,
  nrow(cd8_obj) == nrow(seurat_obj)
)


# ------------------------------------------------------------
# 3.5 Confirm all extracted cells are CD8 T cells
# ------------------------------------------------------------

if (!all(cd8_obj$Transferred_Label == target_cell_type)) {
  stop(
    "Extracted object contains cells that are not labelled CD8_T."
  )
}


# ------------------------------------------------------------
# 3.6 Confirm required provenance metadata survived
# ------------------------------------------------------------

required_cd8_metadata <- c(
  "GSM",
  "Donor",
  "Phase",
  "Transferred_Label",
  "Transfer_Confidence",
  "Neighbour_Agreement",
  "Label_Margin",
  "Low_Confidence"
)

missing_cd8_metadata <- setdiff(
  required_cd8_metadata,
  colnames(cd8_obj@meta.data)
)

if (length(missing_cd8_metadata) > 0) {
  stop(
    "CD8 object is missing required metadata: ",
    paste(missing_cd8_metadata, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 3.7 CD8 clinical-state composition
# ------------------------------------------------------------

cd8_phase_counts <- cd8_obj@meta.data |>
  dplyr::count(Phase, name = "CD8_Cell_Count") |>
  dplyr::arrange(match(Phase, c("NL", "IT", "IA", "AR", "AC")))


# ------------------------------------------------------------
# 3.8 CD8 donor/sample representation
# ------------------------------------------------------------

cd8_donor_counts <- cd8_obj@meta.data |>
  dplyr::count(
    Phase,
    Donor,
    GSM,
    name = "CD8_Cell_Count"
  ) |>
  dplyr::arrange(
    match(Phase, c("NL", "IT", "IA", "AR", "AC")),
    Donor
  )


# ------------------------------------------------------------
# 3.9 Print extraction checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 EXTRACTION CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("Genes: ", nrow(cd8_obj))
message("Clinical states represented: ", dplyr::n_distinct(cd8_obj$Phase))
message("Donors represented: ", dplyr::n_distinct(cd8_obj$Donor))
message("")
message("CD8 cells by clinical state:")
print(cd8_phase_counts)
message("")
message("CD8 cells by donor/sample:")
print(cd8_donor_counts)
message("============================================================")

# ============================================================
# 4. CD8 SUBSET QC CHARACTERIZATION
# ============================================================


# ------------------------------------------------------------
# 4.1 Validate existing QC metadata
# ------------------------------------------------------------

required_qc_metadata <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "RNA_complexity",
  "QC_Flag",
  "QC_Pass",
  "scDblFinder.score",
  "scDblFinder.class",
  "Potential_Doublet"
)

missing_qc_metadata <- setdiff(
  required_qc_metadata,
  colnames(cd8_obj@meta.data)
)

if (length(missing_qc_metadata) > 0) {
  stop(
    "CD8 object is missing required QC metadata: ",
    paste(missing_qc_metadata, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4.2 Confirm CD8 subset passed original QC
# ------------------------------------------------------------

cd8_qc_flag_table <- table(
  cd8_obj$QC_Pass,
  useNA = "ifany"
)

cd8_qc_flag_table


# ------------------------------------------------------------
# 4.3 Confirm doublet status
# ------------------------------------------------------------

cd8_doublet_table <- table(
  cd8_obj$scDblFinder.class,
  useNA = "ifany"
)

cd8_doublet_table


# ------------------------------------------------------------
# 4.4 Summarize continuous QC metrics
# ------------------------------------------------------------

cd8_qc_summary <- cd8_obj@meta.data |>
  dplyr::summarise(
    Cells = dplyr::n(),
    
    nCount_RNA_Median = median(nCount_RNA, na.rm = TRUE),
    nCount_RNA_Mean = mean(nCount_RNA, na.rm = TRUE),
    nCount_RNA_Min = min(nCount_RNA, na.rm = TRUE),
    nCount_RNA_Max = max(nCount_RNA, na.rm = TRUE),
    
    nFeature_RNA_Median = median(nFeature_RNA, na.rm = TRUE),
    nFeature_RNA_Mean = mean(nFeature_RNA, na.rm = TRUE),
    nFeature_RNA_Min = min(nFeature_RNA, na.rm = TRUE),
    nFeature_RNA_Max = max(nFeature_RNA, na.rm = TRUE),
    
    percent_mt_Median = median(percent.mt, na.rm = TRUE),
    percent_mt_Mean = mean(percent.mt, na.rm = TRUE),
    percent_mt_Min = min(percent.mt, na.rm = TRUE),
    percent_mt_Max = max(percent.mt, na.rm = TRUE),
    
    RNA_complexity_Median = median(RNA_complexity, na.rm = TRUE),
    RNA_complexity_Mean = mean(RNA_complexity, na.rm = TRUE),
    RNA_complexity_Min = min(RNA_complexity, na.rm = TRUE),
    RNA_complexity_Max = max(RNA_complexity, na.rm = TRUE)
  )

print(cd8_qc_summary)


# ------------------------------------------------------------
# 4.5 QC summary by clinical state
# ------------------------------------------------------------

cd8_qc_by_phase <- cd8_obj@meta.data |>
  dplyr::group_by(Phase) |>
  dplyr::summarise(
    Cells = dplyr::n(),
    
    Median_nCount_RNA = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_nFeature_RNA = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_percent_mt = median(
      percent.mt,
      na.rm = TRUE
    ),
    
    Median_RNA_complexity = median(
      RNA_complexity,
      na.rm = TRUE
    ),
    
    Doublet_Rate = mean(
      scDblFinder.class == "doublet",
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(
      Phase,
      c("NL", "IT", "IA", "AR", "AC")
    )
  )

print(cd8_qc_by_phase)


# ------------------------------------------------------------
# 4.6 QC summary by donor
# ------------------------------------------------------------

cd8_qc_by_donor <- cd8_obj@meta.data |>
  dplyr::group_by(
    Phase,
    Donor,
    GSM
  ) |>
  dplyr::summarise(
    Cells = dplyr::n(),
    
    Median_nCount_RNA = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_nFeature_RNA = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_percent_mt = median(
      percent.mt,
      na.rm = TRUE
    ),
    
    Median_RNA_complexity = median(
      RNA_complexity,
      na.rm = TRUE
    ),
    
    Doublet_Rate = mean(
      scDblFinder.class == "doublet",
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(
      Phase,
      c("NL", "IT", "IA", "AR", "AC")
    ),
    Donor
  )

print(cd8_qc_by_donor)


# ------------------------------------------------------------
# 4.7 CD8 transfer-confidence summary
# ------------------------------------------------------------

cd8_transfer_summary <- cd8_obj@meta.data |>
  dplyr::count(
    Transfer_Confidence,
    Low_Confidence,
    name = "Cells"
  ) |>
  dplyr::arrange(
    Transfer_Confidence,
    Low_Confidence
  )

print(cd8_transfer_summary)


# ------------------------------------------------------------
# 4.8 Save QC tables
# ------------------------------------------------------------

write.csv(
  cd8_qc_summary,
  file = here(
    "results",
    "tables",
    "phase4_cd8_qc_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_qc_by_phase,
  file = here(
    "results",
    "tables",
    "phase4_cd8_qc_by_phase.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_qc_by_donor,
  file = here(
    "results",
    "tables",
    "phase4_cd8_qc_by_donor.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_transfer_summary,
  file = here(
    "results",
    "tables",
    "phase4_cd8_transfer_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4.9 Checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 QC CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("QC-pass status:")
print(cd8_qc_flag_table)
message("")
message("Doublet status:")
print(cd8_doublet_table)
message("")
message("Overall QC summary:")
print(cd8_qc_summary)
message("")
message("QC summary by clinical state:")
print(cd8_qc_by_phase)
message("============================================================")

# ============================================================
# 5. LOCK CD8 NORMALIZED EXPRESSION REPRESENTATION
# ============================================================


# ------------------------------------------------------------
# 5.1 Confirm normalized RNA data layer
# ------------------------------------------------------------

cd8_rna_layers <- SeuratObject::Layers(
  cd8_obj[["RNA"]]
)

if (!("data" %in% cd8_rna_layers)) {
  stop(
    "Normalized RNA data layer was not found."
  )
}


# ------------------------------------------------------------
# 5.2 Extract normalized CD8 expression
# ------------------------------------------------------------

cd8_data <- SeuratObject::LayerData(
  cd8_obj,
  assay = "RNA",
  layer = "data"
)


# ------------------------------------------------------------
# 5.3 Validate dimensions
# ------------------------------------------------------------

if (
  nrow(cd8_data) != nrow(cd8_obj) ||
  ncol(cd8_data) != ncol(cd8_obj)
) {
  stop(
    "Normalized CD8 expression dimensions do not match ",
    "the CD8 object."
  )
}


# ------------------------------------------------------------
# 5.4 Validate feature and cell order
# ------------------------------------------------------------

if (!identical(
  rownames(cd8_data),
  rownames(cd8_obj)
)) {
  stop(
    "Normalized expression feature order does not match ",
    "the CD8 object."
  )
}

if (!identical(
  colnames(cd8_data),
  colnames(cd8_obj)
)) {
  stop(
    "Normalized expression cell order does not match ",
    "the CD8 object."
  )
}


# ------------------------------------------------------------
# 5.5 Validate finite expression values
# ------------------------------------------------------------

if (any(!is.finite(cd8_data@x))) {
  stop(
    "Non-finite values detected in normalized CD8 expression."
  )
}


# ------------------------------------------------------------
# 5.6 Record expression provenance
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Normalization <- list(
  expression_assay = "RNA",
  expression_layer = "data",
  normalization_status =
    "Existing normalized expression retained",
  raw_count_re_normalization =
    "Not performed because raw integer count matrix is not available in the Phase 2/3 checkpoints",
  rationale =
    "Avoid applying normalization to an already transformed expression representation"
)


# ------------------------------------------------------------
# 5.7 Section 5 checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 EXPRESSION REPRESENTATION CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("Genes: ", nrow(cd8_obj))
message("")
message("RNA layers:")
print(cd8_rna_layers)
message("")
message("Analysis assay: RNA")
message("Analysis layer: data")
message("Expression dimensions: ", nrow(cd8_data), " x ", ncol(cd8_data))
message("Raw-count re-normalization: NOT performed")
message("Existing normalized expression retained: YES")
message("============================================================")

# ============================================================
# 6. CD8-SPECIFIC VARIABLE-FEATURE SELECTION
# ============================================================


# ------------------------------------------------------------
# 6.1 Define CD8 variable-feature parameters
# ------------------------------------------------------------

cd8_nfeatures <- 2000
cd8_hvg_method <- "vst"


# ------------------------------------------------------------
# 6.2 Confirm normalized expression is available
# ------------------------------------------------------------

cd8_rna_layers <- SeuratObject::Layers(
  cd8_obj[["RNA"]]
)

if (!("data" %in% cd8_rna_layers)) {
  stop(
    "Normalized RNA data layer is required for CD8 HVG selection."
  )
}


# ------------------------------------------------------------
# 6.3 Identify variable features
# ------------------------------------------------------------

message("")
message("Running CD8-specific variable-feature selection...")
message(
  "Method: ",
  cd8_hvg_method
)
message(
  "Target features: ",
  cd8_nfeatures
)

cd8_obj <- Seurat::FindVariableFeatures(
  object = cd8_obj,
  assay = "RNA",
  selection.method = cd8_hvg_method,
  nfeatures = cd8_nfeatures,
  verbose = TRUE
)


# ------------------------------------------------------------
# 6.4 Retrieve selected features
# ------------------------------------------------------------

cd8_hvgs <- Seurat::VariableFeatures(
  cd8_obj,
  assay = "RNA"
)


# ------------------------------------------------------------
# 6.5 Validate HVG selection
# ------------------------------------------------------------

if (length(cd8_hvgs) != cd8_nfeatures) {
  stop(
    "Expected ",
    cd8_nfeatures,
    " CD8 variable features but found ",
    length(cd8_hvgs),
    "."
  )
}

if (anyDuplicated(cd8_hvgs) > 0) {
  stop("Duplicated CD8 variable features detected.")
}

if (!all(cd8_hvgs %in% rownames(cd8_obj))) {
  stop(
    "One or more CD8 variable features are absent from the object."
  )
}


# ------------------------------------------------------------
# 6.6 Check for missing feature names
# ------------------------------------------------------------

if (any(is.na(cd8_hvgs)) || any(cd8_hvgs == "")) {
  stop("Missing or empty CD8 variable-feature names detected.")
}


# ------------------------------------------------------------
# 6.7 Record HVG provenance
# ------------------------------------------------------------

cd8_obj@misc$Phase4_HVG <- list(
  method = cd8_hvg_method,
  nfeatures = cd8_nfeatures,
  assay = "RNA",
  layer = "data"
)


# ------------------------------------------------------------
# 6.8 Save CD8 marker-independent HVG table
# ------------------------------------------------------------

cd8_hvg_table <- data.frame(
  Rank = seq_along(cd8_hvgs),
  Gene = cd8_hvgs,
  stringsAsFactors = FALSE
)

write.csv(
  cd8_hvg_table,
  file = here(
    "results",
    "tables",
    "phase4_cd8_variable_features.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6.9 Section 6 checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 VARIABLE-FEATURE CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("HVG method: ", cd8_hvg_method)
message("Requested HVGs: ", cd8_nfeatures)
message("Selected HVGs: ", length(cd8_hvgs))
message("")
message("Top 30 CD8 variable features:")
print(head(cd8_hvgs, 30))
message("============================================================")

# ============================================================
# 7. CHARACTERIZE CD8 VARIABLE-FEATURE COMPOSITION
# ============================================================


# ------------------------------------------------------------
# 7.1 Retrieve CD8 HVGs
# ------------------------------------------------------------

cd8_hvgs <- Seurat::VariableFeatures(
  cd8_obj,
  assay = "RNA"
)

if (length(cd8_hvgs) != 2000) {
  stop(
    "Expected 2,000 CD8 HVGs; found ",
    length(cd8_hvgs),
    "."
  )
}


# ------------------------------------------------------------
# 7.2 Define broad feature categories
# ------------------------------------------------------------

cd8_hvg_categories <- tibble::tibble(
  Gene = cd8_hvgs
) |>
  dplyr::mutate(
    Category = dplyr::case_when(
      
      grepl(
        "^MT-",
        Gene,
        ignore.case = TRUE
      ) ~ "Mitochondrial",
      
      grepl(
        "^RPS|^RPL",
        Gene,
        ignore.case = TRUE
      ) ~ "Ribosomal",
      
      grepl(
        "^TRAV|^TRBV|^TRGV|^TRDV|^TRAJ|^TRBJ|^TRDJ|^TRGJ",
        Gene,
        ignore.case = TRUE
      ) ~ "TCR",
      
      grepl(
        "^IGH|^IGK|^IGL",
        Gene,
        ignore.case = TRUE
      ) ~ "Immunoglobulin",
      
      grepl(
        "MKI67|TOP2A|TYMS|PCNA|STMN1|TUBA1B",
        Gene,
        ignore.case = TRUE
      ) ~ "Cell_Cycle",
      
      grepl(
        "HSPA1A|HSPA1B|FOS|JUN|JUNB|DUSP1|DUSP2|IER2|EGR1|EGR2",
        Gene,
        ignore.case = TRUE
      ) ~ "Stress_Immediate_Early",
      
      TRUE ~ "Other"
    )
  )


# ------------------------------------------------------------
# 7.3 Summarize category composition
# ------------------------------------------------------------

cd8_hvg_category_summary <- cd8_hvg_categories |>
  dplyr::count(
    Category,
    name = "HVG_Count"
  ) |>
  dplyr::mutate(
    Percentage = 100 * HVG_Count / length(cd8_hvgs)
  ) |>
  dplyr::arrange(
    dplyr::desc(HVG_Count)
  )


# ------------------------------------------------------------
# 7.4 Display category summary
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 HVG COMPOSITION")
message("============================================================")

print(cd8_hvg_category_summary)


# ------------------------------------------------------------
# 7.5 Display top HVGs by category
# ------------------------------------------------------------

message("")
message("Top TCR HVGs:")
print(
  cd8_hvg_categories |>
    dplyr::filter(Category == "TCR") |>
    dplyr::slice_head(n = 30)
)

message("")
message("Top stress/immediate-early HVGs:")
print(
  cd8_hvg_categories |>
    dplyr::filter(Category == "Stress_Immediate_Early") |>
    dplyr::slice_head(n = 30)
)


# ------------------------------------------------------------
# 7.6 Save HVG composition table
# ------------------------------------------------------------

write.csv(
  cd8_hvg_categories,
  file = here(
    "results",
    "tables",
    "phase4_cd8_hvg_categories.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_hvg_category_summary,
  file = here(
    "results",
    "tables",
    "phase4_cd8_hvg_category_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7.7 Checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("CD8 HVG COMPOSITION CHECKPOINT")
message("============================================================")
message("Total HVGs: ", length(cd8_hvgs))
message("")
print(cd8_hvg_category_summary)
message("============================================================")

# ============================================================
# SECTION 8 — CD8-specific scaling
# ============================================================

message("============================================================")
message("SECTION 8 — CD8-specific scaling")
message("============================================================")

# ------------------------------------------------------------
# 8.1 Retrieve locked CD8 HVGs
# ------------------------------------------------------------

cd8_hvgs <- Seurat::VariableFeatures(
  cd8_obj,
  assay = "RNA"
)

if (length(cd8_hvgs) != 2000) {
  stop(
    "Expected exactly 2,000 locked CD8 HVGs; found ",
    length(cd8_hvgs), "."
  )
}

if (anyDuplicated(cd8_hvgs) > 0) {
  stop("Locked CD8 HVG list contains duplicated features.")
}

message("Locked CD8 HVGs: ", length(cd8_hvgs))


# ------------------------------------------------------------
# 8.2 Scale the locked CD8 HVGs
# ------------------------------------------------------------

cd8_obj <- Seurat::ScaleData(
  object = cd8_obj,
  assay = "RNA",
  features = cd8_hvgs,
  verbose = TRUE
)


# ------------------------------------------------------------
# 8.3 Confirm scale.data exists
# ------------------------------------------------------------

rna_layers_after_scaling <- SeuratObject::Layers(
  cd8_obj[["RNA"]]
)

message("RNA layers after scaling:")
print(rna_layers_after_scaling)

if (!"scale.data" %in% rna_layers_after_scaling) {
  stop("RNA scale.data layer was not created.")
}


# ------------------------------------------------------------
# 8.4 Retrieve scaled expression
# ------------------------------------------------------------

cd8_scaled_data <- SeuratObject::LayerData(
  cd8_obj,
  assay = "RNA",
  layer = "scale.data"
)

message(
  "Retrieved scale.data dimensions: ",
  nrow(cd8_scaled_data),
  " × ",
  ncol(cd8_scaled_data)
)


# ------------------------------------------------------------
# 8.5 Confirm all locked HVGs are present
# ------------------------------------------------------------

missing_hvgs <- setdiff(
  cd8_hvgs,
  rownames(cd8_scaled_data)
)

if (length(missing_hvgs) > 0) {
  stop(
    "Some locked CD8 HVGs are missing from scale.data: ",
    paste(head(missing_hvgs, 20), collapse = ", ")
  )
}


# ------------------------------------------------------------
# 8.6 Explicitly restore the locked HVG order
# ------------------------------------------------------------

cd8_scaled_data <- cd8_scaled_data[
  cd8_hvgs,
  ,
  drop = FALSE
]

if (!identical(
  rownames(cd8_scaled_data),
  cd8_hvgs
)) {
  stop(
    "Failed to restore the locked CD8 HVG order."
  )
}


# ------------------------------------------------------------
# 8.7 Validate cell order
# ------------------------------------------------------------

if (!identical(
  colnames(cd8_scaled_data),
  colnames(cd8_obj)
)) {
  stop(
    "Scaled expression cell order does not match CD8 object."
  )
}


# ------------------------------------------------------------
# 8.8 Validate dimensions and finite values
# ------------------------------------------------------------

if (!all(
  dim(cd8_scaled_data) ==
  c(length(cd8_hvgs), ncol(cd8_obj))
)) {
  stop(
    "Scaled expression dimensions are incorrect."
  )
}

if (any(!is.finite(cd8_scaled_data))) {
  stop(
    "Scaled expression contains non-finite values."
  )
}


# ------------------------------------------------------------
# 8.9 Scaling diagnostics
# ------------------------------------------------------------

scaled_row_means <- Matrix::rowMeans(cd8_scaled_data)

scaled_row_variances <- apply(
  cd8_scaled_data,
  1,
  var
)

scaling_summary <- data.frame(
  feature = cd8_hvgs,
  mean = scaled_row_means,
  variance = scaled_row_variances
)

message("")
message("Scaling diagnostics:")
message(
  "Maximum absolute row mean: ",
  signif(max(abs(scaled_row_means)), 6)
)

message(
  "Median row variance: ",
  signif(median(scaled_row_variances), 6)
)

message(
  "Variance range: ",
  signif(min(scaled_row_variances), 6),
  " – ",
  signif(max(scaled_row_variances), 6)
)


# ------------------------------------------------------------
# 8.10 Save scaling summary
# ------------------------------------------------------------

# ------------------------------------------------------------
# Restore output directories
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

# Save scaling summary
write.csv(
  scaling_summary,
  file = file.path(
    tables_dir,
    "phase4_cd8_scaling_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  scaling_summary,
  file = file.path(
    tables_dir,
    "phase4_cd8_scaling_summary.csv"
  ),
  row.names = FALSE
)

print(
  file.exists(
    file.path(
      tables_dir,
      "phase4_cd8_scaling_summary.csv"
    )
  )
)


# ------------------------------------------------------------
# 8.11 Checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 8 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("Scaled features: ", nrow(cd8_scaled_data))
message("Expected features: 2000")
message("Feature order locked: YES")
message("Cell order preserved: YES")
message("Finite values: YES")
message("Scaling summary saved.")
message("============================================================")

# ============================================================
# SECTION 9 — CD8 PCA
# ============================================================

message("============================================================")
message("SECTION 9 — CD8 PCA")
message("============================================================")


# ------------------------------------------------------------
# 9.1 Retrieve and validate locked CD8 HVGs
# ------------------------------------------------------------

cd8_hvgs <- Seurat::VariableFeatures(
  cd8_obj,
  assay = "RNA"
)

if (length(cd8_hvgs) != 2000) {
  stop(
    "Expected exactly 2,000 locked CD8 HVGs; found ",
    length(cd8_hvgs), "."
  )
}

if (anyDuplicated(cd8_hvgs) > 0) {
  stop(
    "Locked CD8 HVG list contains duplicated features."
  )
}

message("Locked CD8 HVGs: ", length(cd8_hvgs))


# ------------------------------------------------------------
# 9.2 Confirm scale.data contains all locked HVGs
# ------------------------------------------------------------

rna_layers <- SeuratObject::Layers(
  cd8_obj[["RNA"]]
)

if (!"scale.data" %in% rna_layers) {
  stop(
    "RNA scale.data layer is missing."
  )
}

cd8_scaled_data <- SeuratObject::LayerData(
  cd8_obj,
  assay = "RNA",
  layer = "scale.data"
)

missing_hvgs <- setdiff(
  cd8_hvgs,
  rownames(cd8_scaled_data)
)

if (length(missing_hvgs) > 0) {
  stop(
    "Locked CD8 HVGs are missing from scale.data."
  )
}

if (nrow(cd8_scaled_data) < length(cd8_hvgs)) {
  stop(
    "scale.data contains fewer rows than the locked HVG set."
  )
}

message(
  "scale.data contains all ",
  length(cd8_hvgs),
  " locked CD8 HVGs."
)


# ------------------------------------------------------------
# 9.3 Run PCA
# ------------------------------------------------------------

message("")
message("Running PCA on the 2,000 locked CD8 HVGs...")

cd8_obj <- Seurat::RunPCA(
  object = cd8_obj,
  assay = "RNA",
  features = cd8_hvgs,
  npcs = 50,
  verbose = TRUE
)


# ------------------------------------------------------------
# 9.4 Confirm PCA reduction exists
# ------------------------------------------------------------

if (!"pca" %in% names(cd8_obj@reductions)) {
  stop(
    "PCA reduction was not created."
  )
}

pca_embeddings <- SeuratObject::Embeddings(
  cd8_obj,
  reduction = "pca"
)

message("")
message(
  "PCA embedding dimensions: ",
  nrow(pca_embeddings),
  " × ",
  ncol(pca_embeddings)
)


# ------------------------------------------------------------
# 9.5 Validate PCA cell order
# ------------------------------------------------------------

if (!identical(
  rownames(pca_embeddings),
  colnames(cd8_obj)
)) {
  stop(
    "PCA embedding cell order does not match the CD8 object."
  )
}


# ------------------------------------------------------------
# 9.6 Validate number of PCs
# ------------------------------------------------------------

if (ncol(pca_embeddings) != 50) {
  stop(
    "Expected 50 PCs; found ",
    ncol(pca_embeddings), "."
  )
}


# ------------------------------------------------------------
# 9.7 Validate finite PCA embeddings
# ------------------------------------------------------------

if (any(!is.finite(pca_embeddings))) {
  stop(
    "PCA embeddings contain non-finite values."
  )
}


# ------------------------------------------------------------
# 9.8 Extract PCA variance information
# ------------------------------------------------------------

pca_stdev <- cd8_obj[["pca"]]@stdev

if (length(pca_stdev) != 50) {
  stop(
    "Expected 50 PCA standard deviations; found ",
    length(pca_stdev), "."
  )
}

pca_variance <- pca_stdev^2

pca_variance_summary <- data.frame(
  PC = seq_along(pca_stdev),
  Standard_Deviation = pca_stdev,
  Variance = pca_variance,
  Percent_Variance = 100 * pca_variance / sum(pca_variance)
)

message("")
message("Generating PCA elbow plot...")

elbow_plot <- Seurat::ElbowPlot(
  object = cd8_obj,
  reduction = "pca",
  ndims = 50
)

print(elbow_plot)


# ------------------------------------------------------------
# 9.9 Cumulative variance
# ------------------------------------------------------------

pca_variance_summary$Cumulative_Percent_Variance <-
  cumsum(
    pca_variance_summary$Percent_Variance
  )


# ------------------------------------------------------------
# 9.10 Save PCA variance table
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  pca_variance_summary,
  file = file.path(
    tables_dir,
    "phase4_cd8_pca_variance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9.11 Print PCA variance diagnostics
# ------------------------------------------------------------

message("")
message("First 10 PCs:")
print(
  pca_variance_summary[1:10, ]
)

message("")
message(
  "Cumulative variance after PC10: ",
  round(
    pca_variance_summary$Cumulative_Percent_Variance[10],
    2
  ),
  "%"
)

message(
  "Cumulative variance after PC20: ",
  round(
    pca_variance_summary$Cumulative_Percent_Variance[20],
    2
  ),
  "%"
)

message(
  "Cumulative variance after PC30: ",
  round(
    pca_variance_summary$Cumulative_Percent_Variance[30],
    2
  ),
  "%"
)

message(
  "Cumulative variance after PC50: ",
  round(
    pca_variance_summary$Cumulative_Percent_Variance[50],
    2
  ),
  "%"
)


# ------------------------------------------------------------
# 9.12 PCA checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 9 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", nrow(pca_embeddings))
message("PCA dimensions: ", ncol(pca_embeddings))
message("Expected PCs: 50")
message("HVGs used: 2000")
message("PCA embeddings finite: YES")
message("PCA cell order preserved: YES")
message("PCA variance table saved.")
message("============================================================")

# ============================================================
# SECTION 10 — CD8 nearest-neighbor graph and clustering
# ============================================================

message("============================================================")
message("SECTION 10 — CD8 NEIGHBOR GRAPH + CLUSTERING")
message("============================================================")


# ------------------------------------------------------------
# 10.1 Lock downstream dimensions
# ------------------------------------------------------------

cd8_dims <- 1:20

if (length(cd8_dims) != 20) {
  stop("Expected exactly 20 downstream PCs.")
}

pca_embeddings <- SeuratObject::Embeddings(
  cd8_obj,
  reduction = "pca"
)

if (ncol(pca_embeddings) < max(cd8_dims)) {
  stop(
    "CD8 PCA does not contain enough PCs for the locked ",
    "1:20 dimensionality."
  )
}

message("Downstream PCs locked: 1:20")


# ------------------------------------------------------------
# 10.2 Build nearest-neighbor graph
# ------------------------------------------------------------

message("")
message("Building CD8 nearest-neighbor graph...")

cd8_obj <- Seurat::FindNeighbors(
  object = cd8_obj,
  reduction = "pca",
  dims = cd8_dims,
  k.param = 20,
  compute.SNN = TRUE,
  verbose = TRUE
)


# ------------------------------------------------------------
# 10.3 Confirm neighbor graphs
# ------------------------------------------------------------

graph_names <- names(cd8_obj@graphs)

message("")
message("Graphs present:")
print(graph_names)

if (!"RNA_nn" %in% graph_names) {
  stop("Expected RNA_nn graph was not created.")
}

if (!"RNA_snn" %in% graph_names) {
  stop("Expected RNA_snn graph was not created.")
}


# ------------------------------------------------------------
# 10.4 Validate graph dimensions
# ------------------------------------------------------------

nn_graph <- cd8_obj[["RNA_nn"]]
snn_graph <- cd8_obj[["RNA_snn"]]

expected_cells <- ncol(cd8_obj)

if (nrow(nn_graph) != expected_cells ||
    ncol(nn_graph) != expected_cells) {
  stop("RNA_nn graph dimensions do not match CD8 cell count.")
}

if (nrow(snn_graph) != expected_cells ||
    ncol(snn_graph) != expected_cells) {
  stop("RNA_snn graph dimensions do not match CD8 cell count.")
}

message(
  "RNA_nn dimensions: ",
  nrow(nn_graph),
  " × ",
  ncol(nn_graph)
)

message(
  "RNA_snn dimensions: ",
  nrow(snn_graph),
  " × ",
  ncol(snn_graph)
)


# ------------------------------------------------------------
# 10.5 Cluster CD8 cells
# ------------------------------------------------------------

cd8_resolution <- 0.4

message("")
message(
  "Clustering CD8 cells at resolution ",
  cd8_resolution,
  "..."
)

cd8_obj <- Seurat::FindClusters(
  object = cd8_obj,
  graph.name = "RNA_snn",
  resolution = cd8_resolution,
  algorithm = 1,
  random.seed = 12345,
  verbose = TRUE
)


# ------------------------------------------------------------
# 10.6 Confirm cluster identities
# ------------------------------------------------------------

if (!"seurat_clusters" %in% colnames(cd8_obj@meta.data)) {
  stop(
    "seurat_clusters was not created."
  )
}

cd8_clusters <- cd8_obj$seurat_clusters

if (any(is.na(cd8_clusters))) {
  stop(
    "CD8 clustering produced missing cluster assignments."
  )
}

cluster_counts <- table(cd8_clusters)

message("")
message("CD8 cluster counts:")
print(cluster_counts)

message("")
message(
  "Number of CD8 clusters: ",
  length(cluster_counts)
)


# ------------------------------------------------------------
# 10.7 Confirm every CD8 cell has exactly one cluster
# ------------------------------------------------------------

if (length(cd8_clusters) != ncol(cd8_obj)) {
  stop(
    "Number of cluster assignments does not match CD8 cell count."
  )
}

if (sum(cluster_counts) != ncol(cd8_obj)) {
  stop(
    "Cluster counts do not sum to total CD8 cell count."
  )
}


# ------------------------------------------------------------
# 10.8 Cluster size summary
# ------------------------------------------------------------

cluster_summary <- data.frame(
  cluster = names(cluster_counts),
  n_cells = as.integer(cluster_counts),
  fraction = as.integer(cluster_counts) / ncol(cd8_obj)
)

cluster_summary <- cluster_summary[
  order(
    as.numeric(as.character(cluster_summary$cluster))
  ),
  ,
  drop = FALSE
]

message("")
message("Cluster size summary:")
print(cluster_summary)


# ------------------------------------------------------------
# 10.9 Save cluster table
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  cluster_summary,
  file = file.path(
    tables_dir,
    "phase4_cd8_cluster_sizes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10.10 Record locked clustering parameters
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Clustering <- list(
  reduction = "pca",
  dimensions = "1:20",
  k_param = 20,
  graph = "RNA_snn",
  resolution = cd8_resolution,
  algorithm = 1,
  random_seed = 12345
)


# ------------------------------------------------------------
# 10.11 Save Section 10 checkpoint
# ------------------------------------------------------------

checkpoint_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}

saveRDS(
  cd8_obj,
  file = file.path(
    checkpoint_dir,
    "phase4_cd8_neighbors_clustering.rds"
  )
)


# ------------------------------------------------------------
# 10.12 SECTION 10 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 10 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("Downstream PCs: 1:20")
message("k.param: 20")
message("SNN graph: RNA_snn")
message("Clustering resolution: ", cd8_resolution)
message("Clustering algorithm: 1")
message("Number of clusters: ", length(cluster_counts))
message("Missing cluster assignments: ", sum(is.na(cd8_clusters)))
message("Cluster table saved.")
message("Checkpoint object saved.")
message("============================================================")

# ============================================================
# SECTION 11 — CD8 UMAP
# ============================================================

message("============================================================")
message("SECTION 11 — CD8 UMAP")
message("============================================================")


# ------------------------------------------------------------
# 11.1 Confirm locked dimensionality
# ------------------------------------------------------------

cd8_dims <- 1:20

pca_embeddings <- SeuratObject::Embeddings(
  cd8_obj,
  reduction = "pca"
)

if (ncol(pca_embeddings) < max(cd8_dims)) {
  stop("Required PCs 1:20 are not available.")
}


# ------------------------------------------------------------
# 11.2 Run UMAP
# ------------------------------------------------------------

message("Running CD8 UMAP using PCs 1:20...")

cd8_obj <- Seurat::RunUMAP(
  object = cd8_obj,
  reduction = "pca",
  dims = cd8_dims,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  seed.use = 12345,
  verbose = TRUE
)


# ------------------------------------------------------------
# 11.3 Validate UMAP
# ------------------------------------------------------------

if (!"umap" %in% names(cd8_obj@reductions)) {
  stop("UMAP reduction was not created.")
}

umap_embeddings <- SeuratObject::Embeddings(
  cd8_obj,
  reduction = "umap"
)

message(
  "UMAP dimensions: ",
  nrow(umap_embeddings),
  " × ",
  ncol(umap_embeddings)
)

if (nrow(umap_embeddings) != ncol(cd8_obj)) {
  stop("UMAP cell count does not match CD8 object.")
}

if (ncol(umap_embeddings) != 2) {
  stop("Expected a 2-dimensional UMAP.")
}

if (!identical(
  rownames(umap_embeddings),
  colnames(cd8_obj)
)) {
  stop("UMAP cell order does not match CD8 object.")
}

if (any(!is.finite(umap_embeddings))) {
  stop("UMAP contains non-finite values.")
}


# ------------------------------------------------------------
# 11.4 UMAP — CD8 clusters
# ------------------------------------------------------------

message("")
message("Generating CD8 cluster UMAP...")

p_umap_clusters <- Seurat::DimPlot(
  object = cd8_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::ggtitle(
    "CD8 T-cell subclusters"
  ) +
  ggplot2::theme_classic()


# ------------------------------------------------------------
# 11.5 Save cluster UMAP
# ------------------------------------------------------------

figures_dir <- here::here(
  "results",
  "figures"
)

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "phase4_cd8_umap_clusters.png"
  ),
  plot = p_umap_clusters,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 11.6 UMAP — clinical phase
# ------------------------------------------------------------

message("Generating CD8 clinical-state UMAP...")

p_umap_phase <- Seurat::DimPlot(
  object = cd8_obj,
  reduction = "umap",
  group.by = "Phase",
  raster = TRUE
) +
  ggplot2::ggtitle(
    "CD8 T cells by clinical state"
  ) +
  ggplot2::theme_classic()


# ------------------------------------------------------------
# 11.7 Save clinical-state UMAP
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "phase4_cd8_umap_clinical_phase.png"
  ),
  plot = p_umap_phase,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 11.8 Print plots
# ------------------------------------------------------------

print(p_umap_clusters)
print(p_umap_phase)


# ------------------------------------------------------------
# 11.9 Record UMAP parameters
# ------------------------------------------------------------

cd8_obj@misc$Phase4_UMAP <- list(
  reduction = "pca",
  dimensions = "1:20",
  seed = 12345
)


# ------------------------------------------------------------
# 11.10 Save checkpoint
# ------------------------------------------------------------

checkpoint_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}

saveRDS(
  cd8_obj,
  file = file.path(
    checkpoint_dir,
    "phase4_cd8_umap.rds"
  )
)


# ------------------------------------------------------------
# 11.11 SECTION 11 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 11 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("UMAP dimensions: ", ncol(umap_embeddings))
message("Downstream PCs: 1:20")
message("Clusters visualized: ", length(unique(cd8_obj$seurat_clusters)))
message("UMAP embeddings finite: YES")
message("Cluster UMAP saved.")
message("Clinical-state UMAP saved.")
message("Checkpoint object saved.")
message("============================================================")

# ============================================================
# SECTION 12 — CD8 CLUSTER MARKERS
# ============================================================

message("============================================================")
message("SECTION 12 — CD8 CLUSTER MARKERS")
message("============================================================")


# ------------------------------------------------------------
# 12.1 Confirm CD8 clustering is available
# ------------------------------------------------------------

if (!"seurat_clusters" %in% colnames(cd8_obj@meta.data)) {
  stop(
    "seurat_clusters is not present in the CD8 object."
  )
}

cd8_clusters <- cd8_obj$seurat_clusters

if (any(is.na(cd8_clusters))) {
  stop(
    "Missing CD8 cluster assignments."
  )
}

message(
  "CD8 cells: ",
  ncol(cd8_obj)
)

message(
  "CD8 clusters: ",
  length(unique(cd8_clusters))
)


# ------------------------------------------------------------
# 12.2 Confirm RNA assay and expression layer
# ------------------------------------------------------------

rna_layers <- SeuratObject::Layers(
  cd8_obj[["RNA"]]
)

message("")
message("RNA layers:")
print(rna_layers)

if (!"data" %in% rna_layers) {
  stop(
    "RNA data layer is missing."
  )
}


# ------------------------------------------------------------
# 12.3 Identify positive cluster markers
# ------------------------------------------------------------

message("")
message("Finding positive markers for all CD8 clusters...")
message("Test: Wilcoxon rank-sum")
message("min.pct: 0.25")
message("logfc.threshold: 0.25")

cd8_markers <- Seurat::FindAllMarkers(
  object = cd8_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)


# ------------------------------------------------------------
# 12.4 Validate marker table
# ------------------------------------------------------------

if (!is.data.frame(cd8_markers)) {
  stop(
    "CD8 marker result is not a data.frame."
  )
}

if (nrow(cd8_markers) == 0) {
  stop(
    "No CD8 cluster markers were identified."
  )
}

required_marker_columns <- c(
  "cluster",
  "gene",
  "p_val",
  "p_val_adj",
  "avg_log2FC",
  "pct.1",
  "pct.2"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(cd8_markers)
)

if (length(missing_marker_columns) > 0) {
  stop(
    "Marker table is missing required columns: ",
    paste(
      missing_marker_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 12.5 Restrict to significant markers for interpretation
# ------------------------------------------------------------

cd8_markers_sig <- cd8_markers[
  !is.na(cd8_markers$p_val_adj) &
    cd8_markers$p_val_adj < 0.05,
  ,
  drop = FALSE
]

message("")
message(
  "Total positive marker rows: ",
  nrow(cd8_markers)
)

message(
  "Significant positive marker rows: ",
  nrow(cd8_markers_sig)
)


# ------------------------------------------------------------
# 12.6 Check marker coverage across clusters
# ------------------------------------------------------------

clusters_with_markers <- sort(
  unique(
    as.character(cd8_markers_sig$cluster)
  )
)

all_clusters <- sort(
  unique(
    as.character(cd8_clusters)
  )
)

clusters_without_markers <- setdiff(
  all_clusters,
  clusters_with_markers
)

if (length(clusters_without_markers) > 0) {
  warning(
    "Clusters without significant positive markers: ",
    paste(
      clusters_without_markers,
      collapse = ", "
    )
  )
}

message("")
message("Significant-marker coverage by cluster:")

marker_counts_by_cluster <- table(
  cd8_markers_sig$cluster
)

print(marker_counts_by_cluster)


# ------------------------------------------------------------
# 12.7 Rank markers within each cluster
# ------------------------------------------------------------

cd8_markers_ranked <- cd8_markers_sig[
  order(
    as.numeric(
      as.character(cd8_markers_sig$cluster)
    ),
    -cd8_markers_sig$avg_log2FC,
    cd8_markers_sig$p_val_adj
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12.8 Extract top 10 markers per cluster
# ------------------------------------------------------------

top10_cd8_markers <- cd8_markers_ranked |>
  dplyr::group_by(cluster) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup()

message("")
message("Top 10 markers per CD8 cluster:")

print(
  top10_cd8_markers[
    ,
    c(
      "cluster",
      "gene",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 12.9 Save complete marker table
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  cd8_markers,
  file = file.path(
    tables_dir,
    "phase4_cd8_all_cluster_markers.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd8_markers_sig,
  file = file.path(
    tables_dir,
    "phase4_cd8_significant_cluster_markers.csv"
  ),
  row.names = FALSE
)

write.csv(
  top10_cd8_markers,
  file = file.path(
    tables_dir,
    "phase4_cd8_top10_markers_per_cluster.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12.10 Record marker-analysis parameters
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Markers <- list(
  method = "FindAllMarkers",
  assay = "RNA",
  expression_layer = "data",
  test = "wilcox",
  only_positive = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  significance_threshold =
    "adjusted p-value < 0.05"
)


# ------------------------------------------------------------
# 12.11 Save checkpoint
# ------------------------------------------------------------

checkpoint_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}

saveRDS(
  cd8_obj,
  file = file.path(
    checkpoint_dir,
    "phase4_cd8_markers.rds"
  )
)


# ------------------------------------------------------------
# 12.12 SECTION 12 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 12 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("CD8 clusters: ", length(unique(cd8_clusters)))
message("Total positive marker rows: ", nrow(cd8_markers))
message(
  "Significant positive marker rows: ",
  nrow(cd8_markers_sig)
)
message(
  "Clusters without significant markers: ",
  length(clusters_without_markers)
)
message("Complete marker table saved.")
message("Significant marker table saved.")
message("Top-10 marker table saved.")
message("Checkpoint object saved.")
message("============================================================")

# ============================================================
# SECTION 13 — CD8 MARKER VALIDATION
# ============================================================

message("============================================================")
message("SECTION 13 — CD8 MARKER VALIDATION")
message("============================================================")


# ------------------------------------------------------------
# 13.1 Define biological marker panels
# ------------------------------------------------------------

marker_panels <- list(
  
  Conventional_T_CD8 = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2",
    "CD8A",
    "CD8B"
  ),
  
  Naive_Memory = c(
    "CCR7",
    "LTB",
    "IL7R",
    "MAL",
    "LTB",
    "TCF7",
    "LEF1",
    "MALAT1"
  ),
  
  Cytotoxic = c(
    "NKG7",
    "CCL5",
    "GNLY",
    "GZMB",
    "GZMH",
    "FGFBP2",
    "PRF1",
    "GZMK",
    "CX3CR1"
  ),
  
  NK = c(
    "KLRD1",
    "KLRF1",
    "KLRC1",
    "KLRC2",
    "KLRC3",
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
  
  Dysfunction_Exhaustion = c(
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
  
  CD4_Treg = c(
    "CD4",
    "IL7R",
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "TNFRSF4",
    "TNFRSF18",
    "ICOS"
  ),
  
  Proliferation = c(
    "MKI67",
    "TOP2A",
    "STMN1",
    "TYMS",
    "PCNA"
  )
)


# ------------------------------------------------------------
# 13.2 Keep only genes present in the object
# ------------------------------------------------------------

all_features <- rownames(cd8_obj[["RNA"]])

marker_panels_present <- lapply(
  marker_panels,
  function(x) intersect(x, all_features)
)

marker_panels_present <- marker_panels_present[
  lengths(marker_panels_present) > 0
]

message("")
message("Marker-panel coverage:")

for (panel_name in names(marker_panels_present)) {
  
  message(
    panel_name,
    ": ",
    length(marker_panels_present[[panel_name]]),
    " genes present"
  )
}


# ------------------------------------------------------------
# 13.3 Build ordered marker list
# ------------------------------------------------------------

marker_features <- unique(
  unlist(
    marker_panels_present,
    use.names = FALSE
  )
)

message("")
message(
  "Total unique validation markers present: ",
  length(marker_features)
)


# ------------------------------------------------------------
# 13.4 Generate DotPlot
# ------------------------------------------------------------

message("")
message("Generating CD8 marker-validation DotPlot...")

p_cd8_marker_dotplot <- Seurat::DotPlot(
  object = cd8_obj,
  assay = "RNA",
  features = marker_features,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  ggplot2::ggtitle(
    "CD8 subcluster marker validation"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

print(p_cd8_marker_dotplot)


# ------------------------------------------------------------
# 13.5 Save DotPlot
# ------------------------------------------------------------

figures_dir <- here::here(
  "results",
  "figures"
)

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "phase4_cd8_marker_validation_dotplot.png"
  ),
  plot = p_cd8_marker_dotplot,
  width = 15,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------
# 13.6 Save marker panels
# ------------------------------------------------------------

marker_panel_table <- dplyr::bind_rows(
  lapply(
    names(marker_panels_present),
    function(panel_name) {
      data.frame(
        panel = panel_name,
        gene = marker_panels_present[[panel_name]],
        stringsAsFactors = FALSE
      )
    }
  )
)

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  marker_panel_table,
  file = file.path(
    tables_dir,
    "phase4_cd8_marker_validation_panels.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13.7 SECTION 13 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 13 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("CD8 clusters: ", length(unique(cd8_obj$seurat_clusters)))
message(
  "Validation markers present: ",
  length(marker_features)
)
message("Marker-validation DotPlot saved.")
message("Marker-panel table saved.")
message("============================================================")

# ============================================================
# SECTION 14 — CD8 BIOLOGICAL ANNOTATION
# ============================================================

message("============================================================")
message("SECTION 14 — CD8 BIOLOGICAL ANNOTATION")
message("============================================================")


# ------------------------------------------------------------
# 14.1 Lock cluster-to-state annotations
# ------------------------------------------------------------

cd8_cluster_labels <- c(
  "0"  = "CD8_Memory_GPR183",
  "1"  = "CD8_Naive_Memory",
  "2"  = "CD8_PD1_Dysfunctional",
  "3"  = "NK_like",
  "4"  = "CD8_Activated",
  "5"  = "CD8_Immediate_Early_Stress",
  "6"  = "CD8_Memory_P2RY8",
  "7"  = "GammaDelta_T_like",
  "8"  = "NK_like",
  "9"  = "Ig_Associated_Ambiguous",
  "10" = "Treg_like",
  "11" = "GammaDelta_NK_like_TRM",
  "12" = "Cytotoxic_GammaDelta_like",
  "13" = "GammaDelta_Nonconventional",
  "14" = "GammaDelta_IL7R",
  "15" = "CD8_Cytotoxic_Effector",
  "16" = "NK_like_Nonconventional"
)


# ------------------------------------------------------------
# 14.2 Validate annotation coverage
# ------------------------------------------------------------

observed_clusters <- sort(
  unique(
    as.character(cd8_obj$seurat_clusters)
  )
)

annotated_clusters <- sort(
  names(cd8_cluster_labels)
)

if (!identical(
  observed_clusters,
  annotated_clusters
)) {
  stop(
    "Cluster annotation map does not exactly cover ",
    "the observed CD8 clusters."
  )
}


# ------------------------------------------------------------
# 14.3 Apply annotations
# ------------------------------------------------------------

cd8_obj$CD8_State <- unname(
  cd8_cluster_labels[
    as.character(cd8_obj$seurat_clusters)
  ]
)


# ------------------------------------------------------------
# 14.4 Validate annotations
# ------------------------------------------------------------

if (any(is.na(cd8_obj$CD8_State))) {
  stop(
    "Some CD8 cells received missing biological annotations."
  )
}

state_counts <- table(
  cd8_obj$CD8_State
)

message("")
message("CD8 biological-state counts:")
print(state_counts)


# ------------------------------------------------------------
# 14.5 Save cluster annotation table
# ------------------------------------------------------------

annotation_table <- data.frame(
  cluster = names(cd8_cluster_labels),
  CD8_State = unname(cd8_cluster_labels),
  stringsAsFactors = FALSE
)

annotation_table$n_cells <- as.integer(
  table(
    factor(
      as.character(cd8_obj$seurat_clusters),
      levels = names(cd8_cluster_labels)
    )
  )
)

annotation_table$fraction <- (
  annotation_table$n_cells /
    ncol(cd8_obj)
)


# ------------------------------------------------------------
# 14.6 Save annotation table
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  annotation_table,
  file = file.path(
    tables_dir,
    "phase4_cd8_biological_annotations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14.7 Record annotation framework
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Annotation <- list(
  annotation_basis =
    "Cluster marker expression and marker-panel validation",
  principle =
    "Conservative biological interpretation; no forced canonical CD8 state",
  cluster_labels =
    cd8_cluster_labels,
  note =
    "NK-like, gamma-delta-like, Treg-like and Ig-associated populations retained as observed subpopulations rather than removed."
)


# ------------------------------------------------------------
# 14.8 Save checkpoint
# ------------------------------------------------------------

checkpoint_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}

saveRDS(
  cd8_obj,
  file = file.path(
    checkpoint_dir,
    "phase4_cd8_annotated.rds"
  )
)


# ------------------------------------------------------------
# 14.9 SECTION 14 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 14 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message("Clusters annotated: ", length(cd8_cluster_labels))
message("Missing CD8 states: ", sum(is.na(cd8_obj$CD8_State)))
message("Annotation table saved.")
message("Annotated CD8 object saved.")
message("============================================================")

# ============================================================
# SECTION 15 — CD8 STATE COMPOSITION ACROSS CLINICAL STATES
# ============================================================

message("============================================================")
message("SECTION 15 — CD8 STATE × CLINICAL STATE COMPOSITION")
message("============================================================")


# ------------------------------------------------------------
# 15.1 Validate required metadata
# ------------------------------------------------------------

required_metadata <- c(
  "CD8_State",
  "Phase",
  "Donor"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(cd8_obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 15.2 Cell-level counts
# ------------------------------------------------------------

cell_counts <- cd8_obj@meta.data |>
  dplyr::count(
    Phase,
    CD8_State,
    name = "n_cells"
  )

message("")
message("CD8 state counts by clinical state:")
print(cell_counts)


# ------------------------------------------------------------
# 15.3 Within-phase composition
# ------------------------------------------------------------

phase_totals <- cd8_obj@meta.data |>
  dplyr::count(
    Phase,
    name = "phase_total"
  )

cell_composition <- cell_counts |>
  dplyr::left_join(
    phase_totals,
    by = "Phase"
  ) |>
  dplyr::mutate(
    proportion = n_cells / phase_total,
    percentage = 100 * proportion
  )

message("")
message("Within-clinical-state CD8 composition:")
print(cell_composition)


# ------------------------------------------------------------
# 15.4 Donor-level composition
# ------------------------------------------------------------

donor_state_counts <- cd8_obj@meta.data |>
  dplyr::count(
    Donor,
    Phase,
    CD8_State,
    name = "n_cells"
  )

donor_totals <- cd8_obj@meta.data |>
  dplyr::count(
    Donor,
    name = "donor_total"
  )

donor_composition <- donor_state_counts |>
  dplyr::left_join(
    donor_totals,
    by = "Donor"
  ) |>
  dplyr::mutate(
    proportion = n_cells / donor_total,
    percentage = 100 * proportion
  )

message("")
message(
  "Donor-level composition calculated for ",
  dplyr::n_distinct(donor_composition$Donor),
  " donors."
)


# ------------------------------------------------------------
# 15.5 Donor-level summary by clinical state
# ------------------------------------------------------------

donor_state_summary <- donor_composition |>
  dplyr::group_by(
    Phase,
    CD8_State
  ) |>
  dplyr::summarise(
    n_donors = dplyr::n_distinct(Donor),
    mean_percentage = mean(percentage),
    median_percentage = median(percentage),
    sd_percentage = sd(percentage),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 15.6 Save composition tables
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  cell_counts,
  file = file.path(
    tables_dir,
    "phase4_cd8_state_counts_by_phase.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_composition,
  file = file.path(
    tables_dir,
    "phase4_cd8_state_composition_by_phase.csv"
  ),
  row.names = FALSE
)

write.csv(
  donor_composition,
  file = file.path(
    tables_dir,
    "phase4_cd8_donor_level_composition.csv"
  ),
  row.names = FALSE
)

write.csv(
  donor_state_summary,
  file = file.path(
    tables_dir,
    "phase4_cd8_donor_state_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15.7 Composition heatmap
# ------------------------------------------------------------

composition_heatmap <- cell_composition |>
  dplyr::select(
    Phase,
    CD8_State,
    percentage
  ) |>
  tidyr::pivot_wider(
    names_from = Phase,
    values_from = percentage,
    values_fill = 0
  )

heatmap_matrix <- as.matrix(
  composition_heatmap[
    ,
    setdiff(
      colnames(composition_heatmap),
      "CD8_State"
    ),
    drop = FALSE
  ]
)

rownames(heatmap_matrix) <- composition_heatmap$CD8_State

heatmap_df <- as.data.frame(
  as.table(heatmap_matrix)
)

colnames(heatmap_df) <- c(
  "CD8_State",
  "Phase",
  "Percentage"
)

p_composition_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(
    x = Phase,
    y = CD8_State,
    fill = Percentage
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(
    name = "Percentage"
  ) +
  ggplot2::labs(
    title = "CD8 subpopulation composition across clinical states",
    x = "Clinical state",
    y = "CD8 state"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )


# ------------------------------------------------------------
# 15.8 Save heatmap
# ------------------------------------------------------------

figures_dir <- here::here(
  "results",
  "figures"
)

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "phase4_cd8_state_composition_heatmap.png"
  ),
  plot = p_composition_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

print(p_composition_heatmap)


# ------------------------------------------------------------
# 15.9 Record composition methodology
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Composition <- list(
  cell_level =
    "Counts and within-clinical-state proportions",
  replication_unit =
    "Donor",
  donor_level =
    "Donor-level state proportions summarized by clinical state",
  interpretation =
    "Cell-level composition is descriptive; donor-level summaries preserve biological replication."
)


# ------------------------------------------------------------
# 15.10 Save checkpoint
# ------------------------------------------------------------

checkpoint_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}

saveRDS(
  cd8_obj,
  file = file.path(
    checkpoint_dir,
    "phase4_cd8_composition.rds"
  )
)


# ------------------------------------------------------------
# 15.11 SECTION 15 CHECKPOINT
# ------------------------------------------------------------

message("")
message("============================================================")
message("SECTION 15 CHECKPOINT")
message("============================================================")
message("CD8 cells: ", ncol(cd8_obj))
message(
  "Clinical states: ",
  dplyr::n_distinct(cd8_obj$Phase)
)
message(
  "CD8 states: ",
  dplyr::n_distinct(cd8_obj$CD8_State)
)
message(
  "Donors: ",
  dplyr::n_distinct(cd8_obj$Donor)
)
message("Cell-level composition table saved.")
message("Donor-level composition table saved.")
message("Composition heatmap saved.")
message("Checkpoint object saved.")
message("============================================================")

# ============================================================
# SECTION 16 — FINAL CD8 BIOLOGICAL CHARACTERIZATION
#              AND PATHWAY ENRICHMENT
# ============================================================

message("============================================================")
message("SECTION 16 — FINAL CD8 BIOLOGICAL CHARACTERIZATION")
message("============================================================")


# ------------------------------------------------------------
# 16.1 Required packages
# ------------------------------------------------------------

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Required packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "."
  )
}


# ------------------------------------------------------------
# 16.2 Define biologically interpretable CD8 states
# ------------------------------------------------------------

selected_cd8_states <- c(
  "CD8_Memory_GPR183",
  "CD8_Naive_Memory",
  "CD8_PD1_Dysfunctional",
  "CD8_Activated",
  "CD8_Immediate_Early_Stress",
  "CD8_Memory_P2RY8",
  "CD8_Cytotoxic_Effector"
)

message("")
message("Selected conventional CD8 states:")
print(selected_cd8_states)


# ------------------------------------------------------------
# 16.3 Confirm states exist
# ------------------------------------------------------------

observed_states <- unique(
  as.character(cd8_obj$CD8_State)
)

missing_states <- setdiff(
  selected_cd8_states,
  observed_states
)

if (length(missing_states) > 0) {
  stop(
    "Selected CD8 states are missing: ",
    paste(missing_states, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 16.4 Retrieve significant positive markers
# ------------------------------------------------------------

if (!exists("cd8_markers_sig")) {
  
  cd8_markers_sig <- Seurat::FindAllMarkers(
    object = cd8_obj,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = TRUE
  )
}


# ------------------------------------------------------------
# 16.5 Map selected biological states to their clusters
# ------------------------------------------------------------

state_cluster_map <- data.frame(
  cluster = names(cd8_cluster_labels),
  CD8_State = unname(cd8_cluster_labels),
  stringsAsFactors = FALSE
)

selected_state_clusters <- state_cluster_map[
  state_cluster_map$CD8_State %in% selected_cd8_states,
  ,
  drop = FALSE
]

message("")
message("Selected-state cluster mapping:")
print(selected_state_clusters, row.names = FALSE)


# ------------------------------------------------------------
# 16.6 Extract significant markers for selected states
# ------------------------------------------------------------

selected_markers <- cd8_markers_sig[
  as.character(cd8_markers_sig$cluster) %in%
    selected_state_clusters$cluster,
  ,
  drop = FALSE
]

selected_markers <- selected_markers[
  order(
    as.numeric(
      as.character(selected_markers$cluster)
    ),
    -selected_markers$avg_log2FC,
    selected_markers$p_val_adj
  ),
  ,
  drop = FALSE
]

selected_markers <- dplyr::left_join(
  selected_markers,
  selected_state_clusters,
  by = "cluster"
)

message("")
message(
  "Significant markers across selected CD8 states: ",
  nrow(selected_markers)
)


# ------------------------------------------------------------
# 16.7 Save selected-state marker table
# ------------------------------------------------------------

tables_dir <- here::here(
  "results",
  "tables"
)

if (!dir.exists(tables_dir)) {
  dir.create(
    tables_dir,
    recursive = TRUE
  )
}

write.csv(
  selected_markers,
  file = file.path(
    tables_dir,
    "phase4_cd8_selected_state_markers.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 16.8 Top 15 markers per selected state
# ------------------------------------------------------------

top15_selected_markers <- selected_markers |>
  dplyr::group_by(CD8_State) |>
  dplyr::slice_head(n = 15) |>
  dplyr::ungroup()

write.csv(
  top15_selected_markers,
  file = file.path(
    tables_dir,
    "phase4_cd8_selected_state_top15_markers.csv"
  ),
  row.names = FALSE
)

message("")
message("Top markers by selected CD8 state:")
print(
  top15_selected_markers[
    ,
    c(
      "CD8_State",
      "gene",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj"
    )
  ],
  n = 105,
  row.names = FALSE
)

# ============================================================
# PHASE 4 — FINAL CD8 OBJECT SAVE
# ============================================================

message("============================================================")
message("PHASE 4 — FINAL CD8 OBJECT")
message("============================================================")


# ------------------------------------------------------------
# 1. Final integrity checks
# ------------------------------------------------------------

if (!inherits(cd8_obj, "Seurat")) {
  stop("cd8_obj is not a Seurat object.")
}

if (ncol(cd8_obj) != 44229) {
  stop(
    "Unexpected CD8 cell count: ",
    ncol(cd8_obj)
  )
}

if (!all(cd8_obj$Transferred_Label == "CD8_T")) {
  stop(
    "CD8 object contains cells outside the Phase 3 CD8 definition."
  )
}

if (any(is.na(cd8_obj$CD8_State))) {
  stop(
    "Missing CD8 biological-state annotations."
  )
}

if (any(is.na(cd8_obj$seurat_clusters))) {
  stop(
    "Missing CD8 cluster assignments."
  )
}

if (!"pca" %in% names(cd8_obj@reductions)) {
  stop("PCA reduction is missing.")
}

if (!"umap" %in% names(cd8_obj@reductions)) {
  stop("UMAP reduction is missing.")
}

if (!"RNA_snn" %in% names(cd8_obj@graphs)) {
  stop("RNA_snn graph is missing.")
}

if (length(Seurat::VariableFeatures(
  cd8_obj,
  assay = "RNA"
)) != 2000) {
  stop(
    "Expected exactly 2,000 CD8 HVGs."
  )
}


# ------------------------------------------------------------
# 2. Create final output directory
# ------------------------------------------------------------

final_rds_dir <- here::here(
  "results",
  "rds_objects"
)

if (!dir.exists(final_rds_dir)) {
  dir.create(
    final_rds_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Record final Phase 4 status
# ------------------------------------------------------------

cd8_obj@misc$Phase4_Status <- list(
  
  status =
    "COMPLETE",
  
  analysis_population =
    "CD8_T cells defined by Phase 3 transferred annotation",
  
  n_cells =
    ncol(cd8_obj),
  
  n_genes =
    nrow(cd8_obj),
  
  variable_features =
    length(
      Seurat::VariableFeatures(
        cd8_obj,
        assay = "RNA"
      )
    ),
  
  PCA =
    "50 PCs calculated",
  
  downstream_PCs =
    "1:20",
  
  clustering =
    "RNA_snn, resolution 0.4, algorithm 1",
  
  n_clusters =
    dplyr::n_distinct(
      cd8_obj$seurat_clusters
    ),
  
  UMAP =
    "2-dimensional UMAP using PCs 1:20",
  
  biological_states =
    dplyr::n_distinct(
      cd8_obj$CD8_State
    ),
  
  clinical_states =
    dplyr::n_distinct(
      cd8_obj$Phase
    ),
  
  donors =
    dplyr::n_distinct(
      cd8_obj$Donor
    ),
  
  marker_validation =
    "Completed",
  
  clinical_state_composition =
    "Completed",
  
  pathway_enrichment =
    "NOT performed; reserved for Phase 7",
  
  differential_expression =
    "NOT performed; reserved for Phase 6",
  
  cellchat =
    "NOT performed; reserved for Phase 8",
  
  phase4_endpoint =
    "Section 16.8"
)


# ------------------------------------------------------------
# 4. Save final frozen object
# ------------------------------------------------------------

final_cd8_path <- file.path(
  final_rds_dir,
  "phase4_final_cd8_analysis.rds"
)

saveRDS(
  cd8_obj,
  file = final_cd8_path
)


# ------------------------------------------------------------
# 5. Verify that the saved object can be reloaded
# ------------------------------------------------------------

message("")
message("Verifying saved RDS...")

cd8_final_check <- readRDS(
  final_cd8_path
)

if (!inherits(cd8_final_check, "Seurat")) {
  stop(
    "Saved RDS could not be reloaded as a Seurat object."
  )
}

if (ncol(cd8_final_check) != ncol(cd8_obj)) {
  stop(
    "Reloaded object has a different cell count."
  )
}

if (!identical(
  colnames(cd8_final_check),
  colnames(cd8_obj)
)) {
  stop(
    "Reloaded object cell order differs from the in-memory object."
  )
}

if (!identical(
  cd8_final_check$CD8_State,
  cd8_obj$CD8_State
)) {
  stop(
    "Reloaded biological-state annotations differ."
  )
}


# ------------------------------------------------------------
# 6. Final checkpoint
# ------------------------------------------------------------

message("")
message("============================================================")
message("PHASE 4 FINAL CHECKPOINT")
message("============================================================")
message("Status: COMPLETE")
message("CD8 cells: ", ncol(cd8_final_check))
message("Genes: ", nrow(cd8_final_check))
message(
  "CD8 HVGs: ",
  length(
    Seurat::VariableFeatures(
      cd8_final_check,
      assay = "RNA"
    )
  )
)
message(
  "Computational clusters: ",
  dplyr::n_distinct(
    cd8_final_check$seurat_clusters
  )
)
message(
  "Biological states: ",
  dplyr::n_distinct(
    cd8_final_check$CD8_State
  )
)
message(
  "Clinical states: ",
  dplyr::n_distinct(
    cd8_final_check$Phase
  )
)
message(
  "Donors: ",
  dplyr::n_distinct(
    cd8_final_check$Donor
  )
)
message("PCA: PRESENT")
message("UMAP: PRESENT")
message("RNA_snn: PRESENT")
message("Marker analysis: COMPLETE")
message("Clinical composition: COMPLETE")
message("DE analysis: RESERVED FOR PHASE 6")
message("Pathway enrichment: RESERVED FOR PHASE 7")
message("CellChat: RESERVED FOR PHASE 8")
message("")
message("FINAL RDS:")
message(final_cd8_path)
message("")
message("============================================================")
message("PHASE 4 IS OFFICIALLY WRAPPED. 🎉")
message("============================================================")

