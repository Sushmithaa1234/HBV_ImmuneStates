# ============================================================
# PHASE 6 — DONOR-AWARE DIFFERENTIAL-STATE ANALYSIS
#
# Biological unit:
#   Donor
#
# Inferential compartment:
#   CD8 T-cell populations only
#
# Myeloid:
#   Characterized in Phase 5.
#   Excluded from Phase 6 downstream inference because
#   donor/state representation is insufficient for robust
#   five-state analysis.
#
# Expression:
#   Processed log-CP10K signal.
#
# Aggregation:
#   Donor × CD8 population using aggregateNonCountSignal().
#
# Statistical model:
#   dreamlet donor-level differential-state analysis.
#
# Clinical states:
#   AC, AR, IA, IT, NL
#
# Primary CD8 populations:
#   Determined a priori from donor/state representation:
#   >=10 cells per donor-population
#   AND >=2 retained donors in every clinical state.
#
# Effect:
#   logFC as reported by limma/dreamlet on the modeled
#   expression scale.
#
# FDR:
#   Whole-dreamletResult correction for each contrast,
#   across the eight primary CD8 populations.
# ============================================================


# ============================================================
# SECTION 1 — SETUP
# ============================================================

suppressPackageStartupMessages({
  
  library(Seurat)
  library(SeuratObject)
  
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  
  library(dreamlet)
  library(BiocParallel)
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  
  library(ggplot2)
})

# ------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------

phase6_min_cells <- 10
phase6_primary_min_donors <- 2

phase6_states <- c(
  "AC",
  "AR",
  "IA",
  "IT",
  "NL"
)

phase6_model_states <- c(
  "NL",
  "AC",
  "AR",
  "IA",
  "IT"
)

phase6_contrasts <- c(
  "NL_vs_IT" = "Clinical_StateNL - Clinical_StateIT",
  "NL_vs_IA" = "Clinical_StateNL - Clinical_StateIA",
  "NL_vs_AR" = "Clinical_StateNL - Clinical_StateAR",
  "NL_vs_AC" = "Clinical_StateNL - Clinical_StateAC",
  "IT_vs_IA" = "Clinical_StateIT - Clinical_StateIA",
  "IT_vs_AR" = "Clinical_StateIT - Clinical_StateAR",
  "IT_vs_AC" = "Clinical_StateIT - Clinical_StateAC",
  "IA_vs_AR" = "Clinical_StateIA - Clinical_StateAR",
  "IA_vs_AC" = "Clinical_StateIA - Clinical_StateAC",
  "AR_vs_AC" = "Clinical_StateAR - Clinical_StateAC"
)

# ------------------------------------------------------------
# Input paths
# ------------------------------------------------------------

phase3_path <-
  "results/rds_objects/phase3_final_full_dataset.rds"

phase4_path <-
  "results/rds_objects/phase4_final_cd8_analysis_adjudicated.rds"

phase5_audit_path <-
  "results/rds_objects/phase6_sections1_to5_cd8_audit.rds"

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

out_rds <-
  "results/rds_objects"

out_tables <-
  "results/tables/phase6"

out_figures <-
  "results/figures/phase6"

dir.create(
  out_rds,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out_tables,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out_figures,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n============================================================\n")
cat("PHASE 6 — DONOR-AWARE DIFFERENTIAL-STATE ANALYSIS\n")
cat("============================================================\n")

cat(
  "Seurat:",
  as.character(packageVersion("Seurat")),
  "\n"
)

cat(
  "dreamlet:",
  as.character(packageVersion("dreamlet")),
  "\n"
)


# ============================================================
# SECTION 2 — LOAD AND VALIDATE INPUT OBJECTS
# ============================================================

phase3_obj <-
  readRDS(phase3_path)

cd8_obj <-
  readRDS(phase4_path)

# Required metadata

stopifnot(
  all(
    c(
      "GSM",
      "Donor",
      "Phase"
    ) %in%
      colnames(phase3_obj@meta.data)
  )
)

stopifnot(
  all(
    c(
      "GSM",
      "Donor",
      "Phase",
      "CD8_Cluster"
    ) %in%
      colnames(cd8_obj@meta.data)
  )
)

# Clinical states

observed_states <-
  sort(
    unique(
      as.character(
        phase3_obj$Phase
      )
    )
  )

stopifnot(
  setequal(
    observed_states,
    phase6_states
  )
)

# Each donor belongs to one clinical state

donor_state_check <-
  phase3_obj@meta.data %>%
  dplyr::transmute(
    Donor =
      as.character(Donor),
    Clinical_State =
      as.character(Phase)
  ) %>%
  dplyr::distinct() %>%
  dplyr::count(
    Donor,
    name = "n_states"
  )

stopifnot(
  all(
    donor_state_check$n_states == 1
  )
)


# ============================================================
# SECTION 3 — BUILD CD8 ANALYSIS METADATA
# ============================================================

cd8_metadata <-
  cd8_obj@meta.data %>%
  tibble::rownames_to_column(
    "Cell_ID"
  ) %>%
  dplyr::transmute(
    
    Cell_ID =
      as.character(Cell_ID),
    
    Population_ID =
      paste0(
        "CD8_Cluster_",
        as.character(CD8_Cluster)
      ),
    
    Donor_ID =
      as.character(Donor),
    
    GSM_ID =
      as.character(GSM),
    
    Clinical_State =
      factor(
        as.character(Phase),
        levels = phase6_states
      )
  )

stopifnot(
  nrow(cd8_metadata) ==
    dplyr::n_distinct(
      cd8_metadata$Cell_ID
    ),
  
  dplyr::n_distinct(
    cd8_metadata$Donor_ID
  ) == 23,
  
  dplyr::n_distinct(
    cd8_metadata$Clinical_State
  ) == 5
)

write_csv(
  cd8_metadata,
  file.path(
    out_tables,
    "phase6_cd8_cell_metadata.csv"
  )
)


# ============================================================
# SECTION 4 — DONOR / STATE REPRESENTATION
#
# A donor-population is retained for primary inference when:
#
#   >=10 cells from that population are present in the donor
#
# A population is primary when:
#
#   >=2 retained donors are present in EVERY clinical state.
# ============================================================

donor_population_counts <-
  cd8_metadata %>%
  dplyr::count(
    Population_ID,
    Donor_ID,
    Clinical_State,
    name = "Cells"
  )

state_coverage <-
  donor_population_counts %>%
  dplyr::group_by(
    Population_ID
  ) %>%
  dplyr::summarise(
    
    Donors_10plus =
      dplyr::n_distinct(
        Donor_ID[
          Cells >= phase6_min_cells
        ]
      ),
    
    AC =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "AC" &
            Cells >= phase6_min_cells
        ]
      ),
    
    AR =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "AR" &
            Cells >= phase6_min_cells
        ]
      ),
    
    IA =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "IA" &
            Cells >= phase6_min_cells
        ]
      ),
    
    IT =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "IT" &
            Cells >= phase6_min_cells
        ]
      ),
    
    NL =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "NL" &
            Cells >= phase6_min_cells
        ]
      ),
    
    Total_Cells =
      sum(Cells),
    
    Largest_Donor_Fraction =
      max(Cells) / sum(Cells),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    
    Minimum_State_Donors =
      pmin(
        AC,
        AR,
        IA,
        IT,
        NL
      ),
    
    Primary_Eligible =
      Minimum_State_Donors >=
      phase6_primary_min_donors
  )

write_csv(
  donor_population_counts,
  file.path(
    out_tables,
    "phase6_cd8_donor_population_counts.csv"
  )
)

write_csv(
  state_coverage,
  file.path(
    out_tables,
    "phase6_cd8_state_coverage.csv"
  )
)


# ============================================================
# SECTION 5 — FINAL CD8 POPULATION ADJUDICATION
# ============================================================

cd8_population_eligibility <-
  state_coverage %>%
  
  dplyr::mutate(
    
    Phase6_Role =
      dplyr::case_when(
        
        Primary_Eligible ~
          "PRIMARY",
        
        Minimum_State_Donors == 1 ~
          "SECONDARY",
        
        TRUE ~
          "EXCLUDED"
      ),
    
    Donor_Concentration =
      dplyr::case_when(
        
        Largest_Donor_Fraction >= 0.90 ~
          "Extreme",
        
        Largest_Donor_Fraction >= 0.75 ~
          "Strong",
        
        Largest_Donor_Fraction >= 0.50 ~
          "Moderate",
        
        TRUE ~
          "Low"
      )
  ) %>%
  
  dplyr::arrange(
    Population_ID
  )

write_csv(
  cd8_population_eligibility,
  file.path(
    out_tables,
    "phase6_cd8_population_eligibility.csv"
  )
)

primary_populations <-
  cd8_population_eligibility %>%
  dplyr::filter(
    Phase6_Role == "PRIMARY"
  ) %>%
  dplyr::pull(
    Population_ID
  )

stopifnot(
  length(primary_populations) == 8
)


# ============================================================
# SECTION 6 — BUILD PRIMARY CD8 SCE
# ============================================================

primary_cells <-
  cd8_metadata %>%
  dplyr::filter(
    Population_ID %in%
      primary_populations
  ) %>%
  dplyr::pull(
    Cell_ID
  )

cd8_primary <-
  subset(
    cd8_obj,
    cells = primary_cells
  )

# Join only the processed data layers

cd8_primary <-
  JoinLayers(
    cd8_primary,
    assay = "RNA",
    layers = "data"
  )

data_layers <-
  Layers(
    cd8_primary[["RNA"]],
    search = "^data"
  )

stopifnot(
  length(data_layers) == 1
)

# Extract processed log-CP10K signal

logCP10K_matrix <-
  GetAssayData(
    cd8_primary,
    assay = "RNA",
    layer = "data"
  )

stopifnot(
  all(
    is.finite(
      logCP10K_matrix@x
    )
  ),
  
  all(
    logCP10K_matrix@x >= 0
  )
)


# ============================================================
# SECTION 7 — CONSTRUCT SINGLECELL EXPERIMENT
# ============================================================

sce_cd8 <-
  SingleCellExperiment(
    assays =
      list(
        logCP10K =
          logCP10K_matrix
      ),
    
    colData =
      S4Vectors::DataFrame(
        
        Cell_ID =
          colnames(cd8_primary),
        
        Donor_ID =
          as.character(
            cd8_primary$Donor
          ),
        
        GSM_ID =
          as.character(
            cd8_primary$GSM
          ),
        
        Clinical_State =
          factor(
            as.character(
              cd8_primary$Phase
            ),
            levels =
              phase6_states
          ),
        
        Population_ID =
          paste0(
            "CD8_Cluster_",
            as.character(
              cd8_primary$CD8_Cluster
            )
          )
      )
  )

stopifnot(
  
  ncol(sce_cd8) ==
    ncol(cd8_primary),
  
  nrow(sce_cd8) ==
    nrow(cd8_primary),
  
  dplyr::n_distinct(
    sce_cd8$Donor_ID
  ) == 23,
  
  dplyr::n_distinct(
    sce_cd8$Population_ID
  ) == 8,
  
  dplyr::n_distinct(
    sce_cd8$Clinical_State
  ) == 5
)


# ============================================================
# SECTION 8 — DONOR-LEVEL NON-COUNT AGGREGATION
#
# Averages processed expression within:
#
#   donor × CD8 population
#
# Precision weights are retained by dreamlet.
# No reverse transformation.
# No pseudo-count reconstruction.
# ============================================================

pb_cd8 <-
  aggregateNonCountSignal(
    
    sce_cd8,
    
    assay =
      "logCP10K",
    
    sample_id =
      "Donor_ID",
    
    cluster_id =
      "Population_ID",
    
    min.cells =
      phase6_min_cells,
    
    min.signal =
      0.01,
    
    min.samples =
      4,
    
    min.prop =
      0.4,
    
    verbose =
      TRUE
  )

stopifnot(
  inherits(
    pb_cd8,
    "dreamletProcessedData"
  )
)

pb_path <-
  file.path(
    out_rds,
    "phase6_cd8_donor_population_noncount_aggregated.rds"
  )

saveRDS(
  pb_cd8,
  pb_path
)

# Save aggregation metadata

sce_cell_audit <-
  as.data.frame(
    colData(sce_cd8)
  ) %>%
  dplyr::count(
    Population_ID,
    Donor_ID,
    Clinical_State,
    name = "Cells"
  )

sce_threshold_audit <-
  sce_cell_audit %>%
  dplyr::group_by(
    Population_ID
  ) %>%
  dplyr::summarise(
    
    Donors_10plus =
      dplyr::n_distinct(
        Donor_ID[
          Cells >=
            phase6_min_cells
        ]
      ),
    
    AC =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "AC" &
            Cells >=
            phase6_min_cells
        ]
      ),
    
    AR =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "AR" &
            Cells >=
            phase6_min_cells
        ]
      ),
    
    IA =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "IA" &
            Cells >=
            phase6_min_cells
        ]
      ),
    
    IT =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "IT" &
            Cells >=
            phase6_min_cells
        ]
      ),
    
    NL =
      dplyr::n_distinct(
        Donor_ID[
          Clinical_State == "NL" &
            Cells >=
            phase6_min_cells
        ]
      ),
    
    .groups = "drop"
  )

write_csv(
  sce_cell_audit,
  file.path(
    out_tables,
    "phase6_cd8_sce_donor_population_counts.csv"
  )
)

write_csv(
  sce_threshold_audit,
  file.path(
    out_tables,
    "phase6_cd8_sce_threshold_audit.csv"
  )
)


# ============================================================
# SECTION 9 — FINAL DREAMLET MODEL + CONTRASTS
#
# Cell-means parameterization:
#
#   ~ 0 + Clinical_State
#
# This estimates one state mean per clinical state.
#
# The ten planned pairwise contrasts are then tested directly.
#
# Positive logFC:
#   first state > second state
#
# Negative logFC:
#   first state < second state
# ============================================================

model_data <-
  as.data.frame(
    colData(pb_cd8)
  )

model_data$Clinical_State <-
  factor(
    as.character(
      model_data$Clinical_State
    ),
    levels =
      phase6_model_states
  )

stopifnot(
  !anyNA(
    model_data$Clinical_State
  ),
  
  all(
    table(
      model_data$Clinical_State
    ) > 0
  )
)

# Ensure exact primary assays

cd8_assays <-
  primary_populations

stopifnot(
  identical(
    sort(
      assayNames(pb_cd8)
    ),
    sort(
      cd8_assays
    )
  )
)

# Fit final model

contrast_cd8 <-
  dreamlet(
    
    pb_cd8,
    
    ~ 0 + Clinical_State,
    
    data =
      model_data,
    
    assays =
      cd8_assays,
    
    contrasts =
      phase6_contrasts,
    
    quiet =
      FALSE,
    
    BPPARAM =
      SerialParam()
  )

stopifnot(
  inherits(
    contrast_cd8,
    "dreamletResult"
  )
)

# Model details

model_details <-
  details(
    contrast_cd8
  )

stopifnot(
  all(
    model_details$n_errors == 0
  ),
  
  all(
    !model_details$error_initial
  )
)

write_csv(
  as.data.frame(
    model_details
  ),
  file.path(
    out_tables,
    "phase6_cd8_dreamlet_model_details.csv"
  )
)


# ============================================================
# SECTION 10 — FINAL GENE-LEVEL RESULTS
#
# IMPORTANT:
#
# Results are extracted from the WHOLE dreamletResult.
#
# Therefore adj.P.Val is based on the full set of CD8
# populations tested for each contrast, rather than applying
# FDR independently within each population.
# ============================================================

result_list <-
  vector(
    "list",
    length(phase6_contrasts)
  )

names(result_list) <-
  names(phase6_contrasts)
for (
  contrast_name in
  names(phase6_contrasts)
) {
  
  cat(
    "Extracting:",
    contrast_name,
    "\n"
  )
  
  tab <-
    topTable(
      contrast_cd8,
      coef =
        contrast_name,
      number =
        Inf
    )
  
  tab <-
    as.data.frame(
      tab
    )
  
  tab <-
    tibble::as_tibble(
      tab
    )
  
  tab$Contrast <-
    contrast_name
  
  result_list[[
    contrast_name
  ]] <-
    tab
}

gene_results <-
  dplyr::bind_rows(
    result_list
  )

# Standardize column order

gene_results <-
  gene_results %>%
  dplyr::select(
    assay,
    ID,
    Contrast,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B,
    z.std
  )

gene_results <-
  gene_results %>%
  dplyr::rename(
    Population_ID =
      assay,
    Gene =
      ID
  )


# ============================================================
# SECTION 11 — FINAL RESULT TABLE CHECKS
# ============================================================

expected_rows <-
  sum(
    vapply(
      names(phase6_contrasts),
      function(x) {
        
        nrow(
          topTable(
            contrast_cd8,
            coef = x,
            number = Inf
          )
        )
        
      },
      numeric(1)
    )
  )

stopifnot(
  nrow(gene_results) ==
    expected_rows,
  
  all(
    gene_results$Population_ID
    %in%
      cd8_assays
  ),
  
  all(
    gene_results$Contrast
    %in%
      names(
        phase6_contrasts
      )
  ),
  
  all(
    is.finite(
      gene_results$logFC
    )
  ),
  
  all(
    is.finite(
      gene_results$P.Value
    )
  ),
  
  all(
    is.finite(
      gene_results$adj.P.Val
    )
  ),
  
  all(
    gene_results$adj.P.Val >= 0 &
      gene_results$adj.P.Val <= 1
  )
)

# Number of significant gene × population results

significant_summary <-
  gene_results %>%
  dplyr::group_by(
    Population_ID,
    Contrast
  ) %>%
  dplyr::summarise(
    
    n_genes =
      dplyr::n(),
    
    n_FDR_05 =
      sum(
        adj.P.Val < 0.05
      ),
    
    n_FDR_01 =
      sum(
        adj.P.Val < 0.01
      ),
    
    max_abs_logFC =
      max(
        abs(logFC),
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )

write_csv(
  gene_results,
  file.path(
    out_tables,
    "phase6_cd8_all_gene_level_results.csv"
  )
)

write_csv(
  significant_summary,
  file.path(
    out_tables,
    "phase6_cd8_significant_gene_summary.csv"
  )
)


# ============================================================
# SECTION 12 — FINAL PHASE 6 OBJECT
# ============================================================

phase6_final <-
  list(
    
    dreamlet_result =
      contrast_cd8,
    
    aggregated_data =
      pb_cd8,
    
    model_metadata =
      model_data,
    
    cd8_metadata =
      cd8_metadata,
    
    donor_population_counts =
      donor_population_counts,
    
    state_coverage =
      state_coverage,
    
    population_eligibility =
      cd8_population_eligibility,
    
    parameters =
      list(
        
        min_cells_per_donor =
          phase6_min_cells,
        
        primary_min_donors_per_state =
          phase6_primary_min_donors,
        
        clinical_states =
          phase6_states,
        
        model_state_order =
          phase6_model_states,
        
        primary_populations =
          primary_populations,
        
        contrasts =
          phase6_contrasts,
        
        aggregation =
          "aggregateNonCountSignal",
        
        expression_input =
          "processed log-CP10K",
        
        inference_unit =
          "donor",
        
        downstream_myeloid_analysis =
          FALSE,
        
        effect_measure =
          "dreamlet/limma logFC on modeled expression scale"
      )
  )

final_rds_path <-
  file.path(
    out_rds,
    "phase6_cd8_final_differential_state_analysis.rds"
  )

saveRDS(
  phase6_final,
  final_rds_path
)

# Also save the direct dreamletResult for convenient reuse

saveRDS(
  contrast_cd8,
  file.path(
    out_rds,
    "phase6_cd8_dreamlet_contrasts.rds"
  )
)


# ============================================================
# SECTION 13 — FIGURE 1
# CD8 POPULATION ELIGIBILITY / DONOR COVERAGE
# ============================================================

coverage_long <-
  state_coverage %>%
  
  dplyr::select(
    Population_ID,
    AC,
    AR,
    IA,
    IT,
    NL
  ) %>%
  
  tidyr::pivot_longer(
    cols =
      c(
        AC,
        AR,
        IA,
        IT,
        NL
      ),
    names_to =
      "Clinical_State",
    values_to =
      "Retained_Donors"
  )

fig_coverage <-
  ggplot(
    coverage_long,
    aes(
      x =
        Clinical_State,
      y =
        Population_ID,
      fill =
        Retained_Donors
    )
  ) +
  
  geom_tile(
    color =
      "white"
  ) +
  
  geom_text(
    aes(
      label =
        Retained_Donors
    ),
    size =
      3
  ) +
  
  scale_fill_continuous(
    name =
      "Donors\n>=10 cells"
  ) +
  
  labs(
    title =
      "CD8 donor representation for Phase 6",
    x =
      "Clinical state",
    y =
      "CD8 population"
  ) +
  
  theme_minimal(
    base_size =
      11
  )

ggsave(
  file.path(
    out_figures,
    "phase6_cd8_donor_state_coverage.png"
  ),
  fig_coverage,
  width =
    8,
  height =
    5.5,
  dpi =
    300
)


# ============================================================
# SECTION 14 — FIGURE 2
# NUMBER OF FDR-SIGNIFICANT GENE × POPULATION RESULTS
# ============================================================

fig_fdr <-
  significant_summary %>%
  
  ggplot(
    aes(
      x =
        Contrast,
      y =
        Population_ID,
      fill =
        n_FDR_05
    )
  ) +
  
  geom_tile(
    color =
      "white"
  ) +
  
  geom_text(
    aes(
      label =
        n_FDR_05
    ),
    size =
      2.7
  ) +
  
  scale_fill_continuous(
    name =
      "Genes\nFDR < 0.05"
  ) +
  
  labs(
    title =
      "Differential-state signal across CD8 populations",
    x =
      "Clinical-state contrast",
    y =
      "CD8 population"
  ) +
  
  theme_minimal(
    base_size =
      10
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )

ggsave(
  file.path(
    out_figures,
    "phase6_cd8_fdr05_gene_counts.png"
  ),
  fig_fdr,
  width =
    12,
  height =
    5.5,
  dpi =
    300
)


# ============================================================
# SECTION 15 — FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PHASE 6 COMPLETE\n")
cat("============================================================\n")

cat(
  "\nPrimary CD8 populations:",
  length(primary_populations),
  "\n"
)

cat(
  "Final gene-level result rows:",
  nrow(gene_results),
  "\n"
)

cat(
  "\nFinal RDS:\n",
  final_rds_path,
  "\n"
)

cat(
  "\nFinal gene-level table:\n",
  file.path(
    out_tables,
    "phase6_cd8_all_gene_level_results.csv"
  ),
  "\n"
)

cat(
  "\nSignificant-gene summary:\n",
  file.path(
    out_tables,
    "phase6_cd8_significant_gene_summary.csv"
  ),
  "\n"
)

cat(
  "\nFigures:\n",
  file.path(
    out_figures,
    "phase6_cd8_donor_state_coverage.png"
  ),
  "\n",
  file.path(
    out_figures,
    "phase6_cd8_fdr05_gene_counts.png"
  ),
  "\n"
)

cat("\n============================================================\n")

