# ============================================================
# PHASE 7.01 — DIFFERENTIAL EXPRESSION RESULTS AUDIT
# Audits Phase 6 exported results and checks the prominent
# CD8_Cluster_2 NL-vs-IA signal against donor-level estimates.
#
# This script does NOT refit differential-expression models.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

input_dir <- "results/tables/Phase 6"
table_dir <- "results/tables/Phase 7"
figure_dir <- "results/figures/Phase 7"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- file.path(input_dir, "phase6_cd8_all_gene_level_results.csv")
summary_file <- file.path(input_dir, "phase6_cd8_significant_gene_summary.csv")
model_file <- file.path(input_dir, "phase6_cd8_dreamlet_model_details.csv")
eligibility_file <- file.path(input_dir, "phase6_cd8_population_eligibility.csv")

aggregated_rds_file <- file.path(
  "results/rds_objects",
  "phase6_cd8_donor_population_noncount_aggregated.rds"
)

required_files <- c(
  gene_file, summary_file, model_file,
  eligibility_file, aggregated_rds_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. LOAD PHASE 6 TABLES
# ------------------------------------------------------------

gene_results <- read_csv(gene_file, show_col_types = FALSE)
reported_summary <- read_csv(summary_file, show_col_types = FALSE)
model_details <- read_csv(model_file, show_col_types = FALSE)
population_eligibility <- read_csv(eligibility_file, show_col_types = FALSE)

required_gene_columns <- c(
  "Population_ID", "Gene", "Contrast", "logFC", "AveExpr",
  "t", "P.Value", "adj.P.Val", "B", "z.std"
)
required_summary_columns <- c(
  "Population_ID", "Contrast", "n_genes",
  "n_FDR_05", "n_FDR_01", "max_abs_logFC"
)
required_model_columns <- c(
  "assay", "n_retain", "formula", "n_genes",
  "n_errors", "error_initial"
)
required_eligibility_columns <- c(
  "Population_ID", "AC", "AR", "IA", "IT", "NL",
  "Total_Cells", "Largest_Donor_Fraction",
  "Minimum_State_Donors", "Donors_10plus",
  "Primary_Eligible", "Phase6_Role", "Donor_Concentration"
)

check_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  }
}

check_columns(gene_results, required_gene_columns, "Gene-level results")
check_columns(reported_summary, required_summary_columns, "Reported summary")
check_columns(model_details, required_model_columns, "Model details")
check_columns(
  population_eligibility,
  required_eligibility_columns,
  "Population eligibility"
)

expected_contrasts <- c(
  "NL_vs_IT", "NL_vs_IA", "NL_vs_AR", "NL_vs_AC",
  "IT_vs_IA", "IT_vs_AR", "IT_vs_AC",
  "IA_vs_AR", "IA_vs_AC", "AR_vs_AC"
)
expected_states <- c("AC", "AR", "IA", "IT", "NL")

primary_populations <- population_eligibility %>%
  dplyr::filter(as.logical(Primary_Eligible)) %>%
  dplyr::pull(Population_ID) %>%
  unique() %>%
  sort()

if (length(primary_populations) != 8) {
  stop("Expected eight primary CD8 populations; found ",
       length(primary_populations), ".")
}
if (!setequal(unique(gene_results$Contrast), expected_contrasts)) {
  stop("Gene-level results do not contain the expected ten contrasts.")
}
if (!setequal(unique(gene_results$Population_ID), primary_populations)) {
  stop("Gene-level result populations do not match primary eligibility.")
}

expected_combinations <- tidyr::expand_grid(
  Population_ID = primary_populations,
  Contrast = expected_contrasts
)
observed_combinations <- gene_results %>%
  dplyr::distinct(Population_ID, Contrast)

missing_combinations <- dplyr::anti_join(
  expected_combinations, observed_combinations,
  by = c("Population_ID", "Contrast")
)
unexpected_combinations <- dplyr::anti_join(
  observed_combinations, expected_combinations,
  by = c("Population_ID", "Contrast")
)

if (nrow(missing_combinations) || nrow(unexpected_combinations)) {
  stop("Population-by-contrast completeness check failed.")
}

# ------------------------------------------------------------
# 3. RESULT-TABLE INTEGRITY
# ------------------------------------------------------------

duplicate_gene_keys <- gene_results %>%
  dplyr::count(Population_ID, Contrast, Gene, name = "n") %>%
  dplyr::filter(n > 1)

missing_value_counts <- gene_results %>%
  dplyr::summarise(dplyr::across(
    dplyr::all_of(required_gene_columns),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(
    dplyr::everything(),
    names_to = "Column",
    values_to = "Missing_N"
  )

invalid_statistics <- gene_results %>%
  dplyr::summarise(
    Nonfinite_logFC = sum(!is.finite(logFC)),
    Nonfinite_AveExpr = sum(!is.finite(AveExpr)),
    Nonfinite_t = sum(!is.finite(t)),
    Nonfinite_P_Value = sum(!is.finite(P.Value)),
    Nonfinite_adj_P_Val = sum(!is.finite(adj.P.Val)),
    Nonfinite_B = sum(!is.finite(B)),
    Nonfinite_z_std = sum(!is.finite(z.std)),
    P_out_of_range = sum(P.Value < 0 | P.Value > 1, na.rm = TRUE),
    adjP_out_of_range = sum(adj.P.Val < 0 | adj.P.Val > 1, na.rm = TRUE)
  )

# Recalculate the reported significance summary from gene-level results.
recalculated_summary <- gene_results %>%
  dplyr::group_by(Population_ID, Contrast) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_FDR_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_FDR_01 = sum(adj.P.Val < 0.01, na.rm = TRUE),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  )

summary_comparison <- recalculated_summary %>%
  dplyr::rename(
    n_genes_recalculated = n_genes,
    n_FDR_05_recalculated = n_FDR_05,
    n_FDR_01_recalculated = n_FDR_01,
    max_abs_logFC_recalculated = max_abs_logFC
  ) %>%
  dplyr::full_join(
    reported_summary %>%
      dplyr::select(
        Population_ID, Contrast,
        n_genes_reported = n_genes,
        n_FDR_05_reported = n_FDR_05,
        n_FDR_01_reported = n_FDR_01,
        max_abs_logFC_reported = max_abs_logFC
      ),
    by = c("Population_ID", "Contrast")
  ) %>%
  dplyr::mutate(
    n_genes_match = n_genes_recalculated == n_genes_reported,
    n_FDR_05_match = n_FDR_05_recalculated == n_FDR_05_reported,
    n_FDR_01_match = n_FDR_01_recalculated == n_FDR_01_reported,
    max_abs_logFC_match = dplyr::near(
      max_abs_logFC_recalculated,
      max_abs_logFC_reported,
      tol = 1e-8
    ),
    all_checks_pass = dplyr::coalesce(n_genes_match, FALSE) &
      dplyr::coalesce(n_FDR_05_match, FALSE) &
      dplyr::coalesce(n_FDR_01_match, FALSE) &
      dplyr::coalesce(max_abs_logFC_match, FALSE)
  )

# ------------------------------------------------------------
# 4. MODEL AND DONOR-REPRESENTATION AUDIT
# ------------------------------------------------------------

model_audit <- model_details %>%
  dplyr::mutate(
    error_initial_logical = as.logical(error_initial),
    model_error_free = n_errors == 0 &
      !error_initial_logical,
    formula_includes_clinical_state = grepl(
      "Clinical_State", formula, fixed = TRUE
    )
  )

model_summary <- model_audit %>%
  dplyr::summarise(
    n_models = dplyr::n(),
    n_error_free = sum(model_error_free, na.rm = TRUE),
    min_retained_donors = min(n_retain, na.rm = TRUE),
    max_retained_donors = max(n_retain, na.rm = TRUE),
    min_genes = min(n_genes, na.rm = TRUE),
    max_genes = max(n_genes, na.rm = TRUE)
  )

donor_coverage_long <- population_eligibility %>%
  dplyr::filter(Population_ID %in% primary_populations) %>%
  dplyr::select(Population_ID, dplyr::all_of(expected_states)) %>%
  pivot_longer(
    cols = dplyr::all_of(expected_states),
    names_to = "Clinical_State",
    values_to = "Qualifying_Donors"
  ) %>%
  dplyr::mutate(
    Population_ID = factor(Population_ID, levels = primary_populations),
    Clinical_State = factor(Clinical_State, levels = expected_states)
  )

donor_concentration <- population_eligibility %>%
  dplyr::filter(Population_ID %in% primary_populations) %>%
  dplyr::select(
    Population_ID, Total_Cells, Largest_Donor_Fraction,
    Minimum_State_Donors, Donors_10plus,
    AC, AR, IA, IT, NL, Donor_Concentration
  ) %>%
  dplyr::arrange(dplyr::desc(Largest_Donor_Fraction))

# ------------------------------------------------------------
# 5. OVERALL RESULT SUMMARIES
# ------------------------------------------------------------

contrast_summary <- gene_results %>%
  dplyr::group_by(Contrast) %>%
  dplyr::summarise(
    n_populations = n_distinct(Population_ID),
    n_gene_results = dplyr::n(),
    n_FDR_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_FDR_01 = sum(adj.P.Val < 0.01, na.rm = TRUE),
    proportion_FDR_05 = mean(adj.P.Val < 0.05, na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Contrast = factor(Contrast, levels = expected_contrasts)) %>%
  dplyr::arrange(Contrast)

population_summary <- gene_results %>%
  dplyr::group_by(Population_ID) %>%
  dplyr::summarise(
    n_contrasts = n_distinct(Contrast),
    n_gene_results = dplyr::n(),
    n_FDR_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_FDR_01 = sum(adj.P.Val < 0.01, na.rm = TRUE),
    proportion_FDR_05 = mean(adj.P.Val < 0.05, na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  )

t_vs_zstd <- gene_results %>%
  dplyr::summarise(
    correlation_t_zstd = cor(t, z.std, use = "complete.obs"),
    maximum_absolute_difference = max(abs(t - z.std), na.rm = TRUE),
    identical_within_tolerance = all(abs(t - z.std) < 1e-10)
  )

# ------------------------------------------------------------
# 6. TARGETED DONOR-LEVEL CHECK
#    CD8_Cluster_2, NL_vs_IA; top ten significant genes
# ------------------------------------------------------------

cluster2_target <- gene_results %>%
  dplyr::filter(
    Population_ID == "CD8_Cluster_2",
    Contrast == "NL_vs_IA",
    is.finite(adj.P.Val),
    adj.P.Val < 0.05
  ) %>%
  dplyr::arrange(adj.P.Val, dplyr::desc(abs(t))) %>%
  dplyr::slice_head(n = 10)

if (nrow(cluster2_target) == 0) {
  stop("No significant genes found for CD8_Cluster_2, NL_vs_IA.")
}

pb <- readRDS(aggregated_rds_file)

if (!"CD8_Cluster_2" %in% names(pb)) {
  stop("CD8_Cluster_2 is absent from the aggregated dreamlet object.")
}

cluster2_elist <- pb[["CD8_Cluster_2"]]
E <- cluster2_elist$E
W <- cluster2_elist$weights

if (!identical(dim(E), dim(W)) ||
    !identical(dimnames(E), dimnames(W))) {
  stop("Cluster 2 expression and precision-weight matrices do not align.")
}

if (!all(cluster2_target$Gene %in% rownames(E))) {
  stop("One or more selected genes are absent from the Cluster 2 E matrix.")
}

sample_metadata <- as.data.frame(pb@data)
sample_metadata$Donor_ID <- rownames(sample_metadata)

donors <- colnames(E)
donor_mapping <- sample_metadata[
  match(donors, sample_metadata$Donor_ID),
  ,
  drop = FALSE
]

if (anyNA(donor_mapping$Clinical_State) ||
    anyNA(donor_mapping$GSM_ID)) {
  stop("At least one Cluster 2 donor could not be mapped to sample metadata.")
}

if (anyDuplicated(donors)) {
  stop("Duplicate donor identifiers found in the Cluster 2 expression matrix.")
}

selected_genes <- cluster2_target$Gene

donor_expression <- do.call(rbind, lapply(selected_genes, function(gene) {
  gene_index <- match(gene, rownames(E))
  
  data.frame(
    Gene = gene,
    Donor_ID = donors,
    GSM_ID = donor_mapping$GSM_ID,
    Clinical_State = donor_mapping$Clinical_State,
    Expression = as.numeric(E[gene_index, ]),
    Precision_Weight = as.numeric(W[gene_index, ]),
    stringsAsFactors = FALSE
  )
})) %>%
  dplyr::filter(Clinical_State %in% c("NL", "IA")) %>%
  dplyr::arrange(Gene, Clinical_State, Donor_ID)

donor_level_summary <- do.call(rbind, lapply(selected_genes, function(gene) {
  d <- donor_expression %>% dplyr::filter(Gene == gene)
  nl <- d %>% dplyr::filter(Clinical_State == "NL")
  ia <- d %>% dplyr::filter(Clinical_State == "IA")
  phase6_row <- cluster2_target %>% dplyr::filter(Gene == gene)
  
  nl_mean <- mean(nl$Expression)
  ia_mean <- mean(ia$Expression)
  observed_difference <- nl_mean - ia_mean
  reported_logFC <- phase6_row$logFC[1]
  
  data.frame(
    Gene = gene,
    NL_n = nrow(nl),
    IA_n = nrow(ia),
    NL_mean = nl_mean,
    IA_mean = ia_mean,
    NL_median = median(nl$Expression),
    IA_median = median(ia$Expression),
    NL_minus_IA = observed_difference,
    Reported_logFC = reported_logFC,
    Direction_matches = sign(observed_difference) == sign(reported_logFC),
    NL_max = max(nl$Expression),
    IA_min = min(ia$Expression),
    IA_max = max(ia$Expression),
    NL_min = min(nl$Expression),
    Complete_separation = max(nl$Expression) < min(ia$Expression) ||
      min(nl$Expression) > max(ia$Expression)
  )
}))

# ------------------------------------------------------------
# 7. AUDIT CHECKLIST
# ------------------------------------------------------------

audit_checks <- tibble::tibble(
  Check = c(
    "Eight expected primary populations",
    "Ten expected clinical contrasts",
    "All 80 population-contrast combinations present",
    "No duplicate gene/population/contrast keys",
    "No missing values in required gene-result fields",
    "No nonfinite or out-of-range statistical values",
    "Reported summary reconciles with gene-level results",
    "All eight models are error-free",
    "All model formulas include Clinical_State",
    "Cluster 2 donor IDs map to sample metadata",
    "Top-ten Cluster 2 directions agree with Phase 6 logFC"
  ),
  Passed = c(
    length(primary_populations) == 8,
    setequal(unique(gene_results$Contrast), expected_contrasts),
    nrow(missing_combinations) == 0 &&
      nrow(unexpected_combinations) == 0,
    nrow(duplicate_gene_keys) == 0,
    sum(missing_value_counts$Missing_N) == 0,
    sum(unlist(invalid_statistics)) == 0,
    all(summary_comparison$all_checks_pass),
    all(model_audit$model_error_free),
    all(model_audit$formula_includes_clinical_state),
    !anyNA(donor_mapping$Clinical_State),
    all(donor_level_summary$Direction_matches)
  )
)

overall_pass <- all(audit_checks$Passed)

# ------------------------------------------------------------
# 8. TABLE OUTPUTS
# ------------------------------------------------------------

write_csv(duplicate_gene_keys,
          file.path(table_dir, "phase7_01_duplicate_gene_keys.csv"))
write_csv(missing_value_counts,
          file.path(table_dir, "phase7_01_missing_value_counts.csv"))
write_csv(invalid_statistics,
          file.path(table_dir, "phase7_01_invalid_statistics.csv"))
write_csv(recalculated_summary,
          file.path(table_dir, "phase7_01_recalculated_summary.csv"))
write_csv(summary_comparison,
          file.path(table_dir, "phase7_01_summary_reconciliation.csv"))
write_csv(model_audit,
          file.path(table_dir, "phase7_01_model_diagnostics.csv"))
write_csv(donor_coverage_long,
          file.path(table_dir, "phase7_01_donor_coverage.csv"))
write_csv(donor_concentration,
          file.path(table_dir, "phase7_01_donor_concentration.csv"))
write_csv(contrast_summary,
          file.path(table_dir, "phase7_01_contrast_summary.csv"))
write_csv(population_summary,
          file.path(table_dir, "phase7_01_population_summary.csv"))
write_csv(t_vs_zstd,
          file.path(table_dir, "phase7_01_t_vs_zstd.csv"))
write_csv(cluster2_target,
          file.path(table_dir, "phase7_01_cluster2_NL_vs_IA_top10_genes.csv"))
write_csv(donor_expression,
          file.path(table_dir, "phase7_01_cluster2_NL_vs_IA_donor_expression.csv"))
write_csv(donor_level_summary,
          file.path(table_dir, "phase7_01_cluster2_NL_vs_IA_donor_summary.csv"))
write_csv(audit_checks,
          file.path(table_dir, "phase7_01_audit_checklist.csv"))

# ------------------------------------------------------------
# 9. FIGURE OUTPUTS
# ------------------------------------------------------------

plot_summary <- recalculated_summary %>%
  dplyr::mutate(
    Population_ID = factor(Population_ID, levels = rev(primary_populations)),
    Contrast = factor(Contrast, levels = expected_contrasts)
  )

p_fdr_heatmap <- ggplot(
  plot_summary,
  aes(x = Contrast, y = Population_ID, fill = n_FDR_05)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_FDR_05), size = 2.7) +
  scale_fill_viridis_c(name = "FDR < 0.05") +
  labs(
    title = "FDR-significant genes across CD8 populations",
    x = "Clinical-state contrast",
    y = "CD8 population"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "phase7_01_fdr_heatmap.png"),
  p_fdr_heatmap, width = 12, height = 6, dpi = 300
)

p_donor_coverage <- ggplot(
  donor_coverage_long,
  aes(x = Clinical_State, y = Population_ID, fill = Qualifying_Donors)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = Qualifying_Donors), size = 3) +
  scale_fill_viridis_c(name = "Qualifying donors") +
  labs(
    title = "Donor representation in primary CD8 populations",
    x = "Clinical state",
    y = "CD8 population"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(
  file.path(figure_dir, "phase7_01_donor_coverage.png"),
  p_donor_coverage, width = 9, height = 6, dpi = 300
)

p_fdr_effect <- ggplot(
  gene_results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    alpha = adj.P.Val < 0.05
  )
) +
  geom_point(size = 0.35) +
  facet_wrap(~ Contrast, scales = "free", ncol = 2) +
  scale_alpha_manual(
    values = c(`FALSE` = 0.18, `TRUE` = 0.55),
    guide = "none"
  ) +
  labs(
    title = "Effect estimates and adjusted significance",
    x = "Modelled expression difference (logFC)",
    y = expression(-log[10]("adjusted p-value"))
  ) +
  theme_minimal(base_size = 10)

ggsave(
  file.path(figure_dir, "phase7_01_fdr_vs_logFC.png"),
  p_fdr_effect, width = 12, height = 16, dpi = 300
)

# ------------------------------------------------------------
# 10. TEXT REPORT AND CONSOLE SUMMARY
# ------------------------------------------------------------

n_fdr05 <- sum(gene_results$adj.P.Val < 0.05, na.rm = TRUE)
n_fdr01 <- sum(gene_results$adj.P.Val < 0.01, na.rm = TRUE)

report_lines <- c(
  "PHASE 7.01 — DIFFERENTIAL EXPRESSION RESULTS AUDIT",
  "=================================================",
  "",
  paste("Gene-level result rows:", nrow(gene_results)),
  paste("Primary CD8 populations:", length(primary_populations)),
  paste("Clinical contrasts:", length(expected_contrasts)),
  paste("Population-contrast combinations:", nrow(expected_combinations)),
  "",
  "MODEL DIAGNOSTICS",
  paste("Models:", model_summary$n_models),
  paste("Error-free models:", model_summary$n_error_free),
  paste("Retained donors, minimum:", model_summary$min_retained_donors),
  paste("Retained donors, maximum:", model_summary$max_retained_donors),
  "",
  "RESULT INTEGRITY",
  paste("Duplicate gene keys:", nrow(duplicate_gene_keys)),
  paste("Missing required values:", sum(missing_value_counts$Missing_N)),
  paste("Invalid statistical values:", sum(unlist(invalid_statistics))),
  paste("Summary reconciliation failures:",
        sum(!summary_comparison$all_checks_pass)),
  "",
  "SIGNIFICANCE",
  paste("FDR < 0.05:", n_fdr05),
  paste("FDR < 0.01:", n_fdr01),
  "",
  "CLUSTER 2 NL VS IA DONOR CHECK",
  paste("Genes checked:", nrow(donor_level_summary)),
  paste("Directionally concordant genes:",
        sum(donor_level_summary$Direction_matches), "/", nrow(donor_level_summary)),
  "",
  "AUDIT CHECKLIST",
  paste(audit_checks$Check,
        ifelse(audit_checks$Passed, "PASS", "FAIL"),
        sep = ": "),
  "",
  paste("OVERALL AUDIT:", ifelse(overall_pass, "PASS", "REQUIRES REVIEW")),
  "",
  "LIMITATION:",
  "This script audits exported Phase 6 results and the stored donor-level",
  "aggregated estimates. It does not refit the Phase 6 models or independently",
  "reprocess the original GEO single-cell data."
)

writeLines(
  report_lines,
  file.path(table_dir, "phase7_01_de_results_audit_report.txt")
)

cat("\n============================================================\n")
cat("PHASE 7.01 DIFFERENTIAL EXPRESSION AUDIT COMPLETE\n")
cat("============================================================\n")
cat("Overall audit:", ifelse(overall_pass, "PASS", "REQUIRES REVIEW"), "\n")
cat("Gene-level rows:", nrow(gene_results), "\n")
cat("FDR < 0.05:", n_fdr05, "\n")
cat("FDR < 0.01:", n_fdr01, "\n")
cat(
  "Cluster 2 NL-vs-IA top-gene direction concordance:",
  sum(donor_level_summary$Direction_matches), "/",
  nrow(donor_level_summary), "\n"
)
cat("Tables and report:", table_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat("============================================================\n")

