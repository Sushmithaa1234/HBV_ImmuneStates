# ============================================================
# PHASE 7.02 — PATHWAY ENRICHMENT ANALYSIS
#
# Inputs:
#   Phase 6 donor-aware CD8 differential-state gene-level results
#
# Analyses:
#   1. Rank all tested genes by signed moderated t-statistic.
#   2. Preranked GSEA: Hallmark, Reactome, GO Biological Process.
#   3. Direction-specific ORA for FDR-significant genes, where
#      at least 10 significant genes are available.
#   4. Summarise pathway recurrence across populations/contrasts.
#
# Important:
#   - GSEA uses the complete ranked list for each comparison.
#   - ORA uses the tested genes for that same comparison as universe.
#   - Positive t/NES means enrichment toward the first state in
#     the contrast name; negative means toward the second state.
#   - Clinical states are categorical groups, not a time sequence.
#   - This script does not refit Phase 6 differential-expression models.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. PATHS AND PARAMETERS
# ------------------------------------------------------------

input_dir <- "results/tables/Phase 6"
table_dir <- "results/tables/Phase 7"
figure_dir <- "results/figures/Phase 7"
rds_dir <- "results/rds_objects"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- file.path(input_dir, "phase6_cd8_all_gene_level_results.csv")
eligibility_file <- file.path(input_dir, "phase6_cd8_population_eligibility.csv")

final_rds_file <- file.path(
  rds_dir,
  "phase7_final_enrichment_analysis.rds"
)

required_files <- c(gene_file, eligibility_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

required_packages <- c("msigdbr", "fgsea")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running Phase 7 Script 2: ",
    paste(missing_packages, collapse = ", "),
    "\nCRAN: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))\nBioconductor (if needed): BiocManager::install()"
  )
}

# Analysis thresholds
FDR_CUTOFF <- 0.05
MIN_ORA_GENES <- 10L
GSEA_MIN_SIZE <- 15L
GSEA_MAX_SIZE <- 500L
ORA_MIN_SIZE <- 10L
ORA_MAX_SIZE <- 500L
TOP_PATHWAYS_PER_LIBRARY <- 20L
RANDOM_SEED <- 20261103L

set.seed(RANDOM_SEED)

expected_contrasts <- c(
  "NL_vs_IT", "NL_vs_IA", "NL_vs_AR", "NL_vs_AC",
  "IT_vs_IA", "IT_vs_AR", "IT_vs_AC",
  "IA_vs_AR", "IA_vs_AC", "AR_vs_AC"
)

# ------------------------------------------------------------
# 2. LOAD AND VALIDATE PHASE 6 RESULTS
# ------------------------------------------------------------

gene_results <- readr::read_csv(
  gene_file,
  show_col_types = FALSE
)

population_eligibility <- readr::read_csv(
  eligibility_file,
  show_col_types = FALSE
)

required_gene_columns <- c(
  "Population_ID", "Gene", "Contrast",
  "logFC", "t", "P.Value", "adj.P.Val"
)

missing_columns <- setdiff(required_gene_columns, names(gene_results))
if (length(missing_columns) > 0) {
  stop(
    "Phase 6 gene-level results are missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!all(c("Population_ID", "Primary_Eligible") %in%
         names(population_eligibility))) {
  stop("Population eligibility file lacks required columns.")
}

primary_populations <- population_eligibility %>%
  dplyr::filter(as.logical(Primary_Eligible)) %>%
  dplyr::pull(Population_ID) %>%
  unique() %>%
  sort()

if (length(primary_populations) != 8L) {
  stop("Expected eight Phase 6 primary CD8 populations.")
}

if (!setequal(unique(gene_results$Population_ID), primary_populations)) {
  stop("Gene-level results do not match the eight primary populations.")
}

if (!setequal(unique(gene_results$Contrast), expected_contrasts)) {
  stop("Gene-level results do not contain the expected ten contrasts.")
}

if (anyNA(gene_results$Gene) ||
    anyNA(gene_results$Population_ID) ||
    anyNA(gene_results$Contrast) ||
    any(!is.finite(gene_results$t)) ||
    any(!is.finite(gene_results$adj.P.Val))) {
  stop("Missing or non-finite gene identifiers/statistics detected.")
}

if (any(gene_results$adj.P.Val < 0 | gene_results$adj.P.Val > 1)) {
  stop("Adjusted p-values outside [0, 1] detected.")
}

duplicate_keys <- gene_results %>%
  dplyr::count(Population_ID, Contrast, Gene, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_keys) > 0L) {
  stop(
    "Duplicate Population_ID × Contrast × Gene keys detected. ",
    "Resolve duplicate gene records before enrichment."
  )
}

comparison_index <- gene_results %>%
  dplyr::distinct(Population_ID, Contrast) %>%
  dplyr::arrange(Population_ID, Contrast)

if (nrow(comparison_index) != 80L) {
  stop("Expected 80 population-by-contrast comparisons.")
}

# ------------------------------------------------------------
# 3. RETRIEVE MSIGDB GENE SETS
# ------------------------------------------------------------

# Use the current msigdbr API when available, with compatibility
# for older installed versions.
get_msigdb <- function(collection, subcollection = NULL) {
  msig_formals <- names(formals(msigdbr::msigdbr))
  args <- list()

  if ("db_species" %in% msig_formals) {
    args$db_species <- "HS"
  }

  args$species <- "Homo sapiens"

  if ("collection" %in% msig_formals) {
    args$collection <- collection
    if (!is.null(subcollection)) {
      args$subcollection <- subcollection
    }
  } else {
    args$category <- collection
    if (!is.null(subcollection)) {
      # Older msigdbr versions use the older subcategory names.
      args$subcategory <- if (subcollection == "GO:BP") {
        "BP"
      } else {
        subcollection
      }
    }
  }

  do.call(msigdbr::msigdbr, args)
}

collection_specs <- list(
  Hallmark = list(collection = "H", subcollection = NULL),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  GO_Biological_Process = list(collection = "C5", subcollection = "GO:BP")
)

gene_sets <- list()
gene_set_tables <- list()
gene_set_manifest <- list()

for (library_name in names(collection_specs)) {
  spec <- collection_specs[[library_name]]

  msig_df <- get_msigdb(
    collection = spec$collection,
    subcollection = spec$subcollection
  )

  if (!all(c("gs_name", "gene_symbol") %in% names(msig_df))) {
    stop(
      "Unexpected msigdbr output columns for ",
      library_name,
      ". Expected gs_name and gene_symbol."
    )
  }

  term2gene <- msig_df %>%
    dplyr::transmute(
      Pathway = as.character(gs_name),
      Gene = toupper(as.character(gene_symbol))
    ) %>%
    dplyr::filter(
      !is.na(Pathway), Pathway != "",
      !is.na(Gene), Gene != ""
    ) %>%
    dplyr::distinct()

  pathway_list <- split(term2gene$Gene, term2gene$Pathway)
  pathway_list <- lapply(pathway_list, unique)

  if (length(pathway_list) == 0L) {
    stop("No gene sets retrieved for ", library_name, ".")
  }

  gene_sets[[library_name]] <- pathway_list
  gene_set_tables[[library_name]] <- term2gene

  gene_set_manifest[[library_name]] <- data.frame(
    Library = library_name,
    Collection = spec$collection,
    Subcollection = ifelse(
      is.null(spec$subcollection), "",
      spec$subcollection
    ),
    N_Pathways = length(pathway_list),
    N_Gene_Mappings = nrow(term2gene),
    MSigDB_Version = if ("db_version" %in% names(msig_df)) {
      paste(unique(msig_df$db_version), collapse = ";")
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}

gene_set_manifest <- dplyr::bind_rows(gene_set_manifest)

# Export the gene-set inventory/version used in this analysis.
readr::write_csv(
  gene_set_manifest,
  file.path(table_dir, "phase7_02_gene_set_manifest.csv")
)

# ------------------------------------------------------------
# 4. BUILD RANKED GENE LISTS AND RUN GSEA
# ------------------------------------------------------------

ranked_gene_lists <- list()
ranked_gene_tables <- list()
gsea_result_list <- list()
gsea_error_list <- list()
ranked_gene_coverage_list <- list()

comparison_ids <- seq_len(nrow(comparison_index))

for (i in comparison_ids) {
  pop <- comparison_index$Population_ID[i]
  contrast <- comparison_index$Contrast[i]
  comparison_id <- paste(pop, contrast, sep = "__")

  comparison_data <- gene_results %>%
    dplyr::filter(
      Population_ID == pop,
      Contrast == contrast
    ) %>%
    dplyr::mutate(Gene = toupper(as.character(Gene))) %>%
    dplyr::arrange(dplyr::desc(t), Gene)

  if (anyDuplicated(comparison_data$Gene)) {
    stop("Duplicate gene symbols after case normalisation: ", comparison_id)
  }

  # Signed moderated t-statistic is the ranking metric.
  ranks <- stats::setNames(
    as.numeric(comparison_data$t),
    comparison_data$Gene
  )
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < GSEA_MIN_SIZE) {
    stop("Insufficient ranked genes for GSEA: ", comparison_id)
  }

  ranked_gene_lists[[comparison_id]] <- ranks

  ranked_gene_tables[[comparison_id]] <- data.frame(
    Population_ID = pop,
    Contrast = contrast,
    Rank = seq_along(ranks),
    Gene = names(ranks),
    Ranked_Statistic_t = as.numeric(ranks),
    stringsAsFactors = FALSE
  )

  for (library_name in names(gene_sets)) {
    result_key <- paste(comparison_id, library_name, sep = "__")

    library_genes <- unique(unlist(gene_sets[[library_name]], use.names = FALSE))
    mapped_genes <- intersect(names(ranks), library_genes)

    ranked_gene_coverage_list[[result_key]] <- data.frame(
      Population_ID = pop,
      Contrast = contrast,
      Library = library_name,
      N_Ranked_Genes = length(ranks),
      N_Genes_Mapped = length(mapped_genes),
      Percent_Genes_Mapped = 100 * length(mapped_genes) / length(ranks),
      N_Gene_Sets = length(gene_sets[[library_name]]),
      stringsAsFactors = FALSE
    )

    gsea_out <- tryCatch(
      {
        fgsea::fgseaMultilevel(
          pathways = gene_sets[[library_name]],
          stats = ranks,
          minSize = GSEA_MIN_SIZE,
          maxSize = GSEA_MAX_SIZE,
          eps = 0,
          scoreType = "std",
          nproc = 1
        )
      },
      error = function(e) e
    )

    if (inherits(gsea_out, "error")) {
      gsea_error_list[[result_key]] <- data.frame(
        Population_ID = pop,
        Contrast = contrast,
        Library = library_name,
        Error = conditionMessage(gsea_out),
        stringsAsFactors = FALSE
      )
      next
    }

    gsea_df <- as.data.frame(gsea_out)

    if (nrow(gsea_df) > 0L) {
      gsea_df$LeadingEdge <- vapply(
        gsea_df$leadingEdge,
        function(x) paste(x, collapse = ";"),
        character(1)
      )
      gsea_df$leadingEdge <- NULL

      contrast_parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]

      gsea_df <- gsea_df %>%
        dplyr::mutate(
          Population_ID = pop,
          Contrast = contrast,
          First_State = contrast_parts[1],
          Second_State = contrast_parts[2],
          Library = library_name,
          Enriched_Toward = dplyr::case_when(
            NES > 0 ~ contrast_parts[1],
            NES < 0 ~ contrast_parts[2],
            TRUE ~ NA_character_
          )
        )

      gsea_result_list[[result_key]] <- gsea_df
    }
  }

  if (i %% 10L == 0L || i == length(comparison_ids)) {
    message(
      "GSEA completed for ",
      i, "/", length(comparison_ids),
      " population-contrast comparisons."
    )
  }
}

ranked_gene_table <- dplyr::bind_rows(ranked_gene_tables)
gsea_results <- dplyr::bind_rows(gsea_result_list)
gsea_errors <- dplyr::bind_rows(gsea_error_list)
ranked_gene_coverage <- dplyr::bind_rows(ranked_gene_coverage_list)

if (nrow(gsea_results) == 0L) {
  stop("GSEA returned no results. Inspect the GSEA error log.")
}

# ------------------------------------------------------------
# 5. DIRECTION-SPECIFIC OVER-REPRESENTATION ANALYSIS (ORA)
# ------------------------------------------------------------

# Hypergeometric ORA is calculated directly to keep the tested
# universe explicit and identical to the genes tested in each
# population-contrast model.
run_ora <- function(
  selected_genes,
  universe_genes,
  term2gene,
  min_size = ORA_MIN_SIZE,
  max_size = ORA_MAX_SIZE
) {
  universe_genes <- unique(universe_genes)
  selected_genes <- intersect(unique(selected_genes), universe_genes)

  term_sets <- split(term2gene$Gene, term2gene$Pathway)
  N <- length(universe_genes)
  n <- length(selected_genes)

  if (N == 0L || n == 0L) {
    return(data.frame())
  }

  rows <- lapply(names(term_sets), function(pathway_name) {
    pathway_genes <- intersect(unique(term_sets[[pathway_name]]), universe_genes)
    m <- length(pathway_genes)

    if (m < min_size || m > max_size) {
      return(NULL)
    }

    overlap_genes <- intersect(pathway_genes, selected_genes)
    k <- length(overlap_genes)

    p_value <- stats::phyper(
      q = k - 1L,
      m = m,
      n = N - m,
      k = n,
      lower.tail = FALSE
    )

    data.frame(
      Pathway = pathway_name,
      Gene_Set_Size_in_Universe = m,
      Selected_Gene_Count = n,
      Universe_Gene_Count = N,
      Overlap_Count = k,
      Gene_Ratio = paste0(k, "/", n),
      Background_Ratio = paste0(m, "/", N),
      Fold_Enrichment = if (n > 0 && m > 0) {
        (k / n) / (m / N)
      } else {
        NA_real_
      },
      P_Value = p_value,
      Overlap_Genes = paste(overlap_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  ora_df <- dplyr::bind_rows(rows)

  if (nrow(ora_df) > 0L) {
    ora_df$FDR <- stats::p.adjust(ora_df$P_Value, method = "BH")
    ora_df <- ora_df %>%
      dplyr::filter(Overlap_Count > 0L) %>%
      dplyr::arrange(FDR, P_Value)
  }

  ora_df
}

ora_result_list <- list()
ora_status_list <- list()

for (i in comparison_ids) {
  pop <- comparison_index$Population_ID[i]
  contrast <- comparison_index$Contrast[i]
  comparison_id <- paste(pop, contrast, sep = "__")

  comparison_data <- gene_results %>%
    dplyr::filter(
      Population_ID == pop,
      Contrast == contrast
    ) %>%
    dplyr::mutate(Gene = toupper(as.character(Gene)))

  universe_genes <- unique(comparison_data$Gene)
  significant_genes <- comparison_data %>%
    dplyr::filter(adj.P.Val < FDR_CUTOFF)

  # Positive logFC: higher in the first state.
  # Negative logFC: lower in the first state (higher in the second).
  direction_sets <- list(
    Up_in_first_state = significant_genes %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::pull(Gene),
    Down_in_first_state = significant_genes %>%
      dplyr::filter(logFC < 0) %>%
      dplyr::pull(Gene)
  )

  contrast_parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]

  for (direction_name in names(direction_sets)) {
    selected <- unique(direction_sets[[direction_name]])
    status_key <- paste(comparison_id, direction_name, sep = "__")

    if (length(selected) < MIN_ORA_GENES) {
      for (library_name in names(gene_set_tables)) {
        skip_key <- paste(status_key, library_name, sep = "__")
        ora_status_list[[skip_key]] <- data.frame(
          Population_ID = pop,
          Contrast = contrast,
          Direction = direction_name,
          Library = library_name,
          N_Significant_Genes = length(selected),
          Status = "Skipped: fewer than minimum significant genes",
          N_Pathways_Tested = 0L,
          stringsAsFactors = FALSE
        )
      }
      next
    }

    for (library_name in names(gene_set_tables)) {
      ora_key <- paste(status_key, library_name, sep = "__")

      ora_df <- run_ora(
        selected_genes = selected,
        universe_genes = universe_genes,
        term2gene = gene_set_tables[[library_name]]
      )

      if (nrow(ora_df) > 0L) {
        ora_df <- ora_df %>%
          dplyr::mutate(
            Population_ID = pop,
            Contrast = contrast,
            First_State = contrast_parts[1],
            Second_State = contrast_parts[2],
            Direction = direction_name,
            Library = library_name,
            Enriched_Toward = if (
              direction_name == "Up_in_first_state"
            ) {
              contrast_parts[1]
            } else {
              contrast_parts[2]
            }
          )

        ora_result_list[[ora_key]] <- ora_df
      }

      # Number of eligible pathways tested is calculated using
      # the same comparison-specific universe and size limits.
      term_sets <- split(
        gene_set_tables[[library_name]]$Gene,
        gene_set_tables[[library_name]]$Pathway
      )
      n_tested <- sum(vapply(
        term_sets,
        function(x) {
          size <- length(intersect(unique(x), universe_genes))
          size >= ORA_MIN_SIZE && size <= ORA_MAX_SIZE
        },
        logical(1)
      ))

      ora_status_list[[ora_key]] <- data.frame(
        Population_ID = pop,
        Contrast = contrast,
        Direction = direction_name,
        Library = library_name,
        N_Significant_Genes = length(selected),
        Status = "Run",
        N_Pathways_Tested = n_tested,
        stringsAsFactors = FALSE
      )
    }
  }

  if (i %% 10L == 0L || i == length(comparison_ids)) {
    message(
      "ORA completed for ",
      i, "/", length(comparison_ids),
      " population-contrast comparisons."
    )
  }
}

ora_results <- dplyr::bind_rows(ora_result_list)
ora_status <- dplyr::bind_rows(ora_status_list)

# ------------------------------------------------------------
# 6. PATHWAY COMPARISON TABLES
# ------------------------------------------------------------

gsea_significant <- gsea_results %>%
  dplyr::filter(padj < FDR_CUTOFF) %>%
  dplyr::arrange(Library, Contrast, Population_ID, padj)

gsea_recurrence_by_contrast <- gsea_significant %>%
  dplyr::group_by(Library, Contrast, pathway) %>%
  dplyr::summarise(
    N_Populations_Enriched = dplyr::n_distinct(Population_ID),
    Mean_NES = mean(NES, na.rm = TRUE),
    Min_NES = min(NES, na.rm = TRUE),
    Max_NES = max(NES, na.rm = TRUE),
    N_Positive_NES = sum(NES > 0, na.rm = TRUE),
    N_Negative_NES = sum(NES < 0, na.rm = TRUE),
    Direction_Consistent = N_Positive_NES == 0L || N_Negative_NES == 0L,
    .groups = "drop"
  ) %>%
  dplyr::arrange(Library, Contrast, dplyr::desc(N_Populations_Enriched))

gsea_recurrence_by_population <- gsea_significant %>%
  dplyr::group_by(Library, Population_ID, pathway) %>%
  dplyr::summarise(
    N_Contrasts_Enriched = dplyr::n_distinct(Contrast),
    Mean_NES = mean(NES, na.rm = TRUE),
    N_Positive_NES = sum(NES > 0, na.rm = TRUE),
    N_Negative_NES = sum(NES < 0, na.rm = TRUE),
    Direction_Consistent = N_Positive_NES == 0L || N_Negative_NES == 0L,
    .groups = "drop"
  ) %>%
  dplyr::arrange(Library, Population_ID, dplyr::desc(N_Contrasts_Enriched))

gsea_counts <- gsea_results %>%
  dplyr::group_by(Library, Population_ID, Contrast) %>%
  dplyr::summarise(
    N_Pathways_Tested = dplyr::n(),
    N_FDR_05 = sum(padj < FDR_CUTOFF, na.rm = TRUE),
    N_Positive_NES_FDR_05 = sum(padj < FDR_CUTOFF & NES > 0, na.rm = TRUE),
    N_Negative_NES_FDR_05 = sum(padj < FDR_CUTOFF & NES < 0, na.rm = TRUE),
    .groups = "drop"
  )

ora_significant <- ora_results %>%
  dplyr::filter(FDR < FDR_CUTOFF) %>%
  dplyr::arrange(Library, Contrast, Population_ID, Direction, FDR)

ora_counts <- ora_results %>%
  dplyr::group_by(Library, Population_ID, Contrast, Direction) %>%
  dplyr::summarise(
    N_Pathways_With_Overlap = dplyr::n(),
    N_FDR_05 = sum(FDR < FDR_CUTOFF, na.rm = TRUE),
    .groups = "drop"
  )

ora_recurrence <- ora_significant %>%
  dplyr::group_by(Library, Contrast, Direction, Pathway) %>%
  dplyr::summarise(
    N_Populations_Enriched = dplyr::n_distinct(Population_ID),
    Mean_Fold_Enrichment = mean(Fold_Enrichment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Library, Contrast, Direction,
                 dplyr::desc(N_Populations_Enriched))

# ------------------------------------------------------------
# 7. EXPORT TABLES
# ------------------------------------------------------------

write_table <- function(x, filename) {
  readr::write_csv(x, file.path(table_dir, filename))
}

write_table(
  ranked_gene_table,
  "phase7_02_ranked_gene_lists.csv"
)
write_table(
  gsea_results,
  "phase7_02_GSEA_all_results.csv"
)
write_table(
  gsea_significant,
  "phase7_02_GSEA_FDR05_results.csv"
)
write_table(
  gsea_errors,
  "phase7_02_GSEA_errors.csv"
)
write_table(
  gsea_counts,
  "phase7_02_GSEA_pathway_counts.csv"
)
write_table(
  ranked_gene_coverage,
  "phase7_02_ranked_gene_MSigDB_coverage.csv"
)
write_table(
  gsea_recurrence_by_contrast,
  "phase7_02_GSEA_recurrence_by_contrast.csv"
)
write_table(
  gsea_recurrence_by_population,
  "phase7_02_GSEA_recurrence_by_population.csv"
)
write_table(
  ora_results,
  "phase7_02_ORA_all_results.csv"
)
write_table(
  ora_significant,
  "phase7_02_ORA_FDR05_results.csv"
)
write_table(
  ora_status,
  "phase7_02_ORA_run_status.csv"
)
write_table(
  ora_counts,
  "phase7_02_ORA_pathway_counts.csv"
)
write_table(
  ora_recurrence,
  "phase7_02_ORA_recurrence.csv"
)

# ------------------------------------------------------------
# 8. FIGURES
# ------------------------------------------------------------

# One comparison heatmap per gene-set library. The displayed
# pathways are the 20 most recurrent FDR-significant pathways
# in that library; tiles show NES across populations/contrasts.
for (library_name in names(gene_sets)) {
  library_sig <- gsea_significant %>%
    dplyr::filter(Library == library_name)

  if (nrow(library_sig) == 0L) {
    message("No FDR-significant GSEA pathways for ", library_name,
            "; heatmap skipped.")
    next
  }

  top_pathways <- library_sig %>%
    dplyr::count(pathway, name = "N_Significant_Comparisons") %>%
    dplyr::arrange(
      dplyr::desc(N_Significant_Comparisons),
      pathway
    ) %>%
    dplyr::slice_head(n = TOP_PATHWAYS_PER_LIBRARY) %>%
    dplyr::pull(pathway)

  heatmap_data <- gsea_results %>%
    dplyr::filter(
      Library == library_name,
      pathway %in% top_pathways
    ) %>%
    dplyr::mutate(
      pathway = factor(pathway, levels = rev(top_pathways)),
      Contrast = factor(Contrast, levels = expected_contrasts),
      Significant = padj < FDR_CUTOFF
    )

  p <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(x = Contrast, y = pathway, fill = NES)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::geom_point(
      data = dplyr::filter(heatmap_data, Significant),
      shape = 4,
      size = 1.1,
      stroke = 0.45
    ) +
    ggplot2::facet_wrap(~ Population_ID, ncol = 2) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    ggplot2::labs(
      title = paste0("GSEA pathway comparison — ", library_name),
      subtitle = paste0(
        "Top ", length(top_pathways),
        " recurrent pathways; crosses mark FDR < ",
        FDR_CUTOFF
      ),
      x = "Clinical-state contrast",
      y = "Pathway"
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 6),
      strip.text = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  safe_library_name <- gsub("[^A-Za-z0-9]+", "_", library_name)

  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0("phase7_02_GSEA_NES_heatmap_", safe_library_name, ".png")
    ),
    plot = p,
    width = 15,
    height = 12,
    dpi = 300,
    bg = "white"
  )
}

# GSEA significant pathway counts by population and contrast.
p_gsea_counts <- ggplot2::ggplot(
  gsea_counts,
  ggplot2::aes(x = Contrast, y = N_FDR_05)
) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(Library ~ Population_ID, scales = "free_y") +
  ggplot2::labs(
    title = "FDR-significant GSEA pathways",
    x = "Clinical-state contrast",
    y = "Number of pathways (FDR < 0.05)"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ),
    strip.text = ggplot2::element_text(size = 7)
  )

ggplot2::ggsave(
  file.path(figure_dir, "phase7_02_GSEA_significant_pathway_counts.png"),
  p_gsea_counts,
  width = 17,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ORA significant pathway counts by direction.
if (nrow(ora_counts) > 0L) {
  p_ora_counts <- ggplot2::ggplot(
    ora_counts,
    ggplot2::aes(x = Contrast, y = N_FDR_05, fill = Direction)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(Library ~ Population_ID, scales = "free_y") +
    ggplot2::labs(
      title = "FDR-significant ORA pathways",
      x = "Clinical-state contrast",
      y = "Number of pathways (FDR < 0.05)",
      fill = "Gene direction"
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      ),
      strip.text = ggplot2::element_text(size = 7)
    )

  ggplot2::ggsave(
    file.path(figure_dir, "phase7_02_ORA_significant_pathway_counts.png"),
    p_ora_counts,
    width = 17,
    height = 9,
    dpi = 300,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 9. FINAL PHASE 7 RDS AND RUN REPORT
# ------------------------------------------------------------

analysis_parameters <- list(
  input_gene_file = gene_file,
  input_eligibility_file = eligibility_file,
  FDR_cutoff = FDR_CUTOFF,
  minimum_ORA_significant_genes = MIN_ORA_GENES,
  GSEA_min_size = GSEA_MIN_SIZE,
  GSEA_max_size = GSEA_MAX_SIZE,
  ORA_min_size = ORA_MIN_SIZE,
  ORA_max_size = ORA_MAX_SIZE,
  GSEA_rank_metric = "Signed moderated t-statistic",
  ORA_background = "All genes tested in the same population-contrast",
  ORA_direction = "Up/down gene sets tested separately",
  GSEA_FDR_scope = "Within each population-contrast and gene-set library",
  ORA_FDR_scope = "Within each population-contrast, direction, and gene-set library",
  MSigDB_libraries = names(gene_sets),
  random_seed = RANDOM_SEED
)

package_versions <- c(
  R = as.character(getRversion()),
  msigdbr = as.character(utils::packageVersion("msigdbr")),
  fgsea = as.character(utils::packageVersion("fgsea")),
  dplyr = as.character(utils::packageVersion("dplyr")),
  ggplot2 = as.character(utils::packageVersion("ggplot2"))
)

phase7_results <- list(
  parameters = analysis_parameters,
  package_versions = package_versions,
  gene_set_manifest = gene_set_manifest,
  gene_sets = gene_sets,
  ranked_gene_lists = ranked_gene_lists,
  ranked_gene_table = ranked_gene_table,
  ranked_gene_coverage = ranked_gene_coverage,
  GSEA_all_results = gsea_results,
  GSEA_significant_results = gsea_significant,
  GSEA_errors = gsea_errors,
  GSEA_counts = gsea_counts,
  GSEA_recurrence_by_contrast = gsea_recurrence_by_contrast,
  GSEA_recurrence_by_population = gsea_recurrence_by_population,
  ORA_all_results = ora_results,
  ORA_significant_results = ora_significant,
  ORA_status = ora_status,
  ORA_counts = ora_counts,
  ORA_recurrence = ora_recurrence
)

saveRDS(
  phase7_results,
  final_rds_file,
  compress = "xz"
)

run_report <- c(
  "PHASE 7.02 — PATHWAY ENRICHMENT ANALYSIS",
  "========================================",
  paste0("Input gene-level rows: ", nrow(gene_results)),
  paste0("Primary CD8 populations: ", length(primary_populations)),
  paste0("Population-contrast comparisons: ", nrow(comparison_index)),
  paste0("Ranked gene entries: ", nrow(ranked_gene_table)),
  paste0("Ranked-gene coverage records: ", nrow(ranked_gene_coverage)),
  paste0("GSEA result rows: ", nrow(gsea_results)),
  paste0(
    "FDR-significant GSEA rows: ",
    sum(gsea_results$padj < FDR_CUTOFF, na.rm = TRUE)
  ),
  paste0("GSEA errors: ", nrow(gsea_errors)),
  paste0("ORA result rows with overlap: ", nrow(ora_results)),
  paste0(
    "FDR-significant ORA rows: ",
    sum(ora_results$FDR < FDR_CUTOFF, na.rm = TRUE)
  ),
  paste0(
    "ORA comparison-direction-library runs: ",
    nrow(ora_status)
  ),
  "",
  "GSEA ranking: signed moderated t-statistic.",
  "ORA universe: genes tested in the corresponding population-contrast.",
  "ORA directions: up in first state and down in first state, tested separately.",
  "GSEA FDR: adjusted within each population-contrast-library analysis.",
  "ORA FDR: adjusted within each population-contrast-direction-library analysis.",
  "Positive NES: enrichment toward the first state in the contrast.",
  "Negative NES: enrichment toward the second state in the contrast.",
  "Clinical states are categorical groups; no temporal progression is assumed.",
  "",
  paste0("MSigDB version(s): ",
         paste(unique(gene_set_manifest$MSigDB_Version), collapse = "; ")),
  paste0("Final RDS: ", final_rds_file)
)

writeLines(
  run_report,
  file.path(table_dir, "phase7_02_enrichment_run_report.txt")
)

message("============================================================")
message("PHASE 7.02 — COMPLETE")
message("GSEA results: ", nrow(gsea_results))
message("ORA result rows: ", nrow(ora_results))
message("Final RDS: ", final_rds_file)
message("Tables: ", table_dir)
message("Figures: ", figure_dir)
message("============================================================")