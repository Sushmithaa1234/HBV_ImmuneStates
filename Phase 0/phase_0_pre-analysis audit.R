# ============================================================
# HBV IMMUNE STATES scRNA-seq PROJECT
# PHASE 0 — SAMPLE AUDIT
# ============================================================
#
# Purpose:
#   Establish exactly what samples we have before performing
#   any preprocessing or biological analysis.
#
# This script does NOT:
#   - filter cells
#   - normalize data
#   - perform batch correction
#   - cluster cells
#   - create a Seurat object
#
# ============================================================


# ------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------

library(tidyverse)

setwd("C:/Users/sushm/OneDrive/Desktop/HBV scRNA")

raw_dir <- "Data/Raw"

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


# ------------------------------------------------------------
# 2. FIND LIVER FILES
# ------------------------------------------------------------

liver_files <- list.files(
  raw_dir,
  pattern = "_Liver_.*\\.txt$",
  full.names = TRUE
)

cat("============================================\n")
cat("PHASE 0 — SAMPLE AUDIT\n")
cat("============================================\n\n")

cat("Liver files found:", length(liver_files), "\n\n")

if (length(liver_files) != 23) {
  warning(
    "Expected 23 liver samples, but found ",
    length(liver_files),
    ". Please check the Data/Raw folder."
  )
}


# ------------------------------------------------------------
# 3. OFFICIAL SAMPLE METADATA
# ------------------------------------------------------------
#
# Disease phases:
#
# HC = Healthy control
# IT = Immune tolerant
# IA = Immune active
# AR = Acute resolved
# CR = Chronic resolved
#
# IMPORTANT:
# This mapping is based on the study/GEO metadata.
# Do not infer disease phase from filenames.
# ------------------------------------------------------------

sample_metadata <- tribble(
  
  ~GSM, ~donor, ~phase,
  
  "GSM5519469", "P190604", "HC",
  "GSM5519471", "P190326", "HC",
  "GSM5519472", "P190402", "HC",
  "GSM5519475", "P190716", "HC",
  "GSM5519477", "P190719", "HC",
  "GSM5519483", "Dhc570",  "HC",
  
  "GSM5519484", "D528848", "IT",
  "GSM5519485", "D529074", "IT",
  "GSM5519486", "D529351", "IT",
  "GSM5519487", "D529354", "IT",
  "GSM5519488", "D529409", "IT",
  "GSM5519491", "P190902", "IT",
  
  "GSM5519494", "P190910", "IA",
  "GSM5519495", "P190808", "IA",
  "GSM5519496", "P190801", "IA",
  "GSM5519497", "P190911", "IA",
  "GSM5519499", "P191028", "IA",
  
  "GSM5519502", "P191112", "AR",
  "GSM5519504", "P191008", "AR",
  "GSM5519506", "P191126", "AR",
  
  "GSM5519508", "P191127", "CR",
  "GSM5519510", "P191210", "CR",
  "GSM5519512", "P191217", "CR"
)


# ------------------------------------------------------------
# 4. CHECK METADATA
# ------------------------------------------------------------

cat("Sample metadata:\n\n")
print(sample_metadata)

cat("\nPatients per phase:\n\n")

phase_counts <- sample_metadata %>%
  count(phase, name = "patients")

print(phase_counts)


# ------------------------------------------------------------
# 5. AUDIT EACH FILE
# ------------------------------------------------------------

cat("\n============================================\n")
cat("READING SAMPLE FILES\n")
cat("============================================\n\n")

audit_list <- list()

for (file in liver_files) {
  
  # Extract GSM ID from filename
  gsm <- str_extract(
    basename(file),
    "GSM[0-9]+"
  )
  
  cat("Reading:", gsm, "\n")
  
  # Read expression matrix
  expression_matrix <- read.table(
    file,
    sep = " ",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  # Convert to matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # Number of genes
  n_genes <- nrow(expression_matrix)
  
  # Number of cells
  n_cells <- ncol(expression_matrix)
  
  # Genes detected per cell
  genes_per_cell <- colSums(
    expression_matrix > 0
  )
  
  # Total expression value per cell
  counts_per_cell <- colSums(
    expression_matrix
  )
  
  # Create one-row audit record
  audit_list[[gsm]] <- tibble(
    
    GSM = gsm,
    
    genes_in_matrix = n_genes,
    
    cells_in_matrix = n_cells,
    
    median_genes_per_cell =
      median(genes_per_cell),
    
    mean_genes_per_cell =
      mean(genes_per_cell),
    
    median_expression_per_cell =
      median(counts_per_cell),
    
    mean_expression_per_cell =
      mean(counts_per_cell),
    
    min_genes_per_cell =
      min(genes_per_cell),
    
    max_genes_per_cell =
      max(genes_per_cell)
    
  )
  
  cat(
    "  Genes:", n_genes,
    "| Cells:", n_cells,
    "| Median genes/cell:",
    round(median(genes_per_cell), 1),
    "\n\n"
  )
}


# ------------------------------------------------------------
# 6. COMBINE AUDIT RESULTS
# ------------------------------------------------------------

sample_audit <- bind_rows(audit_list)


# ------------------------------------------------------------
# 7. ADD BIOLOGICAL METADATA
# ------------------------------------------------------------

sample_audit <- sample_metadata %>%
  left_join(
    sample_audit,
    by = "GSM"
  ) %>%
  arrange(
    factor(
      phase,
      levels = c(
        "HC",
        "IT",
        "IA",
        "AR",
        "CR"
      )
    )
  )


# ------------------------------------------------------------
# 8. CHECK FOR MISSING INFORMATION
# ------------------------------------------------------------

cat("\n============================================\n")
cat("METADATA CHECK\n")
cat("============================================\n\n")

cat(
  "Samples in metadata:",
  nrow(sample_metadata),
  "\n"
)

cat(
  "Samples in expression data:",
  nrow(sample_audit),
  "\n"
)

cat(
  "Samples with missing metadata:",
  sum(is.na(sample_audit$phase)),
  "\n"
)

cat(
  "Samples with missing donor ID:",
  sum(is.na(sample_audit$donor)),
  "\n\n"
)


# ------------------------------------------------------------
# 9. PRINT FINAL AUDIT TABLE
# ------------------------------------------------------------

cat("============================================\n")
cat("FINAL SAMPLE AUDIT\n")
cat("============================================\n\n")

print(sample_audit)


# ------------------------------------------------------------
# 10. SAVE AUDIT TABLE
# ------------------------------------------------------------

write_csv(
  sample_audit,
  "results/tables/phase0_sample_audit.csv"
)

cat(
  "\n✓ Saved:\n",
  "results/tables/phase0_sample_audit.csv\n"
)


# ------------------------------------------------------------
# 11. SUMMARY BY DISEASE PHASE
# ------------------------------------------------------------

phase_summary <- sample_audit %>%
  
  group_by(phase) %>%
  
  summarise(
    
    patients = n(),
    
    total_cells =
      sum(cells_in_matrix),
    
    median_cells_per_sample =
      median(cells_in_matrix),
    
    median_genes_per_cell =
      median(median_genes_per_cell),
    
    .groups = "drop"
    
  )


cat("\n============================================\n")
cat("SUMMARY BY DISEASE PHASE\n")
cat("============================================\n\n")

print(phase_summary)


write_csv(
  phase_summary,
  "results/tables/phase0_phase_summary.csv"
)


# ------------------------------------------------------------
# 12. PLOT — CELLS PER SAMPLE
# ------------------------------------------------------------

plot_cells <- ggplot(
  sample_audit,
  aes(
    x = reorder(GSM, cells_in_matrix),
    y = cells_in_matrix,
    fill = phase
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    title = "Number of Cells per Liver Sample",
    x = "Sample",
    y = "Number of Cells"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom"
  )


ggsave(
  "results/figures/phase0_cells_per_sample.png",
  plot_cells,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 13. PLOT — MEDIAN GENES PER CELL
# ------------------------------------------------------------

plot_genes <- ggplot(
  sample_audit,
  aes(
    x = reorder(GSM, median_genes_per_cell),
    y = median_genes_per_cell,
    fill = phase
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    title = "Median Genes Detected per Cell",
    x = "Sample",
    y = "Median Genes per Cell"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "bottom"
  )


ggsave(
  "results/figures/phase0_median_genes_per_cell.png",
  plot_genes,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 14. FINISHED
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PHASE 0 SAMPLE AUDIT COMPLETE\n")
cat("============================================\n\n")

cat(
  "Next step:\n",
  "Review the audit table and plots before making\n",
  "any QC or filtering decisions.\n\n"
)

cat("NO CELLS HAVE BEEN FILTERED.\n")
cat("NO NORMALIZATION HAS BEEN PERFORMED.\n")
cat("NO BATCH CORRECTION HAS BEEN PERFORMED.\n")
cat("NO BIOLOGICAL ANALYSIS HAS BEEN PERFORMED.\n")

colnames(sample_audit)
library(GEOquery)
gse <- getGEO(
  "GSE182159",
  GSEMatrix = TRUE
)
length(gse)
colnames(pData(gse[[1]]))

pdat <- pData(gse[[1]])
unique(pdat$platform_id)
unique(pdat$instrument_model)
unique(pdat$library_strategy)
table(pdat$`Stage:ch1`)
pdat[, c(
  "geo_accession",
  "donor:ch1",
  "Stage:ch1",
  "tissue:ch1",
  "platform_id",
  "instrument_model"
)]

unique(pdat$library_selection)
unique(pdat$library_source)
unique(pdat$`data_processing`)


# ============================================================
# PHASE 0 — FINALIZE VERIFIED SAMPLE METADATA
# ============================================================

# GEO metadata is already loaded as:
# pdat <- pData(gse[[1]])

# ------------------------------------------------------------
# 1. Get the GSM IDs of the 23 liver files we downloaded
# ------------------------------------------------------------

liver_gsm <- str_extract(
  basename(liver_files),
  "GSM[0-9]+"
)

cat("Downloaded liver samples:", length(liver_gsm), "\n")


# ------------------------------------------------------------
# 2. Extract ONLY the liver samples from GEO metadata
# ------------------------------------------------------------

geo_liver_metadata <- pdat %>%
  filter(
    `tissue:ch1` == "Liver",
    geo_accession %in% liver_gsm
  ) %>%
  select(
    GSM = geo_accession,
    donor = `donor:ch1`,
    phase = `Stage:ch1`,
    tissue = `tissue:ch1`,
    platform = platform_id,
    instrument = instrument_model,
    library_selection,
    library_source,
    library_strategy
  )


# ------------------------------------------------------------
# 3. Check that we have exactly 23 samples
# ------------------------------------------------------------

cat(
  "GEO liver samples recovered:",
  nrow(geo_liver_metadata),
  "\n"
)

if (nrow(geo_liver_metadata) != 23) {
  stop(
    "ERROR: Expected 23 liver samples from GEO, but found ",
    nrow(geo_liver_metadata)
  )
}


# ------------------------------------------------------------
# 4. Check for duplicate GSM IDs
# ------------------------------------------------------------

if (anyDuplicated(geo_liver_metadata$GSM) > 0) {
  
  stop(
    "ERROR: Duplicate GSM IDs detected in GEO metadata."
  )
  
} else {
  
  cat("✓ No duplicate GSM IDs\n")
}


# ------------------------------------------------------------
# 5. Check that every downloaded file has GEO metadata
# ------------------------------------------------------------

missing_from_geo <- setdiff(
  liver_gsm,
  geo_liver_metadata$GSM
)

if (length(missing_from_geo) > 0) {
  
  cat(
    "WARNING: These downloaded samples were not found in GEO:\n"
  )
  
  print(missing_from_geo)
  
} else {
  
  cat(
    "✓ Every downloaded liver sample has GEO metadata\n"
  )
}


# ------------------------------------------------------------
# 6. Check tissue
# ------------------------------------------------------------

if (all(geo_liver_metadata$tissue == "Liver")) {
  
  cat("✓ All samples are liver tissue\n")
  
} else {
  
  stop(
    "ERROR: Non-liver samples entered the liver cohort."
  )
}


# ------------------------------------------------------------
# 7. Check official clinical states
# ------------------------------------------------------------

cat("\nClinical states:\n")

print(
  table(geo_liver_metadata$phase)
)


# Expected:
#
# AC = 3
# AR = 3
# IA = 5
# IT = 6
# NL = 6


# ------------------------------------------------------------
# 8. Check donor uniqueness
# ------------------------------------------------------------

cat("\nNumber of unique donors:\n")

print(
  n_distinct(geo_liver_metadata$donor)
)


if (
  n_distinct(geo_liver_metadata$donor) ==
  nrow(geo_liver_metadata)
) {
  
  cat(
    "✓ Each liver sample corresponds to a unique donor\n"
  )
  
} else {
  
  cat(
    "⚠ Some donors have more than one liver sample\n"
  )
}


# ------------------------------------------------------------
# 9. Add the sample-level audit information
# ------------------------------------------------------------

final_metadata <- geo_liver_metadata %>%
  
  left_join(
    sample_audit %>%
      select(
        GSM,
        genes_in_matrix,
        cells_in_matrix,
        median_genes_per_cell,
        mean_genes_per_cell,
        median_expression_per_cell,
        mean_expression_per_cell,
        min_genes_per_cell,
        max_genes_per_cell
      ),
    by = "GSM"
  )


# ------------------------------------------------------------
# 10. Add QC flags
# ------------------------------------------------------------

final_metadata <- final_metadata %>%
  
  mutate(
    
    sample_qc_flag = case_when(
      
      cells_in_matrix < 200 ~
        "MAJOR_OUTLIER",
      
      median_genes_per_cell < 700 ~
        "LOW_COMPLEXITY",
      
      median_genes_per_cell < 900 ~
        "LOWER_COMPLEXITY",
      
      TRUE ~
        "PASS"
    )
    
  )


# ------------------------------------------------------------
# 11. View the final verified metadata
# ------------------------------------------------------------

View(final_metadata)
nrow(final_metadata)
table(final_metadata$phase)
sum(is.na(final_metadata$GSM))
sum(is.na(final_metadata$donor))
sum(is.na(final_metadata$phase))


# ------------------------------------------------------------
# 12. Save it
# ------------------------------------------------------------

write_csv(
  final_metadata,
  "results/tables/final_liver_sample_metadata.csv"
)

cat(
  "\n✓ Final verified metadata saved.\n"
)


# ------------------------------------------------------------
# 13. Final Phase 0 summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PHASE 0 COMPLETE\n")
cat("============================================\n\n")

cat(
  "Liver samples:",
  nrow(final_metadata),
  "\n"
)

cat(
  "Unique donors:",
  n_distinct(final_metadata$donor),
  "\n"
)

cat(
  "Clinical states:\n"
)

print(
  table(final_metadata$phase)
)

cat(
  "\nSample QC flags:\n"
)

print(
  table(final_metadata$sample_qc_flag)
)

cat(
  "\nTechnical metadata:\n"
)

cat(
  "Platform:",
  paste(unique(final_metadata$platform), collapse = ", "),
  "\n"
)

cat(
  "Instrument:",
  paste(unique(final_metadata$instrument), collapse = ", "),
  "\n"
)

cat(
  "Library strategy:",
  paste(unique(final_metadata$library_strategy), collapse = ", "),
  "\n"
)

cat(
  "Library selection:",
  paste(unique(final_metadata$library_selection), collapse = ", "),
  "\n"
)

cat(
  "Library source:",
  paste(unique(final_metadata$library_source), collapse = ", "),
  "\n"
)

cat("\n============================================\n")
cat("READY FOR PHASE 1 — DATA LOADING\n")
cat("============================================\n")
