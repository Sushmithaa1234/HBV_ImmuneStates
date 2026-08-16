# ============================================================
# HBV Immune States scRNA-seq Project
# 00_setup.R
#
# Purpose:
# Install and load all R packages required for the complete HBV immune-state scRNA-seq analysis pipeline.
#
# NOTE:
# renv is NOT used in this project.
# Package versions are recorded at the end of this script.
# ============================================================

# ============================================================
# 1. CRAN PACKAGES
# ============================================================

cran_packages <- c(
  
  # Core data manipulation
  "tidyverse",
  "data.table",
  "here",
  
  # Matrix / graph / general analysis
  "Matrix",
  "igraph",
  
  # General plotting
  "cowplot",
  "gridExtra",
  "ggpubr",
  "scales",
  "viridis",
  
  # Heatmaps
  "pheatmap",
  
  # Survival analysis
  "survival",
  "survminer",
  
  # Molecular signatures
  "msigdbr",
  
  # Installation / GitHub packages
  "remotes"
)

# ============================================================
# 2. BIOCONDUCTOR PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  
  # Single-cell annotation
  "SingleR",
  "celldex",
  
  # Statistical analysis
  "limma",
  "DESeq2",
  
  # Gene-set / pathway analysis
  "AUCell",
  "fgsea",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  
  # Single-cell / matrix infrastructure
  "ComplexHeatmap",
  "BiocNeighbors",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "Biobase",
  
  # GEO data access
  "GEOquery",
  
  # Differential-expression visualization
  "EnhancedVolcano"
)

# ============================================================
# 3. GITHUB PACKAGES
# ============================================================

github_packages <- c(
  "CellChat",
  "monocle3"
)

# ============================================================
# 4. INSTALL CRAN PACKAGES
# ============================================================

cat("\n")
cat("============================================\n")
cat("Installing CRAN packages\n")
cat("============================================\n\n")

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("Installing:", pkg, "\n")
    
    install.packages(pkg)
    
  } else {
    
    cat("Already installed:", pkg, "\n")
  }
}

# ============================================================
# 5. INSTALL BIOCONDUCTOR PACKAGES
# ============================================================

cat("\n")
cat("============================================\n")
cat("Installing Bioconductor packages\n")
cat("============================================\n\n")

for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("Installing:", pkg, "\n")
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  } else {
    
    cat("Already installed:", pkg, "\n")
  }
}

# ============================================================
# 6. INSTALL GITHUB PACKAGES
# ============================================================

cat("\n")
cat("============================================\n")
cat("Installing GitHub packages\n")
cat("============================================\n\n")


# CellChat
if (!requireNamespace("CellChat", quietly = TRUE)) {
  
  cat("Installing CellChat from GitHub...\n")
  
  remotes::install_github(
    "jinworks/CellChat"
  )
  
} else {
  
  cat("Already installed: CellChat\n")
}


# monocle3
if (!requireNamespace("monocle3", quietly = TRUE)) {
  
  cat("Installing monocle3 from GitHub...\n")
  
  remotes::install_github(
    "cole-trapnell-lab/monocle3"
  )
  
} else {
  
  cat("Already installed: monocle3\n")
}

# ============================================================
# 7. VERIFY INSTALLATION
# ============================================================

all_packages <- c(
  cran_packages,
  bioc_packages,
  github_packages
)

cat("\n")
cat("============================================\n")
cat("Verifying package installation\n")
cat("============================================\n\n")

package_status <- data.frame(
  package = all_packages,
  installed = sapply(
    all_packages,
    requireNamespace,
    quietly = TRUE
  ),
  version = sapply(
    all_packages,
    function(x) {
      if (requireNamespace(x, quietly = TRUE)) {
        as.character(packageVersion(x))
      } else {
        NA_character_
      }
    }
  ),
  stringsAsFactors = FALSE
)

print(package_status)


# Stop if anything is missing
missing_packages <- package_status$package[
  !package_status$installed
]

if (length(missing_packages) > 0) {
  
  stop(
    "\nThe following packages could not be installed:\n",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n\nPlease resolve these before running the analysis."
  )
}

# ============================================================
# 8. LOAD ALL PACKAGES
# ============================================================

cat("\n")
cat("============================================\n")
cat("Loading packages\n")
cat("============================================\n\n")

for (pkg in all_packages) {
  
  suppressPackageStartupMessages(
    library(
      pkg,
      character.only = TRUE
    )
  )
  
  cat(
    "Loaded:",
    pkg,
    "—",
    as.character(packageVersion(pkg)),
    "\n"
  )
}

# ============================================================
# 9. COMPUTATIONAL ENVIRONMENT INFORMATION
# ============================================================

cat("\n")
cat("============================================\n")
cat("Computational environment\n")
cat("============================================\n\n")

cat(
  "R version:",
  R.version.string,
  "\n"
)

cat(
  "Platform:",
  R.version$platform,
  "\n"
)

cat(
  "OS:",
  Sys.info()["sysname"],
  "\n"
)

cat(
  "Architecture:",
  Sys.info()["machine"],
  "\n"
)

# ============================================================
# 10. SAVE PACKAGE VERSION RECORD
# ============================================================

dir.create(
  "results/logs",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  package_status,
  "results/logs/package_versions.csv",
  row.names = FALSE
)

# ============================================================
# 11. SAVE SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file = "results/logs/sessionInfo.txt"
)

# ============================================================
# 12. FINAL STATUS
# ============================================================

cat("\n")
cat("============================================\n")
cat("ENVIRONMENT SETUP COMPLETE\n")
cat("============================================\n\n")

cat(
  "All required packages are installed and loaded.\n"
)

cat(
  "Package versions saved to:\n"
)

cat(
  "results/logs/package_versions.csv\n\n"
)

cat(
  "Session information saved to:\n"
)

cat(
  "results/logs/sessionInfo.txt\n\n"
)

cat(
  "The environment is ready for the HBV scRNA-seq pipeline.\n"
)
