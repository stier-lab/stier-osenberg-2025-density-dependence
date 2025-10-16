# ─────────────────────────────────────────────────────────────────────────────
# Reproducibility bootstrap
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
try(renv::restore(), silent = TRUE)
set.seed(20251014)
RNGkind("L'Ecuyer-CMRG")
# ─────────────────────────────────────────────────────────────────────────────

# =============================================================================
# File: 0_libraries.R
# Purpose: Centralized, minimal set of packages for the entire phylogenetic meta-analysis
# Author: Adrian Stier
# Date: 2025-07-09
# =============================================================================



# ----------------------------
# 1. Data wrangling & I/O
# ----------------------------
required_pkgs_tidy <- c(
  "broom",
  "dplyr",
  "forcats",
  "magrittr",
  "purrr",
  "readr",
  "readxl",     # for reading Excel files
  "tidyr",
  "tibble"
)
for (pkg in required_pkgs_tidy) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 2. Project structure & strings
# ----------------------------
for (pkg in c("here", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 3. Phylogenetics
# ----------------------------
for (pkg in c("ape", "phytools", "ggtree", "geiger", "fishtree")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 4. Meta-analysis & modeling
# ----------------------------
required_pkgs_modeling <- c("metafor", "nls2")
for (pkg in required_pkgs_modeling) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 5. Visualization & theming
# ----------------------------
required_pkgs_viz <- c(
  "GGally",
  "RColorBrewer",
  "cowplot",
  "ggbeeswarm",
  "ggplot2",
  "ggtext",
  "patchwork",
  "scales",
  "viridis",
  "wesanderson"
)
for (pkg in required_pkgs_viz) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 6. Spatial data & mapping
# ----------------------------
for (pkg in c("sf", "rnaturalearth", "rnaturalearthdata", "ggspatial")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 7. Tables & reporting
# ----------------------------
required_pkgs_tables <- c("gt", "kableExtra", "knitr")
for (pkg in required_pkgs_tables) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 8. Graphics utilities & layout
# ----------------------------
required_pkgs_graphics <- c("grid", "gridExtra")
for (pkg in required_pkgs_graphics) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 9. Data serialization
# ----------------------------
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required but not installed.")
}
library(jsonlite)

# ----------------------------
# 10. Reproducibility
# ----------------------------
# print session info for reproducibility
message("--- Session information for reproducibility ---")
print(sessionInfo())
