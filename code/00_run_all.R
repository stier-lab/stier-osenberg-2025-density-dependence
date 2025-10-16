###############################################################################
# File: run_all.R
# Purpose: Sequentially run the full analysis pipeline
# Author: Adrian Stier
# Date: 2025-07-09
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. Setup
# ─────────────────────────────────────────────────────────────────────────────
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required but not installed.")
}
here::i_am("code/00_run_all.R")

# Canonical RNG seed to keep any top-level sampling deterministic
set.seed(20250708)

source(file.path("code", "0_libraries.R"))

# create output folders just in case
dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load core libraries
# ─────────────────────────────────────────────────────────────────────────────
message("0_libraries.R sourced: all core packages loaded.")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Data & phylogeny
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "1_data_phylogeny_loading.R"))
message("1_data_phylogeny_loading.R sourced: data and trees available.")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Beta estimation
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "2_beta_estimate.R"))
message("2_beta_estiamte.R sourced: per‐study beta estimates computed.")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Phylogenetic analysis & model selection
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "3_phylogenetic_analysis.R"))
message("3_phylogenetic_analysis.R sourced: phylo signal, model comparisons, variance partitioning, species‐level plots.")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Predation‐paired study analysis
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "4_predators_pairedpredators.R"))
message("4_predators_pairedpredators.R sourced: predator‐presence meta‐analysis & paired comparisons.")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Model selection / AIC‐based reduction
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "5_model_selection_comparison.R"))
message("5_model_selection_comparison.R sourced: model selection via AICc and hierarchical deletion.")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Bivariate predictor vs β plots
# ─────────────────────────────────────────────────────────────────────────────
source(here("code", "6_bivariate_plots_predictors.R"))
message("6_bivariate_plots_predictors.R sourced: bivariate X vs β plots generated.")

# ─────────────────────────────────────────────────────────────────────────────
# 8. Summary
# ─────────────────────────────────────────────────────────────────────────────
message("All scripts executed. Check 'figures/' and 'results/' for outputs.")

# ─────────────────────────────────────────────────────────────────────────────
# Reproducibility footer
if (!dir.exists("results")) dir.create("results", recursive = TRUE, showWarnings = FALSE)

info <- c(
  paste("git_sha:", tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error=function(e) "NA")),
  paste("run_time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  capture.output(sessionInfo())
)
writeLines(info, "results/session-info.txt")

if (dir.exists("figures")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
  manifest <- list(
    git_sha = tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error=function(e) "NA"),
    r_version = R.version.string,
    seed = 20251014,
    produced = list.files("figures", recursive = TRUE)
  )
  jsonlite::write_json(manifest, "results/manifest.json", pretty = TRUE, auto_unbox = TRUE)
}
# ─────────────────────────────────────────────────────────────────────────────
# End of file