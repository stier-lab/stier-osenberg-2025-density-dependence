###############################################################################
# Script: 1_data_phylogeny_loading.R
# Goal:   Assemble the canonical analysis table (`all_dat2`) by merging nonlinear
#         mortality estimates, curated covariates, duration/density summaries,
#         and species-level traits. Also prunes the default phylogeny for reuse.
# Inputs: output/combined_results_2024-09-17.csv,
#         data/covariates-2024-09-30.csv, data/all_studies_looped-2024-09-11.csv,
#         data/manual_densities.csv, data/unique_species_studies.xlsx, phylogenies.
# Outputs: all_dat2 (in memory), pruned_tree (in memory),
#          output/merged-covariates-2024-10-21.csv, output/unique_* CSV/HTML tables.
###############################################################################

# -----------------------------------------------------------------------------
# 1. Dependencies & Preliminary Sources
# -----------------------------------------------------------------------------
# Use the shared loader so every other script can assume the same packages exist.
# This also keeps the order of sourcing consistent across the workflow.
source(here::here("code", "0_libraries.R"))  # Load required packages

# -----------------------------------------------------------------------------
# 2. Load Parameter Estimates & Covariates
# -----------------------------------------------------------------------------
# The goal of this section is to build `all_dat2`, a single table that merges:
#   • Nonlinear mortality fits (α, β) for every substudy
#   • Study-level covariates / metadata from the curation spreadsheet
#   • Duration and density summaries derived from the digitised time-series
# Downstream scripts only touch `all_dat2`, so we do all heavy lifting here.

# 2.1 Primary parameter estimates
# Updated to use 2025-10-15 results with corrected symmetric variance formula
params_primary <- read.csv(here::here("output", "combined_results_2025-10-15.csv")) %>%
  dplyr::select(
    substudy_num,
    alpha,
    beta,
    alpha_variance,
    beta_variance
  ) %>%
  dplyr::rename(
    alpha_nls2        = alpha,
    betanls2_raw      = beta
  ) %>%
  dplyr::mutate(
    alpha_se_nls2      = sqrt(alpha_variance),
    beta_variance_nls2 = beta_variance
  ) %>%
  dplyr::select(
    substudy_num,
    alpha_nls2,
    alpha_se_nls2,
    betanls2_raw,
    beta_variance_nls2
  )

# Some historical CSVs arrived without a trailing newline, which breaks fast readr
# calls on certain systems. `fix_eof()` normalises those files in-place once.
fix_eof <- function(path) {
  txt <- readLines(path, warn = FALSE)
  if (nzchar(txt[length(txt)]) && !grepl("\n$", paste(txt, collapse = "\n"))) {
    writeLines(c(txt, ""), path)
  }
}
fix_eof(here::here("data", "wilson_osenberg_2002.csv"))
fix_eof(here::here("data", "shima_osenberg2003.csv"))

# 2.2 Additional studies without raw NLS fits
schmitt2007_raw <- read.csv(here::here("data", "schmitt_holbrook_2007.csv")) %>%
  dplyr::rename(
    raw_beta    = mean_perm2,  # mean_perm2 holds the beta estimate
    raw_beta_se = se_perm2      # standard error for beta
  ) %>%
  dplyr::transmute(
    substudy_num,
    betanls2_raw      = raw_beta,
    beta_variance_nls2 = raw_beta_se^2
  )

wilson2002_raw <- read.csv(here::here("data", "wilson_osenberg_2002.csv")) %>%
  dplyr::rename(
    raw_beta    = beta_2002,
    raw_beta_se = beta_se_2002
  ) %>%
  dplyr::transmute(
    substudy_num,
    betanls2_raw      = raw_beta,
    beta_variance_nls2 = raw_beta_se^2
  )

shima2003_raw <- read.csv(here::here("data", "shima_osenberg2003.csv")) %>%
  dplyr::rename(
    betanls2_raw       = beta,
    beta_variance_nls2 = beta_variance
  ) %>%
  dplyr::select(substudy_num, betanls2_raw, beta_variance_nls2)

# Bind all parameter sources
params_all <- dplyr::bind_rows(
  params_primary,
  schmitt2007_raw,
  wilson2002_raw,
  shima2003_raw
)

# 2.3 Load study covariates and merge
covariates_raw <- read.csv(here::here("data", "covariates-2024-09-30.csv")) %>%
  dplyr::filter(use_2024 == "yes")
all_dat <- dplyr::right_join(params_all, covariates_raw, by = "substudy_num")

# 2.4 Compute average study duration and attach to only those substudy rows in all_dat
# Duration is re-constructed from the raw time series so that every downstream
# model uses the same notion of “length of experiment” (mean across replicates).
duration_df <- read.csv(here::here("data", "all_studies_looped-2024-09-11.csv")) %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::summarize(duration = mean(t, na.rm = TRUE), .groups = "drop")

#here: left_join keeps only all_dat rows (147), adding duration where available
all_dat2 <- all_dat %>%
  dplyr::left_join(duration_df, by = "substudy_num") %>%
  dplyr::rename(beta_hat = betanls2_raw)

# -----------------------------------------------------------------------------
# 2.5 Estimate densities (mean & median) per substudy
# -----------------------------------------------------------------------------
density_stats <- read.csv(here::here("data", "all_studies_looped-2024-09-11.csv")) %>%
  dplyr::filter(substudy_num %in% all_dat2$substudy_num) %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::summarize(
    mean_density   = mean(n0_m2, na.rm = TRUE),
    median_density = median(n0_m2, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 2.6 Incorporate manual density overrides
# -----------------------------------------------------------------------------
extra_densities <- read.csv(here::here("data", "manual_densities.csv"))

density_merged <- dplyr::full_join(density_stats, extra_densities, by = "substudy_num") %>%
  dplyr::mutate(
    mean_density = dplyr::coalesce(mean_density.x, mean_density.y)
  ) %>%
  dplyr::select(substudy_num, mean_density, median_density)

all_dat2 <- all_dat2 %>%
  dplyr::left_join(density_merged, by = "substudy_num")

# Save merged covariates
readr::write_csv(
  all_dat2,
  here::here("output", "merged-covariates-2024-10-21.csv"),
  na = ""
)

# -----------------------------------------------------------------------------
# 2.7 Derive additional variables
# -----------------------------------------------------------------------------
all_dat2 <- all_dat2 %>%
  dplyr::mutate(
    tropical           = if_else(abs(lat_deci) > 30, "Temperate", "Tropical"),
    betanls2_raw_cm    = beta_hat * 1e4,                  # per m² → per cm²
    betanlsvar_raw_cm  = beta_variance_nls2 * 1e8,        # variance scales by (1e4)^2
    
    # Transform on asinh scale
    betanls2_asinh     = asinh(betanls2_raw_cm),
    
    # Delta-method variance on asinh scale: Var[g(X)] ≈ Var(X) / (1 + X^2)
    betanlsvar_asinh   = betanlsvar_raw_cm / (1 + betanls2_raw_cm^2),
    
    # Standard error on asinh scale
    betanlsse_asinh    = sqrt(betanlsvar_asinh)
  )



message("Unique studies:   ", dplyr::n_distinct(all_dat2$study_num))
message("Unique species:   ", dplyr::n_distinct(all_dat2$g_sp))
message(
  "Duration range:   ",
  paste(range(all_dat2$duration, na.rm = TRUE), collapse = " – ")
)



# -----------------------------------------------------------------------------
# 3. Phylogenetic Trees Loading & Pruning
# -----------------------------------------------------------------------------
# Different analyses need trees with identical tip sets.  We therefore read the
# candidate trees, pick the default (reef fishes), and prune it down to only the
# species present in `all_dat2`.  `pruned_tree` is saved in the parent environment
# so later scripts can reuse it without repeating the pruning logic.

# 3.1 Read in raw phylogenies
reef_fish_tree  <- ape::read.tree(here::here("data", "Reef_fish_all.tacted.newick.tre"))
actino_tree     <- ape::read.tree(here::here("data", "actinopt_12k_treePL.tre"))
custom_tree_12k <- ape::read.tree(here::here("data", "1.newick.tre"))
# fish_phylo_tree <- fishtree_phylogeny()

# 3.2 Choose which tree to prune for downstream analyses
#     (swap in any of the above: e.g. fish_phylo_tree)
phylo_tree <- reef_fish_tree

# 3.3 Ensure species column exists
if (!"g_sp" %in% names(all_dat2)) {
  stop("`all_dat2$g_sp` not found – make sure Section 2 has created all_dat2 with a 'g_sp' column")
}

# 3.4 Determine which tips to keep/drop
keep_sp <- intersect(phylo_tree$tip.label, unique(all_dat2$g_sp))
drop_sp <- setdiff(phylo_tree$tip.label, keep_sp)

# 3.5 Prune away tips not in your data
pruned_tree <- ape::drop.tip(phylo_tree, drop_sp)

message(
  "Pruned phylogeny: kept ", length(keep_sp),
  " species; dropped ", length(drop_sp), " species."
)


# -----------------------------------------------------------------------------
# 4. Load and Merge Species Metadata
# -----------------------------------------------------------------------------
# Species-level trait data (length, weight) are stored in an Excel workbook.
# We join by `(g_sp, study_num, substudy_num)` so replicated species across
# studies inherit the right measurements.
meta_phy_df <- readxl::read_excel(here::here("data", "unique_species_studies.xlsx")) %>%
  dplyr::rename(
    max_length_cm = `max_len (cm)`,
    max_weight_g  = `Max_wt (ga)`
  )
all_dat2 <- all_dat2 %>%
  dplyr::left_join(meta_phy_df, by = c("g_sp", "study_num", "substudy_num")) %>%
  dplyr::mutate(max_length_density = max_length_cm * mean_density)

# -----------------------------------------------------------------------------
# 5. Unique Species-Study Combinations
# -----------------------------------------------------------------------------
species_study_df <- all_dat2 %>%
  dplyr::distinct(g_sp, study_num, substudy_num) %>%
  dplyr::arrange(g_sp, study_num, substudy_num)

readr::write_csv(
  species_study_df,
  here::here("output", "unique_species_studies.csv")
)
message("Saved unique species-study combinations: ", nrow(species_study_df), " rows")

# -----------------------------------------------------------------------------
# 6. Summary Statistics
# -----------------------------------------------------------------------------
message("Number of unique species: ", dplyr::n_distinct(all_dat2$g_sp))
message("Number of unique studies:   ", dplyr::n_distinct(all_dat2$study_num))
message("Number of unique substudies:", dplyr::n_distinct(all_dat2$substudy_num))

# -----------------------------------------------------------------------------
# 7. Publication Information Tables
# -----------------------------------------------------------------------------
study_substudy_info_df <- all_dat2 %>%
  dplyr::select(study_num, substudy_num, Authors, Article.Title, Source.Title, Publication.Year) %>%
  dplyr::distinct()
study_info_df <- all_dat2 %>%
  dplyr::select(study_num, Authors, Article.Title, Source.Title, Publication.Year) %>%
  dplyr::distinct()

readr::write_csv(study_substudy_info_df, here::here("output", "unique_study_substudy_info.csv"))
readr::write_csv(study_info_df, here::here("output", "unique_study_info.csv"))

study_substudy_info_df %>%
  gt::gt() %>%
  gt::tab_header(
    title    = "Study / Substudy → Publication Details",
    subtitle = "Authors • Article Title • Journal • Year"
  ) %>%
  gt::cols_label(
    study_num        = "Study #",
    substudy_num     = "Substudy #",
    Authors          = "Authors",
    Article.Title    = "Article Title",
    Source.Title     = "Journal",
    Publication.Year = "Year"
  ) %>%
  gt::gtsave(filename = here::here("output", "study_substudy_info.html"))

# -----------------------------------------------------------------------------
# 8. Meta-Analysis Citation Table (vectorized + robust author parsing)
# -----------------------------------------------------------------------------

# Helpers ----------------------------------------------------------------------

# split string by first present delimiter among these
.split_first <- function(x, delims = c(";", "\\s+&\\s+", "\\s+and\\s+")) {
  for (d in delims) {
    if (stringr::str_detect(x, d)) {
      return(stringr::str_split(x, d, n = 2, simplify = TRUE)[, 1])
    }
  }
  x
}

# Given a single author token, return "Surname, INITIALS"
.norm_one_author <- function(tok) {
  if (is.null(tok) || is.na(tok)) {
    tok <- ""
  }
  tok <- stringr::str_squish(tok)
  if (tok == "") return(NA_character_)
  # If token contains a comma, assume "Last, First/Initials"
  if (stringr::str_detect(tok, ",")) {
    parts <- stringr::str_split(tok, ",\\s*", n = 2, simplify = TRUE)
    last  <- stringr::str_to_title(parts[,1])
    rhs   <- ifelse(ncol(parts) >= 2, parts[,2], "")
    inits <- stringr::str_replace_all(rhs, "[^A-Za-z]", "")        # keep letters only
    inits <- stringr::str_to_upper(inits)
    inits <- paste(stringr::str_split(inits, "", simplify = TRUE), collapse = "")
    inits <- paste(stringr::str_split(inits, "", simplify = TRUE), collapse = "") # ensure scalar
    if (nzchar(inits)) paste0(last, ", ", inits) else last
  } else {
    # Assume "First Middle Last" → "Last, FM"
    parts <- stringr::str_split(tok, "\\s+", simplify = TRUE)
    if (length(parts) == 0) return(NA_character_)
    last  <- stringr::str_to_title(parts[, ncol(parts), drop = TRUE])
    given <- parts[, seq_len(ncol(parts)-1), drop = TRUE]
    if (length(given) == 0) return(last)
    initials <- stringr::str_to_upper(stringr::str_sub(given, 1, 1))
    initials <- paste0(initials, collapse = "")
    paste0(last, ", ", initials)
  }
}

# Vectorized wrapper: "A; B" / "A & B" / "A and B" → first author normalized
first_author_vec <- function(authors_chr) {
  purrr::map_chr(authors_chr, ~{
    if (is.na(.x) || .x == "") return(NA_character_)
    first_tok <- .split_first(.x)
    .norm_one_author(first_tok)
  })
}

# Build table ------------------------------------------------------------------

meta_table <- all_dat2 %>%
  dplyr::mutate(
    Author = first_author_vec(Authors),
    Year   = Publication.Year,
    Title  = stringr::str_to_sentence(Article.Title),
    Journal= stringr::str_to_title(Source.Title),
    Citation = dplyr::if_else(
      !is.na(Author) & !is.na(Title) & !is.na(Journal) & !is.na(Year),
      paste0(Author, " (", Year, "). ", Title, ". *", Journal, "*."),
      NA_character_
    )
  ) %>%
  dplyr::select(
    Study        = study_num,
    Substudy     = substudy_num,
    Beta         = betanls2_raw_cm,
    Variance     = betanlsvar_raw_cm,
    Author,
    Year,
    Citation,
    DOI,
    g_sp,
    Duration     = duration,
    mean_density
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Year, Author) %>%
  dplyr::mutate(
    Beta     = round(Beta, 3),
    Variance = round(Variance, 5)
  )

# Save CSV + HTML --------------------------------------------------------------

readr::write_csv(
  meta_table,
  here::here("figures", "meta_analysis_table.csv"),
  na = ""
)

meta_table %>%
  gt::gt() %>%
  gt::tab_header(
    title    = gt::md("**Summary of Studies Included in Meta-Analysis**"),
    subtitle = "Effect Size, Variance, and Metadata by Substudy"
  ) %>%
  gt::fmt_number(columns = c("Beta", "Variance"), decimals = 3) %>%
  gt::cols_label(
    Study        = "Study #",
    Substudy     = "Substudy #",
    Beta         = "Effect Size (β)",
    Variance     = "Variance (asinh)",
    Author       = "First Author",
    Year         = "Year",
    DOI          = "DOI",
    g_sp         = "Species",
    Duration     = "Duration (t)",
    mean_density = "Mean Density (m^2)"
  ) %>%
  gt::tab_options(
    table.font.size           = gt::px(13),
    column_labels.font.weight = "bold",
    heading.align             = "center"
  ) %>%
  gt::gtsave(filename = here::here("output", "meta_analysis_table.html"))

# End of 1_data_phylogeny_loading.R
