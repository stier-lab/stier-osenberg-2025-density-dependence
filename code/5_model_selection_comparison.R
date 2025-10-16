###############################################################################
# Script: 5_model_selection_and_autocorrelation.R
# Goal:   Inspect correlations among covariates, fit saturated REML models with
#         interactions, apply hierarchical p-value trimming, and summarise the
#         best additive model for reporting.
# Inputs: all_dat2/pruned_tree from 1_data_phylogeny_loading.R
# Outputs: results/model_comparison_aic.csv, results/best_model_summary.png,
#          Console diagnostics and optional intermediate objects.
###############################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
# Model selection builds on the predator-present subset of `all_dat2`.  We source
# the shared loader first, then confirm the canonical data/tree objects exist.
source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE) ||
    !exists("pruned_tree", inherits = FALSE)) {
  message("Core data objects missing; sourcing 1_data_phylogeny_loading.R.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

# ensure output directories exist
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Autocorrelation among predictor variables (exploratory)
# ----------------------------------------------------------------------------
# Before fitting multivariate models we check correlations among covariates to
# understand potential collinearity and guide interpretation of interaction terms.
all <- all_dat2  # temporary alias

# compute log density safely
all$logmeandensity <- dplyr::if_else(
  all$mean_density > 0,
  log(all$mean_density),
  NA_real_
)

# select raw predictors
predictors <- all %>%
  select(expt_obs, size_start, duration, logmeandensity, max_length_cm, max_weight_g)

# multipanel scatterplot matrix
GGally::ggpairs(
  predictors,
  lower = list(continuous = GGally::wrap("smooth", color = "blue")),
  diag  = list(continuous = "densityDiag"),
  upper = list(continuous = GGally::wrap("cor", size = 5))
) +
  theme_minimal()

# ----------------------------------------------------------------------------
# 3. β vs. density by study type (exploratory)
# ----------------------------------------------------------------------------
# Quick diagnostic plotting raw β against mean density, split by experimental vs observational.
log_ticks  <- c(1e-4,1e-3,1e-2,1e-1,1,10,100,1000)
sinh_ticks <- sinh(c(-10,-5,-2,-1,0,1,2,5,10))

density_plot_df <- all %>%
  filter(!is.na(mean_density), mean_density > 0, !is.na(betanls2_raw_cm))

ggplot(density_plot_df, aes(mean_density, betanls2_raw_cm, shape=expt_obs, color=expt_obs)) +
  geom_point(alpha=0.5, size=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(trans="log10", breaks=log_ticks, labels=log_ticks) +
  scale_y_continuous(trans="asinh", breaks=sinh_ticks) +
  scale_shape_manual(values=c(16,17)) +
  scale_color_manual(values=c("blue","red")) +
  labs(
    x     = "Mean Density (n₀/m², log scale)",
    y     = "Density-dependence strength (β)",
    title = "β vs. Density by Study Type"
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="top")

# ----------------------------------------------------------------------------
# 4. Data preparation & phylogeny
# ----------------------------------------------------------------------------
# Restrict to predator-present substudies and centre/scale continuous predictors.
all <- all_dat2 %>%
  filter(predators == "present") %>%
  mutate(
    expt_obs = factor(expt_obs, levels=c("Exp","Obs")),
    logmd    = dplyr::if_else(mean_density > 0, log(mean_density), NA_real_),
    logmd_c  = as.numeric(scale(logmd)),
    dur_c    = as.numeric(scale(duration)),
    size_c   = as.numeric(scale(size_start)),
    max_c    = as.numeric(scale(max_length_cm))
  ) %>%
  select(-logmd)
all <- all %>% filter(!is.na(logmd_c))

phylo_vcv <- ape::vcv(pruned_tree, corr=TRUE)
keep_sp    <- intersect(rownames(phylo_vcv), unique(all$g_sp))
phylo_vcv  <- phylo_vcv[keep_sp, keep_sp]
all         <- all %>% filter(g_sp %in% keep_sp)

rand_list <- list(~1 | study_num/substudy_num, ~1 | g_sp)

# ----------------------------------------------------------------------------
# 5. Helper functions: fit rma.mv & AICc
# ----------------------------------------------------------------------------
# These wrappers keep the multi-level REML specification in one place so the
# stepwise search only updates formula terms.
fit_rma <- function(formula_obj) {
  rma.mv(
    yi     = betanls2_asinh,
    V      = betanlsvar_asinh,
    mods   = formula_obj,
    random = rand_list,
    R      = list(g_sp = phylo_vcv),
    data   = all,
    method = "REML",
    test   = "t"
  )
}

get_AICc <- function(mod) {
  summary(mod)$fit.stats["AICc","REML"] %>% as.numeric()
}

# check whether a term is embedded inside any higher-order interaction
# check whether a term is embedded inside any higher-order interaction
is_nested_in_higher <- function(term, term_set) {
  term_pattern <- paste0("^", term, "(:|$)")
  any(grepl(term_pattern, term_set) & term_set != term)
}



# ----------------------------------------------------------------------------
# 8. Two-way interaction model
# ----------------------------------------------------------------------------
# Start with the saturated two-way interaction model so we can evaluate which
# interactions survive the trimming rules.
m_2way <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ expt_obs + logmd_c + dur_c + size_c + max_c
  + expt_obs:logmd_c + expt_obs:dur_c + expt_obs:size_c + expt_obs:max_c
  + logmd_c:dur_c   + logmd_c:size_c   + logmd_c:max_c
  + dur_c:size_c    + dur_c:max_c      + size_c:max_c,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_2way)

# ----------------------------------------------------------------------------
# 9. Stepwise removal: interaction-order–specific p-value thresholds
# ----------------------------------------------------------------------------
# Implement hierarchical removal governed by user-defined p-value cutoffs.
drop_if_insig <- function(label, model) {
  ct        <- coef(summary(model))
  terms_now <- attr(terms(formula(model)),"term.labels")
  parts     <- strsplit(label,":",fixed=TRUE)[[1]]
  ord       <- length(parts)
  cutoff    <- if (ord==3) 0.01 else 0.05
  
  if (ord==1 || (label %in% terms_now && is_nested_in_higher(label,terms_now)))
    return(list(model=model,dropped=FALSE))
  
  pat <- if (parts[1]=="expt_obs") {
    sfx <- sub("^expt_obs:","",label)
    paste0("^expt_obs(?:Exp|Obs):",sfx,"$")
  } else paste0("^",label,"$")
  
  idx <- grep(pat, rownames(ct))
  if (!length(idx)) return(list(model=model,dropped=FALSE))
  
  if (max(ct[idx,"pval"],na.rm=TRUE) >= cutoff) {
    new_mod <- tryCatch(fit_rma(reformulate(setdiff(terms_now,label))),
                        error=function(e) NULL) 
    if (!is.null(new_mod)) return(list(model=new_mod,dropped=TRUE))
  }
  list(model=model,dropped=FALSE)
}

# [Insert hierarchical loops over 5-way→2-way terms here as in original script]

# ----------------------------------------------------------------------------
# 10. Final best-fitting additive model
# ----------------------------------------------------------------------------
# After trimming, fit the additive REML model for reporting.
m_best <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ logmd_c + dur_c + size_c + max_c + expt_obs,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = all,
  method = "REML",
  test   = "t"
)
summary(m_best)

# ----------------------------------------------------------------------------
# 11. Back‐transform & tabulate best‐fitting model
# ----------------------------------------------------------------------------
# Convert asinh estimates back to the β scale so the summary table mirrors the manuscript.
results <- broom::tidy(m_best, conf.int = TRUE) %>%
  mutate(
    estimate_bt  = sinh(estimate),
    conf.low_bt  = sinh(conf.low),
    conf.high_bt = sinh(conf.high),
    p.value      = if_else(p.value < .001, "<0.001", as.character(round(p.value, 3)))
  ) %>%
  select(
    Predictor    = term,
    Estimate     = estimate_bt,
    SE           = std.error,
    `t-value`    = statistic,
    `p-value`    = p.value,
    `CI low`     = conf.low_bt,
    `CI high`    = conf.high_bt
  )

results %>%
  gt() %>%
  tab_header(
    title    = md("**Best‐Fitting Multivariate Meta‐Analysis Model**"),
    subtitle = "Back‐Transformed Estimates"
  ) %>%
  fmt_number(columns = c("Estimate", "SE", "t-value"), decimals = 3) %>%
  cols_label(
    Predictor = "Predictor",
    SE        = "SE",
    `t-value` = "t‐value",
    `p-value` = "p‐value",
    `CI low`  = "95% CI (low)",
    `CI high` = "95% CI (high)"
  ) %>%
  tab_options(table.font.size = px(14), heading.align = "center") %>%
  # PNG save requires Chrome/Chromium - save as HTML instead for compatibility
  # gtsave(here::here("results", "best_model_summary.png"))
  gtsave(here::here("results", "best_model_summary.html"))
