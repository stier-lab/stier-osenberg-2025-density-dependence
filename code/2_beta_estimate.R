###############################################################################
# Script: 2_beta_estimate.R
# Goal:   Fit the baseline REML meta-analysis for β on the asinh scale, report
#         pooled estimates/variance components, and produce the combined
#         orchard + Beverton–Holt figure used in the main text.
# Inputs: Requires all_dat2/pruned_tree in memory plus
#         data/all_studies_looped-2024-09-11.csv,
#         output/combined_results_2025-10-15.csv,
#         output/merged-covariates-2024-10-21.csv.
# Outputs: Console summaries, figures/figure1_orchard_bh_plot.{pdf,png}
#          (directories created if missing).
###############################################################################

# ----------------------------------------------------------------------------
# 1. Setup & Preprocessing
# ----------------------------------------------------------------------------
# Load shared packages and, if necessary, create the key data objects.
# `all_dat2` and `pruned_tree` are produced by 1_data_phylogeny_loading.R so
# every downstream script can rely on the same data frame and tree instance.
source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE)) {
  message("all_dat2 not found in environment; sourcing data loader.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

dir.create(here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

#helper
inv_asinh <- function(x) sinh(x)

# Filter out incomplete cases for model fitting
model_data <- all_dat2 %>%
  dplyr::filter(
    !is.na(study_num),
    !is.na(substudy_num),
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  )

# ----------------------------------------------------------------------------
# 2. Fit mixed-effects model
# ----------------------------------------------------------------------------
# The REML model treats substudies nested within studies as random effects.
# We work on the asinh scale for numerical stability, mirroring the later
# publication-bias diagnostics.
meta_model <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  random = list(~1 | study_num/substudy_num),
  data   = model_data,
  method = "REML",
  test   = "t"
)


mu_asinh    <- as.numeric(coef(meta_model))
mu_asinh_ci <- c(meta_model$ci.lb, meta_model$ci.ub)
mu_beta     <- inv_asinh(mu_asinh)
mu_beta_ci  <- inv_asinh(mu_asinh_ci)

QE  <- meta_model$QE; QEp <- meta_model$QEp
sig2 <- meta_model$sigma2; names(sig2) <- c("study","substudy")
k  <- meta_model$k                       # actual number of effects used
p  <- length(coef(meta_model))           # number of fixed effects (1 here)
I2 <- max(0, (QE - (k - p))/QE) * 100


# Helper to format tiny p-values robustly
fmt_p <- function(p, digits = 3) {
  # Use base R's formatter but push the eps very low
  format.pval(p, digits = digits, eps = 1e-300)
}

# For QE, recompute p from the chi-square distribution with log.p to avoid underflow.
p_QE <- {
  df <- meta_model$QEdf
  # try log-scale first; if finite, exponentiate; if -Inf, fall back to a "<" threshold
  lp <- stats::pchisq(QE, df = df, lower.tail = FALSE, log.p = TRUE)
  if (is.finite(lp)) exp(lp) else 0
}

cat("\n[D.1 Robustness & Heterogeneity]\n",
    sprintf("  mu(asinh) = %.3f (%.3f, %.3f)\n", mu_asinh, mu_asinh_ci[1], mu_asinh_ci[2]),
    sprintf("  mu(beta)  = %.2f (%.2f, %.2f)\n",  mu_beta,  mu_beta_ci[1],  mu_beta_ci[2]),
    sprintf("  QE = %.0f; p = %s\n",
            QE,
            if (p_QE > 0) fmt_p(p_QE, digits = 3) else "< 1e-300"),
    sprintf("  sigma^2: study=%.4f, substudy=%.4f\n", sig2["study"], sig2["substudy"]),
    sprintf("  I2 (pseudo, asinh) ≈ %.2f%%\n", I2),
    sep = "")

# ----------------------------------------------------------------------------
# 3) Fixed-effect (intercept) results & back-transform  (robust to version diffs)
# ----------------------------------------------------------------------------
sm <- summary(meta_model)

# Intercept row name is "intrcpt" in metafor
b_mat   <- sm$b
ci_lb   <- sm$ci.lb
ci_ub   <- sm$ci.ub
p_mu    <- sm$pval[1]

stopifnot(identical(rownames(b_mat)[1], "intrcpt"))

mu_asinh    <- as.numeric(b_mat[1, 1])
mu_asinh_ci <- c(ci_lb[1], ci_ub[1])

fixed_effect_results <- tibble::tibble(
  component      = "fixed_effect",
  estimate_asinh = mu_asinh,
  ci_lb_asinh    = mu_asinh_ci[1],
  ci_ub_asinh    = mu_asinh_ci[2],
  estimate       = inv_asinh(mu_asinh),
  ci_lower       = inv_asinh(mu_asinh_ci[1]),
  ci_upper       = inv_asinh(mu_asinh_ci[2]),
  p_value        = sm$pval[1]
)

# convenience values for plotting later
back_transformed_coef       <- fixed_effect_results$estimate
back_transformed_conf_lower <- fixed_effect_results$ci_lower
back_transformed_conf_upper <- fixed_effect_results$ci_upper

# ----------------------------------------------------------------------------
# 4) Variance components (delta-method back-transform)  (preserve correct labels)
# ----------------------------------------------------------------------------
# Use the model’s own random-factor labels to avoid order bugs
vc_vals   <- sm$sigma2
vc_labels <- meta_model$s.names
if (length(vc_labels) == length(vc_vals)) {
  names(vc_vals) <- vc_labels
} else {
  names(vc_vals) <- paste0("sigma2[", seq_along(vc_vals), "]")
}

# Delta-method scale factor for X = sinh(Y), with Y on asinh scale
scale_dm <- cosh(mu_asinh)^2

variance_results <- tibble::enframe(vc_vals, name = "component", value = "var_asinh") |>
  dplyr::mutate(
    estimate = var_asinh * scale_dm,  # approx variance on original (beta) scale
    ci_lower = NA_real_,              # non-trivial to compute; leave NA unless you derive them
    ci_upper = NA_real_
  ) |>
  dplyr::select(component, estimate, ci_lower, ci_upper)
  
# ----------------------------------------------------------------------------
# 5. Compile & (optionally) Save Results
# ----------------------------------------------------------------------------
results_summary <- dplyr::bind_rows(fixed_effect_results, variance_results)
print(results_summary)
# readr::write_csv(results_summary, here::here("results", "beta_meta_summary.csv"))

# ----------------------------------------------------------------------------
# 6. Combined Orchard & Beverton–Holt Figure
# ----------------------------------------------------------------------------
# Figure 1 pairs an orchard plot (meta-analysis overview) with representative
# Beverton–Holt fits.  Everything below reshapes the raw digitised time-series,
# selects four quantile-matched substudies, and overlays the BH curve.

# 6.1 Beverton–Holt prep (kept your formatting/colors/labels)
all_data <- readr::read_csv(
  here("data", "all_studies_looped-2024-09-11.csv"),
  show_col_types = FALSE
)
combined_results <- readr::read_csv(
  here("output", "combined_results_2025-10-15.csv"),
  show_col_types = FALSE
)
covariates <- readr::read_csv(
  here("output", "merged-covariates-2024-10-21.csv"),
  show_col_types = FALSE
)

bh_f <- function(alpha, beta, t_days, x) {
  (x * exp(-alpha * t_days)) / (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
}

usedf <- covariates %>%
  dplyr::filter(use_2024 == "yes", predators == "present") %>%
  dplyr::select(substudy_num, family, predators)

all_data    <- dplyr::inner_join(all_data, usedf, by = "substudy_num")
merged_data <- dplyr::inner_join(all_data, combined_results, by = "substudy_num")

filtered_data <- merged_data %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::filter(n() >= 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(beta_strength = beta)

beta_quantiles <- stats::quantile(combined_results$beta,
                                  probs = c(0.10, 0.5, 0.75, 0.95),
                                  na.rm = TRUE)

closest_values <- sapply(beta_quantiles, function(x) {
  filtered_data$beta[which.min(abs(filtered_data$beta - x))]
})

matching_substudies <- filtered_data %>%
  dplyr::filter(beta %in% closest_values) %>%
  dplyr::distinct(substudy_num)

final_data <- filtered_data %>%
  dplyr::filter(substudy_num %in% matching_substudies$substudy_num) %>%
  dplyr::arrange(beta)

desired_order <- c(212, 249, 76, 256)
final_data <- final_data %>%
  dplyr::filter(substudy_num %in% desired_order) %>%
  dplyr::mutate(substudy_num = factor(substudy_num, levels = desired_order))

predicted_data <- final_data %>%
  dplyr::group_by(substudy_num) %>%
  dplyr::mutate(settler_range = list(seq(0, max(n0_m2, na.rm = TRUE), length.out = 100))) %>%
  tidyr::unnest(cols = c(settler_range)) %>%
  dplyr::mutate(predicted_recruits = bh_f(alpha, beta, t, settler_range))

facet_labels <- final_data %>%
  dplyr::mutate(genus_species = gsub("_", " ", genus_species)) %>%
  dplyr::mutate(
    facet_label = paste0(
      "*", genus_species, "*\n",
      " *β* = ", round(beta * 1e4, 0), " "
    )
  ) %>%
  dplyr::distinct(substudy_num, facet_label)

beverton_holt_plot <- ggplot() +
  geom_point(
    data = final_data, aes(x = n0_m2, y = nt_m2),
    size = 3, alpha = 0.8, color = "#2C3E50",
    position = position_jitter(width = 0.05, height = 0)
  ) +
  geom_line(
    data = predicted_data, aes(x = settler_range, y = predicted_recruits),
    size = 1.5, color = "#E74C3C"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#95A5A6", linewidth = 1) +
  facet_wrap(
    ~substudy_num, ncol = 1, scales = "free",
    labeller = labeller(substudy_num = setNames(facet_labels$facet_label,
                                                facet_labels$substudy_num))
  ) +
  labs(
    x = expression(paste("Initial Density (Settlers per ", m^{-2}, ")")),
    y = expression(paste("Final Density (Recruits per ", m^{-2}, ")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = ggtext::element_markdown(size = 14, face = "bold"),
    axis.text  = element_text(size = 12, color = "#34495E"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white")
  )

# 6.2 Orchard plot (manual) — formatting retained
precision_breaks <- c(1e-12, 1e4, 1e6, 1e8, 1e12)
precision_labels <- c(expression(10^-12), expression(10^4), expression(10^6),
                      expression(10^8), expression(10^12))

global_estimate <- data.frame(
  back_transformed_coef = back_transformed_coef,
  back_transformed_conf_lower = back_transformed_conf_lower,
  back_transformed_conf_upper = back_transformed_conf_upper,
  y_position = "Beta"
)

shape_mapping <- c("212", "249", "76", "256")

beta_all_all <- all_dat2 %>%
  dplyr::select(
    beta_hat,
    beta_variance_nls2,
    study_num,
    substudy_num,
    betanlsvar_raw_cm,
    betanls2_raw_cm,
    betanls2_asinh,
    betanlsvar_asinh
  ) %>%
  dplyr::mutate(
    precision   = 1 / beta_variance_nls2,
    point_shape = ifelse(as.character(substudy_num) %in% shape_mapping, 22, 16)
  )



orchard_manual <- ggplot(beta_all_all, aes(
  x    = betanls2_raw_cm,
  y    = "Beta",
  size = sqrt(precision)
)) +
  ggbeeswarm::geom_quasirandom(
    aes(alpha = precision),
    shape = 21, fill = "#377EB8", color = "#377EB8",
    width = 0.15, alpha = 0.7
  ) +
  geom_segment(
    data = global_estimate,
    aes(x = back_transformed_conf_lower,
        xend = back_transformed_conf_upper,
        y = y_position, yend = y_position),
    color = "black", linewidth = 1.2, inherit.aes = FALSE
  ) +
  geom_point(
    data  = global_estimate,
    aes(x = back_transformed_coef, y = y_position),
    color = "black", size = 5, shape = 21,
    fill = "#377EB8", stroke = 1.5,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 1) +
  geom_vline(xintercept = c(-100, -10, 0, 10, 100, 1000),
             linetype = "dashed", color = "gray70", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", linewidth = 1) +
  labs(
    x    = expression(paste("Strength of density-dependent mortality, ",
                            beta, " (", cm^2, " ", fish^{-1}, " ", day^{-1}, ")")),
    y    = NULL,
    size = "Precision (1/SE)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.y      = element_blank(),
    legend.position   = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    axis.text         = element_text(size = 14),
    axis.title.x      = element_text(size = 16, face = "bold"),
    panel.grid.major  = element_line(color = "gray90"),
    panel.grid.minor  = element_blank(),
    panel.background  = element_rect(fill = "white")
  ) +
  scale_x_continuous(
    trans  = "asinh",
    limits = c(-500, 10000),
    breaks = c(-500, -100, -10, -1, 0, 1, 10, 100, 1000, 2000, 10000)
  ) +
  scale_size_continuous(
    name   = "Precision",
    range  = c(3, 10),
    breaks = sqrt(precision_breaks),
    labels = precision_labels
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        shape  = 16,
        fill   = "#377EB8",
        color  = "#377EB8",
        stroke = 0
      )
    ),
    alpha = "none"
  ) +
  coord_flip()

# 6.3 Combine and save figure (unchanged settings)
bplot2 <- cowplot::plot_grid(
  orchard_manual, beverton_holt_plot,
  labels = c('A', 'B'), label_size = 12, ncol = 2,
  rel_heights = c(1, 2)
)
print(bplot2)

ggplot2::ggsave(
  filename = here::here("figures", "figure1_orchard_bh_plot.pdf"),
  plot     = bplot2,
  width    = 10,
  height   = 11,
  units    = "in",
  dpi      = 300,
  device   = cairo_pdf,   # keep if Cairo is available
  bg       = "white"
)
