###############################################################################
# Script: 8_appendix_E_sensitivity.R
# Goal:   Comprehensive sensitivity analysis for Appendix E including leave-one-out
#         cross-validation, funnel plots, Egger regression, and null-add sensitivity.
#         All analyses conducted on the asinh(β) scale to match the main meta-analysis.
# Inputs: all_dat2 from 1_data_phylogeny_loading.R
# Outputs:
#   - Figure E.1A: Leave-one-out cross-validation
#   - Figure E.1B: Funnel plot with 95%/99% confidence bounds
#   - Figure E.1C: Egger regression test for publication bias
#   - Figure E.1D: Null-add sensitivity analysis (multi-quantile)
#   - Figure E.1 (composite): All panels combined
#   - Table E.1: Summary statistics (tables/Table_E1_sensitivity.csv)
###############################################################################

source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE)) {
  message("all_dat2 not found; sourcing 1_data_phylogeny_loading.R.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

# Publication color for points
pub_blue <- "#2C7BB6"  # nice scientific blue


# ===== Add this utility once near the top (after library(...)) =====
# A gentle, zero-friendly compressing transform for axes (NOT log)

asinh_trans <- function() {
  scales::trans_new(
    name     = "asinh",
    transform = function(x) asinh(x),
    inverse   = function(x) sinh(x)
  )
}

# SI-style numeric labels (compat with modern 'scales')
# Uses label_number() with scale_cut = cut_si("") to mimic label_number_si()
lab_si <- function(accuracy = 1, ...) {
  scales::label_number(accuracy = accuracy, scale_cut = scales::cut_si(""), ...)
}

# ===========================
# 0) Load all_dat2 (canonical)
# ===========================
load_loader <- function() {
  tried <- c(
    here::here("code", "1_data_phylogeny_loading.R"),
    "code/1_data_phylogeny_loading.R",
    "1_data_phylogeny_loading.R"
  )
  for (p in tried) if (file.exists(p)) {
    message(sprintf("[data] sourcing %s", p))
    source(p, local = TRUE)
    return(invisible(TRUE))
  }
  stop("Could not find 1_data_phylogeny_loading.R in: ", paste(tried, collapse=" | "))
}
load_loader()
stopifnot(exists("all_dat2"), is.data.frame(all_dat2))

# ===========================
# Data used in ALL analyses
# ===========================

dat <- all_dat2 %>%
  dplyr::filter(
    !is.na(study_num), !is.na(substudy_num),
    !is.na(betanls2_asinh), !is.na(betanlsvar_asinh),
    !is.na(betanls2_raw_cm), !is.na(betanlsvar_raw_cm)
  ) %>%
  dplyr::transmute(
    study_num    = as.character(study_num),
    substudy_num = as.character(substudy_num),
    
    # ---- asinh scale (analysis/model) ----
    yi_asinh = betanls2_asinh,
    vi_asinh = betanlsvar_asinh,
    se_asinh = sqrt(betanlsvar_asinh),
    
    # ---- raw scale (kept for reference; not for plotting) ----
    beta_raw = betanls2_raw_cm,
    vi_raw   = betanlsvar_raw_cm,
    se_raw   = sqrt(betanlsvar_raw_cm)
  ) %>%
  # Add Egger inputs on the asinh analysis scale (same frame)
  dplyr::mutate(
    precision = 1 / se_asinh,               # 1 / SE(asinh β)
    std_eff   = yi_asinh / se_asinh,        # asinh β / SE(asinh β)
    
    # Aliases used by metafor calls (keeps downstream code working)
    yi = yi_asinh,
    vi = vi_asinh,
    se = se_asinh
  )

# Make sure figure dir exists
fig_dir <- here::here("figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Helpers
inv_asinh <- function(x) sinh(x)

# --- Dynamic asinh-labeled axes for raw-β vectors ---
# Given a vector of raw β values, place major/minor ticks evenly in asinh-space,
# but draw the axis in raw β with labels shown asinh(β).
make_asinh_axis_from_raw_y <- function(raw_vals, n_major = 7, n_minor = 21,
                                       label_digits = 3, expand_mult = c(0.02, 0.03)) {
  rng_asinh <- range(asinh(raw_vals), na.rm = TRUE)
  maj_asinh <- pretty(rng_asinh, n = n_major)
  min_asinh <- pretty(rng_asinh, n = n_minor)
  maj_raw <- sinh(maj_asinh)
  min_raw <- sinh(min_asinh)
  scale_y_continuous(
    breaks       = maj_raw,
    labels       = formatC(maj_asinh, format = "f", digits = label_digits),
    minor_breaks = min_raw,
    expand       = expansion(mult = expand_mult)
  )
}

# --- Axis: asinh-spaced ticks, RAW-β labels (for funnel x) ---
make_asinh_spaced_axis_show_raw_x <- function(
    raw_vals, n_major = 9, n_minor = 45,
    label_func = lab_si(accuracy = 1),
    expand_mult = c(0.02, 0.03),
    limits = NULL
) {
  rng_asinh <- range(asinh(raw_vals), na.rm = TRUE)
  maj_asinh <- pretty(rng_asinh, n = n_major)
  min_asinh <- pretty(rng_asinh, n = n_minor)
  maj_raw <- sinh(maj_asinh)
  min_raw <- sinh(min_asinh)
  scale_x_continuous(
    trans        = asinh_trans(),
    breaks       = maj_raw,
    minor_breaks = min_raw,
    labels       = label_func,
    expand       = expansion(mult = expand_mult),
    limits       = limits
  )
}

# --- Axis: asinh-spaced ticks, RAW-β labels (for y axis) ---
make_asinh_spaced_axis_show_raw_y <- function(raw_vals, n_major = 9, n_minor = 45,
                                              label_func = scales::label_number(accuracy = 0.01),
                                              expand_mult = c(0.02, 0.03)) {
  rng_asinh <- range(asinh(raw_vals), na.rm = TRUE)
  maj_asinh <- pretty(rng_asinh, n = n_major)
  min_asinh <- pretty(rng_asinh, n = n_minor)
  maj_raw <- sinh(maj_asinh)
  min_raw <- sinh(min_asinh)
  scale_y_continuous(
    trans        = asinh_trans(),
    breaks       = maj_raw,
    minor_breaks = min_raw,
    labels       = label_func,
    expand       = expansion(mult = expand_mult)
  )
}

make_asinh_axis_from_raw_x <- function(raw_vals, n_major = 7, n_minor = 21,
                                       label_digits = 3, expand_mult = c(0.02, 0.03)) {
  rng_asinh <- range(asinh(raw_vals), na.rm = TRUE)
  maj_asinh <- pretty(rng_asinh, n = n_major)
  min_asinh <- pretty(rng_asinh, n = n_minor)
  maj_raw <- sinh(maj_asinh)
  min_raw <- sinh(min_asinh)
  scale_x_continuous(
    breaks       = maj_raw,
    labels       = formatC(maj_asinh, format = "f", digits = label_digits),
    minor_breaks = min_raw,
    expand       = expansion(mult = expand_mult)
  )
}


# helper: axis with ticks evenly spaced in asinh-space,
# labels shown as raw values
make_asinh_axis_raw_labels <- function(vals, n_major = 8, n_minor = 30,
                                       which = c("x","y"),
                                       expand_mult = c(0.02,0.03),
                                       label_func = lab_si(accuracy = 1)) {
  which <- match.arg(which)
  rng_asinh <- range(asinh(vals), na.rm = TRUE)
  maj_asinh <- pretty(rng_asinh, n_major)
  min_asinh <- pretty(rng_asinh, n_minor)
  maj_raw <- sinh(maj_asinh)
  min_raw <- sinh(min_asinh)
  
  if (which == "x") {
    scale_x_continuous(
      trans        = asinh_trans(),
      breaks       = maj_raw,
      minor_breaks = min_raw,
      labels       = label_func,
      expand       = expansion(mult = expand_mult)
    )
  } else {
    scale_y_continuous(
      trans        = asinh_trans(),
      breaks       = maj_raw,
      minor_breaks = min_raw,
      labels       = label_func,
      expand       = expansion(mult = expand_mult)
    )
  }
}



# theme
theme_base <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    axis.ticks = element_line(color = "grey60", linewidth = 0.3),
    axis.ticks.length = unit(3, "pt"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )


# ---- Custom Publication Theme (Adrian Style) ----
theme_pub_adrian <- function(base_size = 14, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Grid and panel
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      # Axes
      axis.ticks = element_line(color = "grey60", linewidth = 0.3),
      axis.ticks.length = unit(3, "pt"),
      axis.text = element_text(size = base_size * 0.9, color = "black"),
      axis.title.x = element_text(
        size = base_size * 1.0, margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = base_size * 1.0,
        angle = 90,                      # <--- ROTATES Y AXIS TITLE
        vjust = 0.5,
        margin = margin(r = 10)
      ),
      # Titles
      plot.title = element_text(face = "bold", size = base_size * 1.2, hjust = 0),
      plot.subtitle = element_text(size = base_size * 0.9, color = "grey30"),
      # Legend
      legend.position = "none",
      # Margins for tighter figure alignment
      plot.margin = margin(5, 5, 5, 5)
    )
}

# Optional: set it as the default for the whole script
theme_set(theme_pub_adrian())


# ===========================================================
# D.1 MODEL ROBUSTNESS & HETEROGENEITY (asinh scale)
# ===========================================================
# Everything downstream hangs off this REML fit: the pooled mean, leave-one-out
# diagnostics, funnel/Egger tests, and the null-add perturbations.
meta_REML <- rma.mv(
  yi = yi, V = vi,
  random = list(~ 1 | study_num/substudy_num),
  method = "REML",
  test = "t",
  data = dat
)

mu_asinh    <- as.numeric(coef(meta_REML))
mu_asinh_ci <- c(meta_REML$ci.lb, meta_REML$ci.ub)
mu_beta     <- inv_asinh(mu_asinh)
mu_beta_ci  <- inv_asinh(mu_asinh_ci)

QE  <- meta_REML$QE; QEp <- meta_REML$QEp
sig2 <- meta_REML$sigma2; names(sig2) <- c("study","substudy")
k <- nrow(dat); p <- 1L
I2 <- max(0, (QE - (k - p))/QE) * 100

cat("\n[D.1 Robustness & Heterogeneity]\n",
    sprintf("  mu(asinh) = %.3f (%.3f, %.3f)\n", mu_asinh, mu_asinh_ci[1], mu_asinh_ci[2]),
    sprintf("  mu(beta)  = %.2f (%.2f, %.2f)\n", mu_beta, mu_beta_ci[1], mu_beta_ci[2]),
    sprintf("  QE = %.0f; p < %.2g\n", QE, QEp),
    sprintf("  sigma^2: study=%.4f, substudy=%.4f\n", sig2["study"], sig2["substudy"]),
    sprintf("  I2 (pseudo, asinh) ≈ %.2f%%\n", I2),
    sep = "")

# ---- Figure D-1A: Leave-One-Out (REML, all data) + export ----
# Iterate over substudies, refit the REML model after removing each one, and
# store the shift in the pooled mean.  This mirrors the SI narrative exactly.
loo_fits <- purrr::map(seq_len(k), function(i) {
  fit <- rma.mv(yi, vi, random = ~1|study_num/substudy_num,
                method = "REML", test="t", data = dat[-i, ])
  tibble::tibble(
    left_out_index = i,
    mu_asinh = as.numeric(fit$b),
    ci_lb    = fit$ci.lb,
    ci_ub    = fit$ci.ub
  )
})
loo_tbl <- dplyr::bind_rows(loo_fits) %>%
  dplyr::mutate(mu_beta = sinh(mu_asinh),
                ci_lb_beta = sinh(ci_lb),
                ci_ub_beta = sinh(ci_ub))


loo_rng <- range(loo_tbl$mu_asinh, na.rm = TRUE)

# --- LOO stability sentence (exact values for SI text) ---
# Identify min/max LOO estimates on the asinh scale and the largest shift
loo_min_idx <- which.min(loo_tbl$mu_asinh)
loo_max_idx <- which.max(loo_tbl$mu_asinh)
loo_min_val <- loo_tbl$mu_asinh[loo_min_idx]
loo_max_val <- loo_tbl$mu_asinh[loo_max_idx]

# Full-sample intercept (asinh scale): mu_asinh
shift_min <- abs(mu_asinh - loo_min_val)
shift_max <- abs(mu_asinh - loo_max_val)

if (shift_min >= shift_max) {
  worst_idx <- loo_min_idx
  worst_val <- loo_min_val
  worst_shift <- shift_min
} else {
  worst_idx <- loo_max_idx
  worst_val <- loo_max_val
  worst_shift <- shift_max
}

# Formatted values for direct paste into text
loo_min_fmt <- sprintf("%.3f", loo_min_val)
loo_max_fmt <- sprintf("%.3f", loo_max_val)
mu_asinh_fmt <- sprintf("%.3f", mu_asinh)
worst_shift_fmt_3 <- sprintf("%.3f", worst_shift)
worst_shift_fmt_4 <- sprintf("%.4f", worst_shift)

stability_sentence <- sprintf(
  "Stability: The pooled mean was stable. The LOO estimates ranged from %s to %s (on the asinh scale). Thus, the most influential substudy caused a shift in the estimated mean of %s (its inclusion shifted the estimate from %s to %s).",
  loo_min_fmt, loo_max_fmt, worst_shift_fmt_4, sprintf("%.3f", worst_val), mu_asinh_fmt
)

cat("\n[LOO Stability] ", stability_sentence, "\n", sep = "")

#ranges of values
range(loo_tbl$mu_beta)


# plot: y shows raw β; ticks show asinh(β)
# ---- Figure D-1A (LOO): label β with units -----------------------------------
pA <- ggplot(loo_tbl, aes(left_out_index, mu_beta)) +
  geom_point(size = 2.1, alpha = 0.5, color = pub_blue) +
  geom_hline(yintercept = mu_beta, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    title = "Figure D-1A. Leave-one-out",
    x = "Left-out substudy index",
    y = expression(paste("Pooled mean ", beta, " (cm"^2, " fish"^{-1}, " day"^{-1}, ")"))
  ) +
  theme_pub_adrian()


print(pA)
ggsave(file.path(fig_dir, "Fig_D1A_LOO.png"), pA,
       width = 7.5, height = 4.0, dpi = 600, bg = "white")

print(pA)
ggsave(file.path(fig_dir, "Fig_D1A_LOO.png"), pA, width = 7.5, height = 4.0, dpi = 600, bg = "white")
# export table
readr::write_csv(loo_tbl, file.path(fig_dir, "Fig_D1A_LOO_table.csv"))




# ================================================================
# D.2 PUBLICATION & SMALL-STUDY BIAS (ALL DATA)
# ================================================================

# --- Egger inputs strictly on the asinh scale (for the test) ---
# Trim to positive standard errors and compute the usual precision/std-effect
# quantities for a simple Egger regression.
egger_df <- dat %>%
  dplyr::transmute(
    precision = 1 / se_asinh,   # 1 / SE(asinh β)
    std_eff   = yi_asinh / se_asinh,  # asinh β / SE(asinh β)
    se        = se_asinh
  ) %>%
  dplyr::filter(is.finite(precision), is.finite(std_eff), se > 0)

stopifnot(nrow(egger_df) > 2)

# Egger regression (no transforms inside the model)
fit_egger <- stats::lm(std_eff ~ precision, data = egger_df)
egger_int <- unname(coef(fit_egger)[1])
egger_p   <- summary(fit_egger)$coefficients[1, 4]

# # -------------------------- Funnel plot --------------------------
# # ----- explicit ranges -----
# x_rng <- range(dat$beta_raw, na.rm = TRUE)     # ~[-478, 9495]
# y_min <- 0
# y_max <- 10000 #max(dat$se_raw, na.rm = TRUE)         # ~9.35e9
# 
# # helper to get asinh-spaced breaks but label with raw numbers
# asinh_breaks <- function(vals, n = 9) sinh(pretty(asinh(vals), n = n))
# 
# x_breaks <- asinh_breaks(x_rng, n = 11)
# y_breaks <- asinh_breaks(c(y_min, y_max), n = 7)
# 
# x_labels <- scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)(x_breaks)
# y_labels <- scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)(y_breaks)
# 
# # ----- funnel bounds on the RAW scale (lines are ±1.96*SE around μ on raw scale) -----
# funnel_bounds_raw <- function(mu_beta, se_seq_raw) {
#   tibble::tibble(
#     se_raw  = se_seq_raw,
#     low_raw = mu_beta - 1.96 * se_raw,
#     high_raw= mu_beta + 1.96 * se_raw
#   )
# }
# fb_raw <- funnel_bounds_raw(mu_beta, seq(y_min, y_max, length.out = 1500))
# 
# # ----- PLOT (asinh spacing on both axes, arithmetic labels) -----
# pB <- ggplot(dat, aes(x = beta_raw, y = se_raw)) +
#   geom_point(size = 2.1, alpha = 0.55, color = pub_blue) +
#   geom_line(data = fb_raw, aes(x = low_raw,  y = se_raw, group = 1),
#             linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
#   geom_line(data = fb_raw, aes(x = high_raw, y = se_raw, group = 1),
#             linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
#   geom_vline(xintercept = mu_beta, linetype = "dashed",
#              color = "grey40", linewidth = 0.6) +
#   scale_x_continuous(
#     trans  = asinh_trans(),     # ← spacing, not labeling
#     breaks = x_breaks,
#     labels = x_labels,
#     expand = expansion(mult = c(0.02, 0.04))
#   ) +
#   scale_y_continuous(
#     trans  = asinh_trans(),     # ← spacing, not labeling
#     breaks = y_breaks,
#     labels = y_labels,
#     expand = expansion(mult = c(0.02, 0.04))
#   ) +
#   coord_cartesian(xlim = x_rng, ylim = c(y_min, y_max)) +  # include ALL data
#   labs(
#     title = "Figure D-1B. Funnel",
#     x = expression(paste("Effect size ", beta, " (cm"^2, " fish"^{-1}, " day"^{-1}, ")")),
#     y = expression(paste("Standard error, SE(", beta, ")"))
#   ) +
#   theme_pub_adrian()
# 
# print(pB)
# 
# 
# 
# #####
# #pB Flippled
# #####
# #####
# # pB Flipped
# #####
# 
# # ----- explicit ranges -----
# x_rng <- range(dat$beta_raw, na.rm = TRUE)
# y_min <- 0
# y_max <- 10000
# 
# # X-axis breaks (keep asinh_breaks for x-axis)
# x_breaks <- asinh_breaks(x_rng, n = 11)
# x_labels <- scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)(x_breaks)
# 
# # Y-axis breaks - MANUAL
# y_breaks <- c(0, 10, 50, 100, 500, 1000, 5000, 10000)
# y_labels <- scales::label_number(scale_cut = scales::cut_si(""), accuracy = 1)(y_breaks)
# 
# # ----- Recalculate funnel bounds with NEW y_max -----
# fb_raw <- funnel_bounds_raw(mu_beta, seq(y_min, y_max, length.out = 1500))
# 
# # ----- PLOT -----
# pB <- ggplot(dat, aes(x = beta_raw, y = se_raw)) +
#   geom_point(size = 2.1, alpha = 0.55, color = pub_blue) +
#   geom_line(data = fb_raw, aes(x = low_raw,  y = se_raw, group = 1),
#             linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
#   geom_line(data = fb_raw, aes(x = high_raw, y = se_raw, group = 1),
#             linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
#   geom_vline(xintercept = mu_beta, linetype = "dashed",
#              color = "grey40", linewidth = 0.6) +
#   scale_x_continuous(
#     trans  = asinh_trans(),
#     breaks = x_breaks,
#     labels = x_labels,
#     expand = expansion(mult = c(0.02, 0.04))
#   ) +
#   scale_y_continuous(
#     trans  = scales::trans_new(
#       "asinh_reverse",
#       transform = function(x) -asinh(x),
#       inverse = function(x) sinh(-x)
#     ),
#     breaks = y_breaks,
#     labels = y_labels,
#     limits = c(y_max, y_min),
#     expand = expansion(mult = c(0.01, 0.01)),  # <--- Small expansion instead of c(0,0)
#     oob = scales::squish
#   ) +
#   coord_cartesian(xlim = x_rng, clip = "on") +  # <--- Added clip = "on"
#   labs(
#     title = "Figure D-1B. Funnel",
#     x = expression(paste("Effect size ", beta, " (cm"^2, " fish"^{-1}, " day"^{-1}, ")")),
#     y = expression(paste("Standard error, SE(", beta, ")"))
#   ) +
#   theme_pub_adrian()
# 
# pB




#####
# pB-asinh: Funnel plot in transformed (asinh) space
# This shows the data as it's analyzed (straight funnel lines)
#####


# ----- Ranges in asinh space - USE ACTUAL DATA RANGE -----
x_rng_asinh <- range(dat$yi_asinh, na.rm = TRUE)
y_min_asinh <- 0
y_max_asinh <- max(dat$se_asinh, na.rm = TRUE)  

# Add a tiny bit of padding
y_max_asinh <- y_max_asinh * 1.05  # 5% padding

# Create breaks in asinh space
x_breaks_asinh <- pretty(x_rng_asinh, n = 11)
y_breaks_asinh <- pretty(c(y_min_asinh, y_max_asinh), n = 9)

# Labels show asinh values with nice formatting
x_labels_asinh <- formatC(x_breaks_asinh, format = "f", digits = 2)
y_labels_asinh <- formatC(y_breaks_asinh, format = "f", digits = 2)

# ----- Funnel bounds in ASINH space (these will be STRAIGHT lines) -----
funnel_bounds_asinh <- function(mu_asinh, se_seq_asinh) {
  tibble::tibble(
    se_asinh  = se_seq_asinh,
    low_asinh = mu_asinh - 1.96 * se_seq_asinh,
    high_asinh= mu_asinh + 1.96 * se_seq_asinh
  )
}

# Generate funnel bounds
fb_asinh <- funnel_bounds_asinh(mu_asinh, seq(y_min_asinh, y_max_asinh, length.out = 1500))

# ----- PLOT in pure asinh space - ZOOMED -----
pB_asinh <- ggplot(dat, aes(x = yi_asinh, y = se_asinh)) +
  geom_point(size = 2.1, alpha = 0.55, color = pub_blue) +
  # Funnel bounds (perfectly straight!)
  geom_line(data = fb_asinh, aes(x = low_asinh,  y = se_asinh, group = 1),
            linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
  geom_line(data = fb_asinh, aes(x = high_asinh, y = se_asinh, group = 1),
            linewidth = 0.5, color = "grey50", inherit.aes = FALSE) +
  # Pooled mean line
  geom_vline(xintercept = mu_asinh, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  scale_x_continuous(
    breaks = x_breaks_asinh,
    labels = x_labels_asinh,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_reverse(  # Reverse so small SE is at top
    breaks = y_breaks_asinh,
    labels = y_labels_asinh,
    limits = c(y_max_asinh, y_min_asinh),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(xlim = x_rng_asinh, clip = "on") +
  labs(
    title = "Figure D-1B (asinh). Funnel plot in transformed space",
    x = expression(paste("Effect size asinh(", beta, ")")),
    y = expression(paste("Standard error, SE(asinh(", beta, "))"))
  ) +
  theme_pub_adrian()

print(pB_asinh)

# Save asinh funnel plot
ggsave(file.path(fig_dir, "Fig_D1B_Funnel_asinh.png"), pB_asinh,
       width = 7.5, height = 4.0, dpi = 600, bg = "white")

# Quick comparison plot (if pB exists from earlier code)
if (exists("pB")) {
  comparison <- pB / pB_asinh + 
  plot_annotation(
    title = "Funnel plot comparison",
    subtitle = "Top: Raw scale with asinh transforms (curved lines)\nBottom: Pure asinh space (straight lines)"
  )

  print(comparison)

  ggsave(file.path(fig_dir, "Fig_D1B_Comparison.png"), comparison,
         width = 7.5, height = 8, dpi = 600, bg = "white")
}


# ---- Figure D-1C: Egger

#subset values so that it's clean 


# 5 precision points are >1500
# 5 std eff are >11000
# ---- Limits from your note ----
xlim_vals <- c(0, 1500)
ylim_vals <- c(-1e2, 11000)

# ---- Plot-window counts (for the scatter) ----
n_total <- nrow(dat)
n_dropped_view <- sum(
  dat$precision < xlim_vals[1] |
    dat$precision > xlim_vals[2] |
    dat$std_eff  < ylim_vals[1] |
    dat$std_eff  > ylim_vals[2],
  na.rm = TRUE
)
n_kept_view <- n_total - n_dropped_view

cat(sprintf("Plot limits: x=[%.0f, %.0f], y=[%.0f, %.0f]\n",
            xlim_vals[1], xlim_vals[2], ylim_vals[1], ylim_vals[2]))
cat(sprintf("Scatter keeps %d / %d points (drops %d by view limits)\n",
            n_kept_view, n_total, n_dropped_view))

# ---- Regression subset: drop the 10 extremes you specified ----
reg_df <- dat %>%
  dplyr::filter(
    is.finite(precision), is.finite(std_eff),
    precision <= 1500,          # drop 5 precision outliers
    std_eff   <= 11000          # drop 5 std-eff outliers
  )

n_reg <- nrow(reg_df)
n_drop_reg <- n_total - n_reg
cat(sprintf("Regression uses %d / %d points (dropped %d: precision>1500 or std_eff>11000)\n",
            n_reg, n_total, n_drop_reg))

# ---- Fit regression on the subset ----
fit <- lm(std_eff ~ precision, data = reg_df)
co <- coef(fit)
egger_int <- unname(co[1]); egger_slope <- unname(co[2])
egger_p <- summary(fit)$coefficients[2,4]

# ---- Prediction grid over the plotted x-range ----
pred_df <- data.frame(precision = seq(xlim_vals[1], xlim_vals[2], length.out = 400))
pred_df$std_eff_hat <- predict(fit, newdata = pred_df)

# ---- Plot: scatter (all data, clipped by coord) + regression line (subset fit) ----
pC <- ggplot(dat, aes(x = precision, y = std_eff)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_line(data = pred_df, aes(y = std_eff_hat), color = "#D97706", linewidth = 1.2) +
  coord_cartesian(xlim = xlim_vals, ylim = ylim_vals) +
  labs(
    title = "Figure D-1C. Egger",
    x = "Precision (1 / SE of asinh β)",
    y = "Standardized effect (asinh β / SE)"
  )+
  theme_pub_adrian()



# ==========================================
# D.3 NULL-ADD SENSITIVITY (ALL data baseline)
# ==========================================
null_add_curve <- function(N_seq = 0:10, se_ref = quantile(dat$se, 0.95, na.rm = TRUE)) {
  stopifnot(is.numeric(se_ref), is.finite(se_ref), se_ref > 0)
  vi_ref <- se_ref^2
  
  make_aug <- function(n_add) {
    if (n_add == 0) return(dat)
    add_df <- tibble(
      study_num    = paste0("null_study_",    seq_len(n_add)),
      substudy_num = paste0("null_substudy_", seq_len(n_add)),
      yi = 0, vi = vi_ref, se = se_ref
    ) %>% mutate(beta_raw = sinh(yi))
    bind_rows(dat, add_df)
  }
  
  map_dfr(N_seq, function(n_add) {
    dat_aug <- make_aug(n_add)
    fit <- rma.mv(yi, vi, random = ~ 1 | study_num/substudy_num,
                  method = "REML", test = "t", data = dat_aug)
    tibble(n_add = n_add,
           mu_asinh = as.numeric(fit$b),
           mu_beta  = inv_asinh(as.numeric(fit$b)))
  })
}

null_path <- null_add_curve(0:10)
delta_5 <- null_path %>% filter(n_add == 5) %>% pull(mu_asinh) - mu_asinh
cat(sprintf("[D.3 Null-add] Δμ(asinh) at N=5 = %.4f\n", delta_5))

range(null_path$mu_beta)

# Plot: y shows raw β mean; ticks show asinh
pD3 <- ggplot(null_path, aes(n_add, mu_beta)) +
  geom_line(linewidth = 0.7, color = "grey35") +
  geom_point(size = 2.1, alpha = 0.5, color = pub_blue) +
  geom_hline(yintercept = mu_beta, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  labs(
    title = "Figure D-1D. Null-add sensitivity",
    x = "Number of added null substudies",
    y = expression(paste("Pooled mean ", beta, " (cm"^2, " fish"^{-1}, " day"^{-1}, ")"))
  ) +
  theme_pub_adrian()

print(pD3)

ggsave(file.path(fig_dir, "Fig_D1D_Null.png"), pD3, width = 7.5, height = 4.0, dpi = 600, bg = "white")

# ==========================================
# Composite (everything in one figure)
# ==========================================
fig_all <- (pA + pB) / (pC + pD3) + patchwork::plot_layout(guides = "collect")
print(fig_all)
ggsave(file.path(fig_dir, "Fig_D1_ALL.png"), fig_all,
       width = 12, height = 10, dpi = 600, bg = "white")

# ==========================================
# Summary table (console + CSV)
# ==========================================
safe_num <- function(x, digits = 3) ifelse(is.na(x), NA, signif(as.numeric(x), digits))

summary_tbl <- tibble::tibble(
  section   = c(
    "Model (REML)", "Model (REML)",
    "Heterogeneity", "Heterogeneity", "Heterogeneity",
    "LOO", "LOO", "LOO", "LOO",
    "Egger", "Egger",
    "Begg", "Begg",
    "Null-add"
  ),
  metric    = c(
    "μ (asinh)", "μ (β, raw)",
    "Q_E", "I² (%)", "σ² (study/substudy)",
    "μ_min (asinh)", "μ_max (asinh)", "β_min (raw)", "β_max (raw)",
    "Intercept", "p-value",
    "Kendall τ", "p-value",
    "Δμ(asinh) at N=5"
  ),
  estimate  = c(
    safe_num(mu_asinh), safe_num(mu_beta),
    safe_num(QE), safe_num(I2), sprintf("%.4f / %.4f", sig2["study"], sig2["substudy"]),
    safe_num(loo_rng[1]), safe_num(loo_rng[2]), safe_num(sinh(loo_rng[1])), safe_num(sinh(loo_rng[2])),
    safe_num(egger_int), safe_num(egger_p),
    safe_num(unname(begg$estimate)), safe_num(unname(begg$p.value)),
    safe_num(delta_5)
  ),
  ci_lower  = c(
    safe_num(mu_asinh_ci[1]), safe_num(mu_beta_ci[1]),
    NA, NA, NA,
    NA, NA, NA, NA,
    NA, NA,
    NA, NA,
    NA
  ),
  ci_upper  = c(
    safe_num(mu_asinh_ci[2]), safe_num(mu_beta_ci[2]),
    NA, NA, NA,
    NA, NA, NA, NA,
    NA, NA,
    NA, NA,
    NA
  ),
  notes     = c(
    sprintf("k=%d; random = ~1|study/substudy; method=REML", k), "Back-transform via sinh()",
    sprintf("p-value < %.2g", QEp), "Pseudo-I² from Q_E", "Variance components on asinh scale",
    "From leave-one-out fits", "From leave-one-out fits", "Back-transformed via sinh()", "Back-transformed via sinh()",
    "Egger on asinh-scale variables", "", "Begg rank correlation on asinh scale", "", 
    "Null effects at SE = 95th percentile"
  )
)

# Print a clean console header
cat("\n================ Appendix D: Summary (key statistics) ================\n")
cat(sprintf("REML μ (asinh) = %.3f  CI [%.3f, %.3f]\n", mu_asinh, mu_asinh_ci[1], mu_asinh_ci[2]))
cat(sprintf("REML μ (β raw) = %.2f  CI [%.2f, %.2f]\n", mu_beta,  mu_beta_ci[1],  mu_beta_ci[2]))
cat(sprintf("Heterogeneity: Q_E = %.0f (p < %.2g); I² ≈ %.2f%%; σ²(study/substudy) = %.4f / %.4f\n",
            QE, QEp, I2, sig2['study'], sig2['substudy']))
cat(sprintf("LOO (asinh μ): min = %.3f, max = %.3f  |  LOO (β raw): min = %.2f, max = %.2f\n",
            loo_rng[1], loo_rng[2], sinh(loo_rng[1]), sinh(loo_rng[2])))
cat(sprintf("Egger: intercept = %.3f, p = %.3g  |  Begg: τ = %.3f, p = %.3g\n",
            egger_int, egger_p, unname(begg$estimate), unname(begg$p.value)))
cat(sprintf("Null-add: Δμ(asinh) at N=5 = %.4f\n", delta_5))
cat("=====================================================================\n\n")

# Pretty console table if knitr is available
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(summary_tbl, align = "llrrr", digits = 3,
                     col.names = c("Section", "Metric", "Estimate", "CI Lower", "CI Upper", "Notes")))
} else {
  print(summary_tbl)
}

# Save CSV for supplement
out_csv <- file.path(fig_dir, "AppendixD_summary_table.csv")
readr::write_csv(summary_tbl, out_csv)
message(sprintf("[saved] %s", out_csv))
