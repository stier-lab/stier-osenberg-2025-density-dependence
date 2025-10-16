###############################################################################
# Script: 4_predators_pairedpredators.R
# Goal:   Quantify how predator presence modifies β using phylogenetically aware
#         REML models, interrogate paired experimental designs, and generate the
#         predator-related figures/tables for the manuscript.
# Inputs: all_dat2/pruned_tree from 1_data_phylogeny_loading.R
# Outputs: results/predator_* CSV/PNG files, figures/predators_vs_beta_boxplot.png,
#          figures/Fig4_paired_vs_unpaired.png, console summaries.
###############################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
# Predation analyses reuse the shared `all_dat2` table and pruned tree.
# Sourcing the loaders defensively keeps the script usable on its own.
source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE) ||
    !exists("pruned_tree", inherits = FALSE)) {
  message("Core data objects missing; sourcing 1_data_phylogeny_loading.R.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

# Create output directories if needed
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)


# ──────────────────────────────────────────────────────────────────────────────
# 1. Predator‐presence summary
# ──────────────────────────────────────────────────────────────────────────────
# Quick checks to understand the balance of predator-present vs predator-absent
# substudies.  These counts inform later interpretation.
predator_counts <- all_dat2 %>%
  count(predators) %>%
  mutate(proportion = n / sum(n))

print(predator_counts)

cat("Fraction with predators present: ",
    round(predator_counts$proportion[predator_counts$predators == "present"], 3),
    "\n")
cat("Fraction with predators absent:  ",
    round(predator_counts$proportion[predator_counts$predators == "absent"], 3),
    "\n")


# ──────────────────────────────────────────────────────────────────────────────
# 2. Build & trim phylogenetic VCV to match data
# ──────────────────────────────────────────────────────────────────────────────
# We prepare a correlation matrix for the species random effect so that each
# model can include phylogenetic structure when requested.
phylo_vcv_full <- ape::vcv(pruned_tree, corr = TRUE)

good_sp    <- intersect(all_dat2$g_sp, rownames(phylo_vcv_full))
model_data <- all_dat2 %>%
  filter(
    g_sp %in% good_sp,
    !is.na(betanls2_asinh),
    !is.na(betanlsvar_asinh)
  )
phylo_vcv <- phylo_vcv_full[good_sp, good_sp]

missing_sp <- setdiff(unique(all_dat2$g_sp), good_sp)
if (length(missing_sp) > 0) {
  warning("Dropped species not in tree: ",
          paste(missing_sp, collapse = ", "))
}


# ──────────────────────────────────────────────────────────────────────────────
# 3. Mixed‐effects meta‐analysis (null vs predators)
# ──────────────────────────────────────────────────────────────────────────────
# Baseline REML model (no predator effect) vs model with predator factor.
# Both include study/substudy nesting and a phylogenetic species effect.
m_null <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ 1,
  random = list(~1 | study_num/substudy_num,
                ~1 | g_sp),
  R      = list(g_sp = phylo_vcv),
  data   = model_data,
  method = "REML"
)

m_pred <- update(m_null, mods = ~ predators, test = "t")
print(AIC(m_null, m_pred))


# ──────────────────────────────────────────────────────────────────────────────
# 4. Extract & back‐transform predator‐model coefficients
# ──────────────────────────────────────────────────────────────────────────────
# Everything is estimated on the asinh scale, so we back-transform for reporting.
cfs        <- coef(summary(m_pred))
β_abs_log  <- cfs["intrcpt",         "estimate"]
β_pre_log  <- cfs["predatorspresent","estimate"]
ci_abs_log <- c(m_pred$ci.lb[1], m_pred$ci.ub[1])
ci_pre_log <- c(m_pred$ci.lb[2], m_pred$ci.ub[2])

β_abs  <- sinh(β_abs_log)
β_pre  <- sinh(β_abs_log + β_pre_log)
CI_abs <- sinh(ci_abs_log)
CI_pre <- sinh(β_abs_log + ci_pre_log)
pct_ch <- 100 * (β_pre - β_abs) / abs(β_abs)

cat(sprintf(
  "Meta β no predators:    %.3f [%.3f, %.3f]\n",
  β_abs, CI_abs[1], CI_abs[2]
))
cat(sprintf(
  "Meta β with predators:  %.3f [%.3f, %.3f]\n",
  β_pre, CI_pre[1], CI_pre[2]
))
cat(sprintf("Percent change:         %.1f%%\n", pct_ch))


# ──────────────────────────────────────────────────────────────────────────────
# 5. Summary table (gt)
# ──────────────────────────────────────────────────────────────────────────────
# A concise table for the supplement capturing counts and effect sizes.
res_tbl <- tibble::tribble(
  ~Metric,                  ~Estimate,                                                       ~`95% CI`,
  "Studies w/ predators",   sprintf("%d", predator_counts$n[predator_counts$predators == "present"]),    "",
  "Studies w/o predators",  sprintf("%d", predator_counts$n[predator_counts$predators == "absent"]),     "",
  "Fraction present",       sprintf("%.3f", predator_counts$proportion[predator_counts$predators == "present"]), "",
  "Fraction absent",        sprintf("%.3f", predator_counts$proportion[predator_counts$predators == "absent"]),   "",
  "β (no predators)",       sprintf("%.3f", β_abs),                                            sprintf("[%.3f, %.3f]", CI_abs[1], CI_abs[2]),
  "β (with predators)",     sprintf("%.3f", β_pre),                                            sprintf("[%.3f, %.3f]", CI_pre[1], CI_pre[2]),
  "Pct change",             sprintf("%.1f%%", pct_ch),                                         ""
)

res_tbl %>%
  gt::gt() %>%
  gt::tab_header(
    title    = md("**Predator-Presence Meta-Analysis**"),
    subtitle = "Effect sizes & study counts"
  ) %>%
  gt::cols_label(
    Metric   = "Metric",
    Estimate = "Estimate",
    `95% CI` = "95% Confidence Interval"
  ) %>%
  gt::fmt_missing(columns = everything(), missing_text = "") %>%
  gt::tab_options(
    table.font.size           = px(14),
    heading.align             = "center",
    column_labels.font.weight = "bold"
  ) %>%
  # PNG save requires Chrome/Chromium - commented out for compatibility
  # gt::gtsave(here::here("results", "predator_presence_summary.png"))
  gt::gtsave(here::here("results", "predator_presence_summary.html"))


# ──────────────────────────────────────────────────────────────────────────────
# 6. All‐vs‐paired studies comparison
# ──────────────────────────────────────────────────────────────────────────────
# Compare results across all studies and the subset with explicit predator pairs.
summarize_pred_effect <- function(mod) {
  cf    <- coef(summary(mod))
  i_log <- cf["intrcpt",         "estimate"]
  m_log <- cf["predatorspresent","estimate"]
  b0    <- sinh(i_log)
  b1    <- sinh(i_log + m_log)
  ci0   <- sinh(c(mod$ci.lb[1], mod$ci.ub[1]))
  ci1   <- sinh(i_log + c(mod$ci.lb[2], mod$ci.ub[2]))
  
  tibble(
    set      = NA_character_,
    beta0    = b0,
    ci0_lo   = ci0[1],
    ci0_hi   = ci0[2],
    beta1    = b1,
    ci1_lo   = ci1[1],
    ci1_hi   = ci1[2],
    pct_mean = 100 * (b1 - b0) / abs(b0)
  )
}

m_all      <- m_pred
paired_df  <- model_data %>% filter(paired_pred == "paired")
m_paired   <- update(m_all, data = paired_df)

sum_all    <- summarize_pred_effect(m_all)    %>% mutate(set = "All studies")
sum_paired <- summarize_pred_effect(m_paired) %>% mutate(set = "Paired only")

comp_tbl <- bind_rows(sum_all, sum_paired) %>%
  transmute(
    `Study set`   = set,
    `β no pred`   = sprintf("%.2f [%.2f, %.2f]", beta0, ci0_lo, ci0_hi),
    `β with pred` = sprintf("%.2f [%.2f, %.2f]", beta1, ci1_lo, ci1_hi),
    `% inc (mean)` = sprintf("%.1f%%", pct_mean)
  )

comp_tbl %>%
  gt::gt() %>%
  gt::tab_header(
    title    = md("**Predator‐Presence: All vs Paired**"),
    subtitle = "Meta‐analytic β estimates"
  ) %>%
  gt::cols_label(
    `Study set`   = "Study set",
    `β no pred`   = md("**β** (no predators)"),
    `β with pred` = md("**β** (with predators)")
  ) %>%
  gt::tab_options(heading.align = "center", table.font.size = px(14)) %>%
  # PNG save requires Chrome/Chromium - commented out for compatibility
  # gt::gtsave(here::here("results", "predator_all_vs_paired.png"))
  gt::gtsave(here::here("results", "predator_all_vs_paired.html"))



# ──────────────────────────────────────────────────────────────────────────────
# 6b. Paired experiments: within-pair meta-analytic test of predators
# ──────────────────────────────────────────────────────────────────────────────
# Treat paired manipulations as blocking factors so we can isolate within-pair contrast.

# Defensive checks
if (!all(c("paired_pred","paired_substudy_num","predators") %in% names(model_data))) {
  stop("Expected columns 'paired_pred', 'paired_substudy_num', and 'predators' are missing in model_data.")
}

# Keep only experimental, paired rows; standardize predator labels
paired_df <- model_data %>%
  dplyr::filter(expt_obs == "Exp", paired_pred == "paired", !is.na(paired_substudy_num)) %>%
  dplyr::mutate(
    predators = tolower(trimws(predators)),
    predators = dplyr::case_when(
      predators %in% c("present","predator_present","with","yes") ~ "present",
      predators %in% c("absent","control","without","no")         ~ "absent",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(predators)) %>%
  dplyr::mutate(predators = factor(predators, levels = c("absent","present"))) %>%
  droplevels()

# Keep only those pairs that truly have BOTH predator levels
paired_ids <- paired_df %>%
  dplyr::distinct(paired_substudy_num, predators) %>%
  dplyr::count(paired_substudy_num, name = "n_levels") %>%
  dplyr::filter(n_levels == 2) %>%
  dplyr::pull(paired_substudy_num)

paired_df <- paired_df %>% dplyr::filter(paired_substudy_num %in% paired_ids)

cat("# paired groups with BOTH levels: ", dplyr::n_distinct(paired_df$paired_substudy_num), "\n")

# --- Fit within-pair meta-regression (predators effect is a within-pair contrast)
m_paired_only <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ predators,                          # test predatorspresent
  random = list(
    ~ 1 | paired_substudy_num,                   # block by pair
    ~ 1 | study_num/substudy_num,                # retain study/substudy structure
    ~ 1 | g_sp                                   # species (can be phylo-correlated)
  ),
  R      = list(g_sp = phylo_vcv[levels(factor(paired_df$g_sp)), 
                                 levels(factor(paired_df$g_sp))]),
  data   = paired_df,
  method = "REML",
  test   = "t"
)

sm <- coef(summary(m_paired_only))
est_log <- sm["predatorspresent","estimate"]
se_log  <- sm["predatorspresent","se"]
t_val   <- sm["predatorspresent","tval"]
df_val  <- sm["predatorspresent","df"]
p_val   <- sm["predatorspresent","pval"]

# Back-transform predicted means for absent/present (as in Sec. 4)
new_pred_6b <- tibble::tibble(
  predators = factor(c("absent","present"), levels = c("absent","present"))
)
X_new_6b <- model.matrix(~ predators, data = new_pred_6b)
preds_6b <- predict(m_paired_only, newmods = X_new_6b[, "predatorspresent", drop = FALSE])

bt <- function(x) sinh(x)
paired_abs_hat  <- bt(preds_6b$pred[1])
paired_abs_lo   <- bt(preds_6b$ci.lb[1])
paired_abs_hi   <- bt(preds_6b$ci.ub[1])
paired_pres_hat <- bt(preds_6b$pred[2])
paired_pres_lo  <- bt(preds_6b$ci.lb[2])
paired_pres_hi  <- bt(preds_6b$ci.ub[2])
pct_change_6b   <- 100 * (paired_pres_hat - paired_abs_hat) / abs(paired_abs_hat)

cat(sprintf(
  "Paired meta-analytic test (present vs absent): estimate=%.3f (asinh), SE=%.3f, t=%.2f, df=%.1f, p=%.4g\n",
  est_log, se_log, t_val, df_val, p_val
))
cat(sprintf(
  "Paired β (absent):  %.3f [%.3f, %.3f]\n", paired_abs_hat,  paired_abs_lo,  paired_abs_hi
))
cat(sprintf(
  "Paired β (present): %.3f [%.3f, %.3f]\n", paired_pres_hat, paired_pres_lo, paired_pres_hi
))
cat(sprintf("Percent change (present vs absent): %.1f%%\n", pct_change_6b))

# --- Optional: ML LRT within pairs (model comparison)
m0_pair_ML <- update(m_paired_only, mods = ~ 1, method = "ML")
m1_pair_ML <- update(m_paired_only,              method = "ML")
LRT_6b  <- 2 * (as.numeric(logLik(m1_pair_ML)) - as.numeric(logLik(m0_pair_ML)))
p_LR_6b <- pchisq(LRT_6b, df = 1, lower.tail = FALSE)
cat(sprintf("Paired ML LRT: LRT = %.3f, p = %.4g\n", LRT_6b, p_LR_6b))

# --- Save a tidy results row
res_6b <- tibble::tibble(
  n_pairs                 = dplyr::n_distinct(paired_df$paired_substudy_num),
  predators_coef_asinh    = est_log,
  predators_se_asinh      = se_log,
  predators_t             = t_val,
  predators_df            = df_val,
  predators_p             = p_val,
  beta_absent_hat         = paired_abs_hat,
  beta_absent_ci_lo       = paired_abs_lo,
  beta_absent_ci_hi       = paired_abs_hi,
  beta_present_hat        = paired_pres_hat,
  beta_present_ci_lo      = paired_pres_lo,
  beta_present_ci_hi      = paired_pres_hi,
  percent_change_present  = pct_change_6b,
  LRT_within_pairs        = LRT_6b,
  LRT_p                   = p_LR_6b
)

readr::write_csv(res_6b, here::here("results", "predator_paired_withinpair_test.csv"))

# ──────────────────────────────────────────────────────────────────────────────
# 7. Boxplot + meta‐predictions
# ──────────────────────────────────────────────────────────────────────────────
# Visualise raw data alongside model-predicted means and confidence intervals.
new_pred <- tibble(
  predators = factor(c("absent", "present"), levels = c("absent", "present"))
)
X_new <- model.matrix(~ predators, data = new_pred)
preds <- predict(m_pred, newmods = X_new[, "predatorspresent", drop = FALSE])

new_pred <- new_pred %>%
  mutate(
    beta_hat = sinh(preds$pred),
    ci_lo    = sinh(preds$ci.lb),
    ci_hi    = sinh(preds$ci.ub)
  )

pred_cols <- c(absent = "#1E90FF", present = "#FF4500")

p1 <- ggplot(model_data, aes(predators, betanls2_raw_cm)) +
  geom_boxplot(aes(fill = predators), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = predators), width = 0.2, size = 2) +
  geom_point(
    data = new_pred,
    aes(predators, beta_hat),
    shape = 21, size = 5, color = "black", stroke = 1
  ) +
  geom_errorbar(
    data = new_pred,
    aes(predators, beta_hat, ymin = ci_lo, ymax = ci_hi),
    width = 0.1, size = 1
  ) +
  scale_fill_manual(values = pred_cols) +
  scale_color_manual(values = pred_cols) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(x = "Predator Presence", y = "β (density‐dependence)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

ggsave(
  here::here("figures", "predators_vs_beta_boxplot.png"),
  plot = p1, width = 5, height = 5, dpi = 300
)


# ──────────────────────────────────────────────────────────────────────
# 8. Paired‐vs‐unpaired study plot (revised styling)
# ──────────────────────────────────────────────────────────────────────
# Overlay paired lines with unpaired jittered points to illustrate heterogeneity.
# ---- Fix paired vs unpaired data (do NOT use distinct() here) ----
# 1) Standardize predator labels on the modeled data
md_pred <- model_data %>%
  mutate(
    predators = tolower(trimws(predators)),
    predators = case_when(
      predators %in% c("present","predator_present","with","yes") ~ "present",
      predators %in% c("absent","control","without","no")         ~ "absent",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(predators)) %>%
  mutate(predators = factor(predators, levels = c("absent","present"))) %>%
  droplevels()

# 2) Keep only paired_substudy_num that truly have BOTH levels
paired_ids <- md_pred %>%
  filter(!is.na(paired_substudy_num)) %>%
  distinct(paired_substudy_num, predators) %>%
  count(paired_substudy_num, name = "n_levels") %>%
  filter(n_levels == 2) %>%
  pull(paired_substudy_num)

# 3) Build full data frames for paired and unpaired
paired_full <- md_pred %>%
  filter(!is.na(paired_substudy_num), paired_substudy_num %in% paired_ids)

unpaired_studies <- md_pred %>%
  filter(is.na(paired_substudy_num) | !(paired_substudy_num %in% paired_ids))

cat("Paired IDs with both levels: ", length(paired_ids), "\n")
cat("Unpaired rows: ", nrow(unpaired_studies), "\n")

# 4) Collapse each pair to ONE value per condition (median on raw scale)
paired_pts <- paired_full %>%
  mutate(y_raw = sinh(betanls2_asinh)) %>%
  group_by(paired_substudy_num, predators) %>%
  summarise(y = median(y_raw, na.rm = TRUE), .groups = "drop")

paired_segments <- paired_pts %>%
  pivot_wider(names_from = predators, values_from = y) %>%
  filter(is.finite(absent), is.finite(present))

cat("Complete Absent–Present pairs for plotting: ", nrow(paired_segments), "\n")

# 5) Colors for paired groups: Brewer "Paired" without yellow
raw_paired_pal <- brewer.pal(12, "Paired")
base_qual      <- raw_paired_pal[ raw_paired_pal != "#FFFF99" ]

pair_ids <- sort(unique(paired_pts$paired_substudy_num))
if (length(pair_ids) > length(base_qual)) {
  base_qual <- rep(base_qual, length.out = length(pair_ids))
}
paired_cols <- setNames(base_qual[seq_along(pair_ids)], pair_ids)
paired_fill <- setNames(alpha(paired_cols, 0.7), names(paired_cols))

# 6) Plot: unpaired gray points; paired colored segments + endpoints; meta preds
p2 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
  
  geom_jitter(
    data   = unpaired_studies,
    aes(x = predators, y = sinh(betanls2_asinh)),
    shape  = 21, width = 0.1, size = 2.5,
    fill   = "gray85", color = "gray50", stroke = 0.6, alpha = 0.8
  ) +
  
  geom_segment(
    data = paired_segments,
    aes(x = "absent", xend = "present",
        y = absent,   yend = present,
        color = factor(paired_substudy_num)),
    linewidth = 1
  ) +
  geom_point(
    data  = paired_pts %>% filter(predators == "absent"),
    aes(x = predators, y = y, fill = factor(paired_substudy_num)),
    shape = 21, size = 2, color = "black", stroke = 0.8
  ) +
  geom_point(
    data  = paired_pts %>% filter(predators == "present"),
    aes(x = predators, y = y, fill = factor(paired_substudy_num)),
    shape = 21, size = 2, color = "black", stroke = 0.8
  ) +
  
  geom_errorbar(
    data = new_pred,
    aes(x = predators, y = beta_hat, ymin = ci_lo, ymax = ci_hi),
    width = 0.15, linewidth = 1, color = "black", lineend = "round"
  ) +
  geom_point(
    data = new_pred,
    aes(x = predators, y = beta_hat),
    shape = 23, size = 4, fill = "#1E3A5F", color = "black", stroke = 1.2
  ) +
  
  scale_color_manual(values = paired_cols, guide = "none") +
  scale_fill_manual(values  = paired_fill,  guide = "none") +
  scale_x_discrete(labels = c(absent = "Absent", present = "Present")) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  ) +
  labs(
    x = "Predator Presence",
    y = expression(
      paste("Strength of density-dependent mortality, ",
            beta, " (", cm^2, ~ fish^-1, ~ day^-1, ")")
    )
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

print(p2)

ggsave(
  here::here("figures", "Fig4_paired_vs_unpaired.png"),
  plot = p2, width = 6, height = 6, dpi = 300, bg = "white"
)
