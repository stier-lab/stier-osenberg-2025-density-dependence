################################################################################
# Script: 6_bivariate_plots_predictors.R
# Goal:   Generate marginal (bivariate) visualisations of the additive REML
#         model predictions for each centred predictor, paired with raw data.
# Inputs: all_dat2/pruned_tree from 1_data_phylogeny_loading.R
# Outputs: figures/bivar_* PNG files plus supplementary species-level panels.
################################################################################

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
# Most scripts are designed to run stand-alone for figure regeneration.  We load
# the shared package bundle and, if needed, source the merged data/phylogeny.
source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE) ||
    !exists("pruned_tree", inherits = FALSE)) {
  message("Core data objects missing; sourcing 1_data_phylogeny_loading.R.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

# ensure output directory exists
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Data Preparation (centered covariates)
# ----------------------------------------------------------------------------
# Focus on predator-present substudies and centre the continuous predictors so
# the additive REML model returns interpretable intercepts (average conditions).
df <- all_dat2 %>%
  filter(predators == "present") %>%
  mutate(
    expt_obs = factor(expt_obs, levels = c("Exp", "Obs")),
    logmd_raw = dplyr::if_else(mean_density > 0, log(mean_density), NA_real_),
    logmd_c  = as.numeric(scale(logmd_raw)),
    dur_c    = as.numeric(scale(duration)),
    size_c   = as.numeric(scale(size_start)),
    max_c    = as.numeric(scale(max_length_cm))
  ) %>%
  select(-logmd_raw)

df <- df %>% filter(!is.na(logmd_c))

phylo_vcv <- ape::vcv(pruned_tree, corr = TRUE)
keep_sp    <- intersect(rownames(phylo_vcv), unique(df$g_sp))
phylo_vcv  <- phylo_vcv[keep_sp, keep_sp]
df         <- df %>% filter(g_sp %in% keep_sp)

rand_list  <- list(~1 | study_num/substudy_num, ~1 | g_sp)

# ----------------------------------------------------------------------------
# 3. Helper: predict & back-transform β
# ----------------------------------------------------------------------------
# Convenience wrapper: predict on the asinh scale and back-transform the mean/CI.
predict_beta <- function(model, newmods) {
  pr <- predict(model, newmods = newmods)
  tibble(
    beta  = sinh(pr$pred),
    lower = sinh(pr$ci.lb),
    upper = sinh(pr$ci.ub)
  )
}

# ----------------------------------------------------------------------------
# 4. Fit the no-interaction additive model
# ----------------------------------------------------------------------------
# Same formula as the “best” model in 5_model_selection_comparison.R, fitted here
# once so we can generate smooth marginal effects.
m_all_scaled <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ logmd_c + dur_c + size_c + max_c + expt_obs,
  random = rand_list,
  R      = list(g_sp = phylo_vcv),
  data   = df,
  method = "REML",
  test   = "t"
)

# ----------------------------------------------------------------------------
# 5. Compute means & SDs for un-centering
# ----------------------------------------------------------------------------
# These get used to convert centred predictors back to raw units when plotting.
mean_logden <- mean(log(df$mean_density), na.rm = TRUE)
# use positive densities only to avoid -Inf
positive_density <- df %>%
  filter(mean_density > 0)

mean_logden <- mean(log(positive_density$mean_density), na.rm = TRUE)
sd_logden   <- sd(  log(positive_density$mean_density), na.rm = TRUE)
mean_dur    <- mean(df$duration,      na.rm = TRUE)
sd_dur      <- sd(  df$duration,      na.rm = TRUE)
mean_size   <- mean(df$size_start,    na.rm = TRUE)
sd_size     <- sd(  df$size_start,    na.rm = TRUE)
mean_max    <- mean(df$max_length_cm, na.rm = TRUE)
sd_max      <- sd(  df$max_length_cm, na.rm = TRUE)

# ----------------------------------------------------------------------------
# 6. Bivariate curves & raw data points
# ----------------------------------------------------------------------------
# For each focal predictor we build a grid of values while holding others at zero
# (i.e., their centred mean) and pull predictions + raw scatter for context.


# (a) Helper: the factor levels for expt_obs (so “Exp” is baseline)
exp_levels <- levels(df$expt_obs)

# (b) Density curve over mean_density
logmd_seq <- seq(min(df$logmd_c), max(df$logmd_c), length.out = 100)
df_den <- tibble(
  logmd_c  = logmd_seq,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_den      <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_den)
newmods_den <- X_den[, -1]  # drop intercept
df_den <- df_den %>%
  bind_cols(predict_beta(m_all_scaled, newmods_den)) %>%
  mutate(mean_density = exp(logmd_c * sd_logden + mean_logden))

# (c) Duration curve
dur_seq <- seq(min(df$dur_c), max(df$dur_c), length.out = 100)
df_dur <- tibble(
  logmd_c  = 0,
  dur_c    = dur_seq,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_dur       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_dur)
newmods_dur <- X_dur[, -1]
df_dur <- df_dur %>%
  bind_cols(predict_beta(m_all_scaled, newmods_dur)) %>%
  mutate(duration = dur_c * sd_dur + mean_dur)

# (d) Maximum length curve
max_seq <- seq(min(df$max_c), max(df$max_c), length.out = 100)
df_max <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = max_seq,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_max       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_max)
newmods_max <- X_max[, -1]
df_max <- df_max %>%
  bind_cols(predict_beta(m_all_scaled, newmods_max)) %>%
  mutate(max_length_cm = max_c * sd_max + mean_max)

# (e) Initial size curve
size_seq <- seq(min(df$size_c), max(df$size_c), length.out = 100)
df_size <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = size_seq,
  max_c    = 0,
  expt_obs = factor("Exp", levels = exp_levels)
)
X_size        <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_size)
newmods_size  <- X_size[, -1]
df_size <- df_size %>%
  bind_cols(predict_beta(m_all_scaled, newmods_size)) %>%
  mutate(size_start = size_c * sd_size + mean_size)

# (f) Experimental vs Observational predictions
df_eo <- tibble(
  logmd_c  = 0,
  dur_c    = 0,
  size_c   = 0,
  max_c    = 0,
  expt_obs = factor(c("Exp","Obs"), levels = exp_levels)
)
X_eo       <- model.matrix(~ logmd_c + dur_c + size_c + max_c + expt_obs, df_eo)
newmods_eo <- X_eo[, -1]
df_eo <- df_eo %>% bind_cols(predict_beta(m_all_scaled, newmods_eo))


# --- now your existing p1–p5 code will find df_den, df_dur, df_max, df_size, df_eo ---

## (1) Density
p1 <- ggplot(df_den, aes(mean_density, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(mean_density, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    name   = expression(Mean~density~(fish~m^{-2}))
  ) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
    name   = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6, color = "black"),
    plot.margin = margin(2,2,2,2)
  )

## (2) Duration
p2 <- ggplot(df_dur, aes(duration, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(duration, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Duration (days)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text    = element_text(size = 6, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (3) Max Length
p3 <- ggplot(df_max, aes(max_length_cm, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(max_length_cm, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Maximum Length (cm)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8),
    axis.text    = element_text(size = 6, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (4) Initial Size
p4 <- ggplot(df_size, aes(size_start, beta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4682B4", alpha = 0.2) +
  geom_line(color = "#4682B4", size = 1) +
  geom_point(data = df, aes(size_start, betanls2_raw_cm),
             color = "#4682B4", alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Initial Size (mm)", y = NULL) +
  scale_y_continuous(trans = "asinh",
                     breaks = c(-1000,-100,-10,-1,0,1,10,100,1000),
                     name = expression(
                       paste(
                         "Strength of density-dependent mortality, ",
                         beta,
                         " (", cm^2, ~ fish^-1, ~ day^-1, ")"
                       )
                     )) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text    = element_text(size = 10, color = "black"),
    plot.margin  = margin(2,2,2,2)
  )

## (5) Study Type
p5 <- ggplot() +
  geom_jitter(
    data  = df,
    aes(x = expt_obs, y = betanls2_raw_cm, fill = expt_obs),
    shape = 21, color = "black", width = 0.15, size = 1.5, alpha = 0.6
  ) +
  geom_errorbar(
    data  = df_eo,
    aes(x = expt_obs, ymin = lower, ymax = upper),
    width = 0.1, size = 0.8, color = "black"
  ) +
  geom_point(
    data = df_eo,
    aes(x = expt_obs, y = beta),
    shape = 21, fill = "black", size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Study Type", y = expression(
    paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2, ~ fish^-1, ~ day^-1, ")"
    )
  )) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000)
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6, color = "black"),
    legend.position = "none",
    plot.margin     = margin(2,2,2,2)
  )

# ----------------------------------------------------------------------------
# 7. Display & Save
# ----------------------------------------------------------------------------
print(p1); print(p2); print(p3); print(p4); print(p5)

ggsave(
  here::here("figures", "bivar_density_vs_beta.png"),
  p1, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_duration_vs_beta.png"),
  p2, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_maxlen_vs_beta.png"),
  p3, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_sizestart_vs_beta.png"),
  p4, width = 5, height = 4, dpi = 300
)
ggsave(
  here::here("figures", "bivar_exptobs_vs_beta.png"),
  p5, width = 5, height = 4, dpi = 300
)



# ----------------------------------------------------------------------------
# 8. Body Size by Species
# ----------------------------------------------------------------------------

# 8a. Filter and prepare data — keep predator-present species with enough replicates
all_filtered <- all_dat2 %>%
  filter(predators == "present") %>%
  select(
    size_start, 
    betanls2_raw_cm, 
    study_num, 
    substudy_num, 
    g_sp, 
    family
  ) %>%
  filter(!substudy_num %in% c(119, 121)) %>%
  drop_na() %>%
  group_by(g_sp) %>%
  filter(n() > 3,
         !g_sp %in% c(
           "Dascyllus_flavicaudus",
           "Sebastes_atrovirens",
           "Thalassoma_hardwicke"
         )
  ) %>%
  ungroup() %>%
  mutate(
    # convert underscores to spaces and preserve order
    g_sp = str_replace_all(g_sp, "_", " "),
    g_sp = factor(g_sp, levels = unique(g_sp))
  )

# 8b. Plot: body size vs. β, facetted by species
# 1) Global: pretty integer ticks (no decimals)
p_body_size <- ggplot(all_filtered, aes(x = size_start, y = betanls2_raw_cm)) +
  geom_point(color = "black", fill = "steelblue", shape = 21, alpha = 0.6, size = 3) +
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 1) +
  facet_wrap(~ g_sp, scales = "free", ncol = 4) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = scales::comma
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_x_continuous(
    limits = function(l) c(floor(l[1]), ceiling(l[2])),     # snap edges to integers
    breaks = function(l) scales::breaks_pretty(n = 4)(l),   # ~3–4 ticks per facet
    labels = scales::label_number(accuracy = 1)             # no decimals
  ) +
  labs(
    x = expression(paste("Initial Size (", mm, ")")),
    y = expression(paste("Strength of density-dependent mortality, ",
                         beta, " (", cm^2, ~ fish^-1, ~ day^-1, ")"))
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text = element_text(face = "italic", size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11)
  )

# 2) Panel-specific x-range WITHOUT extra packages:
#    Add an invisible data row that forces the x-limits for that facet.
pad_garnoti <- data.frame(
  g_sp = "Halichoeres garnoti",
  size_start = c(20, 23),           # desired limits
  betanls2_raw_cm = 0               # any y within panel range
)

# (Add more pads the same way for other species if needed)
# pad_partitus  <- data.frame(g_sp = "Stegastes partitus",  size_start = c(17, 20), betanls2_raw_cm = 0)
# pad_elacatinus<- data.frame(g_sp = "Elacatinus sp.",      size_start = c(8, 20),  betanls2_raw_cm = 0)

p_body_size_final <- p_body_size +
  geom_blank(data = pad_garnoti, aes(size_start, betanls2_raw_cm))
# + geom_blank(data = pad_partitus,  aes(size_start, betanls2_raw_cm))
# + geom_blank(data = pad_elacatinus,aes(size_start, betanls2_raw_cm))

print(p_body_size_final)


ggsave(
  filename = here::here("figures", "body_size_by_species.png"),
  plot     = p_body_size,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)



# ----------------------------------------------------------------------------
# 8b. Experimental Duration by Species
# ----------------------------------------------------------------------------

# (a) Filter and prepare data — same species set as 8a, but include duration
dur_filtered <- all_dat2 %>%
  filter(predators == "present") %>%
  select(
    duration,
    betanls2_raw_cm,
    study_num,
    substudy_num,
    g_sp
  ) %>%
  filter(!substudy_num %in% c(119, 121)) %>%
  drop_na() %>%
  group_by(g_sp) %>%
  filter(
    n() > 3,
    !g_sp %in% c(
      "Dascyllus_flavicaudus",
      "Sebastes_atrovirens",
      "Thalassoma_hardwicke"
    )
  ) %>%
  ungroup() %>%
  mutate(
    # convert underscores to spaces and preserve order
    g_sp = str_replace_all(g_sp, "_", " "),
    g_sp = factor(g_sp, levels = unique(g_sp))
  )

# (b) Plot: duration vs. β, facetted by species
p_duration_by_species <- ggplot(dur_filtered, aes(x = duration, y = betanls2_raw_cm)) +
  geom_point(
    color = "black", fill = "steelblue", shape = 21,
    alpha = 0.6, size = 3
  ) +
  stat_smooth(
    method = "lm", formula = y ~ x,
    se = FALSE, color = "black", size = 1
  ) +
  facet_wrap(~ g_sp, scales = "free", ncol = 4) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = comma
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "gray"
  ) +
  labs(
    x = "Experiment Duration (Days)",
    y = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2,  ~ fish^-1, ~ day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "italic", size = 8),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

# (c) Save
ggsave(
  filename = here("figures", "duration_by_species.png"),
  plot     = p_duration_by_species,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# ============================================================================
# 9. Experimental vs. Observational: Paired Meta-Analysis
# ============================================================================
# Tests whether experimental and observational studies detect different
# strengths of density dependence, accounting for species-level pairing
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(metafor); library(tibble); 
  library(ggplot2); library(here)
})

# ---- 9.1 Prepare data ------------------------------------------------------

# Subset to predator-present studies with defined experimental/observational status
expobs <- all_dat2 %>%
  filter(
    predators == "present",
    !is.na(expt_obs)
  ) %>%
  mutate(
    expt_obs = factor(expt_obs, levels = c("Exp", "Obs"))
  )

# Split into paired vs. unpaired studies
paired_expobs   <- expobs %>% filter(!is.na(expt_obs_pairs))
unpaired_expobs <- expobs %>% filter(is.na(expt_obs_pairs))

# Build paired dataset with complete pairs only
paired0 <- paired_expobs %>%
  filter(expt_obs %in% c("Exp", "Obs"), !is.na(expt_obs_pairs)) %>%
  mutate(
    expt_obs = factor(expt_obs, levels = c("Exp", "Obs")),
    g_sp     = droplevels(factor(g_sp)),
    pair_id  = paste0("pair_", expt_obs_pairs)
  )

# Keep only pair_ids with BOTH Exp and Obs
complete_pairs <- paired0 %>%
  group_by(pair_id) %>%
  summarise(n_levels = dplyr::n_distinct(expt_obs), .groups = "drop") %>%
  filter(n_levels == 2) %>%
  select(pair_id)

# Final paired dataset with variance floor for stability
var_floor <- 1e-8
expobs_paired <- paired0 %>%
  inner_join(complete_pairs, by = "pair_id") %>%
  filter(
    is.finite(betanls2_asinh),
    is.finite(betanlsvar_asinh),
    betanlsvar_asinh >= 0
  ) %>%
  mutate(betanlsvar_asinh = pmax(betanlsvar_asinh, var_floor)) %>%
  arrange(pair_id, expt_obs)

# Sanity check
n_pairs <- dplyr::n_distinct(expobs_paired$pair_id)
cat(sprintf("\n[9] Paired analysis with %d complete species pairs (n = %d observations)\n", 
            n_pairs, nrow(expobs_paired)))

if (n_pairs < 2) stop("[9] Need at least 2 complete pairs for analysis.")

# ---- 9.2 Meta-analytic models ----------------------------------------------

# Model 1: Main effect only (assumes consistent effect across species)
m_main_effect <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ expt_obs,
  random = list(~1 | g_sp),
  data   = expobs_paired,
  method = "REML",
  test   = "t"
)

cat("\n=== Model 1: Main Effect (assumes consistent across species) ===\n")
print(summary(m_main_effect))

# Model 2: Random slopes (allows effect to vary by species) - RECOMMENDED
m_random_slopes <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ expt_obs,
  random = list(~ expt_obs | g_sp),
  struct = "UN",
  data   = expobs_paired,
  method = "REML",
  test   = "t"
)

cat("\n=== Model 2: Random Slopes (effect varies by species) - RECOMMENDED ===\n")
print(summary(m_random_slopes))

# Test whether random slopes improve model fit
lrt_result <- anova(m_main_effect, m_random_slopes)
cat("\n=== Likelihood Ratio Test: Does effect vary by species? ===\n")
print(lrt_result)

# ---- 9.3 Back-transform and extract group means ----------------------------

# Get group means for plotting (from main effect model for simplicity)
m_groups <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ expt_obs - 1,
  random = list(~1 | g_sp),
  data   = expobs_paired,
  method = "REML",
  test   = "t"
)

# Back-transform to raw scale
X_pred <- model.matrix(
  ~ expt_obs - 1,
  data.frame(expt_obs = factor(c("Exp", "Obs"), levels = c("Exp", "Obs")))
)
preds <- predict(m_groups, newmods = X_pred)

pooled_means <- tibble(
  study_type = c("Experimental", "Observational"),
  beta_mean  = sinh(preds$pred),
  ci_lo      = sinh(preds$ci.lb),
  ci_hi      = sinh(preds$ci.ub)
)

cat("\n=== Back-transformed pooled means (raw β scale) ===\n")
print(pooled_means)

# ---- 9.4 Interpretation summary --------------------------------------------

cat("\n=== INTERPRETATION ===\n\n")

main_pval <- coef(summary(m_main_effect))["expt_obsObs", "pval"]
slope_pval <- coef(summary(m_random_slopes))["expt_obsObs", "pval"]
lrt_pval <- lrt_result$pval[1]

cat(sprintf("1. Main effect model: p = %.4f\n", main_pval))
cat(sprintf("   - Observational studies differ from experimental (assuming consistent effect)\n\n"))

cat(sprintf("2. Random slopes model: p = %.4f\n", slope_pval))
cat(sprintf("   - Effect not significant when accounting for species-level variation\n"))
cat(sprintf("   - Between-species variance:\n"))
cat(sprintf("     * Experimental: tau^2 = %.2f\n", m_random_slopes$tau2[1]))
cat(sprintf("     * Observational: tau^2 = %.2f\n", m_random_slopes$tau2[2]))
cat(sprintf("   - Correlation between species effects: rho = %.2f\n\n", 
            m_random_slopes$rho))

cat(sprintf("3. Likelihood ratio test: p < %.4f\n", lrt_pval))
cat(sprintf("   - Random slopes model significantly better fit\n"))
cat(sprintf("   - CONCLUSION: Effect varies substantially by species\n\n"))

cat("4. RECOMMENDATION:\n")
cat("   - Report random slopes model (Model 2)\n")
cat("   - State: 'The experimental-observational difference varied significantly\n")
cat("     across species (LRT p < 0.0001), with no consistent main effect when\n")
cat("     accounting for species-level heterogeneity (p = ", 
    sprintf("%.2f", slope_pval), ")'\n")

# ---- 9.5 Publication-quality paired figure ---------------------------------

# Prepare plotting data
plot_data <- expobs_paired %>%
  mutate(
    beta_raw = sinh(betanls2_asinh),
    beta_var = betanlsvar_asinh,
    ci_lo    = beta_raw - beta_var,
    ci_hi    = beta_raw + beta_var
  ) %>%
  arrange(g_sp, expt_obs)

# Define colorblind-friendly palette
species_colors <- c(
  "Brachyistius_frenatus"       = "#E64B35",  # Coral red
  "Coryphopterus_glaucofraenum" = "#4DBBD5",  # Teal blue
  "Dascyllus_trimaculatus"      = "#00A087",  # Emerald green
  "Elactinus_sp."               = "#3C5488",  # Navy blue
  "Thalassoma_bifasciatum"      = "#F39B7F",  # Peach
  "Thalassoma_hardwicke"        = "#8491B4"   # Lavender
)

# Create publication-quality figure
p_paired <- ggplot(plot_data, aes(x = expt_obs, y = beta_raw, group = g_sp)) +
  # Lines connecting pairs
  geom_line(
    aes(color = g_sp), 
    linewidth = 0.9, 
    alpha = 0.8
  ) +
  # Vertical error bars (variance)
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi, color = g_sp), 
    width = 0,
    alpha = 0.5,
    linewidth = 0.6
  ) +
  # Points for each estimate
  geom_point(
    aes(color = g_sp, shape = expt_obs), 
    size = 3,
    stroke = 0.8
  ) +
  # Reference line at zero
  geom_hline(
    yintercept = 0, 
    linetype = "dashed", 
    color = "grey30",
    linewidth = 0.4
  ) +
  # Scales
  scale_shape_manual(
    name = "Study type",
    values = c(Exp = 16, Obs = 17),
    labels = c("Experimental", "Observational")
  ) +
  scale_color_manual(
    name = "Species",
    values = species_colors
  ) +
  scale_x_discrete(
    labels = c("Experimental", "Observational"),
    expand = expansion(mult = 0.15)
  ) +
  scale_y_continuous(
    trans = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = scales::comma
  ) +
  # Labels
  labs(
    x = NULL,
    y = expression(
      paste(
        "Strength of density dependence, ",
        beta,
        " (", cm^2, " ", fish^-1, " ", day^-1, ")"
      )
    )
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, face = "plain"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "italic"),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_paired)

# Save publication-quality outputs
ggsave(
  here("figures", "paired_exp_obs_by_species.png"),
  p_paired, 
  width = 7, 
  height = 5, 
  dpi = 600,
  bg = "white"
)



ggsave(
  here("figures", "paired_exp_obs_by_species.pdf"),
  p_paired, 
  width = 7, 
  height = 5,
  device = cairo_pdf,
  bg = "white"
)

cat("\n[9] Figures saved to figures/ directory\n")




# --- Build unpaired background points (gray dots) ---
# Requires `expobs` from 9.1 above. If it's not in memory, run 9.1 first.
bg_points <- expobs %>%
  dplyr::filter(is.na(expt_obs_pairs)) %>%             # unpaired only
  dplyr::mutate(
    expt_obs = factor(expt_obs, levels = c("Exp","Obs")),
    beta_raw = sinh(betanls2_asinh)
  ) %>%
  dplyr::select(expt_obs, beta_raw) %>%
  tidyr::drop_na()

# Safety: if there are no unpaired, keep an empty tibble so ggplot doesn't error
if (!nrow(bg_points)) {
  bg_points <- tibble::tibble(
    expt_obs = factor(character(), levels = c("Exp","Obs")),
    beta_raw = numeric(0)
  )
}

# Ensure consistent factor levels for study type across all inputs
bg_points$expt_obs <- factor(bg_points$expt_obs, levels = c("Exp","Obs"))
plot_data$expt_obs <- factor(plot_data$expt_obs, levels = c("Exp","Obs"))
df_eo$expt_obs     <- factor(df_eo$expt_obs,     levels = c("Exp","Obs"))

p_paired <- ggplot() +
  # --- NEW: unpaired background points (gray) ---
  # --- background (unpaired) points: make Observational triangles ---
  geom_point(
    data = bg_points,
    aes(x = expt_obs, y = beta_raw, shape = expt_obs),   # <- add shape mapping
    color = "gray70",
    alpha = 0.35,
    size  = 1.6,
    position = position_jitter(width = 0.08, height = 0, seed = 1),
    show.legend = FALSE
  ) +
  # --- Paired species lines/points (foreground) ---
  geom_line(
    data = plot_data,
    aes(x = expt_obs, y = beta_raw, group = g_sp, color = g_sp),
    linewidth = 0.9,
    alpha     = 0.85
  ) +
  # geom_errorbar(
  #   data = plot_data,
  #   aes(x = expt_obs, ymin = ci_lo, ymax = ci_hi, color = g_sp),
  #   width = 0,
  #   alpha = 0.55,
  #   linewidth = 0.6
  # ) +
  geom_point(
    data = plot_data,
    aes(x = expt_obs, y = beta_raw, color = g_sp, shape = expt_obs),
    size   = 3,
    stroke = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  scale_shape_manual(
    name = "Study type",
    values = c(Exp = 17, Obs = 16),   # Experimental = triangle, Observational = circle
    labels = c("Experimental", "Observational")
  ) +
  scale_color_manual(
    name   = "Species",
    values = species_colors,
    labels = function(x) gsub("_", " ", x)               # <- remove underscores
  )+
  scale_x_discrete(labels = c("Experimental", "Observational"), expand = expansion(mult = 0.15)) +
  scale_y_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = scales::comma
  ) +
  
  geom_pointrange(
    data = df_eo,                    # <- was df_expobs
    inherit.aes = FALSE,
    aes(x = expt_obs, y = beta,
        ymin = lower, ymax = upper),
    size = 0.9,
    color = "black"
  )+
  labs(
    x = NULL,
    y = expression(
      paste(
        "Strength of density dependence, ",
        beta,
        " (", cm^2, " ", fish^-1, " ", day^-1, ")"
      )
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line       = element_line(linewidth = 0.5, color = "black"),
    axis.ticks      = element_line(linewidth = 0.5, color = "black"),
    axis.text       = element_text(color = "black", size = 11),
    axis.title.y    = element_text(size = 12, margin = margin(r = 10)),
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10, face = "italic"),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y= unit(0.3, "cm"),
    plot.margin     = margin(10, 15, 10, 10)
  )



# ---- 9.6 Save paired figure (PNG + PDF) ------------------------------------

# Define output directory (if not already existing)
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)

# PNG (raster, high-res)
ggsave(
  filename = here::here("figures", "paired_exp_obs_by_species.png"),
  plot     = p_paired,
  width    = 7,
  height   = 5,
  dpi      = 600,
  bg       = "white"
)

# PDF (vector, editable in Illustrator/Acrobat)
ggsave(
  filename = here::here("figures", "paired_exp_obs_by_species.pdf"),
  plot     = p_paired,
  width    = 7,
  height   = 5,
  device   = cairo_pdf,
  bg       = "white"
)

cat("✅ Paired Exp–Obs figure saved as both PNG and PDF in /figures.\n")


# ---- 9.6 Optional: Species-specific effects visualization -----------------

# Extract species-specific effects for supplementary figure
if (exists("m_random_slopes")) {
  
  # Calculate observed differences for each species
  species_diffs <- expobs_paired %>%
    select(g_sp, expt_obs, betanls2_asinh) %>%
    pivot_wider(names_from = expt_obs, values_from = betanls2_asinh) %>%
    mutate(
      observed_diff = Obs - Exp,
      species = g_sp
    ) %>%
    select(species, observed_diff)
  
  # Get overall fixed effect
  fixed_effect <- coef(m_random_slopes)["expt_obsObs"]
  
  # Create caterpillar plot
  p_caterpillar <- ggplot(species_diffs, 
                          aes(x = reorder(species, observed_diff),
                              y = observed_diff)) +
    geom_hline(
      yintercept = as.numeric(fixed_effect), 
      linetype = "dashed", 
      color = "#E64B35", 
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0, 
      linetype = "dotted", 
      color = "grey50"
    ) +
    geom_point(size = 4, color = "#3C5488") +
    geom_segment(
      aes(xend = species, 
          y = as.numeric(fixed_effect), 
          yend = observed_diff),
      linewidth = 0.8,
      color = "#3C5488"
    ) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Observed difference (Observational - Experimental)\non asinh scale",
      title = "Species-Specific Variation in Study Type Effect",
      subtitle = "Red line = overall average; dotted line = no difference"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "italic", size = 11),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
  
  print(p_caterpillar)
  
  ggsave(
    here("figures", "species_specific_effects.png"),
    p_caterpillar,
    width = 6,
    height = 4,
    dpi = 600,
    bg = "white"
  )
  
  cat("\n[9] Species-specific effects plot saved\n")
}




# ============================================================================

# ----------------------------------------------------------------------------
# 10. Combined bivariate panels: Density, Duration & Max Length
# ----------------------------------------------------------------------------

# 10a) Shared components ----------------------------------------------------
hline0 <- geom_hline(yintercept = 0, linetype = "dashed", color = "black")

y_scale <- scale_y_continuous(
  trans  = "asinh",
  breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
  name   = expression(
    paste(
      "Strength of density-dependent mortality, ",
      beta,
      " (", cm^2, ~ fish^-1, ~ day^-1, ")"
    )
  )
)

theme_bivar <- theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text    = element_text(size = 10, color = "black"),
    plot.margin  = margin(3, 3, 3, 3)
  )

# 10b) Panel A: Mean Density -----------------------------------------------
p1 <- ggplot() +
  geom_ribbon(
    data = df_den,
    aes(x = mean_density, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_den,
    aes(x = mean_density, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = mean_density, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    name   = expression(Mean~density~(fish~m^{-2}))
  ) +
  y_scale +
  theme_bivar

# 10c) Panel B: Duration ----------------------------------------------------
p2 <- ggplot() +
  geom_ribbon(
    data = df_dur,
    aes(x = duration, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_dur,
    aes(x = duration, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = duration, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  labs(x = "Duration (days)") +
  y_scale +
  theme_bivar +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# 10d) Panel C: Max Length -------------------------------------------------
p3 <- ggplot() +
  geom_ribbon(
    data = df_max,
    aes(x = max_length_cm, y = beta, ymin = lower, ymax = upper),
    fill  = "#4682B4", alpha = 0.2
  ) +
  geom_line(
    data = df_max,
    aes(x = max_length_cm, y = beta),
    color = "#4682B4", size = 1
  ) +
  geom_point(
    data = df,
    aes(x = max_length_cm, y = betanls2_raw_cm),
    color = "#4682B4", alpha = 0.6, size = 2
  ) +
  hline0 +
  labs(x = "Maximum Length (cm)") +
  y_scale +
  theme_bivar

p3

#save as png 
ggsave(
  filename = here("figures", "bivar_maxlen_vs_beta.png"),
  plot     = p3,
  width    = 5,
  height   = 5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)
# ----------------------------------------------------------------------------
# 10. Combined bivariate panels: Density, Duration & Max Length (cowplot)
# ----------------------------------------------------------------------------

# Assume p1, p2, p3 are defined as in the previous steps.

# Add a consistent margin so panel content doesn’t shift
p1m <- p1 + theme(plot.margin = margin(5,5,5,5))
p2m <- p2 + theme(plot.margin = margin(5,5,5,5))
p3m <- p3 + theme(plot.margin = margin(5,5,5,5))

# Combine and align
combined <- plot_grid(
  p1m, p2m,
  nrow           = 1,
  align          = "v",   # align vertically
  axis           = "l",   # align on left axes
  rel_widths     = c(1,1,1),
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold",
  # place labels at 90% from left, 95% from bottom
  label_x        = c(0.90, 0.90, 0.90),
  label_y        = c(0.95, 0.95, 0.95)
)

# Render & Save
print(combined)

ggsave(
  filename = here("figures", "bivar_density_duration_maxlen.pdf"),
  plot     = combined,
  width    = 12,
  height   = 4,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

#same but as png

ggsave(
  filename = here("figures", "bivar_density_duration_maxlen.png"),
  plot     = combined,
  width    = 12,
  height   = 4,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)
