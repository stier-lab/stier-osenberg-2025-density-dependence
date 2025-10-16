# Script: 3_phylo_species_analysis.R
# Goal:   Estimate species-level β effects, map them onto the pruned phylogeny,
#         quantify phylogenetic signal, compare phylogenetic vs. non-phylogenetic
#         REML models, partition variance (wrasse example), and export publication
#         figures/tables for species and families.
# Inputs: all_dat2 and pruned_tree from 1_data_phylogeny_loading.R
# Outputs: results/species_* CSVs, results/variance_* files,
#          figures/phylo_beta_circular.png, figures/Figure3_* images, console logs.

# ----------------------------------------------------------------------------
# 1. Dependencies & Data Load
# ----------------------------------------------------------------------------
# These scripts all run in sequence via 00_run_all.R, but we guard against
# interactive use by sourcing the shared loader and checking the key objects.
source(here::here("code", "0_libraries.R"))

if (!exists("all_dat2", inherits = FALSE) ||
    !exists("pruned_tree", inherits = FALSE)) {
  message("Core data objects missing; sourcing 1_data_phylogeny_loading.R.")
  source(here::here("code", "1_data_phylogeny_loading.R"))
}

# Ensure output directories exist
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("results"), showWarnings = FALSE, recursive = TRUE)

# Define color palette
pal <- wesanderson::wes_palette("Zissou1", n = 100, type = "continuous")
# ----------------------------------------------------------------------------
# 2. Prepare species-level dataset
# ----------------------------------------------------------------------------
# Collapse to one row per substudy/species, keep the asinh-transformed effect
# and any traits used later.  The name standardisation step matches tree tip
# labels and stabilises joins across multiple data sources.
species_df <- all_dat2 %>%
  dplyr::select(
    g_sp,
    betanls2_raw_cm,
    betanls2_asinh,
    betanlsvar_asinh,
    study_num,
    substudy_num,
    family,
    predators,
    paired_pred,
    expt_obs_pairs,
    duration,
    mean_density,
    max_length_cm,
    max_weight_g
  ) %>%
  # Standardize a few species names
  dplyr::mutate(
    g_sp = stringr::str_replace_all(g_sp, c(
      "Sebastes_spp_NA"           = "Sebastes_mystinus",
      "Fundulus_heteroclitus"     = "Fundulus_heteroclitus_heteroclitus",
      "Diplodus_sargus"           = "Diplodus_sargus_sargus",
      "Elactinus_sp\\."           = "Elacatinus_evelynae"
    ))
  ) %>%
  dplyr::filter(!is.na(betanls2_asinh), !is.na(betanlsvar_asinh))

# ----------------------------------------------------------------------------
# 3. Species-level meta-analysis (no phylogeny)
# ----------------------------------------------------------------------------
# First, estimate species means ignoring phylogenetic structure so the results
# can be compared against the phylogenetic models below and plotted directly.
sp_model <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ g_sp - 1,                   # species-specific estimates
  random = list(~1 | study_num/substudy_num),
  data   = species_df,
  method = "REML"
)

sp_coefs <- coef(summary(sp_model)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  dplyr::rename(
    est_asinh   = estimate,
    ci_lb_asinh = ci.lb,
    ci_ub_asinh = ci.ub
  ) %>%
  dplyr::mutate(
    g_sp     = stringr::str_remove(term, "^g_sp"),
    beta_est = sinh(est_asinh),
    ci.lb    = sinh(ci_lb_asinh),
    ci.ub    = sinh(ci_ub_asinh)
  ) %>%
  dplyr::select(g_sp, beta_est, ci.lb, ci.ub)

readr::write_csv(sp_coefs, here::here("results", "species_effects_no_phylo.csv"))

# ----------------------------------------------------------------------------
# 4. Map species estimates onto phylogeny
# ----------------------------------------------------------------------------
# Create data frame for plotting
beta_plot_df <- sp_coefs %>%
  dplyr::rename(mean = beta_est)

# Simple circular tree with colored points
p_phylo_simple <- ggtree::ggtree(pruned_tree, layout = "circular") %<+% beta_plot_df +
  ggtree::geom_tippoint(aes(color = mean), size = 3) +
  ggplot2::scale_color_viridis_c(trans = "asinh", name = "Beta") +
  ggtree::geom_tiplab(size = 2, offset = 1.5) +
  ggplot2::theme(legend.position = "right")

# Detailed circular tree: black outline + white-to-red fill
p_phylo <- ggtree::ggtree(pruned_tree, layout = "circular") %<+% beta_plot_df +
  ggtree::geom_tippoint(
    aes(fill = mean),
    shape = 21, size = 5, color = "black", stroke = 1
  ) +
  ggplot2::scale_fill_gradientn(
    colours = pal,
    name    = "Beta",
    trans   = "asinh",
    limits  = c(min(beta_plot_df$mean), max(beta_plot_df$mean)),
    breaks  = c(-10, 0, 10, 100, 1000, 10000),
    labels  = c("-10", "0", "10", "100", "1,000", "10,000"),
    guide   = guide_colorbar(
      frame.colour    = "black",
      frame.linewidth = 0.5,
      ticks.colour    = "black",
      ticks.linewidth = 0.5
    )
  ) +
  ggtree::geom_tiplab(size = 2, offset = 2) +
  ggtree::theme_tree2() +
  ggplot2::theme(
    axis.line   = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

# Save the detailed plot
ggplot2::ggsave(
  filename = here::here("figures", "phylo_beta_circular.png"),
  plot     = p_phylo,
  width    = 8, height = 8, dpi = 300
)


# ----------------------------------------------------------------------------
# 5. Phylogenetic signal tests (with parametric bootstrap CI for observed K)
# ----------------------------------------------------------------------------
# We quantify the strength of phylogenetic signal in β using Pagel's λ and
# Blomberg's K.  A parametric bootstrap gives uncertainty for K on the BM null.

# 5.1 Prepare trait vector named by species
betavec <- beta_plot_df$mean
names(betavec) <- beta_plot_df$g_sp

# 5.2 Keep only species present in both tree and data
common_sp <- intersect(pruned_tree$tip.label, names(betavec))
betavec   <- betavec[common_sp]
pruned_tree <- ape::drop.tip(
  pruned_tree,
  setdiff(pruned_tree$tip.label, common_sp)
)

# 5.3 Pagel's λ test (likelihood‐ratio against λ = 0)
lambda_res <- phytools::phylosig(
  pruned_tree, betavec,
  method = "lambda", test = TRUE
)

# 5.4 Observed Blomberg's K + permutation p-value
k_test <- phytools::phylosig(
  pruned_tree, betavec,
  method = "K", test = TRUE
)

k_obs  <- k_test$K      # your observed K
k_pval <- k_test$P      # permutation p-value

# ----------------------------------------------------------------------------
# 5.5 Parametric bootstrap under a BM model (σ² estimated from your data)
# ----------------------------------------------------------------------------
# 1) Fit BM to get σ²̂
bm_fit     <- geiger::fitContinuous(pruned_tree, betavec, model = "BM")
sigma2_hat <- bm_fit$opt$sigsq
# sigma2_hat is our MLE of σ² (rate of trait evolution)

# 2) Bootstrap K under BM(σ²̂)
nsim   <- 500
boot_K <- numeric(nsim)

for (i in seq_len(nsim)) {
  sim_trait <- phytools::fastBM(pruned_tree, sig2 = sigma2_hat, internal = FALSE)
  boot_K[i] <- phytools::phylosig(
    pruned_tree, sim_trait,
    method = "K", test = FALSE
  )
}

# 3) 95% CI around the observed K
ci_K_obs <- quantile(boot_K, probs = c(0.025, 0.975))

# ----------------------------------------------------------------------------
# 5.6 Report results
# ----------------------------------------------------------------------------
cat("Pagel's λ:\n")
print(lambda_res)

cat("\nBlomberg's K:\n",
    "  Observed K     =", round(k_obs, 3), "\n",
    "  Permutation p =", signif(k_pval, 3), "\n")

cat("\n95% CI for K (parametric bootstrap):\n")
print(ci_K_obs)




# ----------------------------------------------------------------------------
# 6. Compare models: no-phylo vs. phylo-random (AIC + LRT)
# ----------------------------------------------------------------------------
# Fit REML models with and without a phylogenetic random effect for the species
# intercept.  The boundary-corrected LRT follows standard meta-analytic practice.

# Build phylogenetic VCV matrix
phylo_vcv <- ape::vcv(pruned_tree, corr = TRUE)

# Restrict to species present in the VCV matrix
species_phylo <- species_df %>%
  dplyr::filter(g_sp %in% rownames(phylo_vcv)) %>%
  droplevels()

# (a) Null model: no species-level term
m_nosp <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ 1,
  random = list(~ 1 | study_num/substudy_num),
  data   = species_phylo,
  method = "REML"
)

# (b) Alternative: species-level term with phylogenetic correlation
m_gsp <- metafor::rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ 1,
  random = list(~ 1 | study_num/substudy_num, ~ 1 | g_sp),
  R      = list(g_sp = phylo_vcv),
  data   = species_phylo,
  method = "REML"
)

# ---- AIC ----
comp_aic <- AIC(m_nosp, m_gsp)
print(comp_aic)
readr::write_csv(as.data.frame(comp_aic),
                 here::here("results", "model_comparison_aic.csv"))

# ---- REML LRT (boundary-corrected) ----
ll0 <- as.numeric(logLik(m_nosp))   # store once, reuse if needed later
ll1 <- as.numeric(logLik(m_gsp))
LRT <- 2 * (ll1 - ll0)
p_chi     <- pchisq(LRT, df = 1, lower.tail = FALSE)
p_mixture <- 0.5 * p_chi   # 0.5*χ²0 + 0.5*χ²1

lrt_tab <- data.frame(
  logLik_no_phylo    = ll0,
  logLik_phylo       = ll1,
  LRT                = LRT,
  df                 = 1L,
  p_chisq_df1        = p_chi,
  p_boundary_mixture = p_mixture
)

print(lrt_tab)
readr::write_csv(lrt_tab,
                 here::here("results", "species_variance_LRT.csv"))

# ----------------------------------------------------------------------------
# 7. Inference on species variance (phylo-random effect): CIs (and optional bootstrap)
# ----------------------------------------------------------------------------

# Profile-likelihood CI for the species variance component (component #3)
ci_sigma3_obj <- confint(m_gsp, sigma2 = 3, level = 0.95)

# Robust coercion helper (handles matrix/data.frame or list with $random/$beta)
to_ci_df <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
      df <- tibble::rownames_to_column(df, "parameter")
    } else if (is.null(df$parameter)) {
      df$parameter <- rownames(x)
    }
    rownames(df) <- NULL
    return(df)
  }
  if (!is.null(x$random)) {
    df <- as.data.frame(x$random, stringsAsFactors = FALSE)
    df <- tibble::rownames_to_column(df, "parameter")
    rownames(df) <- NULL
    return(df)
  }
  if (!is.null(x$beta)) {
    df <- as.data.frame(x$beta, stringsAsFactors = FALSE)
    df <- tibble::rownames_to_column(df, "parameter")
    rownames(df) <- NULL
    return(df)
  }
  stop("Unhandled confint() return structure; consider updating 'metafor'.")
}

ci_sigma3_df <- to_ci_df(ci_sigma3_obj)
print(ci_sigma3_df)
readr::write_csv(ci_sigma3_df, here::here("results", "species_variance_CI.csv"))

# (Optional) Parametric bootstrap LRT for extra rigor (run only if needed)
# set.seed(1)
# B <- 1000
# LRT_sim <- numeric(B)
# for (b in 1:B) {
#   sim <- metafor::simulate_rma_mv(m_nosp)  # or simulate_rma.mv() on older versions
#   fit0 <- update(m_nosp, yi = sim$yi, V = sim$V)
#   fit1 <- update(m_gsp,  yi = sim$yi, V = sim$V)
#   LRT_sim[b] <- 2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0)))
# }
# p_boot <- mean(LRT_sim >= LRT)
# print(tibble::tibble(LRT = LRT, p_boundary_mixture = p_mixture, p_boot = p_boot))
# readr::write_csv(tibble::tibble(LRT = LRT, p_boundary_mixture = p_mixture, p_boot = p_boot),
#                  here::here("results", "species_variance_LRT_boot.csv"))


# ----------------------------------------------------------------------------
# 8. Wrasse example: within vs. among variance
# ----------------------------------------------------------------------------
# Demonstration of variance partitioning for a focal clade (wrasses).  Helpful
# for narrative text and to double-check that among-species variance is material.
wrasse_data <- all_dat2 %>%
  filter(g_sp %in% c("Thalassoma_bifasciatum", "Halichoeres_garnoti")) %>%
  select(g_sp, study_num, substudy_num, betanls2_raw_cm)


# 1) Compute species means and overall mean
species_means <- wrasse_data %>%
  group_by(g_sp) %>%
  summarise(mu_sp = mean(betanls2_raw_cm), n_sp = n(), .groups="drop")

grand_mu <- with(species_means, sum(mu_sp * n_sp) / sum(n_sp))

# 2) Among-species variance
among_var <- sum( (species_means$mu_sp - grand_mu)^2 * species_means$n_sp ) /
  (sum(species_means$n_sp) - 1)

# 3) Within-species variance (weighted average)
within_var <- wrasse_data %>%
  left_join(species_means, by="g_sp") %>%
  mutate(dev = betanls2_raw_cm - mu_sp) %>%
  summarise(w = sum(dev^2) / (n() - n_distinct(g_sp))) %>%
  pull(w)

# 4) Percent partition
partition <- tibble(
  component = c("within_species", "among_species"),
  variance  = c(within_var, among_var)
) %>%
  mutate(percent = variance / sum(variance) * 100)

partition

# ----------------------------------------------------------------------------
# 9. Plot species by family ordered β
# ----------------------------------------------------------------------------
# Several flavours of species-level plots appear in the supplement.  The code
# below keeps the wrangling steps explicit so the same ordered factors can feed
# multiple styling choices (Darjeeling palette, phylogenetic ordering, etc.).


# ──────────────────────────────────────────────────────────────────────
# 1. Prepare the data exactly as for the model (incl. study IDs!)
# ──────────────────────────────────────────────────────────────────────
all2 <- all_dat2 %>%
  filter(!is.na(betanls2_asinh),
         !is.na(betanlsvar_asinh),
         predators == "present") %>%     # correct level
  rename(beta_raw = betanls2_raw_cm) %>%
  select(g_sp, study_num, substudy_num, beta_raw, betanls2_asinh, betanlsvar_asinh) %>%
  droplevels()

# ──────────────────────────────────────────────────────────────────────
# 2. Fit the random‐effects species model (no intercept)
# ──────────────────────────────────────────────────────────────────────
m_sp_ni <- rma.mv(
  yi     = betanls2_asinh,
  V      = betanlsvar_asinh,
  mods   = ~ g_sp - 1,
  random = list(~1 | study_num/substudy_num),
  data   = all2,
  method = "REML",
  test   = "t"
)

# ──────────────────────────────────────────────────────────────────────
# 3. Pull out which species actually went into that model
# ──────────────────────────────────────────────────────────────────────
model_species <- sub("^g_sp","", rownames(coef(summary(m_sp_ni))))

# ──────────────────────────────────────────────────────────────────────
# 4. Extract & back‐transform estimates + CIs
# ──────────────────────────────────────────────────────────────────────
sp_coefs <- coef(summary(m_sp_ni)) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  rename(
    est_asinh   = estimate,
    ci.lb_asinh = ci.lb,
    ci.ub_asinh = ci.ub
  ) %>%
  mutate(
    g_sp     = sub("^g_sp","",term),
    beta_est = sinh(est_asinh),
    ci.lb    = sinh(ci.lb_asinh),
    ci.ub    = sinh(ci.ub_asinh)
  ) %>%
  select(g_sp, beta_est, ci.lb, ci.ub)

# ──────────────────────────────────────────────────────────────────────
# 5. Restrict both to the modeled species
# ──────────────────────────────────────────────────────────────────────
all2     <- all2    %>% filter(g_sp %in% model_species) %>% droplevels()
sp_coefs <- sp_coefs %>% filter(g_sp %in% model_species)

# ──────────────────────────────────────────────────────────────────────
# 6. Add fish‐family for coloring
# ──────────────────────────────────────────────────────────────────────
family_df <- all_dat2 %>% distinct(g_sp, family)
all2     <- left_join(all2,     family_df, by = "g_sp")
sp_coefs <- left_join(sp_coefs, family_df, by = "g_sp")

# ──────────────────────────────────────────────────────────────────────
# 7. Order families by mean β, then species by β within each family
# ──────────────────────────────────────────────────────────────────────
family_order <- sp_coefs %>%
  group_by(family) %>%
  summarise(mean_beta = mean(beta_est)) %>%
  arrange(desc(mean_beta)) %>%
  pull(family)

species_order <- sp_coefs %>%
  mutate(family = factor(family, levels = family_order)) %>%
  arrange(family, desc(beta_est)) %>%
  pull(g_sp)

all2$g_sp     <- factor(all2$g_sp,     levels = species_order)
sp_coefs$g_sp <- factor(sp_coefs$g_sp, levels = species_order)
# ──────────────────────────────────────────────────────────────────────
# After step 7 (you already have `family_order` and `species_order`)
# ──────────────────────────────────────────────────────────────────────

# 7b) make sure `family` is a factor in the same order
all2$family    <- factor(all2$family,    levels = family_order)
sp_coefs$family <- factor(sp_coefs$family, levels = family_order)

# 8a) drop any non‐finite β’s so geom_point/geom_jitter won’t silently remove them
all2     <- all2    %>% filter(is.finite(beta_raw))
sp_coefs <- sp_coefs %>% filter(is.finite(beta_est), is.finite(ci.lb), is.finite(ci.ub))

# ──────────────────────────────────────────────────────────────────────
# After you’ve built all2, sp_coefs, family_order & species_order…
# ──────────────────────────────────────────────────────────────────────

# 1) Build Darjeeling1 ×9
palette_darjeeling9 <- wes_palette("Darjeeling1",
                                   type = "continuous",
                                   n    = 9)

# 2) Re-draw your plot with that palette
p_sp_darjeeling <- ggplot() +
  # species‐level point estimates
  geom_point(
    data    = sp_coefs,
    aes(x = beta_est, y = g_sp, colour = family),
    shape = 17, size = 4
  ) +
  # species‐level CIs
  geom_errorbarh(
    data    = sp_coefs,
    aes(y = g_sp, xmin = ci.lb, xmax = ci.ub, colour = family),
    height = 0, size = 0.8
  ) +
  # raw sub‐study dots
  geom_jitter(
    data    = all2,
    aes(x = beta_raw, y = g_sp, colour = family),
    width   = 0.2, height = 0, size = 2, alpha = 0.5
  ) +
  
  
  # vertical zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  # asinh x‐axis
  scale_x_continuous(
    trans  = "asinh",
    breaks = c(-1000,-100,-10,-1,0,1,10,100,1000,10000),
    labels = scales::comma
  ) +
  # italic y‐axis labels (no underscores)
  scale_y_discrete(
    labels = function(x) {
      txt <- gsub("_","~", x)
      parse(text = paste0("italic(", txt, ")"))
    }
  ) +
  # Darjeeling1×9 manual palette, legend in same top→bottom order
  scale_colour_manual(
    values = palette_darjeeling9,
    breaks = family_order,
    name   = "Fish Family",
    guide  = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", 
        cm^2,   # cm²
        ~ fish^-1, 
        ~ day^-1,
        ")"
      )
    ),
    y = "",
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted")
  )

# 3) Render & save
print(p_sp_darjeeling)
ggsave("figures/Figure3_species_by_family_order_beta_darjeeling9.png",
       p_sp_darjeeling,
       width = 10, height = 8,
       units = "in", dpi = 300,bg = "white")



#alternative ordering 

# ──────────────────────────────────────────────────────────────────────
# After you’ve built all2, sp_coefs, family_order & species_order…
# Phylogenetic reordering + continuous Darjeeling1 palette
# ──────────────────────────────────────────────────────────────────────
# 1) Define the phylogenetic family order
family_phylo_order <- c(
  "Embiotocidae",  # perches
  "Plesiopidae",   # longfins
  "Pomacentridae", # damselfishes
  "Gobiidae",      # gobies
  "Serranidae",    # groupers
  "Sparidae",      # porgies
  "Cottidae",      # sculpins
  "Sebastidae",    # rockfishes
  "Labridae"       # wrasses
)

# 2) Re‐level your factor in both data frames
sp_coefs <- sp_coefs %>%
  mutate(family = factor(family, levels = family_phylo_order))

all2 <- all2 %>%
  mutate(family = factor(family, levels = family_phylo_order))

# 3) Grab a continuous 9‐colour Darjeeling1 ramp
palette_darjeeling9 <- wes_palette(
  name = "Darjeeling1",
  type = "continuous",
  n    = length(family_phylo_order)
)

# ──────────────────────────────────────────────────────────────────────
# After you’ve built all2, sp_coefs, family_order & species_order…
# Use RColorBrewer Paired palette with yellow removed and legend reversed
# ──────────────────────────────────────────────────────────────────────
# 1) Lock your species factor so the y‐axis follows species_order
sp_coefs <- sp_coefs %>%
  mutate(g_sp = factor(g_sp, levels = species_order))
all2 <- all2 %>%
  mutate(g_sp = factor(g_sp, levels = species_order))

# 2) Extract the families in the order they appear on the y‐axis
legend_family_order <- sp_coefs %>%
  arrange(g_sp) %>%
  pull(family) %>%
  unique()

# 3) Grab the Paired palette (12 colours), drop the yellow, then take as many as you need
raw_paired <- brewer.pal(12, "Paired")
# "#FFFF99" is the yellow in that palette
paired_no_yellow <- raw_paired[ raw_paired != "#FFFF99" ]
qual_cols <- paired_no_yellow[1:length(legend_family_order)]
names(qual_cols) <- legend_family_order

# 4) Build the plot: raw dots → CIs → means
p_sp_clean <- ggplot() +
  # a) raw sub‐study points behind
  geom_jitter(
    data    = all2,
    aes(x = beta_raw, y = g_sp, colour = family),
    width   = 0.2, height = 0, size = 2, alpha = 0.5
  ) +
  # b) species‐level confidence intervals
  geom_errorbarh(
    data    = sp_coefs,
    aes(y = g_sp, xmin = ci.lb, xmax = ci.ub, colour = family),
    height = 0, size = 0.8
  ) +
  # c) species‐level mean triangles on top
  geom_point(
    data    = sp_coefs,
    aes(x = beta_est, y = g_sp, colour = family),
    shape = 17, size = 4
  ) +
  # vertical zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  # asinh x‐axis
  scale_x_continuous(
    trans  = "asinh",
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  # italic species names on y‐axis
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic(", gsub("_","~",x), ")"))
  ) +
  # manual Paired palette (no yellow) + enforce legend order
  scale_colour_manual(
    values = qual_cols,
    limits = legend_family_order,
    name   = "Fish Family"
  ) +
  # reverse legend so it matches top→bottom of the plot
  guides(colour = guide_legend(reverse = TRUE)) +
  labs(
    x = expression(
      paste(
        "Strength of density-dependent mortality, ",
        beta,
        " (", cm^2, ~ fish^-1, ~ day^-1, ")"
      )
    ),
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted")
  )

# 5) Render & save
print(p_sp_clean)

ggsave(
  "figures/Figure3_species_by_family_order_beta_clean.png",
  plot  = p_sp_clean,
  width = 10, height = 8,
  units = "in", dpi = 300,
  bg    = "white"
)
