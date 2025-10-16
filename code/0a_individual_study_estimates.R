###############################################################################
# Script: 0a_individual_study_estimates.R
# Goal:   Estimate alpha and beta parameters for each substudy using either
#         binomial likelihood (manual grid search) or NLS2 (for "Other_No_Counts").
#         Combines results to produce combined_results_2025-10-15.csv for main pipeline.
#
# IMPORTANT FLAGS:
#   - SKIP_IF_EXISTS (line 22): Set to TRUE to skip if combined_results_2025-10-15.csv exists
#   - USE_FAST_LOADER (line ~505): Set to FALSE for proper NLS2 estimation (slow, hours)
#                                   Set to TRUE to load existing Other_No_Counts results (fast, 2 min)
#
# SMART CACHING:
#   If combined_results_2025-10-15.csv exists, the script will automatically:
#   - Load existing NLS2 results for substudies already present in the file
#   - Only run NLS2 fitting for substudies that are missing
#   This allows for incremental updates without re-running expensive computations
#
# Inputs: data/all_studies_looped-2024-09-11.csv (raw mortality time-series),
#         data/covariates-2024-09-30.csv (study metadata & filters),
#         data/range8132024.csv (parameter search ranges for grid search)
# Outputs: output/combined_results_2025-10-15.csv (used by 1_data_phylogeny_loading.R),
#          output/substudy_plots_2025-10-15.pdf (4-panel diagnostic plots),
#          output/substudy_other_nocount_plots.pdf (NLS2 fits),
#          output/substudy_results_2025-10-15.csv (detailed Binomial results)
###############################################################################

# ----------------------------------------------------------------------------
# 0. Skip-if-exists flag (for pipeline integration)
# ----------------------------------------------------------------------------
SKIP_IF_EXISTS <- TRUE  # Set to FALSE to force regeneration

if (!SKIP_IF_EXISTS || !file.exists("output/combined_results_2025-10-15.csv")) {
  # Only run if we're not skipping or if the file doesn't exist
  
  cat("\n================================================================================\n")
  cat("           RUNNING: Generating combined_results_2025-10-15.csv                  \n")
  cat("================================================================================\n\n")
  
  # ----------------------------------------------------------------------------
  # 1. Setup & Dependencies
  # ----------------------------------------------------------------------------
  # All packages loaded via 0_libraries.R (sourced by 00_run_all.R)
  # Packages used: tidyverse, gridExtra, kableExtra, nls2, ggtext

  # Create output directories
  dir.create("output/individual_plots/", showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # 2. Load Data
  # ----------------------------------------------------------------------------
  # Load raw mortality time-series data
  all_data <- read.csv("data/all_studies_looped-2024-09-11.csv")
  names(all_data) <- c("study_num", "substudy_num", "genus_species", "settlers", "recruits",
                       "t", "area", "no_pts", "other_factor", "pdf_num", "notes")
  
  # Load covariates (UPDATED to Sept 30 version)
  covariates <- read.csv("data/covariates-2024-09-30.csv")
  
  # Load parameter search ranges for grid search
  likdf2 <- read.csv("data/range8132024.csv")
  
  # Validate that required columns exist in range file
  required_cols <- c("substudy_num", "alpha_min", "alpha_max", "beta_min", "beta_max")
  if (!all(required_cols %in% names(likdf2))) {
    stop("Missing required columns in range8132024.csv: ",
         paste(setdiff(required_cols, names(likdf2)), collapse = ", "))
  }
  
  # Filter for substudies to use
  usedf <- covariates %>%
    select(use_2024, substudy_num, estimation_approach)
  
  # Join and filter data
  all_data <- inner_join(all_data, usedf, by = "substudy_num")
  all_data <- filter(all_data, use_2024 == "yes")
  
  # Validate we have data after filtering
  if (nrow(all_data) == 0) {
    stop("No data remaining after filtering for use_2024 == 'yes'")
  }
  
  # Calculate weights for likelihood estimation
  all_data$weight <- all_data$area * all_data$no_pts * all_data$other_factor
  
  # Check estimation approaches
  unique_estimation_approach <- unique(all_data$estimation_approach)
  cat("\n=== Data Summary ===\n")
  cat("Estimation approaches in data:", paste(unique_estimation_approach, collapse = ", "), "\n")
  
  count_per_estimation_approach <- table(all_data$estimation_approach)
  cat("Count per estimation approach:\n")
  print(count_per_estimation_approach)
  cat("\n")
  
  # ----------------------------------------------------------------------------
  # 3. Core Functions
  # ----------------------------------------------------------------------------
  
  # Beverton-Holt recruitment function
  bh_f <- function(alpha, beta, t_days, x) {
    predicted_recruits <- (x * exp(-alpha * t_days)) /
      (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
    return(predicted_recruits)
  }
  
  # Beverton-Holt function with safeguards
  bh_f_safeguarded <- function(alpha, beta, t_days, x) {
    predicted_recruits <- (x * exp(-alpha * t_days)) /
      (1 + (beta * x * (1 - exp(-alpha * t_days)) / alpha))
    # Ensure recruits are within [0, settlers]
    predicted_recruits[predicted_recruits > x] <- x[predicted_recruits > x]
    return(predicted_recruits)
  }
  
  # Negative log-likelihood function (binomial)
  neg_loglik_safeguarded <- function(params, t_days, settlers, recruits, weight) {
    alpha <- params[1]
    beta <- params[2]
    y_hat <- bh_f(alpha, beta, t_days, settlers)
    probs <- pmax(1e-10, pmin(y_hat / settlers, 1 - 1e-10))
    ll <- sum(log(probs) * recruits * weight +
                log(1 - probs) * (settlers * weight - recruits * weight))
    
    if (is.na(ll) || is.infinite(ll)) {
      return(.Machine$double.xmax)  # Return very large value if likelihood is NA or Inf
    }
    
    return(-ll)
  }
  
  # Function to compute second derivative (curvature) of negative log-likelihood
  second_derivative <- function(f, param, index, h = 1e-5, ...) {
    base_params <- param
    base_params[index] <- base_params[index] - h
    f1 <- f(base_params, ...)
    
    base_params[index] <- base_params[index] + 2 * h
    f2 <- f(base_params, ...)
    
    base_params[index] <- param[index]  # Reset to original
    f0 <- f(base_params, ...)
    
    return((f1 + f2 - 2 * f0) / (h^2))
  }
  
  # ----------------------------------------------------------------------------
  # 4. Manual Grid Search Function (for Binomial approach)
  # ----------------------------------------------------------------------------
  manual_neg_loglik <- function(substudy, all_data, results_df, likdf2) {
    # Check if the substudy is in the valid list
    if (!(substudy %in% unique(all_data$substudy_num))) {
      return(NULL)
    }
    
    data_sub <- filter(all_data, substudy_num == substudy)
    t_days <- data_sub$t
    settlers <- data_sub$settlers
    recruits <- data_sub$recruits
    weight <- data_sub$weight
    
    # Fetch the refined ranges from likdf2 for the given substudy
    prev_alphabeta <- filter(likdf2, substudy_num == substudy)
    refined_alpha_range <- c(prev_alphabeta$alpha_min, prev_alphabeta$alpha_max)
    refined_beta_range <- c(prev_alphabeta$beta_min, prev_alphabeta$beta_max)
    
    # Initialize variables to store minimum negative log likelihood
    min_neg_loglik <- Inf
    alpha_manual <- NA
    beta_manual <- NA
    
    # Loop through alpha and beta values within the refined ranges
    for (alpha in seq(refined_alpha_range[1], refined_alpha_range[2], length.out = 100)) {
      for (beta in seq(refined_beta_range[1], refined_beta_range[2], length.out = 100)) {
        current_neg_loglik <- neg_loglik_safeguarded(c(alpha, beta), t_days, settlers, recruits, weight)
        
        if (current_neg_loglik < min_neg_loglik) {
          min_neg_loglik <- current_neg_loglik
          alpha_manual <- alpha
          beta_manual <- beta
        }
      }
    }
    
    # Compute SEs using second derivative method
    alpha_se <- 1 / sqrt(abs(second_derivative(neg_loglik_safeguarded,
                                               c(alpha_manual, beta_manual), 1,
                                               t_days = t_days, settlers = settlers,
                                               recruits = recruits, weight = weight)))
    beta_se <- 1 / sqrt(abs(second_derivative(neg_loglik_safeguarded,
                                              c(alpha_manual, beta_manual), 2,
                                              t_days = t_days, settlers = settlers,
                                              recruits = recruits, weight = weight)))
    
    # Create scatter plot with model fit
    scatter_plot2 <- ggplot(data_sub, aes(x = settlers, y = recruits)) +
      geom_point(alpha = 0.5) +
      stat_function(fun = bh_f,
                    args = list(alpha = alpha_manual, beta = beta_manual, t_days = mean(data_sub$t)),
                    color = "red", aes(linetype = "Manual")) +
      labs(title = paste("Substudy", substudy),
           x = "Settlers (n0_m2)",
           y = "Recruits (nt_m2)") +
      xlim(0, max(data_sub$settlers)) +
      geom_abline(slope = 1, linetype = "dashed", color = "gray") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      scale_linetype_manual(values = c("Manual" = "dashed"))
    
    # 1D profile for alpha
    alpha_range <- seq(refined_alpha_range[1], refined_alpha_range[2], length.out = 100)
    ll_alpha <- sapply(alpha_range, function(a)
      neg_loglik_safeguarded(c(a, beta_manual), t_days, settlers, recruits, weight))
    df_alpha <- data.frame(alpha = alpha_range, likelihood = ll_alpha)
    min_alpha_nll <- min(ll_alpha)
    
    # Find the alpha value at the minimum NLL
    min_alpha_index <- which.min(ll_alpha)
    alpha_at_min_nll <- alpha_range[min_alpha_index]
    
    # Calculate confidence interval bounds for alpha
    alpha_ci_lower <- if (any(ll_alpha >= min_alpha_nll + 2 & alpha_range < alpha_at_min_nll)) {
      max(alpha_range[ll_alpha >= min_alpha_nll + 2 & alpha_range < alpha_at_min_nll])
    } else {
      NA
    }
    
    alpha_ci_upper <- if (any(ll_alpha >= min_alpha_nll + 2 & alpha_range > alpha_at_min_nll)) {
      min(alpha_range[ll_alpha >= min_alpha_nll + 2 & alpha_range > alpha_at_min_nll])
    } else {
      NA
    }
    
    # Calculate alpha variance using symmetric formula (matches manuscript)
    alpha_variance <- if (!is.na(alpha_ci_lower) & !is.na(alpha_ci_upper)) {
      ((alpha_ci_upper - alpha_ci_lower) / (2 * 1.96))^2
    } else {
      NA
    }
    
    # 1D profile for beta
    beta_range <- seq(refined_beta_range[1], refined_beta_range[2], length.out = 100)
    ll_beta <- sapply(beta_range, function(b)
      neg_loglik_safeguarded(c(alpha_manual, b), t_days, settlers, recruits, weight))
    df_beta <- data.frame(beta = beta_range, likelihood = ll_beta)
    min_beta_nll <- min(ll_beta)
    min_beta_index <- which.min(ll_beta)
    beta_at_min_nll <- beta_range[min_beta_index]
    
    # Calculate confidence interval bounds for beta
    beta_ci_lower <- if (any(ll_beta >= min_beta_nll + 2 & beta_range < beta_at_min_nll)) {
      max(beta_range[ll_beta >= min_beta_nll + 2 & beta_range < beta_at_min_nll])
    } else {
      NA
    }
    
    beta_ci_upper <- if (any(ll_beta >= min_beta_nll + 2 & beta_range > beta_at_min_nll)) {
      min(beta_range[ll_beta >= min_beta_nll + 2 & beta_range > beta_at_min_nll])
    } else {
      NA
    }
    
    # Calculate beta variance using symmetric formula (matches manuscript)
    beta_variance <- if (!is.na(beta_ci_lower) & !is.na(beta_ci_upper)) {
      ((beta_ci_upper - beta_ci_lower) / (2 * 1.96))^2
    } else {
      NA
    }
    
    # Generate alpha profile plot
    p_alpha <- ggplot(df_alpha, aes(x = alpha, y = likelihood)) +
      geom_line(linewidth = 0.8) +
      geom_vline(aes(xintercept = alpha_manual), color = "red", linetype = "dashed", linewidth = 0.8)
    
    if (!is.na(alpha_ci_lower)) {
      p_alpha <- p_alpha +
        geom_vline(aes(xintercept = alpha_ci_lower), color = "blue", linetype = "dotted") +
        annotate("text", x = alpha_ci_lower, y = max(ll_alpha) * 0.8,
                 label = paste("Alpha CI Lower:", round(alpha_ci_lower, 4)), hjust = 1, color = "blue")
    }
    
    if (!is.na(alpha_ci_upper)) {
      p_alpha <- p_alpha +
        geom_vline(aes(xintercept = alpha_ci_upper), color = "blue", linetype = "dotted") +
        annotate("text", x = alpha_ci_upper, y = max(ll_alpha) * 0.9,
                 label = paste("Alpha CI Upper:", round(alpha_ci_upper, 4)), hjust = 1, color = "blue")
    }
    
    p_alpha <- p_alpha +
      annotate("text", x = alpha_manual, y = max(ll_alpha),
               label = paste("Alpha Manual:", round(alpha_manual, 4)), hjust = 1, color = "red") +
      labs(title = paste("1D Negative Log-Likelihood Profile for Alpha - Substudy", substudy),
           x = "Alpha", y = "Negative Log-Likelihood") +
      theme_minimal()
    
    # Generate beta profile plot
    p_beta <- ggplot(df_beta, aes(x = beta, y = likelihood)) +
      geom_line(linewidth = 0.8) +
      geom_vline(aes(xintercept = beta_manual), color = "red", linetype = "dashed", linewidth = 0.8)
    
    if (!is.na(beta_ci_lower)) {
      p_beta <- p_beta +
        geom_vline(aes(xintercept = beta_ci_lower), color = "blue", linetype = "dotted") +
        annotate("text", x = beta_ci_lower, y = max(ll_beta) * 0.8,
                 label = paste("Beta CI Lower:", round(beta_ci_lower, 4)), hjust = 1, color = "blue")
    }
    
    if (!is.na(beta_ci_upper)) {
      p_beta <- p_beta +
        geom_vline(aes(xintercept = beta_ci_upper), color = "blue", linetype = "dotted") +
        annotate("text", x = beta_ci_upper, y = max(ll_beta) * 0.9,
                 label = paste("Beta CI Upper:", round(beta_ci_upper, 4)), hjust = 1, color = "blue")
    }
    
    p_beta <- p_beta +
      annotate("text", x = beta_manual, y = max(ll_beta),
               label = paste("Beta Manual:", round(beta_manual, 4)), hjust = 1, color = "red") +
      labs(title = paste("1D Negative Log-Likelihood Profile for Beta - Substudy", substudy),
           x = "Beta", y = "Negative Log-Likelihood") +
      theme_minimal()
    
    # Compute 2D likelihood profile
    alpha_range <- seq(refined_alpha_range[1], refined_alpha_range[2], length.out = 100)
    beta_range <- seq(refined_beta_range[1], refined_beta_range[2], length.out = 100)
    
    ll_grid <- outer(alpha_range, beta_range, Vectorize(function(a, b)
      neg_loglik_safeguarded(c(a, b), t_days, settlers, recruits, weight)))
    
    # Convert matrix to data frame for plotting
    df_2d <- expand.grid(alpha = alpha_range, beta = beta_range)
    df_2d$likelihood <- as.vector(ll_grid)
    
    # Calculate threshold for 95% CI
    threshold_95 <- min(df_2d$likelihood) + qchisq(0.95, 2) / 2
    
    # Create 2D heatmap
    p_2d_heatmap <- ggplot(df_2d, aes(x = alpha, y = beta)) +
      geom_tile(aes(fill = likelihood)) +
      scale_fill_viridis_c(direction = -1, name = "Neg. Log-Lik") +
      geom_contour(aes(z = likelihood), breaks = threshold_95, color = "white", linewidth = 0.8) +
      geom_point(data = data.frame(alpha = alpha_manual, beta = beta_manual),
                 aes(x = alpha, y = beta), color = "red", size = 3, shape = 4, stroke = 2) +
      labs(title = paste("2D Negative Log-Likelihood Heatmap with 95% CI - Substudy", substudy),
           x = "Alpha", y = "Beta") +
      theme_minimal()
    
    # Save 1D profile data frames to CSV files
    output_dir <- "output/individual_plots/"
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.csv(df_alpha, paste0(output_dir, "alpha_profile_substudy_", substudy, ".csv"), row.names = FALSE)
    write.csv(df_beta, paste0(output_dir, "beta_profile_substudy_", substudy, ".csv"), row.names = FALSE)
    
    # Save minimum NLL values and confidence intervals
    min_nll_values <- data.frame(
      substudy = substudy,
      min_alpha_nll = min_alpha_nll,
      alpha_ci_lower = alpha_ci_lower,
      alpha_ci_upper = alpha_ci_upper,
      alpha_variance = alpha_variance,
      min_beta_nll = min_beta_nll,
      beta_ci_lower = beta_ci_lower,
      beta_ci_upper = beta_ci_upper,
      beta_variance = beta_variance
    )
    write.csv(min_nll_values, paste0(output_dir, "min_nll_values_substudy_", substudy, ".csv"), row.names = FALSE)
    
    # Return plots and results
    return(list(
      scatter_plot = scatter_plot2,
      alpha_manual = alpha_manual,
      beta_manual = beta_manual,
      p_2d = p_2d_heatmap,
      p_alpha = p_alpha,
      p_beta = p_beta,
      alpha_se = alpha_se,
      beta_se = beta_se,
      df_beta = df_beta,
      df_alpha = df_alpha,
      min_alpha_nll = min_alpha_nll,
      min_beta_nll = min_beta_nll,
      alpha_variance = alpha_variance,
      beta_variance = beta_variance,
      alpha_ci_lower = alpha_ci_lower,
      alpha_ci_upper = alpha_ci_upper,
      beta_ci_lower = beta_ci_lower,
      beta_ci_upper = beta_ci_upper
    ))
  }
  
  # ----------------------------------------------------------------------------
  # 5. Test Manual Grid Search Function (Single Substudy Example)
  # ----------------------------------------------------------------------------
  cat("\n=== Testing manual grid search on single substudy ===\n")
  
  # Filter for Binomial approach
  all_data_binomial <- filter(all_data, estimation_approach == "Binomial")
  valid_substudies <- sort(unique(all_data_binomial$substudy_num))
  
  # Run model for a single substudy (example)
  substudy_selected <- valid_substudies[1]
  cat("Testing substudy:", substudy_selected, "\n")
  
  results_df2 <- NULL  # Not used in function
  results <- manual_neg_loglik(substudy = substudy_selected,
                               all_data = all_data_binomial,
                               results_df = results_df2,
                               likdf2 = likdf2)
  
  if (!is.null(results)) {
    grid.arrange(results$scatter_plot, results$p_alpha, results$p_beta, results$p_2d,
                 ncol = 2, nrow = 2)
    cat("Test successful!\n")
  }
  
  # ----------------------------------------------------------------------------
  # 6. Run Manual Grid Search for All Binomial Substudies
  # ----------------------------------------------------------------------------
  cat("\n=== Running manual grid search for all Binomial substudies ===\n")
  
  # Reload and filter data for clean run
  all_data <- read.csv("data/all_studies_looped-2024-09-11.csv")
  names(all_data) <- c("study_num", "substudy_num", "genus_species", "settlers", "recruits",
                       "t", "area", "no_pts", "other_factor", "pdf_num", "notes")
  
  covariates <- read.csv("data/covariates-2024-09-30.csv")
  
  usedf <- covariates %>%
    select(use_2024, substudy_num, estimation_approach)
  
  all_data <- inner_join(all_data, usedf, by = "substudy_num")
  all_data$weight <- all_data$area * all_data$no_pts * all_data$other_factor
  all_data <- filter(all_data, use_2024 == "yes", estimation_approach == "Binomial")
  
  valid_substudies <- sort(unique(all_data$substudy_num))
  
  # Create PDF to save all plots
  pdf("output/substudy_plots_2025-10-15.pdf", width = 11, height = 8.5)
  
  # Initialize results data frame (symmetric variance formula matches manuscript)
  results_df <- data.frame(
    substudy_num = integer(),
    alpha_hat = numeric(),
    beta_hat = numeric(),
    alpha_se = numeric(),
    beta_se = numeric(),
    alpha_ci_lower = numeric(),
    alpha_ci_upper = numeric(),
    alpha_variance = numeric(),
    beta_ci_lower = numeric(),
    beta_ci_upper = numeric(),
    beta_variance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through all valid substudies
  for (substudy in valid_substudies) {
    cat("\n=== Processing substudy:", substudy, "(", which(valid_substudies == substudy),
        "of", length(valid_substudies), ") ===\n")
    results <- manual_neg_loglik(substudy = substudy,
                                 all_data = all_data,
                                 results_df = NULL,
                                 likdf2 = likdf2)
    
    if (!is.null(results)) {
      # Plot combined grid
      grid.arrange(results$scatter_plot, results$p_alpha, results$p_beta, results$p_2d,
                   ncol = 2, nrow = 2)
      
      # Append results to data frame (symmetric variance formula matches manuscript)
      results_df <- rbind(results_df, data.frame(
        substudy_num = substudy,
        alpha_hat = results$alpha_manual,
        beta_hat = results$beta_manual,
        alpha_se = results$alpha_se,
        beta_se = results$beta_se,
        alpha_ci_lower = results$alpha_ci_lower,
        alpha_ci_upper = results$alpha_ci_upper,
        alpha_variance = results$alpha_variance,
        beta_ci_lower = results$beta_ci_lower,
        beta_ci_upper = results$beta_ci_upper,
        beta_variance = results$beta_variance
      ))
    }
  }
  
  # Close PDF device
  dev.off()
  
  # Save results to CSV
  output_file <- "output/substudy_results_2025-10-15.csv"
  write.csv(results_df, output_file, row.names = FALSE)
  cat("\n=== Manual grid search completed ===\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total substudies processed:", nrow(results_df), "\n")
  
  # ----------------------------------------------------------------------------
  # 7. NLS2 Fitting for "Other_No_Counts" Substudies
  # ----------------------------------------------------------------------------
  # NOTE: NLS2 brute-force with 1000×1000 grid is extremely slow (hours per substudy).
  # Set USE_FAST_LOADER = TRUE to skip NLS2 and load existing results instead.
  # If combined_results_2025-10-15.csv exists, this script will automatically load
  # existing results for substudies already present and only run NLS2 for missing ones.
  USE_FAST_LOADER <- FALSE  # Set to TRUE for quick runs, FALSE for proper estimation
  
  # Check if combined_results file exists and load existing results
  existing_results_file <- "output/combined_results_2025-10-15.csv"
  existing_results <- NULL
  existing_substudies <- c()
  
  if (file.exists(existing_results_file)) {
    cat("\n=== Found existing combined_results file ===\n")
    cat("Loading:", existing_results_file, "\n")
    existing_results <- read.csv(existing_results_file)
    existing_substudies <- existing_results$substudy_num
    cat("Found", length(existing_substudies), "substudies with existing results\n")
  }
  
  if (USE_FAST_LOADER) {
    cat("\n=== FAST MODE: Loading existing Other_No_Counts results ===\n")
    
    # Check for covariates to identify Other_No_Counts substudies
    covariates <- read.csv("data/covariates-2024-09-30.csv")
    usedf <- covariates %>%
      filter(use_2024 == "yes", estimation_approach == "Other_No_Counts") %>%
      select(substudy_num, estimation_approach)
    
    cat("Found", nrow(usedf), "substudies with Other_No_Counts approach\n")
    
    # Load existing combined_results file to extract Other_No_Counts estimates
    # First try current file, then fall back to 2024-09-17 file
    if (!is.null(existing_results)) {
      results_to_load <- existing_results
      cat("Using existing results from:", existing_results_file, "\n")
    } else if (file.exists("output/combined_results_2024-09-17.csv")) {
      results_to_load <- read.csv("output/combined_results_2024-09-17.csv")
      cat("Using existing results from: output/combined_results_2024-09-17.csv\n")
    } else {
      stop("No existing results file found. Set USE_FAST_LOADER = FALSE to run NLS2 fitting.")
    }
    
    # Filter for Other_No_Counts substudies
    param_table <- results_to_load %>%
      filter(substudy_num %in% usedf$substudy_num) %>%
      select(substudy_num, alpha, beta, alpha_variance, beta_variance)
    
    names(param_table) <- c("substudy", "alpha", "beta", "alpha_variance", "beta_variance")
    
    cat("Total Other_No_Counts substudies loaded:", nrow(param_table), "\n")
    
  } else {
    cat("\n=== Running NLS2 fitting for Other_No_Counts substudies ===\n")
    if (!is.null(existing_results)) {
      cat("NOTE: Will skip substudies already present in", existing_results_file, "\n")
    }
    cat("WARNING: This uses 1000×1000 brute-force grid and may take hours!\n\n")
    
    # Reload and filter data
    all_data <- read.csv("data/all_studies_looped-2024-09-11.csv")
    names(all_data) <- c("study_num", "substudy_num", "genus_species", "settlers", "recruits",
                         "t", "area", "no_pts", "other_factor", "pdf_num", "notes")
    
    covariates <- read.csv("data/covariates-2024-09-30.csv")
    
    usedf <- covariates %>%
      select(use_2024, substudy_num, estimation_approach)
    
    all_data <- inner_join(all_data, usedf, by = "substudy_num")
    all_data <- filter(all_data, use_2024 == "yes")
    all_data_nc <- filter(all_data, estimation_approach == "Other_No_Counts")
    
    substudy_list <- unique(all_data_nc$substudy_num)
    
    if (length(substudy_list) == 0) {
      cat("No substudies found with 'Other_No_Counts' approach.\n")
    } else {
      cat("Found", length(substudy_list), "substudies with Other_No_Counts approach\n")
    }
    
    # Separate substudies into those with existing results and those needing NLS2
    if (!is.null(existing_results)) {
      substudies_with_results <- substudy_list[substudy_list %in% existing_substudies]
      substudies_needing_nls2 <- substudy_list[!substudy_list %in% existing_substudies]
      
      cat("  - Substudies with existing results:", length(substudies_with_results), "\n")
      cat("  - Substudies needing NLS2 fitting:", length(substudies_needing_nls2), "\n")
      
      # Load existing results for substudies that already have them
      if (length(substudies_with_results) > 0) {
        param_table_existing <- existing_results %>%
          filter(substudy_num %in% substudies_with_results) %>%
          select(substudy_num, alpha, beta, alpha_variance, beta_variance)
        names(param_table_existing) <- c("substudy", "alpha", "beta", "alpha_variance", "beta_variance")
      } else {
        param_table_existing <- data.frame(
          substudy = numeric(),
          alpha = numeric(),
          beta = numeric(),
          alpha_variance = numeric(),
          beta_variance = numeric()
        )
      }
    } else {
      substudies_needing_nls2 <- substudy_list
      cat("  - All", length(substudies_needing_nls2), "substudies need NLS2 fitting\n")
      param_table_existing <- data.frame(
        substudy = numeric(),
        alpha = numeric(),
        beta = numeric(),
        alpha_variance = numeric(),
        beta_variance = numeric()
      )
    }
    
    # Initialize parameter table for new NLS2 fits
    param_table_new <- data.frame(
      substudy = numeric(),
      alpha = numeric(),
      beta = numeric(),
      alpha_variance = numeric(),
      beta_variance = numeric()
    )
    
    # Define Beverton-Holt function for NLS2
    bh_f_safeguarded_nls2 <- function(alpha, beta, t_days, settlers) {
      predicted_recruits <- (settlers * exp(-alpha * t_days)) /
        (1 + (beta * settlers * (1 - exp(-alpha * t_days)) / alpha))
      return(predicted_recruits)
    }
    
    # Create PDF for Other_No_Counts plots (only if there are substudies needing NLS2)
    if (length(substudies_needing_nls2) > 0) {
      pdf("output/substudy_other_nocount_plots.pdf", width = 5, height = 5)
    }
    
    # Loop through each Other_No_Counts substudy that needs NLS2 fitting
    for (substudy in substudies_needing_nls2) {
      cat("\n=== Processing Other_No_Counts substudy:", substudy, "(", which(substudies_needing_nls2 == substudy),
          "of", length(substudies_needing_nls2), ") ===\n")
      
      # Filter data for current substudy
      data_sub <- filter(all_data_nc, substudy_num == substudy)
      t_days <- data_sub$t
      settlers <- data_sub$settlers
      recruits <- data_sub$recruits
      
      # Define grid for brute-force search in nls2
      st1 <- expand.grid(
        alpha = seq(1E-10, 0.1, len = 1000),
        beta = seq(0, 0.05, len = 1000)
      )
      
      # Fit model using nls2
      nls2fit <- nls2(
        recruits ~ (settlers * exp(-alpha * t_days)) /
          (1 + (beta * settlers * (1 - exp(-alpha * t_days)) / alpha)),
        data = data_sub,
        start = st1,
        algorithm = "brute-force",
        control = nls.control(maxiter = 50000, tol = 1e-10)
      )
      
      # Extract coefficients
      exact_params <- coef(nls2fit)
      exact_alpha <- exact_params["alpha"]
      exact_beta <- exact_params["beta"]
      
      # Extract standard errors
      summary_fit <- summary(nls2fit)
      alpha_std_error <- summary_fit$coefficients["alpha", "Std. Error"]
      beta_std_error <- summary_fit$coefficients["beta", "Std. Error"]
      
      alpha_variance <- alpha_std_error^2
      beta_variance <- beta_std_error^2
      
      # Store parameters in table for NEW fits
      param_table_new <- rbind(param_table_new, data.frame(
        substudy = substudy,
        alpha = exact_alpha,
        beta = exact_beta,
        alpha_variance = alpha_variance,
        beta_variance = beta_variance
      ))
      
      # Generate smooth predictions for plotting
      settlers_smooth <- seq(0, max(settlers), length.out = 1000)
      y_hat_smoothnls2 <- bh_f_safeguarded_nls2(exact_alpha, exact_beta, data_sub$t[1], settlers_smooth)
      
      # Plot data with model fit
      p <- ggplot(data = data_sub, aes(x = settlers, y = recruits)) +
        geom_point() +
        geom_line(data = data.frame(settlers_smooth, y_hat_smoothnls2),
                  aes(x = settlers_smooth, y = y_hat_smoothnls2),
                  color = "red", linewidth = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        labs(title = paste("Substudy", substudy, "- Fit with NLS2"),
             x = "Settlers", y = "Recruits") +
        theme_minimal()
      
      print(p)
    }
    
    # Close PDF device (only if it was opened)
    if (length(substudies_needing_nls2) > 0) {
      dev.off()
    }
    
    # Combine existing results with newly fitted results
    param_table <- bind_rows(param_table_existing, param_table_new)
    
    # Display parameter table
    cat("\n=== NLS2 fitting completed ===\n")
    if (!is.null(existing_results)) {
      cat("Total Other_No_Counts substudies:\n")
      cat("  - Loaded from existing results:", nrow(param_table_existing), "\n")
      cat("  - Newly fitted with NLS2:", nrow(param_table_new), "\n")
      cat("  - Combined total:", nrow(param_table), "\n\n")
    } else {
      cat("Total Other_No_Counts substudies processed:", nrow(param_table), "\n\n")
    }
    print(param_table)
    
    # Generate fancy HTML table
    if (nrow(param_table) > 0) {
      fancy_table <- param_table %>%
        kable("html", caption = "Model Parameter Estimates for Other_No_Counts Substudies") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
      
      print(fancy_table)
    }
  }
  
  # ----------------------------------------------------------------------------
  # 8. Combine Results from Both Approaches
  # ----------------------------------------------------------------------------
  cat("\n=== Combining results from both approaches ===\n")
  
  # Use the NEWLY generated Binomial results from Section 6
  # Select relevant columns and rename to match
  binomial_results <- results_df %>%
    select(substudy_num, alpha_hat, beta_hat, alpha_variance, beta_variance)
  
  names(binomial_results) <- c("substudy_num", "alpha", "beta", "alpha_variance", "beta_variance")
  
  # Rename param_table columns to match
  names(param_table) <- c("substudy_num", "alpha", "beta", "alpha_variance", "beta_variance")
  
  # Combine the two tables: param_table (Other_No_Counts) + binomial_results (Binomial)
  combined_results <- bind_rows(param_table, binomial_results)
  
  # Remove row names and sort by substudy_num
  rownames(combined_results) <- NULL
  combined_results_sorted <- combined_results %>%
    arrange(substudy_num)
  
  # Check for duplicates
  duplicates <- combined_results_sorted %>%
    filter(duplicated(substudy_num) | duplicated(substudy_num, fromLast = TRUE))
  
  if (nrow(duplicates) > 0) {
    cat("WARNING: Found duplicate substudy numbers:\n")
    print(duplicates)
  } else {
    cat("No duplicate substudy numbers found.\n")
  }
  
  # Define output path (hardcoded for reproducibility)
  output_path <- "output/combined_results_2025-10-15.csv"
  
  # Write combined results to CSV
  write.csv(combined_results_sorted, file = output_path, row.names = FALSE)
  
  cat("\n================================================================================\n")
  cat("                    SCRIPT COMPLETED SUCCESSFULLY!                              \n")
  cat("================================================================================\n\n")
  cat("Combined results written to:", output_path, "\n")
  cat("  - Binomial substudies:", nrow(binomial_results), "\n")
  cat("  - Other_No_Counts substudies:", nrow(param_table), "\n")
  cat("  - Total substudies:", nrow(combined_results_sorted), "\n\n")
  cat("This file is used by 1_data_phylogeny_loading.R for the main meta-analysis.\n")
  cat("\nAll outputs saved to output/ directory:\n")
  cat("  - substudy_plots_2025-10-15.pdf (132 Binomial substudies, 4 plots each)\n")
  cat("  - substudy_results_2025-10-15.csv (detailed Binomial results with variances)\n")
  cat("  - substudy_other_nocount_plots.pdf (Other_No_Counts fits)\n")
  cat("  - combined_results_2025-10-15.csv (all substudies for meta-analysis)\n")
  cat("  - individual_plots/ (CSV files for each substudy profile)\n\n")
  
  # Close the skip-if-exists conditional block
} else {
  cat("\n================================================================================\n")
  cat("                  SKIPPING: combined_results_2025-10-15.csv already exists      \n")
  cat("================================================================================\n")
  cat("Set SKIP_IF_EXISTS <- FALSE in 0a_individual_study_estimates.R to force regeneration.\n\n")
}
