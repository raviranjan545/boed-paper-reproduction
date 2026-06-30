# Declare present directory
here::i_am("1_scripts/temperature_5_comparing_designs_with_scoring_rules.R")

# Load libraries
library(here)
library(scoringRules)
library(MASS)
library(dplyr)
library(Matrix)
library(conflicted)
conflicts_prefer(dplyr::select())

# Load Prior
source(here("1_scripts",
            "support_file_temperature_prior.R"))

### 1. Set up and generate global scaling factors

# Set up the number of simulations
n_sims <- 1000
n_points_list <- c(5, 7, 10, 15)

# Load the "params" data frame containing true values
params <- read.csv(here("2_designs_and_other_simulation_inputs",
                        "1000_simulated_norberg_curve_parameter_combinations_original_scale.csv")) %>%
    select(c(a, b, tmax, tmin))

# Centring temperature 
centre_constant <- 20.88513513513514 

# Generate reference scales
set.seed(87293)
prior_raw <- norberg_prior(100000) # Draw 20k samples from prior (a, b, tmax, tmin)
prior_raw <- prior_raw[, colnames(prior_raw) != "sig2"]

# A. Raw Scale (a, b, tmax, tmin)
sds_raw <- apply(prior_raw, 2, sd)

# B. Log Scale (log_a, log_b, tmax, tmin)
prior_log <- prior_raw
prior_log[,"a"] <- log(prior_raw[,"a"])
prior_log[,"b"] <- log(prior_raw[,"b"])
colnames(prior_log) <- c("loga", "logb", "tmax", "tmin")
sds_log <- apply(prior_log, 2, sd)

# B. Corrected Scale (loga_corrected, log_b, tmax, tmin)
prior_corr <- prior_log
prior_corr[,"loga"] <- prior_log[,"loga"] + (prior_raw[,"b"] * centre_constant)
colnames(prior_corr) <- c("loga_corr", "logb", "tmax", "tmin")
sds_corr <- apply(prior_corr, 2, sd)

### 2. Scoring function
score_simulation <- function(csv_path, 
                             true_row, 
                             global_sds_raw, 
                             global_sds_log, 
                             global_sds_corr) {
    
    # 1. Read and prep data
    raw_df <- read.csv(csv_path, check.names = FALSE)
    
    # Extract log scale parameters (as estimated by brms)
    post_base <- raw_df %>%
        select(loga = `b_loga_Intercept`,
               logb = `b_logb_Intercept`,
               tmax = `b_tmax_Intercept`,
               tmin = `b_tmin_Intercept`) %>%
        as.matrix()
    
    # Prepare Truth (Raw Scale)
    true_base <- as.numeric(true_row[c("a", "b", "tmax", "tmin")])
    names(true_base) <- c("a", "b", "tmax", "tmin")
    
    # Helper: ES calculation
    calc_es <- function(post_mat, true_vec, global_sds) {
        means <- colMeans(post_mat)
        post_centered <- sweep(post_mat, 2, means, "-")
        post_scaled   <- sweep(post_centered, 2, global_sds, "/")
        true_scaled   <- (true_vec - means) / global_sds
        return(es_sample(y = true_scaled, dat = t(post_scaled)))
    }
    
    # Helper: CRPS calculation
    calc_crps <- function(post_col, true_val, name_prefix) {
        val <- crps_sample(y = true_val, dat = post_col)
        names(val) <- name_prefix
        return(val)
    }
    
    # 1. RAW SCALE (a, b, tmax, tmin)
    # Transform Posterior to Raw
    post_raw <- post_base
    post_raw[, "loga"] <- exp(post_base[, "loga"]) 
    post_raw[, "logb"] <- exp(post_base[, "logb"]) 
    colnames(post_raw) <- c("a", "b", "tmax", "tmin")
    
    # Scores
    es_raw <- calc_es(post_raw, true_base, global_sds_raw)
    names(es_raw) <- "ES_raw"
    
    # Calculate ALL CRPS here (since this is the natural scale)
    c_a    <- calc_crps(post_raw[,"a"],    true_base["a"],    "CRPS_a")
    c_b    <- calc_crps(post_raw[,"b"],    true_base["b"],    "CRPS_b")
    c_tmax <- calc_crps(post_raw[,"tmax"], true_base["tmax"], "CRPS_tmax")
    c_tmin <- calc_crps(post_raw[,"tmin"], true_base["tmin"], "CRPS_tmin")
    
    scores_1 <- c(es_raw, c_a, c_b, c_tmax, c_tmin)
    
    # 2. LOG SCALE (loga, logb, tmax, tmin)
    # Posterior is already loga, logb
    post_log <- post_base
    colnames(post_log) <- c("loga", "logb", "tmax", "tmin")
    
    # Transform Truth
    true_log <- true_base
    true_log["a"] <- log(true_base["a"])
    true_log["b"] <- log(true_base["b"])
    names(true_log) <- c("loga", "logb", "tmax", "tmin")
    
    # Scores
    es_log <- calc_es(post_log, true_log, global_sds_log)
    names(es_log) <- "ES_log"
    
    # Only calculate CRPS for the CHANGED variables (loga, logb)
    c_loga <- calc_crps(post_log[,"loga"], true_log["loga"], "CRPS_loga")
    c_logb <- calc_crps(post_log[,"logb"], true_log["logb"], "CRPS_logb")
    
    scores_2 <- c(es_log, c_loga, c_logb)
    
    # 3. CORRECTED SCALE (loga_corr, logb, tmax, tmin)
    # Transform Posterior: loga_corr = loga + exp(logb) * C
    post_corr <- post_log
    post_corr[, "loga"] <- post_base[, "loga"] + (exp(post_base[, "logb"]) * centre_constant)
    colnames(post_corr) <- c("loga_corr", "logb", "tmax", "tmin")
    
    # Transform Truth
    true_corr <- true_log
    true_corr["loga"] <- true_log["loga"] + (true_base["b"] * centre_constant)
    names(true_corr) <- c("loga_corr", "logb", "tmax", "tmin")
    
    # Scores
    es_corr <- calc_es(post_corr, true_corr, global_sds_corr)
    names(es_corr) <- "ES_corr"
    
    # Only calculate CRPS for the CHANGED variable (loga_corr)
    c_loga_corr <- calc_crps(post_corr[,"loga_corr"], true_corr["loga_corr"], "CRPS_loga_corr")
    
    scores_3 <- c(es_corr, c_loga_corr)
    
    return(c(scores_1, scores_2, scores_3))
}

### 3. Main loop

# Store results per sample size
results_by_n <- vector("list", length(n_points_list))
names(results_by_n) <- as.character(n_points_list)
base_dir <- here("temperature_curves_posteriors") 


for (i in 1:length(n_points_list)) {
    
    # Select number of points 
    n_points <- n_points_list[i]
    
    # Define lists of saved model files
    fit_opt_list <- list.files(
        base_dir, 
        pattern = paste0("^norberg_posterior_opt_", n_points, ".*\\.csv$"), 
        full.names = TRUE
    )
    fit_unif_list <- list.files(
        base_dir, 
        pattern = paste0("^norberg_posterior_unif_", n_points, ".*\\.csv$"), 
        full.names = TRUE
    )
    
    results_opt <- vector("list", n_sims)
    results_unif <- vector("list", n_sims)
    
    for (j in 1:n_sims) {
        
        print(c(i, j))
        true_row <- params[j, ]
        
        # Process Optimal
        if (j <= length(fit_opt_list)) {
            results_opt[[j]] <- score_simulation(fit_opt_list[[j]], 
                                                 true_row, sds_raw, sds_log, sds_corr)
        }
        
        # Process Uniform
        if (j <= length(fit_unif_list)) {
            results_unif[[j]] <- score_simulation(fit_unif_list[[j]], 
                                                  true_row, sds_raw, sds_log, sds_corr)
        }
    }
    
    # Bind and Label
    df_opt <- as.data.frame(do.call(rbind, results_opt))
    df_unif <- as.data.frame(do.call(rbind, results_unif))

    df_opt$design  <- "opt"
    df_unif$design <- "unif"
    df_opt$simulation <- seq_len(n_sims)
    df_unif$simulation <- seq_len(n_sims)
    
    df_this <- rbind(df_opt, df_unif)
    df_this$n_points <- n_points
    
    results_by_n[[as.character(n_points)]] <- df_this
}

# Combine everything
all_scores <- do.call(rbind, results_by_n)

# Save results
write.csv(all_scores, 
          here("3_simulation_result_summaries",
               "norberg_fit_opt_unif_scoring_rules_comparison.csv"), 
          row.names = FALSE)
