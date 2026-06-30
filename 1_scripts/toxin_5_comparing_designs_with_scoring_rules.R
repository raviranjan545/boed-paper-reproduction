# Declare present directory
here::i_am("1_scripts/toxin_5_comparing_designs_with_scoring_rules.R")

# Load libraries
library(here)
library(dplyr)
library(scoringRules) 

# Set up the number of simulations
n_sims <- 1000

# Load the "params" data frame containing true values
params <- read.csv(here("2_designs_and_other_simulation_inputs", 
                        "1000_simulated_toxin_curve_parameter_combinations_original_scale.csv")) %>%
    select(c(mu_max, ec50, slope))

# Define number of points
n_points_list <- c(5, 7, 10, 15)

# Calculate SDs of priors for standardisation
sd_lognormal <- function(meanlog, sdlog) {
    sqrt((exp(sdlog^2) - 1) * exp(2*meanlog + sdlog^2))
    # equivalently:
    # exp(meanlog + 0.5*sdlog^2) * sqrt(exp(sdlog^2) - 1)
}

global_sd_mumax <- sd_lognormal(meanlog = 0.01, sdlog = 0.4)
global_sd_ec50  <- sd_lognormal(meanlog = 4, sdlog = 1)
global_sd_slope  <- sd_lognormal(meanlog = 1, sdlog = 0.5)


# Define function to score a single simulation (one RDS path + one row of truth)
score_one_sim <- function(csv_path,
                          true_vals,   
                          cols = c("b_mumax_Intercept", "b_ec50_Intercept", "b_slope_Intercept")) {
    
    draws_df <- read.csv(csv_path, check.names = FALSE)
    
    # Extract the posterior columns
    post_samples <- as.matrix(draws_df[, cols])
    colnames(post_samples) <- c("mu_max", "ec50", "slope")
    
    means <- colMeans(post_samples)
    post_centered <- sweep(post_samples, 2, means, "-")
    
    # 3. Scale (GLOBAL FIXED SDs)
    fixed_sds <- c(global_sd_mumax, global_sd_ec50, global_sd_slope)
    
    post_scaled <- sweep(post_centered, 2, fixed_sds, "/")
    
    # 4. Scale true values by global SD)
    true_scaled <- (true_vals - means) / fixed_sds
    
    es_joint <- es_sample(y = true_scaled, dat = t(post_scaled))
    crps_mu  <- crps_sample(y = true_vals[1], dat = post_samples[, "mu_max"])
    crps_ec50   <- crps_sample(y = true_vals[2], dat = post_samples[, "ec50"])
    crps_slope   <- crps_sample(y = true_vals[3], dat = post_samples[, "slope"])
    
    c(ES = es_joint, CRPS_mumax = crps_mu, CRPS_ec50 = crps_ec50, CRPS_slope = crps_slope)
}

# Store results per sample size
results_by_n <- vector("list", length(n_points_list))
names(results_by_n) <- as.character(n_points_list)

base_dir <- here("toxin_curves_posteriors")



# Loop through all saved model files and summarise
for (i in 1:length(n_points_list)){
    
    # Select number of points 
    n_points <- n_points_list[i]
    
    # Define lists of saved model files
    fit_opt_list <- list.files(
        base_dir, 
        pattern = paste0("^toxin_posterior_opt_", n_points, ".*\\.csv$"), 
        full.names = TRUE
    )
    fit_unif_list <- list.files(
        base_dir, 
        pattern = paste0("^toxin_posterior_unif_", n_points, ".*\\.csv$"), 
        full.names = TRUE
    )
    
    # Initialize lists to store the results for both designs
    results_opt <- vector("list", n_sims)
    results_unif <- vector("list", n_sims)
    
    # Loop over the 1000 model fits for each design and summarise model fits per simulation
    for (j in 1:n_sims) {
        
        print(c(i, j))
        true_vals <- c(params$mu_max[j], params$ec50[j], params$slope[j])
        
        # Calculate the scores
        results_opt[[j]] <- score_one_sim(
            csv_path = fit_opt_list[[j]],
            true_vals = true_vals
        )
        results_unif[[j]] <- score_one_sim(
            csv_path = fit_unif_list[[j]],
            true_vals = true_vals
        )
    }
    
    # Convert to data frames for this n_points
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
               "toxin_fit_opt_unif_scoring_rules_comparison.csv"), 
          row.names = FALSE
)

