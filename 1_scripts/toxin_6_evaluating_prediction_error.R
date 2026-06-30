# Declare present directory
here::i_am("1_scripts/toxin_6_evaluating_prediction_error.R")

# Load libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(roahd)

# Toxin function from Ritz (2010)
loglogis <- function(toxin, mu_max, ec50, slope){
    gr <- mu_max / (1 + ((toxin/ec50)^slope))
    gr
}

# Load the "params" data frame containing true values
params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "1000_simulated_toxin_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(simulation = row_number())

# Load experimental designs
opt_tox5_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                     "toxin_SIG_optimal_design_5_points.rds")),
                        sort(final.d[[besti]]))
opt_tox7_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                     "toxin_SIG_optimal_design_7_points.rds")),
                        sort(final.d[[besti]]))
opt_tox10_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "toxin_SIG_optimal_design_10_points.rds")),
                         sort(final.d[[besti]]))
opt_tox15_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "toxin_SIG_optimal_design_15_points.rds")),
                         sort(final.d[[besti]]))

# Number of simulations and design points
n_sims <- 1000
n_points_list <- c(5, 7, 10, 15)

# Design space (log-spaced concentrations plus zero) and dense evaluation grid.
env_grid <- c(0, 10^seq(log10(0.1), log10(1000), length.out = 101))

# Function to identify points not in a design
not_in_design <- function(x, design_points, tol = 1e-8) {
    sapply(x, function(xx) all(abs(xx - design_points) > tol))
}

# Function to calculate prediction errors
calc_prediction_error <- function(mu_true, ec50_true, slope_true,
                                  mu_pred, ec50_pred, slope_pred,
                                  eval_grid) {
    
    true_curve <- loglogis(eval_grid, mu_max = mu_true,
                           ec50 = ec50_true, slope = slope_true)
    pred_curve <- loglogis(eval_grid, mu_max = mu_pred,
                           ec50 = ec50_pred, slope = slope_pred)
    
    data.frame(
        mae = mean(abs(pred_curve - true_curve)),
        rmse = sqrt(mean((pred_curve - true_curve)^2))
    )
}

# Function to extract parameter estimates for a central curve from a posterior
get_estimates <- function(csv_path,
                          eval_grid,
                          cols = c("b_mumax_Intercept",
                                   "b_ec50_Intercept",
                                   "b_slope_Intercept")) {
    
    draws_df <- read.csv(csv_path, check.names = FALSE)
    post_samples <- as.matrix(draws_df[, cols])
    colnames(post_samples) <- c("mu_max", "ec50", "slope")
    
    # Optional subsampling (uncomment if MBD becomes too slow on full posterior)
    # n_sub <- 2000
    # if (nrow(post_samples) > n_sub) {
    #   set.seed(123)
    #   sub_idx <- sample.int(nrow(post_samples), n_sub)
    #   post_samples <- post_samples[sub_idx, ]
    # }
    
    n_draws <- nrow(post_samples)
    
    # Build the matrix of posterior curves over eval_grid
    # Each row is f(eval_grid; theta^(s)). Used for the MBD calculation
    curve_matrix <- matrix(NA, nrow = n_draws, ncol = length(eval_grid))
    for (s in seq_len(n_draws)) {
        curve_matrix[s, ] <- loglogis(
            toxin  = eval_grid,
            mu_max = post_samples[s, "mu_max"],
            ec50   = post_samples[s, "ec50"],
            slope  = post_samples[s, "slope"]
        )
    }
    
    # Identify central curve using Modified Band Depth (functional median by containment)
    f_data  <- fData(grid = eval_grid, values = curve_matrix)
    depths  <- MBD(f_data)
    mbd_idx <- which.max(depths)
    
    mu_mbd    <- post_samples[mbd_idx, "mu_max"]
    ec50_mbd  <- post_samples[mbd_idx, "ec50"]
    slope_mbd <- post_samples[mbd_idx, "slope"]
    
    data.frame(
        mu_max_mbd = mu_mbd,
        ec50_mbd   = ec50_mbd,
        slope_mbd  = slope_mbd
    )
}

# Helper: compute MAE and RMSE for one (mu, ec50, slope) set on both grids
errors_on_both_grids <- function(mu_true, ec50_true, slope_true,
                                 mu_pred, ec50_pred, slope_pred,
                                 grid_full, grid_offdesign,
                                 suffix) {
    
    err_full <- calc_prediction_error(mu_true, ec50_true, slope_true,
                                      mu_pred, ec50_pred, slope_pred,
                                      grid_full)
    err_off  <- calc_prediction_error(mu_true, ec50_true, slope_true,
                                      mu_pred, ec50_pred, slope_pred,
                                      grid_offdesign)
    
    out <- data.frame(
        err_full$mae,
        err_full$rmse,
        err_off$mae,
        err_off$rmse
    )
    names(out) <- paste0(
        c("mae_full_", "rmse_full_", "mae_offdesign_", "rmse_offdesign_"),
        suffix
    )
    out
}

# Storage
results_by_n <- vector("list", length(n_points_list))
names(results_by_n) <- as.character(n_points_list)

base_dir <- here("toxin_curves_posteriors")

# Loop through sample sizes
for (i in seq_along(n_points_list)) {
    
    n_points <- n_points_list[i]
    print(n_points)
    
    # Uniform design (evenly spaced indices in env_grid, matching the
    # simulation script's construction)
    unif_design_indices <- round(seq(1, 102, length.out = n_points), 0)
    unif_levels <- env_grid[unif_design_indices]
    
    # Optimal design
    opt_levels <- switch(
        as.character(n_points),
        "5"  = opt_tox5_levels,
        "7"  = opt_tox7_levels,
        "10" = opt_tox10_levels,
        "15" = opt_tox15_levels
    )
    
    # Off-design grid excludes union of both designs
    combined_design_levels <- sort(unique(c(unif_levels, opt_levels)))
    grid_full <- env_grid
    grid_offdesign <- env_grid[
        not_in_design(env_grid, combined_design_levels, tol = 1e-8)
    ]
    
    # Posterior CSV file paths
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
    
    # Storage for both designs at this sample size
    results_opt  <- vector("list", n_sims)
    results_unif <- vector("list", n_sims)
    
    # Loop over the 1000 fits per design
    for (j in seq_len(n_sims)) {
        
        print(c(i, j))
        
        mu_true    <- params$mu_max[j]
        ec50_true  <- params$ec50[j]
        slope_true <- params$slope[j]
        
        # Extract the four parameter sets (and SDs) from each posterior
        ests_opt <- get_estimates(
            csv_path  = fit_opt_list[[j]],
            eval_grid = env_grid
        )
        ests_unif <- get_estimates(
            csv_path  = fit_unif_list[[j]],
            eval_grid = env_grid
        )
        
        # Errors for the two grids (per design)
        errs_opt <- errors_on_both_grids(mu_true, ec50_true, slope_true,
                                         ests_opt$mu_max_mbd, ests_opt$ec50_mbd,
                                         ests_opt$slope_mbd,
                                         grid_full, grid_offdesign, "mbd")
        errs_unif <- errors_on_both_grids(mu_true, ec50_true, slope_true,
                                          ests_unif$mu_max_mbd, ests_unif$ec50_mbd,
                                          ests_unif$slope_mbd,
                                          grid_full, grid_offdesign, "mbd")
        
        results_opt[[j]] <- data.frame(
            simulation = j,
            design     = "opt",
            n_points   = n_points,
            mu_max_true = mu_true,
            ec50_true   = ec50_true,
            slope_true  = slope_true,
            ests_opt,
            errs_opt
        )
        results_unif[[j]] <- data.frame(
            simulation = j,
            design     = "unif",
            n_points   = n_points,
            mu_max_true = mu_true,
            ec50_true   = ec50_true,
            slope_true  = slope_true,
            ests_unif,
            errs_unif
        )
    }
    
    df_this <- rbind(
        do.call(rbind, results_opt),
        do.call(rbind, results_unif)
    )
    
    results_by_n[[as.character(n_points)]] <- df_this
}

# Combine across sample sizes
all_prediction_errors <- bind_rows(results_by_n)

write.csv(
    all_prediction_errors,
    here(
        "3_simulation_result_summaries",
        "toxin_fit_opt_unif_prediction_error_comparison.csv"),
    row.names = FALSE
)

# Summary table for checking
summary_by_n_and_design <- all_prediction_errors %>%
    group_by(n_points, design) %>%
    summarise(
        mean_mae_full_mbd       = mean(mae_full_mbd, na.rm = TRUE),
        mean_rmse_full_mbd      = mean(rmse_full_mbd, na.rm = TRUE),
        mean_mae_offdesign_mbd  = mean(mae_offdesign_mbd, na.rm = TRUE),
        mean_rmse_offdesign_mbd = mean(rmse_offdesign_mbd, na.rm = TRUE),
        .groups = "drop"
    )

print(summary_by_n_and_design)

# Pointwise prediction error across the design space, using the MBD curve
# For every n_points x design x grid point, the absolute and squared error of
# the MBD curve relative to the true curve, computed per simulation and then
# summarised across simulations into a pointwise MAE / RMSE. 

# Function to calculate pointwise prediction errors
calc_prediction_error_pointwise <- function(mu_true, ec50_true, slope_true,
                                            mu_pred, ec50_pred, slope_pred,
                                            eval_grid) {
    
    true_curve <- loglogis(eval_grid, mu_max = mu_true,
                           ec50 = ec50_true, slope = slope_true)
    pred_curve <- loglogis(eval_grid, mu_max = mu_pred,
                           ec50 = ec50_pred, slope = slope_pred)
    
    data.frame(
        env = eval_grid,
        abs_error = abs(pred_curve - true_curve),
        sq_error = (pred_curve - true_curve)^2
    )
}

# Calculate pointwise MBD errors rowwise over all simulations, using the
# already-computed MBD parameters and true values
all_pointwise_prediction_errors <- all_prediction_errors %>%
    rowwise() %>%
    mutate(
        pointwise_err = list(
            calc_prediction_error_pointwise(
                mu_true    = mu_max_true,
                ec50_true  = ec50_true,
                slope_true = slope_true,
                mu_pred    = mu_max_mbd,
                ec50_pred  = ec50_mbd,
                slope_pred = slope_mbd,
                eval_grid  = env_grid
            )
        )
    ) %>%
    ungroup() %>%
    select(simulation, design, n_points, pointwise_err) %>%
    tidyr::unnest(pointwise_err)

write.csv(
    all_pointwise_prediction_errors,
    here(
        "3_simulation_result_summaries",
        "toxin_fit_opt_unif_pointwise_prediction_error_mbd.csv"),
    row.names = FALSE
)
