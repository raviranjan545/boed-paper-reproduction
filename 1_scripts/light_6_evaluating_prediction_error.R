# Declare present directory
here::i_am("1_scripts/light_6_evaluating_prediction_error.R")

# Load libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(roahd)

# Eilers-Peeters function
eilers_peeters <- function(light, mu_max, alpha, i_opt){
    gr <- ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) +
                                   ((1 - ((2 * mu_max)/(alpha * i_opt))) * light) +
                                   (mu_max / alpha)))
    gr
}

# Load the "params" data frame containing true values
params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "1000_simulated_eilerspeeters_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(simulation = row_number())

# Load experimental designs
opt_light5_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs",
                                       "light_SIG_optimal_design_5_points.rds")),
                          sort(final.d[[besti]]))
opt_light7_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs",
                                       "light_SIG_optimal_design_7_points.rds")),
                          sort(final.d[[besti]]))
opt_light10_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs",
                                        "light_SIG_optimal_design_10_points.rds")),
                           sort(final.d[[besti]]))
opt_light15_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs",
                                        "light_SIG_optimal_design_15_points.rds")),
                           sort(final.d[[besti]]))

# Number of simulations and design points
n_sims <- 1000
n_points_list <- c(5, 7, 10, 15)

# Dense evaluation grid
env_grid <- seq(10, 1000, by = 10)

# Function to identify points not in a design
not_in_design <- function(x, design_points, tol = 1e-8) {
    sapply(x, function(xx) all(abs(xx - design_points) > tol))
}

# Function to calculate prediction errors
calc_prediction_error <- function(mu_true, alpha_true, iopt_true,
                                  mu_pred, alpha_pred, iopt_pred,
                                  eval_grid) {
    
    true_curve <- eilers_peeters(eval_grid, mu_max = mu_true,
                                 alpha = alpha_true, i_opt = iopt_true)
    pred_curve <- eilers_peeters(eval_grid, mu_max = mu_pred,
                                 alpha = alpha_pred, i_opt = iopt_pred)
    
    data.frame(
        mae = mean(abs(pred_curve - true_curve)),
        rmse = sqrt(mean((pred_curve - true_curve)^2))
    )
}

# Function to extract parameter estimates for a central curve from a posterior
get_estimates <- function(csv_path,
                          eval_grid,
                          cols = c("b_mumax_Intercept",
                                   "b_alpha_Intercept",
                                   "b_iopt_Intercept")) {
    
    draws_df <- read.csv(csv_path, check.names = FALSE)
    post_samples <- as.matrix(draws_df[, cols])
    colnames(post_samples) <- c("mu_max", "alpha", "i_opt")
    
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
        curve_matrix[s, ] <- eilers_peeters(
            light  = eval_grid,
            mu_max = post_samples[s, "mu_max"],
            alpha  = post_samples[s, "alpha"],
            i_opt  = post_samples[s, "i_opt"]
        )
    }
    
    # Identify central curve using Modified Band Depth (functional median by containment)
    f_data  <- fData(grid = eval_grid, values = curve_matrix)
    depths  <- MBD(f_data)
    mbd_idx <- which.max(depths)
    
    mu_mbd    <- post_samples[mbd_idx, "mu_max"]
    alpha_mbd <- post_samples[mbd_idx, "alpha"]
    iopt_mbd  <- post_samples[mbd_idx, "i_opt"]
    
    data.frame(
        mu_max_mbd  = mu_mbd,
        alpha_mbd   = alpha_mbd,
        iopt_mbd    = iopt_mbd
    )
}

# Helper: compute MAE and RMSE for one (mu, alpha, iopt) triple on both grids
errors_on_both_grids <- function(mu_true, alpha_true, iopt_true,
                                 mu_pred, alpha_pred, iopt_pred,
                                 grid_full, grid_offdesign,
                                 suffix) {
    
    err_full <- calc_prediction_error(mu_true, alpha_true, iopt_true,
                                      mu_pred, alpha_pred, iopt_pred,
                                      grid_full)
    err_off  <- calc_prediction_error(mu_true, alpha_true, iopt_true,
                                      mu_pred, alpha_pred, iopt_pred,
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

base_dir <- here("light_curves_posteriors")

# Loop through sample sizes
for (i in seq_along(n_points_list)) {
    
    n_points <- n_points_list[i]
    print(n_points)
    
    # Uniform design
    unif_levels <- round(seq(10, 1000, length.out = n_points), -1)
    
    # Optimal design
    opt_levels <- switch(
        as.character(n_points),
        "5"  = opt_light5_levels,
        "7"  = opt_light7_levels,
        "10" = opt_light10_levels,
        "15" = opt_light15_levels
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
        pattern = paste0("^eilerspeeters_posterior_opt_", n_points, ".*\\.csv$"),
        full.names = TRUE
    )
    fit_unif_list <- list.files(
        base_dir,
        pattern = paste0("^eilerspeeters_posterior_unif_", n_points, ".*\\.csv$"),
        full.names = TRUE
    )
    
    # Storage for both designs at this sample size
    results_opt  <- vector("list", n_sims)
    results_unif <- vector("list", n_sims)
    
    # Loop over the 1000 fits per design
    for (j in seq_len(n_sims)) {
        
        print(c(i, j))
        
        mu_true    <- params$mu_max[j]
        alpha_true <- params$alpha[j]
        iopt_true  <- params$i_opt[j]
        
        # Extract the represntative parameters from each posterior
        ests_opt <- get_estimates(
            csv_path  = fit_opt_list[[j]],
            eval_grid = env_grid
        )
        ests_unif <- get_estimates(
            csv_path  = fit_unif_list[[j]],
            eval_grid = env_grid
        )
        
        # Errors for the two grids (per design)
        errs_opt <- errors_on_both_grids(
            mu_true, alpha_true, iopt_true,
            ests_opt$mu_max_mbd,
            ests_opt$alpha_mbd,
            ests_opt$iopt_mbd,
            grid_full, grid_offdesign, "mbd")
        
        errs_unif <- errors_on_both_grids(
            mu_true, alpha_true, iopt_true,
            ests_unif$mu_max_mbd,
            ests_unif$alpha_mbd,
            ests_unif$iopt_mbd,
            grid_full, grid_offdesign, "mbd")
        
        results_opt[[j]] <- data.frame(
            simulation  = j,
            design      = "opt",
            n_points    = n_points,
            mu_max_true = mu_true,
            alpha_true  = alpha_true,
            iopt_true   = iopt_true,
            ests_opt,
            errs_opt
        )
        results_unif[[j]] <- data.frame(
            simulation  = j,
            design      = "unif",
            n_points    = n_points,
            mu_max_true = mu_true,
            alpha_true  = alpha_true,
            iopt_true   = iopt_true,
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
        "eilerspeeters_fit_opt_unif_prediction_error_comparison.csv"),
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
calc_prediction_error_pointwise <- function(mu_true, alpha_true, iopt_true,
                                            mu_pred, alpha_pred, iopt_pred,
                                            eval_grid) {
    
    true_curve <- eilers_peeters(eval_grid, mu_max = mu_true,
                                 alpha = alpha_true, i_opt = iopt_true)
    pred_curve <- eilers_peeters(eval_grid, mu_max = mu_pred,
                                 alpha = alpha_pred, i_opt = iopt_pred)
    
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
                alpha_true = alpha_true,
                iopt_true  = iopt_true,
                mu_pred    = mu_max_mbd,
                alpha_pred = alpha_mbd,
                iopt_pred  = iopt_mbd,
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
        "eilerspeeters_fit_opt_unif_pointwise_prediction_error_mbd.csv"),
    row.names = FALSE
)
