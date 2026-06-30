# Declare present directory
here::i_am("1_scripts/temperature_6_evaluating_prediction_error.R")

# Load libraries
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(roahd)

# Temperature function (Norberg)
norberg <- function(temp, a, b, tmax, tmin){
    gr <- (a*exp(b*temp))*(tmax - temp)*(temp - tmin)
    gr
}

# Load the params" data frame containing true values
params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "1000_simulated_norberg_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(simulation = row_number())

# Load experimental designs
opt_temp5_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "temperature_SIG_optimal_design_5_points.rds")),
                         sort(final.d[[besti]]))
opt_temp7_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "temperature_SIG_optimal_design_7_points.rds")),
                         sort(final.d[[besti]]))
opt_temp10_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                       "temperature_SIG_optimal_design_10_points.rds")),
                          sort(final.d[[besti]]))
opt_temp15_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                       "temperature_SIG_optimal_design_15_points.rds")),
                          sort(final.d[[besti]]))


# Number of simulations and design points
n_sims <- 1000
n_points_list <- c(5, 7, 10, 15)

# Dense evaluation grid
env_grid <- seq(5, 35, by = 0.5)

# Function to identify points not in a design
not_in_design <- function(x, design_points, tol = 1e-8) {
    sapply(x, function(xx) all(abs(xx - design_points) > tol))
}

# Function to calculate prediction errors
calc_prediction_error <- function(a_true, b_true, tmax_true, tmin_true,
                                  a_pred, b_pred, tmax_pred, tmin_pred,
                                  eval_grid) {
    
    true_curve <- norberg(eval_grid, a = a_true, b = b_true,
                          tmax = tmax_true, tmin = tmin_true)
    pred_curve <- norberg(eval_grid, a = a_pred, b = b_pred,
                          tmax = tmax_pred, tmin = tmin_pred)
    
    data.frame(
        mae = mean(abs(pred_curve - true_curve)),
        rmse = sqrt(mean((pred_curve - true_curve)^2))
    )
}

# Function to extract parameter estimates for a central curve from a posterior
# 
# NOTE on parameterisation: brms estimates a and b on the LOG scale
# (b_loga_Intercept, b_logb_Intercept), with tmax, tmin on the natural scale.
get_estimates <- function(csv_path,
                          eval_grid,
                          cols = c("b_loga_Intercept",
                                   "b_logb_Intercept",
                                   "b_tmax_Intercept",
                                   "b_tmin_Intercept")) {
    
    draws_df <- read.csv(csv_path, check.names = FALSE)
    post_log <- as.matrix(draws_df[, cols])
    colnames(post_log) <- c("loga", "logb", "tmax", "tmin")
    
    # Natural-scale draws (a, b exponentiated; tmax, tmin unchanged)
    post_nat <- post_log
    post_nat[, "loga"] <- exp(post_log[, "loga"])
    post_nat[, "logb"] <- exp(post_log[, "logb"])
    colnames(post_nat) <- c("a", "b", "tmax", "tmin")
    
    # Optional subsampling (uncomment if MBD becomes too slow on full posterior)
    # n_sub <- 2000
    # if (nrow(post_log) > n_sub) {
    #   set.seed(123)
    #   sub_idx <- sample.int(nrow(post_log), n_sub)
    #   post_log <- post_log[sub_idx, ]
    #   post_nat <- post_nat[sub_idx, ]
    # }
    
    n_draws <- nrow(post_nat)
    
    # Build the matrix of posterior curves over eval_grid
    # Each row is f(eval_grid; theta^(s)). Used for the MBD calculation
    curve_matrix <- matrix(NA, nrow = n_draws, ncol = length(eval_grid))
    for (s in seq_len(n_draws)) {
        curve_matrix[s, ] <- norberg(
            temp = eval_grid,
            a    = post_nat[s, "a"],
            b    = post_nat[s, "b"],
            tmax = post_nat[s, "tmax"],
            tmin = post_nat[s, "tmin"]
        )
    }
    
    # Identify central curve using Modified Band Depth (functional median by containment)
    f_data  <- fData(grid = eval_grid, values = curve_matrix)
    depths  <- MBD(f_data)
    mbd_idx <- which.max(depths)
    
    a_mbd    <- post_nat[mbd_idx, "a"]
    b_mbd    <- post_nat[mbd_idx, "b"]
    tmax_mbd <- post_nat[mbd_idx, "tmax"]
    tmin_mbd <- post_nat[mbd_idx, "tmin"]
    
    data.frame(
        a_mbd    = a_mbd,
        b_mbd    = b_mbd,
        tmax_mbd = tmax_mbd,
        tmin_mbd = tmin_mbd
    )
}

# Helper: compute MAE and RMSE for one (a, b, tmax, tmin) set on both grids
errors_on_both_grids <- function(a_true, b_true, tmax_true, tmin_true,
                                 a_pred, b_pred, tmax_pred, tmin_pred,
                                 grid_full, grid_offdesign,
                                 suffix) {
    
    err_full <- calc_prediction_error(a_true, b_true, tmax_true, tmin_true,
                                      a_pred, b_pred, tmax_pred, tmin_pred,
                                      grid_full)
    err_off  <- calc_prediction_error(a_true, b_true, tmax_true, tmin_true,
                                      a_pred, b_pred, tmax_pred, tmin_pred,
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

base_dir <- here("temperature_curves_posteriors")

# Loop through sample sizes
for (i in seq_along(n_points_list)) {
    
    n_points <- n_points_list[i]
    print(n_points)
    
    # Uniform design
    unif_levels <- round(seq(5, 35, length.out = n_points) * 2) / 2
    
    # Optimal design
    opt_levels <- switch(
        as.character(n_points),
        "5"  = opt_temp5_levels,
        "7"  = opt_temp7_levels,
        "10" = opt_temp10_levels,
        "15" = opt_temp15_levels
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
        pattern = paste0("^norberg_posterior_opt_", n_points, ".*\\.csv$"),
        full.names = TRUE
    )
    fit_unif_list <- list.files(
        base_dir,
        pattern = paste0("^norberg_posterior_unif_", n_points, ".*\\.csv$"),
        full.names = TRUE
    )
    
    # Storage for both designs at this sample size
    results_opt  <- vector("list", n_sims)
    results_unif <- vector("list", n_sims)
    
    # Loop over the 1000 fits per design
    for (j in seq_len(n_sims)) {
        
        print(c(i, j))
        
        a_true    <- params$a[j]
        b_true    <- params$b[j]
        tmax_true <- params$tmax[j]
        tmin_true <- params$tmin[j]
        
        # Extract the four parameter sets (and SDs) from each posterior
        ests_opt <- get_estimates(
            csv_path  = fit_opt_list[[j]],
            eval_grid = env_grid
        )
        ests_unif <- get_estimates(
            csv_path  = fit_unif_list[[j]],
            eval_grid = env_grid
        )
        
        # Errors for the four methods x two grids (per design)
        errs_opt <- errors_on_both_grids(
            a_true, b_true, tmax_true, tmin_true,
            ests_opt$a_mbd, ests_opt$b_mbd,
            ests_opt$tmax_mbd, ests_opt$tmin_mbd,
            grid_full, grid_offdesign, "mbd")
        
        errs_unif <- errors_on_both_grids(
            a_true, b_true, tmax_true, tmin_true,
            ests_unif$a_mbd, ests_unif$b_mbd,
            ests_unif$tmax_mbd, ests_unif$tmin_mbd,
            grid_full, grid_offdesign, "mbd")
        
        results_opt[[j]] <- data.frame(
            simulation = j,
            design     = "opt",
            n_points   = n_points,
            a_true     = a_true,
            b_true     = b_true,
            tmax_true  = tmax_true,
            tmin_true  = tmin_true,
            ests_opt,
            errs_opt
        )
        results_unif[[j]] <- data.frame(
            simulation = j,
            design     = "unif",
            n_points   = n_points,
            a_true     = a_true,
            b_true     = b_true,
            tmax_true  = tmax_true,
            tmin_true  = tmin_true,
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
        "norberg_fit_opt_unif_prediction_error_comparison.csv"),
    row.names = FALSE
)

# Summary table
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
calc_prediction_error_pointwise <- function(a_true, b_true, tmax_true, tmin_true,
                                            a_pred, b_pred, tmax_pred, tmin_pred,
                                            eval_grid) {
    
    true_curve <- norberg(eval_grid, a = a_true, b = b_true,
                          tmax = tmax_true, tmin = tmin_true)
    pred_curve <- norberg(eval_grid, a = a_pred, b = b_pred,
                          tmax = tmax_pred, tmin = tmin_pred)
    
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
                a_true    = a_true,
                b_true    = b_true,
                tmax_true = tmax_true,
                tmin_true = tmin_true,
                a_pred    = a_mbd,
                b_pred    = b_mbd,
                tmax_pred = tmax_mbd,
                tmin_pred = tmin_mbd,
                eval_grid = env_grid
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
        "norberg_fit_opt_unif_pointwise_prediction_error_mbd.csv"
    ),
    row.names = FALSE
)
