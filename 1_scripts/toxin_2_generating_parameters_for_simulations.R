# Declare present directory
here::i_am("1_scripts/toxin_2_generating_parameters_for_simulations.R")

# Load libraries
library(here)
library(roahd)

# Set number of simulations
n_sims <- 1000

# Toxin function from Ritz (2010)
loglogis <- function(toxin, mu_max, ec50, slope){
    gr <- mu_max / (1 + ((toxin/ec50)^slope))
    gr
}

## Generate true parameter values for simulations 

# Set seed for reproducible generation of parameter values
set.seed(47157891)

# Generate true parameter values
params <- data.frame(
    mu_max = rlnorm(n = n_sims, meanlog = 0.01, sdlog = 0.4),
    ec50 = rlnorm(n = n_sims, meanlog = 4, sdlog = 1),
    slope = rlnorm(n = n_sims, meanlog = 1, sdlog = 0.5),
    sig = rlnorm(n = n_sims, meanlog = -2.3, sdlog = 0.1)
)

# Save parameter values
write.csv(params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "1000_simulated_toxin_curve_parameter_combinations_original_scale.csv"),
          row.names = FALSE)



## Identify the central curve by Modified Band Depth (functional median by containment)

# Set driver values 
driver_values <- c(0, 10^seq(log10(0.1), log10(1000), length.out = 101)) 

# Choose central curve by Modified Band Depth (functional median by containment)
post_mat <- as.matrix(params[, c("mu_max", "ec50", "slope")])
curve_matrix <- matrix(NA, nrow = 1000, ncol = length(driver_values))

for (s in 1:1000) {
    curve_matrix[s, ] <- loglogis(
        toxin  = driver_values,
        mu_max = post_mat[s, "mu_max"],
        ec50   = post_mat[s, "ec50"],
        slope  = post_mat[s, "slope"]
    )
}

f_data  <- fData(grid = driver_values, values = curve_matrix)
depths  <- MBD(f_data)
mbd_idx <- which.max(depths)

mbd_curve_params <- data.frame(
    mu_max = post_mat[mbd_idx, "mu_max"],
    ec50 = post_mat[mbd_idx, "ec50"],
    slope = post_mat[mbd_idx, "slope"]
)


# Save central curve parameter values
write.csv(mbd_curve_params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "toxin_central_curve_mbd_parameters.csv"),
          row.names = FALSE)

