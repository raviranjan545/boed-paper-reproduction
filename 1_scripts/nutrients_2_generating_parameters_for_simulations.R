# Declare present directory
here::i_am("1_scripts/nutrients_2_generating_parameters_for_simulations.R")

# Load libraries
library(here)
library(roahd)

# Set number of simulations
n_sims <- 1000

# Monod function 
monod_halfsat <- function(resource, mu_max, k){
    gr <- (mu_max * resource) / (resource + k)
    gr
}

## Generate true parameter values for simulations 

# Set seed for reproducible generation of parameter values
set.seed(198725)

# Generate true parameter values
params <- data.frame(mu_max = rlnorm(n = n_sims, meanlog = 0.01, sdlog = 0.4),
                     k = rlnorm(n = n_sims, meanlog = 0.3, sdlog = 0.5),
                     sig = rlnorm(n = n_sims, meanlog = -2.3, sdlog = 0.1))

# Save parameter values
write.csv(params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "1000_simulated_monod_curve_parameter_combinations_original_scale.csv"),
          row.names = FALSE)


## Identify the central curve by Modified Band Depth (functional median by containment)

# Set driver values 
driver_values <- seq(0, 25, 0.01)

# Identify central curve by Modified Band Depth (functional median by containment)
post_mat <- as.matrix(params[, c("mu_max", "k")])
curve_matrix <- matrix(NA, nrow = 1000, ncol = length(driver_values))

for (s in 1:1000) {
    curve_matrix[s, ] <- monod_halfsat(
        resource = driver_values,
        mu_max   = post_mat[s, "mu_max"],
        k        = post_mat[s, "k"]
    )
}

f_data  <- fData(grid = driver_values, values = curve_matrix)
depths  <- MBD(f_data)
mbd_idx <- which.max(depths)

mbd_curve_params <- data.frame(
    mu_max = post_mat[mbd_idx, "mu_max"],
    k = post_mat[mbd_idx, "k"]
)

# Save central curve parameter values
write.csv(mbd_curve_params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "monod_central_curve_mbd_parameters.csv"),
          row.names = FALSE)

