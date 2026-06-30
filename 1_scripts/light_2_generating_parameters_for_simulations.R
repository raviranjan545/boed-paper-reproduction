# Declare present directory
setwd(this.path::this.dir())   # move into 1_scripts, where this file lives
here::i_am("1_scripts/light_2_generating_parameters_for_simulations.R")

# Load libraries
library(here)
library(roahd)

# Set number of simulations
n_sims <- 1000

# Eilers-Peeters function 
eilers_peeters <- function(light, mu_max, alpha, i_opt){
    gr <- ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) +
                                   ((1 - ((2 * mu_max)/(alpha * i_opt))) * light) +
                                   (mu_max / alpha)))
    gr
}


## Generate true parameter values for simulations 

# Set seed for reproducible generation of parameter values
set.seed(128751)

# Generate parameter values
params <- data.frame(
    mu_max = rlnorm(n = n_sims, meanlog = 0.01, sdlog = 0.4),
    alpha = rlnorm(n = n_sims, meanlog = -3, sdlog = 0.8),
    i_opt = rlnorm(n = n_sims, meanlog = 5.5, sdlog = 0.3),
    sig = rlnorm(n = n_sims, meanlog = -2.3, sdlog = 0.1)
)

# Save parameter values
write.csv(params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "1000_simulated_eilerspeeters_curve_parameter_combinations_original_scale.csv"),
          row.names = FALSE)


## Identify the central curve by Modified Band Depth (functional median by containment)

# Set driver values 
driver_values <- seq(0, 1000, 1)

post_mat <- as.matrix(params[, c("mu_max", "alpha", "i_opt")])
curve_matrix <- matrix(NA, nrow = 1000, ncol = length(driver_values))

for (s in 1:1000) {
    curve_matrix[s, ] <- eilers_peeters(
        light  = driver_values,
        mu_max = post_mat[s, "mu_max"],
        alpha  = post_mat[s, "alpha"],
        i_opt  = post_mat[s, "i_opt"]
    )
}

f_data  <- fData(grid = driver_values, values = curve_matrix)
depths  <- MBD(f_data)
mbd_idx <- which.max(depths)

mbd_curve_params <- data.frame(
    mu_max = post_mat[mbd_idx, "mu_max"],
    alpha = post_mat[mbd_idx, "alpha"],
    i_opt = post_mat[mbd_idx, "i_opt"]
)

# Save central curve parameter values
write.csv(mbd_curve_params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "eilerspeeters_central_curve_mbd_parameters.csv"),
          row.names = FALSE)

