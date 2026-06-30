# Declare present directory
here::i_am("1_scripts/temperature_2_generating_parameters_for_simulations.R")

# Load libraries
library(here)
library(roahd)
library(Matrix)
library(MASS)

# Set number of simulations
n_sims <- 1000

# Temperature function (Norberg)
norberg <- function(temp, a, b, tmax, tmin){
    gr <- (a*exp(b*temp))*(tmax - temp)*(temp - tmin)
    gr
}

## Generate true parameter values - with correlations

# Target marginals on ORIGINAL scale
mean_a <- 0.002; sd_a <- 0.0002
mean_b <- 0.025; sd_b <- 0.005
mean_tmax <- 35;  sd_tmax <- 2
mean_tmin <- 2;   sd_tmin <- 3

# Convert to lognormal params for a,b (sdlog first, then meanlog)
sdlog_a   <- sqrt(log(1 + (sd_a/mean_a)^2))
meanlog_a <- log(mean_a) - 0.5 * sdlog_a^2
sdlog_b   <- sqrt(log(1 + (sd_b/mean_b)^2))
meanlog_b <- log(mean_b) - 0.5 * sdlog_b^2

# Correlations on (log a, log b, tmax, tmin)
rho_ab     <- -0.6
rho_a_tmax <- -0.5
rho_a_tmin <-  0.3

C <- diag(4)
dimnames(C) <- list(c("loga","logb","tmax","tmin"), c("loga","logb","tmax","tmin"))
C["loga","logb"] <- C["logb","loga"] <- rho_ab
C["loga","tmax"] <- C["tmax","loga"] <- rho_a_tmax
C["loga","tmin"] <- C["tmin","loga"] <- rho_a_tmin
C <- as.matrix(nearPD(C, corr = TRUE)$mat) 

# MVN on (log a, log b, tmax, tmin)
mu_vec <- c(meanlog_a, meanlog_b, mean_tmax, mean_tmin)
sd_vec <- c(sdlog_a,   sdlog_b,   sd_tmax,  sd_tmin)
Sigma  <- diag(sd_vec) %*% C %*% diag(sd_vec)
L      <- t(chol(Sigma))

# Simulate and save parameters
set.seed(93785)
draws <- mvrnorm(n_sims, mu = mu_vec, Sigma = Sigma)

set.seed(20350297)
params <- data.frame(
    a    = exp(draws[,1]), # convert back from logged to original scale
    b    = exp(draws[,2]), # convert back from logged to original scale
    tmax = draws[,3],
    tmin = draws[,4],
    sig = rlnorm(n = n_sims, meanlog = -2.3, sdlog = 0.1)
)

# Save parameter values
write.csv(params, 
          here(
              "2_designs_and_other_simulation_inputs", 
              "1000_simulated_norberg_curve_parameter_combinations_original_scale.csv"),
          row.names = FALSE)



## Identify the central curve by Modified Band Depth (functional median by containment)

# Set driver values 
driver_values <- seq(5, 35, 0.5)

post_mat <- as.matrix(params[, c("a", "b", "tmax", "tmin")])
curve_matrix <- matrix(NA, nrow = 1000, ncol = length(driver_values))

for (s in 1:1000) {
    curve_matrix[s, ] <- norberg(
        temp = driver_values,
        a    = post_mat[s, "a"],
        b    = post_mat[s, "b"],
        tmax = post_mat[s, "tmax"],
        tmin = post_mat[s, "tmin"]
    )
}

f_data  <- fData(grid = driver_values, values = curve_matrix)
depths  <- MBD(f_data)
mbd_idx <- which.max(depths)

mbd_curve_params <- data.frame(
    a = post_mat[mbd_idx, "a"],
    b = post_mat[mbd_idx, "b"],
    tmax = post_mat[mbd_idx, "tmax"],
    tmin = post_mat[mbd_idx, "tmin"]
)

# Save central curve parameter values
write.csv(mbd_curve_params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "norberg_central_curve_mbd_parameters.csv"),
          row.names = FALSE)

