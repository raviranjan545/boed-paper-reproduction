# Declare present directory
here::i_am("1_scripts/temperature_3_simulation_and_bayesian_fitting.R")

# Load libraries
library(here)
library(brms)
library(dplyr)
library(Matrix)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set number of simulations
n_sims <- 1000

# Save model fits to this subfolder
temperature_folder <- "temperature_curves_simulations_and_fits"

# Temperature function (norberg)
norberg <- function(temp, a, b, tmax, tmin){
    gr <- (a*exp(b*temp))*(tmax - temp)*(temp - tmin)
    gr
}

# Read in true parameter values
params <- read.csv(here(
    "2_designs_and_other_simulation_inputs", 
    "1000_simulated_norberg_curve_parameter_combinations_original_scale.csv")
)

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

# Define single T0 value
T0_universal <- mean(
    c(
        ((round(c(seq(5, 35, length.out = 5)) * 2)) / 2),
        ((round(c(seq(5, 35, length.out = 7)) * 2)) / 2),
        ((round(c(seq(5, 35, length.out = 10)) * 2)) / 2),
        ((round(c(seq(5, 35, length.out = 15)) * 2)) / 2),
        opt_temp5_levels, 
        opt_temp7_levels, 
        opt_temp10_levels, 
        opt_temp15_levels
    )
)


# Define centred equation
norberg_formula_centered <- bf(
    growth_rate ~
        # a * exp(b*T0) * exp(b*temp_c)
        (exp(loga) * exp(exp(logb) * T0) * exp(exp(logb) * temp_c)) *
        # (tmax - T) and (T - tmin) with T = temp_c + T0
        ((tmax - T0) - temp_c) * (temp_c - (tmin - T0)),
    loga + logb + tmax + tmin ~ 1, nl = TRUE
)

## Define prior with correlations
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

# Define the multivariate normal prior block for (loga, logb, tmax, tmin) 
stan_code <- "
{
  vector[4] theta;
  theta[1] = b_loga[1];
  theta[2] = b_logb[1];
  theta[3] = b_tmax[1];
  theta[4] = b_tmin[1];
  target += multi_normal_cholesky_lpdf(theta | mu_theta, L_theta);
}
"
stanvars <- stanvar(mu_vec, name = "mu_theta") +
    stanvar(L, name = "L_theta") +
    stanvar(scode = stan_code, block = "model")

# Specify sigma prior separately
prior_sigma <- prior(lognormal(-2.3, 0.1), class = "sigma")

# Define function for generating initial guesses
inits_fun <- function() list(
    loga = meanlog_a + rnorm(1, 0, 0.01),
    logb = meanlog_b + rnorm(1, 0, 0.02),
    tmax = mean_tmax + rnorm(1, 0, 0.2),
    tmin = mean_tmin + rnorm(1, 0, 0.3)
)


# Simulate a dataset for compiling model
growthdat_opt <- data.frame(temperature = c(10, 15, 20, 25, 30))

growthdat_opt <- growthdat_opt %>%
    mutate(growth_rate = norberg(temperature, 
                                 a = params$a[1],
                                 b = params$b[1],
                                 tmax = params$tmax[1],
                                 tmin = params$tmin[1]) +
               rnorm(growthdat_opt$temperature, 0, params$sig[1]),
           T0 = T0_universal,
           temp_c = temperature - T0)

# Compile the model once outside the loop (will throw a warning)
compiled_model <- brm(
    formula = norberg_formula_centered,
    data    = growthdat_opt,
    prior   = prior_sigma,
    stanvars = stanvars,
    chains  = 1, 
    iter = 1
)

# Remove data no longer needed
rm(growthdat_opt)

# Define number of points to simulate
n_points_list <- c(5, 7, 10, 15)


# Loop through all sample sizes, then loop through 1000 parameter combinations and fit
for (i in 1:length(n_points_list)){
    
    # Select number of data points
    n_points <- n_points_list[i]
    
    # Define uniform design
    growthdat_unif <- (round(data.frame(temperature = seq(5, 35, length.out = n_points)) * 2)) / 2
    
    # Define optimal design
    if (n_points == 5){
        growthdat_opt <- data.frame(temperature = opt_temp5_levels)
    } else if (n_points == 7){
        growthdat_opt <- data.frame(temperature = opt_temp7_levels)
    } else if (n_points == 10){
        growthdat_opt <- data.frame(temperature = opt_temp10_levels)
    } else if (n_points == 15){
        growthdat_opt <- data.frame(temperature = opt_temp15_levels)
    }
    
    # Centring temperature for better fitting
    growthdat_unif$T0 <- T0_universal
    growthdat_unif$temp_c <- growthdat_unif$temperature - growthdat_unif$T0
    
    growthdat_opt$T0 <- T0_universal
    growthdat_opt$temp_c <- growthdat_opt$temperature - growthdat_opt$T0
    
    # Loop through all parameter combinations, simulate data, fit model, and save fitted model
    for (j in 1:n_sims){
        
        print(c(i, j))
        
        # Optimal design - generate dataset with random noise
        set.seed(4091825 * n_points + j)
        growthdat_opt$growth_rate <- norberg(
            growthdat_opt$temperature,
            a = params$a[j],
            b = params$b[j],
            tmax = params$tmax[j],
            tmin = params$tmin[j]
        ) +
            rnorm(n_points, 0, params$sig[j])
        
        # Uniform design - generate dataset with random noise
        set.seed(4091825 * n_points + j)
        growthdat_unif$growth_rate <- norberg(
            growthdat_unif$temperature,
            a = params$a[j],
            b = params$b[j],
            tmax = params$tmax[j],
            tmin = params$tmin[j]
        ) +
            rnorm(n_points, 0, params$sig[j])
        
        # Define the filename format
        opt_base <- sprintf("norberg_fit_opt_%spts_%04d", n_points, j)
        unif_base <- sprintf("norberg_fit_unif_%spts_%04d", n_points, j)
        
        # Optimal design - fit Norberg function
        fit_opt <- update(
            compiled_model,
            newdata = growthdat_opt,
            prior   = prior_sigma,
            stanvars = stanvars,
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 107582,   
            file = here(temperature_folder, paste0(opt_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 20
            )
        )
        
        # Uniform design - fit Norberg function
        fit_unif <- update(
            compiled_model,
            newdata = growthdat_unif,
            prior   = prior_sigma,
            stanvars = stanvars,
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 107582, 
            file = here(temperature_folder, paste0(unif_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 20
            )
        )
    }
}



