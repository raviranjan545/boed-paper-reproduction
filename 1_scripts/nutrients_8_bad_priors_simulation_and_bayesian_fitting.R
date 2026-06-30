# Declare present directory
here::i_am("1_scripts/nutrients_8_bad_priors_simulation_and_bayesian_fitting.R")

# Load libraries
library(here)
library(brms)
library(dplyr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set number of simulations
n_sims <- 1000

# Save model fits to this subfolder
nutrient_bad_priors_folder <- "nutrient_curves_bad_priors_simulations_and_fits"

# Monod function 
monod_halfsat <- function(resource, mu_max, k){
    gr <- (mu_max * resource) / (resource + k)
    gr
}

# Set seed for reproducible generation of parameter values
set.seed(10795)

# Generate NEW true parameter values for bad prior test
params <- data.frame(mu_max = rlnorm(n = n_sims, meanlog = 0.01 + log(0.7), sdlog = 0.4),
                     k = rlnorm(n = n_sims, meanlog = 0.3 + log(0.7), sdlog = 0.5),
                     sig = rlnorm(n = n_sims, meanlog = -2.3, sdlog = 0.1))

# Save parameter values
write.csv(params, 
          here(
              "2_designs_and_other_simulation_inputs",
              "bad_priors_1000_simulated_monod_curve_parameter_combinations_original_scale.csv"),
          row.names = FALSE)

# Load experimental designs
opt_nut5_levels <- with(readRDS(here(
    "2_designs_and_other_simulation_inputs",
    "nutrients_SIG_optimal_design_5_points.rds")),
    sort(final.d[[besti]]))
opt_nut7_levels <- with(readRDS(here(
    "2_designs_and_other_simulation_inputs",
    "nutrients_SIG_optimal_design_7_points.rds")),
    sort(final.d[[besti]]))
opt_nut10_levels <- with(readRDS(here(
    "2_designs_and_other_simulation_inputs",
    "nutrients_SIG_optimal_design_10_points.rds")),
    sort(final.d[[besti]]))
opt_nut15_levels <- with(readRDS(here(
    "2_designs_and_other_simulation_inputs",
    "nutrients_SIG_optimal_design_15_points.rds")),
    sort(final.d[[besti]]))

# Define rescaled equation
monod_formula_rescaled <- bf(
    growth_rate ~ ((mumax*exp(0.01)) * resource / ((k*exp(0.3)) + resource)),
    mumax + k ~ 1, nl = TRUE
)

# Define bad priors for rescaled equation
prior_monod_rescaled <- 
    prior(lognormal(0, 0.4), nlpar = "mumax", lb = 0) +
    prior(lognormal(0, 0.5), nlpar = "k", lb = 0) +
    prior(lognormal(-2.3, 0.1), class = "sigma")

# Simulate a dataset for compiling model
growthdat_opt <- data.frame(resource = c(1, 3, 5, 10, 20))
growthdat_opt$growth_rate <- monod_halfsat(growthdat_opt$resource, 
                                           mu_max = params$mu_max[1],
                                           k = params$k[1]) + 
    rnorm(length(growthdat_opt$resource), 0, params$sig[1])

# Define function for generating initial guesses - RStan backend
inits_fun <- function() list(
    mumax = 1 + rnorm(1, 0, 0.05),
    k = 1 + rnorm(1, 0, 0.05)
)

# Compile the model once outside the loop (will throw a warning)
compiled_model <- brm(
    formula = monod_formula_rescaled,
    data = growthdat_opt,  
    prior = prior_monod_rescaled, 
    iter = 1, 
    chains = 1,
    cores = 1,
    sample_prior = "no",
    control = list(
        adapt_delta = 0.95,
        max_treedepth = 15
    )
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
    growthdat_unif <- round(data.frame(resource = seq(0.2, 25, length.out = n_points)), 1)
    
    # Define optimal design
    if (n_points == 5){
        growthdat_opt <- data.frame(resource = opt_nut5_levels)
    } else if (n_points == 7){
        growthdat_opt <- data.frame(resource = opt_nut7_levels)
    } else if (n_points == 10){
        growthdat_opt <- data.frame(resource = opt_nut10_levels)
    } else if (n_points == 15){
        growthdat_opt <- data.frame(resource = opt_nut15_levels)
    }
    
    # Loop through all parameter combinations, simulate data, fit model, and save fitted model
    for (j in 1:n_sims){
        
        print(c(i, j))
        
        # Optimal design - generate dataset with random noise
        set.seed(2935 * n_points + j)
        growthdat_opt$growth_rate <- monod_halfsat(growthdat_opt$resource, 
                                                   mu_max = params$mu_max[j],
                                                   k = params$k[j]) + 
            rnorm(n_points, 0, params$sig[j])
        
        # Uniform design - generate dataset with random noise
        set.seed(2935 * n_points + j)
        growthdat_unif$growth_rate <- monod_halfsat(growthdat_unif$resource, 
                                                    mu_max = params$mu_max[j],
                                                    k = params$k[j]) + 
            rnorm(n_points, 0, params$sig[j])
        
        # Define the filename format
        opt_base <- sprintf("monod_badprior_fit_opt_%spts_%04d", n_points, j)
        unif_base <- sprintf("monod_badprior_fit_unif_%spts_%04d", n_points, j)
        
        # Optimal design - fit Monod function
        fit_opt <- update(
            compiled_model,
            newdata = growthdat_opt,
            prior = prior_monod_rescaled, 
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4, 
            cores = 4, 
            seed = 8715,
            file = here(nutrient_bad_priors_folder, paste0(opt_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.95,
                max_treedepth = 15
            )
        )
        
        # Uniform design - fit Monod function
        fit_unif <- update(
            compiled_model,
            newdata = growthdat_unif, 
            prior = prior_monod_rescaled,
            init = inits_fun,
            init_r = 0, 
            iter = 5000,
            chains = 4, 
            cores = 4, 
            seed = 8715,
            file = here(nutrient_bad_priors_folder, paste0(unif_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.95,
                max_treedepth = 15
            )
        )
    }
}
