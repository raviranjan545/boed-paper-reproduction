# Declare present directory
here::i_am("1_scripts/toxin_3_simulation_and_bayesian_fitting.R")

# Load libraries
library(here)
library(brms)
library(dplyr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set number of simulations
n_sims <- 1000

# Save model fits to this subfolder
toxin_folder <- "toxin_curves_simulations_and_fits"

# Toxin function from Ritz (2010)
loglogis <- function(toxin, mu_max, ec50, slope){
    gr <- mu_max / (1 + ((toxin/ec50)^slope))
    gr
}


# Set seed for reproducible generation of parameter values
set.seed(47157891)

# Read in true parameter values
params <- read.csv(here(
    "2_designs_and_other_simulation_inputs",
    "1000_simulated_toxin_curve_parameter_combinations_original_scale.csv")
)

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

# Define rescaled equation
toxin_formula_rescaled <- bf(
    growth_rate ~ ((mumax*exp(0.01)) / (1 + (toxin/(ec50*exp(4)))^(slope*exp(1)))),
    mumax + ec50 + slope ~ 1, nl = TRUE
)

# Define priors for rescaled equation
prior_toxin_rescaled <- 
    prior(lognormal(0, 0.4), nlpar = "mumax", lb = 0) +
    prior(lognormal(0, 1), nlpar = "ec50", lb = 0) +
    prior(lognormal(0, 0.5), nlpar = "slope", lb = 0) +
    prior(lognormal(-2.3, 0.1), class = "sigma")

# Simulate a dataset for compiling model
growthdat_opt <- data.frame(toxin = c(0, 1, 10, 100, 1000))
growthdat_opt$growth_rate <- loglogis(growthdat_opt$toxin, 
                                      mu_max = params$mu_max[1],
                                      ec50 = params$ec50[1],
                                      slope = params$slope[1]
) + 
    rnorm(length(growthdat_opt$toxin), 0, params$sig[1])

# Define function for generating initial guesses - RStan backend
inits_fun <- function() list(
    mumax = 1 + rnorm(1, 0, 0.05),
    ec50 = 1 + rnorm(1, 0, 0.05),
    slope = 1 + rnorm(1, 0, 0.05)
)

# Compile the model once outside the loop (will throw a warning)
compiled_model <- brm(
    formula = toxin_formula_rescaled,
    data = growthdat_opt,  
    prior = prior_toxin_rescaled, 
    iter = 1, 
    chains = 1
)

# Remove data no longer needed
rm(growthdat_opt)

# Define design space
design_space <- c(0, 10^seq(log10(0.1), log10(1000), length.out = 101)) 

# Define number of points to simulate
n_points_list <- c(5, 7, 10, 15)

# Loop through all sample sizes, then loop through 1000 parameter combinations and fit
for (i in 1:length(n_points_list)){
    
    # Select number of data points
    n_points <- n_points_list[i]
    
    # Define uniform design differently for toxin function
    # Levels are uniformly spaced on a log scale, after which a zero level is added
    unif_design_indices <- round(seq(1, 102, length.out = n_points), 0)
    
    # Define uniform design
    growthdat_unif <- data.frame(toxin = design_space[unif_design_indices])
    
    # Define optimal design
    if (n_points == 5){
        growthdat_opt <- data.frame(toxin = opt_tox5_levels)
    } else if (n_points == 7){
        growthdat_opt <- data.frame(toxin = opt_tox7_levels)
    } else if (n_points == 10){
        growthdat_opt <- data.frame(toxin = opt_tox10_levels)
    } else if (n_points == 15){
        growthdat_opt <- data.frame(toxin = opt_tox15_levels)
    }
    
    # Loop through all parameter combinations, simulate data, fit model, and save fitted model
    for (j in 1:n_sims){
        
        print(c(i, j))
        
        # Optimal design - generate dataset with random noise
        set.seed(60194 * n_points + j)
        growthdat_opt$growth_rate <- loglogis(growthdat_opt$toxin,
                                              mu_max = params$mu_max[j],
                                              ec50 = params$ec50[j],
                                              slope = params$slope[j]) +
            rnorm(n_points, 0, params$sig[j])

        # Uniform design - generate dataset with random noise
        set.seed(60194 * n_points + j)
        growthdat_unif$growth_rate <- loglogis(growthdat_unif$toxin,
                                               mu_max = params$mu_max[j],
                                               ec50 = params$ec50[j],
                                               slope = params$slope[j]) +
            rnorm(n_points, 0, params$sig[j])

        # Define the filename format
        opt_base <- sprintf("toxin_fit_opt_%spts_%04d", n_points, j)
        unif_base <- sprintf("toxin_fit_unif_%spts_%04d", n_points, j)
        
        # Optimal design - fit log-logistic function
        fit_opt <- update(
            compiled_model,
            newdata = growthdat_opt,
            prior = prior_toxin_rescaled,
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 9814,
            file = here(toxin_folder, paste0(opt_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 15
            )
        )
        
        # Uniform design - fit log-logistic function
        fit_unif <- update(
            compiled_model,
            newdata = growthdat_unif,
            prior = prior_toxin_rescaled,
            init = inits_fun,
            init_r = 0, 
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 9814,
            file = here(toxin_folder, paste0(unif_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 15
            )
        )
    }
}
