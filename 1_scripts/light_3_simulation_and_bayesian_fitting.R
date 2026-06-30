# Declare present directory
here::i_am("1_scripts/light_3_simulation_and_bayesian_fitting.R")

# Load libraries
library(here)
library(brms)
library(dplyr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set number of simulations
n_sims <- 1000

# Save model fits to this subfolder
light_folder <- "light_curves_simulations_and_fits"

# Eilers-Peeters function 
eilers_peeters <- function(light, mu_max, alpha, i_opt){
    gr <- ((mu_max * light) / (((mu_max / (alpha * i_opt^2)) * (light^2)) +
                                   ((1 - ((2 * mu_max)/(alpha * i_opt))) * light) +
                                   (mu_max / alpha)))
    gr
}

# Read in true parameter values
params <- read.csv(here(
    "2_designs_and_other_simulation_inputs",
    "1000_simulated_eilerspeeters_curve_parameter_combinations_original_scale.csv")
)

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

# Define rescaled equation
eilers_peeters_formula_rescaled <- bf(
    growth_rate ~ (((mumax*exp(0.01)) * light) / 
                       (((mumax*exp(0.01) / ((alpha*exp(-3)) * (iopt*exp(5.5))^2)) * (light^2)) +
                            ((1 - ((2 * mumax*exp(0.01))/((alpha*exp(-3)) * (iopt*exp(5.5))))) * light) +
                            (mumax*exp(0.01) / (alpha*exp(-3))))),
    mumax + alpha + iopt ~ 1, nl = TRUE
)

# Define priors for rescaled equation
prior_eilers_peeters_rescaled <- 
    prior(lognormal(0, 0.4), nlpar = "mumax", lb = 0) +
    prior(lognormal(0, 0.8), nlpar = "alpha", lb = 0) +
    prior(lognormal(0, 0.3), nlpar = "iopt", lb = 0) +
    prior(lognormal(-2.3, 0.1), class = "sigma")

# Simulate a dataset for compiling model
growthdat_opt <- data.frame(light = c(10, 50, 100, 200, 1000))
growthdat_opt$growth_rate <- eilers_peeters(growthdat_opt$light, 
                                            mu_max = params$mu_max[1],
                                            alpha = params$alpha[1],
                                            i_opt = params$i_opt[1]) + 
    rnorm(length(growthdat_opt$light), 0, params$sig[1])

# Define function for generating initial guesses
inits_fun <- function() list(
    mumax = 1 + rnorm(1, 0, 0.05),
    alpha = 1 + rnorm(1, 0, 0.05),
    iopt = 1 + rnorm(1, 0, 0.05)
)

# Compile the model once outside the loop (will throw a warning)
compiled_model <- brm(
    formula = eilers_peeters_formula_rescaled,
    data = growthdat_opt,  
    prior = prior_eilers_peeters_rescaled, 
    iter = 1, 
    chains = 1
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
    growthdat_unif <- round(data.frame(light = seq(10, 1000, length.out = n_points)), -1)
    
    if (n_points == 5){
        growthdat_opt <- data.frame(light = opt_light5_levels)
    } else if (n_points == 7){
        growthdat_opt <- data.frame(light = opt_light7_levels)
    } else if (n_points == 10){
        growthdat_opt <- data.frame(light = opt_light10_levels)
    } else if (n_points == 15){
        growthdat_opt <- data.frame(light = opt_light15_levels)
    }
    
    # Loop through all parameter combinations, simulate data, fit model, and save fitted model
    for (j in 1:n_sims){
        
        print(c(i, j))
        
        # Optimal design - generate dataset with random noise
        set.seed(12915 * n_points + j)
        growthdat_opt$growth_rate <- eilers_peeters(growthdat_opt$light,
                                                    mu_max = params$mu_max[j],
                                                    alpha = params$alpha[j],
                                                    i_opt = params$i_opt[j]) +
            rnorm(n_points, 0, params$sig[j])
        
        # Uniform design - generate dataset with random noise
        set.seed(12915 * n_points + j)
        growthdat_unif$growth_rate <- eilers_peeters(growthdat_unif$light,
                                                     mu_max = params$mu_max[j],
                                                     alpha = params$alpha[j],
                                                     i_opt = params$i_opt[j]) +
            rnorm(n_points, 0, params$sig[j])
        
        # Define the filename format
        opt_base <- sprintf("eilerspeeters_fit_opt_%spts_%04d", n_points, j)
        unif_base <- sprintf("eilerspeeters_fit_unif_%spts_%04d", n_points, j)
        
        # Optimal design - fit Eilers-Peeters function
        fit_opt <- update(
            compiled_model,
            newdata = growthdat_opt,
            prior = prior_eilers_peeters_rescaled,
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 1926,
            file = here(light_folder, paste0(opt_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 20
            )
        )
        
        # Uniform design - fit Eilers-Peeters function
        fit_unif <- update(
            compiled_model,
            newdata = growthdat_unif,
            prior = prior_eilers_peeters_rescaled,
            init = inits_fun,
            init_r = 0,
            iter = 5000,
            chains = 4,
            cores = 4,
            seed = 1926,
            file = here(light_folder, paste0(unif_base, ".rds")),
            file_refit = "always",
            refresh = 0,
            control = list(
                adapt_delta = 0.99,
                max_treedepth = 20
            )
        )
    }
}



