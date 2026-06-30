# Declare present directory
here::i_am("1_scripts/nutrients_4_checking_and_saving_posteriors.R")

# Load libraries
library(here)
library(brms)
library(posterior)
library(dplyr)

# Factors for conversion
mu_max_adjustment <- exp(0.01)
k_adjustment <- exp(0.3)

# Define number of points
n_points_list <- c(5, 7, 10, 15)

# Create function to identify fits that have not converged (Rhat > 1.01)
check_rhat <- function(fit, simulation, n_points, design, file) {
    
    rh <- rhat(fit)
    
    if (any(rh > 1.01, na.rm = TRUE)) {
        tibble(
            simulation = simulation,
            n_points   = n_points,
            design     = design,
            file       = basename(file),
            max_rhat   = max(rh, na.rm = TRUE)
        )
    }
}

# Define list to store problematic fits for examination and possible re-fitting
fits_to_examine <- list()

# Loop through all saved model files and summarise
for (i in 1:length(n_points_list)){
    
    # Select number of points 
    n_points <- n_points_list[i]
    
    # Define lists of saved model files
    fit_opt_list <- list.files(here("nutrient_curves_simulations_and_fits"), 
                               pattern = paste0("^monod_fit_opt_", n_points), 
                               full.names = TRUE)
    fit_unif_list <- list.files(here("nutrient_curves_simulations_and_fits"), 
                                pattern = paste0("^monod_fit_unif_", n_points), 
                                full.names = TRUE)
    
    # Loop over the 1000 model fits for each design and summarise model fits per simulation
    for (j in 1:length(fit_unif_list)) {
        
        print(c(i, j))
        
        # Read in the fitted models
        fit_opt <- readRDS(fit_opt_list[j])
        fit_unif <- readRDS(fit_unif_list[j])
        
        # Check Rhat values and store simulation information if any >1.01         
        opt_check <- check_rhat(fit_opt, j, n_points, "opt", fit_opt_list[j])
        if (!is.null(opt_check)) {
            fits_to_examine[[length(fits_to_examine) + 1L]] <- opt_check
        }
        
        unif_check <- check_rhat(fit_unif, j, n_points, "unif", fit_unif_list[j])
        if (!is.null(unif_check)) {
            fits_to_examine[[length(fits_to_examine) + 1L]] <- unif_check
        }
        
        # Extract posterior draws as a draws_df
        draws_opt <- as_tibble(fit_opt) %>%
            select(c("b_mumax_Intercept", "b_k_Intercept")) %>%
            mutate(b_mumax_Intercept = b_mumax_Intercept * mu_max_adjustment,
                   b_k_Intercept = b_k_Intercept * k_adjustment
            )
        
        draws_unif <- as_tibble(fit_unif) %>%
            select(c("b_mumax_Intercept", "b_k_Intercept")) %>%
            mutate(b_mumax_Intercept = b_mumax_Intercept * mu_max_adjustment,
                   b_k_Intercept = b_k_Intercept * k_adjustment
            )
        
        # Save posteriors as CSVs for subsequent analysis
        write.csv(draws_opt, 
                  here("nutrient_curves_posteriors/", 
                       paste0(sprintf("monod_posterior_opt_%spts_%04d", n_points, j), 
                              ".csv")), 
                  row.names = FALSE)
        write.csv(draws_unif, 
                  here("nutrient_curves_posteriors/", 
                       paste0(sprintf("monod_posterior_unif_%spts_%04d", n_points, j), 
                              ".csv")), 
                  row.names = FALSE)
        
    }
}

# Save fits that need examination, if any
# Check and re-fit if needed
fits_to_examine <- bind_rows(fits_to_examine)

print(fits_to_examine)

if(nrow(fits_to_examine) > 0){
    write.csv(fits_to_examine,
              here("3_simulation_result_summaries",
                   paste0("monod_problem_fits_to_check_", 
                          n_points, 
                          ".csv")
              ), 
              row.names = FALSE)
}
