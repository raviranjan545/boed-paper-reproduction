# Declare present directory and create directory structure if not done
here::i_am("1_scripts/temperature_1_SIG_optimal_design_calculation.R")
source(here::here("1_scripts", "00_support_file_create_directories_and_install_packages.R"))

# Load necessary packages
library(here)
library(parallel)
library(acebayes)
library(MASS)
library(Matrix)

# Source the new stable C++ SIG calculation function
source(here(
    "1_scripts",
    "00_support_file_stable_SIG_CPP_function.R"))

# Call the amended utility function
source(here(
    "1_scripts",
    "00_support_file_stable_utility_function.R"))

# Load the prior function
source(here(
    "1_scripts",
    "00_support_file_temperature_prior.R"))


# to assure that the new utility function will be able to call 
# other hidden functions from the package.
environment(utilitynlmTemp) <- asNamespace("acebayes")

# to assure that other functions from the package will call 
# the updated version of the function.
assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")

## Define the number of starting conditions
C <- 40

## IMPORTANT: Number of cores to be used
# We recommend running this code on a machine with at least 40 cores. 
# With fewer cores, it might take more than a week to finish with current settings.
# If you want to run it on a laptop, we recommend reducing:
# the number of starting conditions C and 
# these settings in pacenlm below: N1 = 40, N2 = 200, B = c(40000, 2000)
# More details about the settings in Section 2, SI
cores <- detectCores() - 1

## Set seed for reproducibility.
set.seed(1)

# Lower and upper limits of the support of the design space
low<-5
upp<-35

## Prior function is in the temperature_prior.R script
## Create a discretized grid.
## Steps of 0.5 Celsius
limits <- function(d,i,j){
    grid <- seq(from = low, to = upp,length.out=61)
    grid
}

# ########## n = 5 ###########################
# ## Specify the sample size (number of experimental units).
n<-5

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 61),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("temp")))
}


#calculate the design
tempResc5 <- pacenlm(formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
                     start.d = start.d, prior = norberg_prior,
                     N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                     upper = upp, limits = limits,
                     method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",tempResc5$time))

# Save design
saveRDS(tempResc5,
        here(
            "2_designs_and_other_simulation_inputs",
            "temperature_SIG_optimal_design_5_points.rds")
)

########## n = 7 ###########################
## Specify the sample size (number of experimental units).
n<-7

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 61),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("temp")))
}

#calculate the design
tempResc7 <- pacenlm(formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
                     start.d = start.d, prior = norberg_prior,
                     N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                     upper = upp, limits = limits,
                     method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",tempResc7$time))

# Save design
saveRDS(tempResc7,
        here(
            "2_designs_and_other_simulation_inputs",
            "temperature_SIG_optimal_design_7_points.rds")
)

######################################################################################################
## Specify the sample size (number of experimental units).
n<-10

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 61),n), 
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("temp")))
}

#calculate the design
tempResc10 <- pacenlm(formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
                      start.d = start.d, prior = norberg_prior, 
                      N1 = 40, N2 = 200, B = c(40000, 2000), lower = low, 
                      upper = upp, limits = limits,
                      method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",tempResc10$time))

# Save design
saveRDS(tempResc10,
        here(
            "2_designs_and_other_simulation_inputs",
            "temperature_SIG_optimal_design_10_points.rds")
)

######################################################################################################
## Specify the sample size (number of experimental units).
n<-15

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 61),n), 
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("temp")))
}

#calculate the design
tempResc15 <- pacenlm(formula = ~ (a*exp(b*temp))*(tmax - temp)*(temp - tmin),
                      start.d = start.d, prior = norberg_prior, 
                      N1 = 40, N2 = 200, B = c(40000, 2000), lower = low, 
                      upper = upp, limits = limits,
                      method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",tempResc15$time))

# Save design
saveRDS(tempResc15,
        here(
            "2_designs_and_other_simulation_inputs",
            "temperature_SIG_optimal_design_15_points.rds")
)
