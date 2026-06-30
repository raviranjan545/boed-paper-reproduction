# Declare present directory and create directory structure if not done
here::i_am("1_scripts/nutrients_1_SIG_optimal_design_calculation.R")
source(here::here("1_scripts", "00_support_file_create_directories_and_install_packages.R"))

# Load necessary packages
library(here)
library(parallel)
library(acebayes)

# Source the new stable C++ SIG calculation function
source(here(
    "1_scripts",
    "00_support_file_stable_SIG_CPP_function.R"))

# Call the amended utility function
source(here(
    "1_scripts",
    "00_support_file_stable_utility_function.R"))

#to assure that the new utility function will be able to call other hidden functions from the package.
environment(utilitynlmTemp) <- asNamespace("acebayes")

#to assure that other functions from the package will call the updated version of the function.
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
low<-0.2
upp<-25

#Define the prior function
monod_prior <- function(B){
    mumax <-rlnorm(n = B, meanlog = 0, sdlog = 0.4)
    k <- rlnorm(n = B, meanlog = 0, sdlog = 0.5)
    sig2 <- rlnorm(n = B, meanlog = -2.3*2, sdlog = 0.1*2)
    out <- cbind(mumax, k, sig2)
    colnames(out) <- c("mumax", "k","sig2")
    out}

## Create a discretised grid using a function. 
limits <- function(d,i,j){
    grid <- seq(from = low, to = upp,length.out=249)
    grid
}

########## n = 5 ###########################
## Specify the sample size (number of experimental units).
n<-5

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 249),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("x")))
}

# Calculate the design
nutResc5 <- pacenlm(formula = ~ ((mumax*exp(0.01)) * x /((k*exp(0.3)) + x)),
                    start.d = start.d, prior = monod_prior,
                    N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                    upper = upp, limits = limits,
                    method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",nutResc5$time))

# Save design
saveRDS(nutResc5,        
        here(
            "2_designs_and_other_simulation_inputs",
            "nutrients_SIG_optimal_design_5_points.rds")
)

########## n = 7 ###########################
## Specify the sample size (number of experimental units).
n<-7

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 249),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("x")))
}

# Calculate the design
nutResc7 <- pacenlm(formula = ~ ((mumax*exp(0.01)) * x /((k*exp(0.3)) + x)),
                    start.d = start.d, prior = monod_prior,
                    N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                    upper = upp, limits = limits,
                    method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",nutResc7$time))

# Save design
saveRDS(nutResc7,        
        here(
            "2_designs_and_other_simulation_inputs",
            "nutrients_SIG_optimal_design_7_points.rds")
)

# ###### n = 10 ######
## Specify the sample size (number of experimental units).
n<-10

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 249),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("x")))
}

# Calculate the design
nutResc10 <- pacenlm(formula = ~ ((mumax*exp(0.01)) * x /((k*exp(0.3)) + x)),
                     start.d = start.d, prior = monod_prior,
                     N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                     upper = upp, limits = limits,
                     method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",nutResc10$time))

# Save design
saveRDS(nutResc10,
        here(
            "2_designs_and_other_simulation_inputs",
            "nutrients_SIG_optimal_design_10_points.rds")
)

############n = 15 #######
## Specify the sample size (number of experimental units).
n<-15

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 249),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("x")))
}

# Calculate the design
nutResc15 <- pacenlm(formula = ~ ((mumax*exp(0.01)) * x /((k*exp(0.3)) + x)),
                     start.d = start.d, prior = monod_prior,
                     N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                     upper = upp, limits = limits,
                     method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",nutResc15$time))

# Save design
saveRDS(nutResc15,
        here(
            "2_designs_and_other_simulation_inputs",
            "nutrients_SIG_optimal_design_15_points.rds")
)
