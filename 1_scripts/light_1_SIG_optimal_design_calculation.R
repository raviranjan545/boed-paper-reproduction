# Declare present directory and create directory structure if not done
here::i_am("1_scripts/light_1_SIG_optimal_design_calculation.R")
source(here::here("1_scripts", "00_support_file_create_directories_and_install_packages.R"))

# Load libraries
library(here)
library(parallel)
library(acebayes)

#Source the new stable C++ SIG calculation function
source(here(
    "1_scripts",
    "00_support_file_stable_SIG_CPP_function.R"))

#Call the amended utility function
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
low<-10
upp<-1000

#Define the prior function
eilers_peeters_prior <- function(B){
    mu_max <-rlnorm(n = B, meanlog = 0, sdlog = 0.4)
    alpha <- rlnorm(n = B, meanlog = 0, sdlog = 0.8)
    i_opt <- rlnorm(n = B, meanlog = 0, sdlog = 0.3)
    sig2 <- rlnorm(n = B, meanlog = -2.3*2, sdlog = 0.1*2)
    out <- cbind(mu_max,alpha,i_opt,sig2)
    colnames(out) <- c("mu_max","alpha","i_opt","sig2")
    out}

## Create a discretised grid using a function. 
limits <- function(d,i,j){
    grid <- seq(from = low, to = upp,length.out=100)
    grid
}

########## n = 5 ###########################
## Specify the sample size (number of experimental units).
n<-5

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 100),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("light")))
}

# Calculate the design
lightResc5 <- pacenlm(formula = ~ (((mu_max*exp(0.01)) * light) / (((mu_max*exp(0.01) / ((alpha*exp(-3)) * (i_opt*exp(5.5))^2)) * (light^2)) +
                                                                       ((1 - ((2 * mu_max*exp(0.01))/((alpha*exp(-3)) * (i_opt*exp(5.5))))) * light) +
                                                                       (mu_max*exp(0.01) / (alpha*exp(-3))))),
                      start.d = start.d, prior = eilers_peeters_prior,
                      N1 = 40, N2 = 200, B = c(40000, 2000),
                      lower = low, upper = upp, limits = limits,
                      method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",lightResc5$time))

# Save design
saveRDS(lightResc5,
        here(
            "2_designs_and_other_simulation_inputs",
            "light_SIG_optimal_design_5_points.rds")
)

########## n = 7 ###########################
## Specify the sample size (number of experimental units).
n<-7

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 100),n),
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("light")))
}

# Calculate the design
lightResc7 <- pacenlm(formula = ~ (((mu_max*exp(0.01)) * light) / (((mu_max*exp(0.01) / ((alpha*exp(-3)) * (i_opt*exp(5.5))^2)) * (light^2)) +
                                                                       ((1 - ((2 * mu_max*exp(0.01))/((alpha*exp(-3)) * (i_opt*exp(5.5))))) * light) +
                                                                       (mu_max*exp(0.01) / (alpha*exp(-3))))),
                      start.d = start.d, prior = eilers_peeters_prior,
                      N1 = 40, N2 = 200, B = c(40000, 2000), lower = low,
                      upper = upp, limits = limits,
                      method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",lightResc7$time))

# Save design
saveRDS(lightResc7,
        here(
            "2_designs_and_other_simulation_inputs",
            "light_SIG_optimal_design_7_points.rds")
)

######################################################################################################
## Specify the sample size (number of experimental units).
n<-10

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 100),n), 
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("light")))
}

# Calculate the design
lightResc10 <- pacenlm(formula = ~ (((mu_max*exp(0.01)) * light) / (((mu_max*exp(0.01) / ((alpha*exp(-3)) * (i_opt*exp(5.5))^2)) * (light^2)) +
                                                                        ((1 - ((2 * mu_max*exp(0.01))/((alpha*exp(-3)) * (i_opt*exp(5.5))))) * light) +
                                                                        (mu_max*exp(0.01) / (alpha*exp(-3))))),
                       start.d = start.d, prior = eilers_peeters_prior, 
                       N1 = 40, N2 = 200, B = c(40000, 2000), lower = low, 
                       upper = upp, limits = limits,
                       method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",lightResc10$time))

# Save design
saveRDS(lightResc10,
        here(
            "2_designs_and_other_simulation_inputs",
            "light_SIG_optimal_design_10_points.rds")
)

######################################################################################################
## Specify the sample size (number of experimental units).
n<-15

# Store all starting conditions in a list
start.d <- list()
for(i in 1:C){
    start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 100),n), 
                           nrow = n, ncol = 1,
                           dimnames = list(as.character(1:n), c("light")))
}

# Calculate the design
lightResc15 <- pacenlm(formula = ~ (((mu_max*exp(0.01)) * light) / (((mu_max*exp(0.01) / ((alpha*exp(-3)) * (i_opt*exp(5.5))^2)) * (light^2)) +
                                                                        ((1 - ((2 * mu_max*exp(0.01))/((alpha*exp(-3)) * (i_opt*exp(5.5))))) * light) +
                                                                        (mu_max*exp(0.01) / (alpha*exp(-3))))),
                       start.d = start.d, prior = eilers_peeters_prior, 
                       N1 = 40, N2 = 200, B = c(40000, 2000), lower = low, 
                       upper = upp, limits = limits,
                       method = "MC",criterion="SIG",mc.cores=cores)

print(paste("n=",n,"finished with time=",lightResc15$time))

# Save design
saveRDS(lightResc15,
        here(
            "2_designs_and_other_simulation_inputs",
            "light_SIG_optimal_design_15_points.rds")
)
