# ==============================================================================
# Exploratory Script: 5-Point SIG Optimal Design for Monod Kinetics
# ==============================================================================
# This script is designed for interactive exploration. It calculates a 5-point 
# optimal experimental design using the acebayes package, plots the parameter 
# priors and visualizes the resulting design. We strongly recommend reading the
# overview of ACE in the Supplementary Information of the paper before
# exploring designs with this script. 
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Package Loading
# ------------------------------------------------------------------------------
# We use the R package acebayes to do the design calculation
# We load acebayes for the design calculation, and the tidyverse/patchwork 
# suite for design plotting.

# Check for required packages and install if missing
local({
  pkgs <- c(
    "acebayes",
    "dplyr",
    "ggplot2",
    "parallel",
    "patchwork",
    "tidyr"
  )
  
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  
  if (length(missing)) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
})

# Load libraries
library(acebayes)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
theme_set(theme_minimal(base_size = 14))

# NOTE: The original script sourced custom C++ functions for stability. 
# library(here)
# here::i_am("1_scripts/nutrients_1_SIG_optimal_design_calculation.R")

# The original package's code occasionally results in Inf or -Inf SIG values 
# This happens due to numerical issues: underflow in case of 
# large negative log-likelihood values.
# We have fixed this issue and are sourcing the modified code below
# The files are in the Github repo.
source(here("1_scripts", "support_file_stable_SIG_CPP_function.R"))
source(here("1_scripts", "support_file_stable_utility_function.R"))
environment(utilitynlmTemp) <- asNamespace("acebayes")
assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")

# Set seed for reproducibility of the starting designs during exploration
set.seed(42)

# ------------------------------------------------------------------------------
# 2. Define Prior Distribution & Visualize
# ------------------------------------------------------------------------------
# The Monod model parameters (mumax, k, and measurement error variance sig2) 
# are drawn from lognormal distributions.

monod_prior <- function(B){
  mumax <- rlnorm(n = B, meanlog = 0.01, sdlog = 0.4)
  k <- rlnorm(n = B, meanlog = 0.3, sdlog = 0.5)
  sig2  <- rlnorm(n = B, meanlog = -2.3*2, sdlog = 0.1*2)
  
  out <- cbind(mumax, k, sig2)
  colnames(out) <- c("mumax", "k","sig2")
  return(out)
}

# --- Plot the Priors ---
# Draw 10,000 samples to visualize what the prior space looks like
prior_samples <- as.data.frame(monod_prior(10000))

p_mu <- ggplot(prior_samples, aes(x = mumax)) +
  geom_density(fill = "#56B4E9", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = expression("Prior for" ~ mu[max]), x = expression(mu[max]), y = "Density")

p_k <- ggplot(prior_samples, aes(x = k)) +
  geom_density(fill = "#009E73", alpha = 0.6, color = NA) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Prior for K", x = "K", y = "Density")

# Show prior plots side-by-side using patchwork
print(p_mu | p_k)

# ------------------------------------------------------------------------------
# 3. Setup Design Space and ACE Algorithm parameters
# ------------------------------------------------------------------------------
#Range of nutrient concentrations at which an experimental unit can be placed
low <- 0.2
upp <- 25

#Number of experimental units in the design
n   <- 5   # We are strictly looking for a 5-point design

# Experimental units can't typically be placed at any arbitrary value
# Limits function to constrain the acebayes search space to a grid
# Right now, this is set such that the step size is 0.1
# This can be adjusted to reflect experimental settings
limits <- function(d, i, j){
  seq(from = low, to = upp, length.out = 249)
}

# ------------------------------------------------------------------------------
# 4. Calculate Optimal Design
# ------------------------------------------------------------------------------
#First we try a single starting design and look at the final design we get
#Define a starting design by randomly picking 5 points from the grid
start.d <- matrix(sample(seq(from = low, to = upp, length.out = 249), n),
                  nrow = n, ncol = 1,
                  dimnames = list(as.character(1:n), c("x")))

# acenlm finds the SIG-optimal design (criterion = "SIG"),
# averaging over priors using Monte Carlo (method = "MC").
# N1 is the number of iterations in Phase 1.
# N2 is the number of iterations in Phase 2.
# B[1] is the sample size to use in the Bayesian comparison of two designs' SIG values 
# B[2] is the Monte Carlo sample size to use for integration 
# (N1, N2, B) are reduced for rapid exploration as compared to the manuscript
# Increase them later for more precision. 
# progress = TRUE prints out the current iteration
# This takes about 2 minutes to run on a Mac M3
# Running times might vary depending on the computer involved.
# Changing any of the settings will change the run time as well
nut_5_single_start <- acenlm(
  formula = ~ (mumax * x / (k + x)),
  start.d = start.d, 
  prior   = monod_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000), # Reduced from 40, 200, c(40000, 2000)
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", progress = TRUE
)


#Look at the summary
nut_5_single_start

# Look at the design at the end of Phase 1
# This is the design without any attempts to replicate experimental units
nut_5_single_start$phase1.d

# Look at whether the algorithm converged
# This is a plot of the estimated SIG value of the design at every iteration of Phase I
# Reminder that each SIG value is an estimate with error
# Therefore, the focus should not be on the highest value
# But whether the SIG estimate increases initially and
# consistently stays around a high value by the end.
plot(nut_5_single_start$phase1.trace, 
     type = "b",           # "b" means plot BOTH points and lines
     pch = 19,             # pch = 19 gives solid circle points
     col = "blue",         # Line and point color
     xlab = "Iteration (Phase I)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

# Look at the design at the end of Phase 2
# This is the design after attempts to replicate experimental units
# This is also the final design.
nut_5_single_start$phase2.d

# Look at whether Phase II converged
plot(nut_5_single_start$phase2.trace, 
     type = "b",           # "b" means plot BOTH points and lines
     pch = 19,             # pch = 19 gives solid circle points
     col = "blue",         # Line and point color
     xlab = "Iteration (Phase II)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

####### Mac/Linux version: Designs from several starting designs ########
# If you are on a Windows machine skip ahead to the Windows block


# So far we have calculated a design from only one starting design
# ACE is a stochastic algorithm, so to make sure it doesn't get stuck in a local optima, 
# it needs to be started from multiple starting designs. 
# The best design is then chosen out of all the final designs.
# For exploration, we start with only a few starting designs (C) 
# Note that start.d is a list of starting designs now
C <- 4 
start.d <- list()
for(i in 1:C){
  start.d[[i]] <- matrix(sample(seq(from = low, to = upp, length.out = 249), n),
                         nrow = n, ncol = 1,
                         dimnames = list(as.character(1:n), c("x")))
}

# pacenlm is a wrapper over acenlm 
# that runs acenlm starting from several starting designs in parallel.
# pacenlm uses mclapply for parallelization.
# Since mclapply only works on Mac and Linux, pacenlm will not work on Windows
# However, Windows users can use the 'parallel' function to parallelize acenlm
# Number of cores to use (leaves 1 free to keep the computer usable)
cores <- max(1, detectCores() - 1)

# As with acenlm above, (N1, N2, B) are reduced for rapid exploration.
# mc.cores specifies how many cores pacenlm will use
# This takes about 2 minutes and 15 seconds to run on a Mac M3
# Running times might vary depending on the computer involved.
# Changing any of the settings will change the run time as well
nut_5_multiple_starts <- pacenlm(
  formula = ~ (mumax * x / (k + x)),
  start.d = start.d, 
  prior   = monod_prior,
  N1 = 15, N2 = 50, B = c(5000, 1000),
  lower = low, upper = upp, limits = limits,
  method = "MC", criterion = "SIG", mc.cores = cores
)

#Look at the summary, including how long it took.
nut_5_multiple_starts

# This is the best final design.
nut_5_multiple_starts$d

# These are ALL the final designs.
nut_5_multiple_starts$final.d

# For the best final design, you can look at the design at the end of Phase 1
# Reminder that this is the design without any attempts to replicate experimental units
nut_5_multiple_starts$phase1.d

# For the best final design,look at whether Phase I converged
plot(nut_5_multiple_starts$phase1.trace, 
     type = "b",           # "b" means plot BOTH points and lines
     pch = 19,             # pch = 19 gives solid circle points
     col = "blue",         # Line and point color
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

# Look at the design at the end of Phase 2
# This is the design after attempts to replicate experimental units
# This is also the final design.
nut_5_single_start$phase2.d

# For the best final design,look at whether Phase II converged
plot(nut_5_multiple_starts$phase2.trace, 
     type = "b",           # "b" means plot BOTH points and lines
     pch = 19,             # pch = 19 gives solid circle points
     col = "blue",         # Line and point color
     xlab = "Iteration (Phase I Stage)", 
     ylab = "Approximate Expected SIG",
     main = "ACE Algorithm Convergence Trace")

####### Windows version: Designs from several starting designs ########
# pacenlm uses mclapply (fork-based), which does not work on Windows.
# This block reproduces pacenlm by running acenlm on each starting design
# across a PSOCK cluster (parLapply), then selecting the best final design
# the same way pacenlm does: re-evaluate each design n.assess times and
# compare the MEAN estimated SIG (each SIG value is a noisy MC estimate).
#
# Mac/Linux users: use the pacenlm block above instead.
# Windows users: comment out the pacenlm block and use this.
# ==============================================================================

# cores <- max(1, detectCores() - 1)
# 
# # 1. Create a PSOCK cluster (works on Windows, Mac, and Linux)
# cl <- makeCluster(cores)
# 
# # 2. PSOCK workers start empty: load packages and ship needed objects to each.
# #    (Forking would inherit these automatically; PSOCK does not.)
# clusterEvalQ(cl, library(acebayes))
# clusterExport(cl,
#               varlist = c("monod_prior", "low", "upp", "n", "limits"),
#               envir   = environment())
# 
# # 3. Source the stable versions of utility calculations.
# clusterEvalQ(cl, {
#   source(here::here("1_scripts", "support_file_stable_SIG_CPP_function.R"))
#   source(here::here("1_scripts", "support_file_stable_utility_function.R"))
#   environment(utilitynlmTemp) <- asNamespace("acebayes")
#   assignInNamespace("utilitynlm", utilitynlmTemp, ns = "acebayes")
# })
# 
# 
# # 4. Run acenlm once per starting design, in parallel across the cluster.
# run_one_start <- function(sd) {
#   acenlm(
#     formula = ~ ((mumax*exp(0.01)) * x / ((k*exp(0.3)) + x)),
#     start.d = sd,
#     prior   = monod_prior,
#     N1 = 15, N2 = 50, B = c(5000, 1000),
#     lower = low, upper = upp, limits = limits,
#     method = "MC", criterion = "SIG", progress = FALSE
#   )
# }
# 
# ace_runs <- parLapply(cl, start.d, run_one_start)
# 
# # 5. Select the best final design the way pacenlm does.
# #    You could simply pick the design with the highest SIG estimate.
# #    But each design's SIG is a noisy MC estimate, so we re-evaluate each final
# #    design n.assess times and compare the MEAN (not a single value).
# n.assess <- 20
# B_assess <- 5000   # = B[1] used above
# 
# assess_one <- function(run, n.assess, B_assess) {
#   u <- run$utility                 # the design-scale utility acenlm returns
#   d <- run$phase2.d                # the final (Phase II) design for this start
#   ev <- numeric(n.assess)
#   for (k in seq_len(n.assess)) ev[k] <- mean(u(d = d, B = B_assess))
#   ev
# }
# 
# # Ship the assessment helper + args, then evaluate in parallel.
# clusterExport(cl, varlist = c("assess_one", "n.assess", "B_assess"),
#               envir = environment())
# eval_list <- parLapply(cl, ace_runs, assess_one, n.assess, B_assess)
# 
# stopCluster(cl)   # always shut the cluster down when done
# 
# # 7. Best start = highest mean estimated SIG, and rebuild a pacenlm-like object
# #    so the downstream code (which reads $d, $final.d, $phase1.d, the traces)
# #    works unchanged.
# besti <- which.max(vapply(eval_list, mean, numeric(1)))
# 
# #This returns an equivalent object back
# nut_5_multiple_starts <- list(
#   d            = ace_runs[[besti]]$phase2.d,                 # best final design
#   final.d      = lapply(ace_runs, function(r) r$phase2.d),   # all final designs
#   phase1.d     = ace_runs[[besti]]$phase1.d,
#   phase2.d     = ace_runs[[besti]]$phase2.d,
#   phase1.trace = ace_runs[[besti]]$phase1.trace,
#   phase2.trace = ace_runs[[besti]]$phase2.trace,
#   eval         = eval_list[[besti]],
#   besti        = besti
# )
# 
# ------------------------------------------------------------------------------
# 5. Extract and Visualize the Resulting Design
# ------------------------------------------------------------------------------
# Extract the best design points from the acebayes object
SIG_opt_design <- nut_5_multiple_starts$d

# Prepare design to be plotted
design_df <- data.frame(S = as.numeric(SIG_opt_design), Design = "Optimal (n=5)")

# Aggregate counts for each unique level
designs_agg <- design_df %>% 
  group_by(S, Design) %>% 
  summarise(count = n(), .groups = 'drop')

# Uncount and map to a Y-axis stack position
# Adjust the vertical distance between two points stacked on top of each other
stack_step <- 0.02
y_base <- 0

dot_df <- designs_agg %>%
  uncount(count, .id = "stack_id") %>%
  mutate(y = y_base + (stack_id - 1) * stack_step)

# --- Plot the Design ---
p_design <- ggplot() +
  # Draw a baseline to anchor the points
  geom_hline(yintercept = -0.05, color = "black", linewidth = 0.5) +
  geom_point(data = dot_df, aes(x = S, y = y),
             fill = "#0072B2", shape = 21, alpha = 0.85, 
             size = 6, colour = "white", stroke = 0.5) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 25, 5)) +
  coord_cartesian(ylim = c(-0.1, max(dot_df$y) + 0.2)) +
  labs(title = "5-point SIG-optimal design for Monod",
       subtitle = "Each point represents one experimental unit",
       x = expression(Nutrient~concentration~(mu*M)), 
       y = NULL) +
  theme(plot.title = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

# Display the final design plot
print(p_design)

