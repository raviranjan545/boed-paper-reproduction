# Load necessary libraries
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr) 

# --- 1. Define Parameters and Setup ---
theme_set(theme_minimal(base_size = 14))
cols <- c(Custom = "#0072B2", Uniform = "orange") 

set.seed(42) # For reproducible simulations

#True values of parameters
mu_max_true <- 2
K_true <- 6

#Nutrient range
S_range <- c(0, 25)

#Grid of nutrient values
n_grid <- 1000
S_seq <- seq(S_range[1], S_range[2], length.out = n_grid)

#Monod function
monod_func <- function(S, mu, K) { (mu * S) / (K + S) }

# Calculate gradient of the Monod function
monod_gradient <- function(S, mu, K) {
  d_mu <- S / (K + S)
  d_K <- (-mu * S) / ((K + S)^2)
  return(c(d_mu, d_K))
}

# --- 2. Calculate Generalized Leverage ---
X <- matrix(NA, nrow = n_grid, ncol = 2)
for (i in 1:n_grid) {
  X[i, ] <- monod_gradient(S_seq[i], mu_max_true, K_true)
}

XtX_inv <- solve(t(X) %*% X)
h_ii <- rowSums((X %*% XtX_inv) * X)
gen_leverage <- n_grid * h_ii

# Find the S value that maximizes leverage
S_max_lev <- S_seq[which.max(gen_leverage)]
S_opt_point <- round(S_max_lev, 1)

plot_data <- data.frame(
  S = S_seq,
  Response = monod_func(S_seq, mu_max_true, K_true),
  GenLeverage = gen_leverage
)

# --- 3. Define the designs and styling ---
design_unif <- data.frame(S = seq(S_range[1], S_range[2], length.out = 5), Design = "Uniform")
design_opt <- data.frame(S = c(rep(S_opt_point, 3), rep(25, 2)), Design = "Custom")

designs_agg <- rbind(
  design_unif %>% group_by(S, Design) %>% summarise(count = n(), .groups = 'drop'),
  design_opt %>% group_by(S, Design) %>% summarise(count = n(), .groups = 'drop')
)

# Dynamically scale the stacked dots to sit below the leverage curve
max_lev <- max(plot_data$GenLeverage)
stack_step <- max_lev * 0.04
y_opt_base <- -max_lev * 0.08
y_uni_base <- -max_lev * 0.15

dot_df <- designs_agg %>%
  uncount(count, .id = "stack_id") %>%
  mutate(
    y = ifelse(Design == "Custom",
               y_opt_base + (stack_id - 1) * stack_step,
               y_uni_base + (stack_id - 1) * stack_step),
    Design = factor(Design, levels = c("Custom", "Uniform"))
  )

# --- 4. Simulation ---
n_sims <- 10000 
noise_sd <- 0.1 
sim_results <- data.frame()

for (i in 1:n_sims) {
  # --- Uniform design simulation ---
  y_unif <- monod_func(design_unif$S, mu_max_true, K_true) + rnorm(5, 0, noise_sd)
  fit_unif <- try(nls(y ~ (mu * S) / (K + S), 
                      data = list(S = design_unif$S, y = y_unif), 
                      start = list(mu = 2.5, K = 5)), silent = TRUE)
  if (!inherits(fit_unif, "try-error")) {
    sim_results <- rbind(sim_results, data.frame(Sim = i, Design = "Uniform", 
                                                 K_est = coef(fit_unif)["K"],
                                                 mu_est = coef(fit_unif)["mu"]))
  }
  
  # --- Custom design simulation ---
  y_opt <- monod_func(design_opt$S, mu_max_true, K_true) + rnorm(5, 0, noise_sd)
  fit_opt <- try(nls(y ~ (mu * S) / (K + S), 
                     data = list(S = design_opt$S, y = y_opt), 
                     start = list(mu = 2.5, K = 5)), silent = TRUE)
  if (!inherits(fit_opt, "try-error")) {
    sim_results <- rbind(sim_results, data.frame(Sim = i, Design = "Custom", 
                                                 K_est = coef(fit_opt)["K"],
                                                 mu_est = coef(fit_opt)["mu"]))
  }
}

sim_results$Design <- factor(sim_results$Design, levels = c("Custom", "Uniform"))

mu_mean_unif <- mean(sim_results$mu_est[sim_results$Design == "Uniform"])
mu_mean_opt  <- mean(sim_results$mu_est[sim_results$Design == "Custom"])
K_mean_unif  <- mean(sim_results$K_est[sim_results$Design == "Uniform"])
K_mean_opt   <- mean(sim_results$K_est[sim_results$Design == "Custom"])

# --- 5. Create the plots: Fig. 1 ---

# Panel a: Monod Function
p1 <- ggplot(plot_data, aes(x = S, y = Response)) +
  geom_line(color = "black", linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 25.5), breaks = seq(0, 25, 5), expand = expansion(add = c(0.8, 0.5))) + 
labs(title = "a) Monod function", 
       y = expression(Growth~rate~(day^{-1})), 
       x = expression(Nutrient~concentration~(mu*M))) +
  theme(plot.title = element_text(size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

dummy <- data.frame(Design = factor(c("Custom", "Uniform"), levels = c("Custom", "Uniform")))

# Panel b: Generalized leverage with overlaid designs
p2 <- ggplot() +
  # The leverage curve
  geom_line(data = plot_data, aes(x = S, y = GenLeverage), color = "black", linewidth = 1.2) +
  # The stacked dots 
  geom_point(data = dot_df, aes(x = S, y = y, fill = Design),
             shape = 21, alpha = 0.75, size = 5, colour = "white", stroke = 0.12, show.legend = FALSE) +
  geom_rect(data = dummy, aes(fill = Design), inherit.aes = FALSE,
            xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(0, 25.5), breaks = seq(0, 25, 5), expand = expansion(add = c(0.8, 0.5))) +
  labs(title = "b) Generalized leverage & designs", 
       y = "Generalized Leverage", 
       x = expression(Nutrient~concentration~(mu*M))) +
  guides(fill = guide_legend(direction = "horizontal",
                             keywidth  = unit(3, "lines"),
                             keyheight = unit(1, "lines"))) +
  coord_cartesian(ylim = c(-max_lev * 0.15, max_lev), clip = "off") +
  theme(plot.title = element_text(size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# Panel c: mu_max estimation distributions
p3 <- ggplot(sim_results, aes(x = mu_est, fill = Design)) +
  geom_density(alpha = 0.6, color = NA, show.legend = FALSE) +
  geom_vline(xintercept = mu_max_true, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = mu_mean_unif, color = "orange", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = mu_mean_opt, color = cols["Custom"], linetype = "solid", linewidth = 1) +  scale_fill_manual(values = cols) +
  coord_cartesian(xlim = c(1.5, 2.8)) + 
  labs(title = expression(paste("c) ", mu[max], " estimates ")), 
       y = "Density", 
       x = expression(Estimated~mu[max]~(day^{-1}))) +
  theme(plot.title = element_text(size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Panel d: K estimation distributions
p4 <- ggplot(sim_results, aes(x = K_est, fill = Design)) +
  geom_density(alpha = 0.6, color = NA, show.legend = FALSE) +
  geom_vline(xintercept = K_true, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = K_mean_unif, color = "orange", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = K_mean_opt, color = cols["Custom"], linetype = "solid", linewidth = 1) +  scale_fill_manual(values = cols) +
  coord_cartesian(xlim = c(0, 12)) + 
  labs(title = "d) K estimates", 
       y = "Density", 
       x = expression(Estimated~K~(mu*M))) +
  theme(plot.title = element_text(size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


# --- 6. Combine with patchwork ---
final_plot <- (p1 | p2) / guide_area() / (p3 | p4) + 
  plot_layout(guides = "collect", heights = c(1, 0.1, 1)) & 
  theme(legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0))

# Print to view
print(final_plot)

# Save the figure: NEEDS TO BE UPDATED
ggsave("/Users/rranjan/Documents/SideProjects/SideProjects/BayesianExperimentalDesigns/Manuscript/Plots/main_text/monod_Custom_vs_uniform_design_updated.pdf", plot = final_plot, width = 10, height = 8, units = "in")
