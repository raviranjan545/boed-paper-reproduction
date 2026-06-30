# Declare present directory
here::i_am("1_scripts/nutrients_7_visualising_design_comparison.R")

# Load libraries
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggdist)
library(patchwork)
library(ggfortify)
library(ggridges)
library(scales)
theme_set(theme_minimal(base_size = 14))

###### Read in files

## Source plotting function
source(here(
    "1_scripts", 
    "00_support_file_ridgeline_plot_function.R")
)

## True parameter values
params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "1000_simulated_monod_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(curve_id = row_number())

## Central curve parameter values
central_curve_params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "monod_central_curve_mbd_parameters.csv")
) 


## Experimental designs

# Optimal
opt_nut5_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                     "nutrients_SIG_optimal_design_5_points.rds")),
                        sort(final.d[[besti]]))
opt_nut7_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                     "nutrients_SIG_optimal_design_7_points.rds")),
                        sort(final.d[[besti]]))
opt_nut10_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "nutrients_SIG_optimal_design_10_points.rds")),
                         sort(final.d[[besti]]))
opt_nut15_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "nutrients_SIG_optimal_design_15_points.rds")),
                         sort(final.d[[besti]]))

# Uniform
unif_nut5_levels <- round(seq(0.2, 25, length.out = 5), 1)
unif_nut7_levels <- round(seq(0.2, 25, length.out = 7), 1)
unif_nut10_levels <- round(seq(0.2, 25, length.out = 10), 1)
unif_nut15_levels <- round(seq(0.2, 25, length.out = 15), 1)

# Joining optimal and uniform designs into one data frame
designs <- data.frame(levels = c(opt_nut5_levels, opt_nut7_levels, 
                                 opt_nut10_levels, opt_nut15_levels,
                                 unif_nut5_levels, unif_nut7_levels, 
                                 unif_nut10_levels, unif_nut15_levels)) %>%
    mutate(design = c(rep("optimal", length(levels)/2), 
                      rep("uniform", length(levels)/2)),
           n_points = c(rep(5,5), rep(7,7), rep(10,10), rep(15,15),
                        rep(5,5), rep(7,7), rep(10,10), rep(15,15)))

# Removing objects no longer needed
rm(opt_nut5_levels, opt_nut7_levels, opt_nut10_levels, opt_nut15_levels,
   unif_nut5_levels, unif_nut7_levels, unif_nut10_levels, unif_nut15_levels)

## CRPS and ES scores
all_scores <- read.csv(here("3_simulation_result_summaries", 
                            "monod_fit_opt_unif_scoring_rules_comparison.csv")) %>%
    mutate(design = factor(design, levels = c("unif", "opt")))

# Calculate percentage of simulations where optimal design is better for plotting
pct_opt_better <- all_scores %>%
    select(n_points, simulation, design, CRPS_mumax, CRPS_k, ES) %>%
    pivot_wider(
        names_from = design,
        values_from = c(CRPS_mumax, CRPS_k, ES)
    ) %>%
    group_by(n_points) %>%
    summarise(
        mumax_opt_better = 100 * mean(CRPS_mumax_opt < CRPS_mumax_unif, na.rm = TRUE),
        k_opt_better = 100 * mean(CRPS_k_opt < CRPS_k_unif, na.rm = TRUE),
        ES_opt_better = 100 * mean(ES_opt < ES_unif, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        label_mumax = paste0(round(mumax_opt_better), "%"),
        label_k = paste0(round(k_opt_better), "%"),
        label_ES = paste0(round(ES_opt_better), "%")
    )

## Aggregate prediction error distributions
agg_err <- read.csv(here("3_simulation_result_summaries", 
                         "monod_fit_opt_unif_prediction_error_comparison.csv")) %>%
    select(design, n_points, simulation, mae_offdesign_mbd) %>%
    rename(agg_MAE = mae_offdesign_mbd) %>%
    mutate(design = factor(design, levels = c("unif", "opt"))
    )

# Calculate percentage of simulations where optimal design is better for plotting
agg_opt_better <- agg_err %>%
    pivot_wider(
        names_from = design,
        values_from = agg_MAE
    ) %>%
    group_by(n_points) %>%
    summarise(
        perc_opt_better = mean(opt < unif) * 100,
        .groups = "drop"
    ) %>%
    mutate(
        label_perc = paste0(round(perc_opt_better), "%")
    )

## Pointwise prediction error distributions
point_err <- read.csv(here("3_simulation_result_summaries", 
                           "monod_fit_opt_unif_pointwise_prediction_error_mbd.csv")) %>%
    group_by(design, n_points, env) %>%
    summarise(point_MAE = mean(abs_error),
              .groups = "drop_last") %>%
    mutate(design = factor(design, levels = c("unif", "opt")))


###### PLOTS 

#### 1. Curves

# Monod function 
monod_halfsat <- function(resource, mu_max, k){
    gr <- (mu_max * resource) / (resource + k)
    gr
}

# Set driver values 
driver_values <- seq(0, 25, 0.01)

# Generate central curve predictions
dat_mbd_curve <- tibble(
    driver = driver_values,
    mu_max = central_curve_params$mu_max,
    k = central_curve_params$k
) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, k), 
                               monod_halfsat))

# Generate curve dataset for plotting
dat_plot <- params %>%
    filter(curve_id < 11) %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, k), monod_halfsat))

# Plot legend info
cols <- c(
    optimal = "purple",
    uniform = "orange"
)
dummy <- data.frame(Design = names(cols))  

# Plot
a1 <- ggplot(dat_plot, aes(x = driver, y = y, group = curve_id)) +
    geom_line(alpha = 0.1) +
    geom_line(
        data = dat_mbd_curve,
        aes(x = driver, y = y),
        inherit.aes = FALSE,
        linewidth = 1.2
    ) +
    scale_x_continuous(
        limits = c(0, 25.5),
        breaks = c(0, 5, 10, 15, 20, 25),
        expand = expansion(mult = c(0.005, 0.005))
    ) +
    labs(
        title = "a) Nutrient curves",
        x = expression(Nutrient~concentration~(mu*M)),
        y = expression(Growth~rate~(day^{-1}))
    ) +
    scale_fill_manual(values = cols, name = NULL,
                      labels = c("Optimal", "Uniform")) +
    geom_rect(data = dummy, aes(fill = Design), inherit.aes = FALSE,
              xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf) +  
    guides(fill = guide_legend(direction = "horizontal",
                               keywidth  = unit(3, "lines"),
                               keyheight = unit(1, "lines"))) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 15),
        legend.margin = margin(t = 1, b = -1)
    )
a1

#### 2. Designs

# Plot calculations and parameters
n_order <- c(5, 7, 10, 15)

point_size <- 3.5
stack_step <- 0.007 * point_size
row_gap    <- 0.08

lane_offset <- c(
    optimal = 0.01,
    uniform = -0.03
)

dot_counts <- designs %>%
    mutate(
        level_plot = round(levels, 2),
        design = factor(design, levels = c("optimal", "uniform")),
        n_points = factor(n_points, levels = n_order)
    ) %>%
    count(n_points, design, level_plot, name = "count") %>%
    mutate(
        design_chr = as.character(design),
        rel_ymin = lane_offset[design_chr],
        rel_ymax = rel_ymin + (count - 1) * stack_step
    )

row_layout <- dot_counts %>%
    group_by(n_points) %>%
    summarise(
        rel_min = min(rel_ymin),
        rel_max = max(rel_ymax),
        row_height = rel_max - rel_min,
        .groups = "drop"
    ) %>%
    arrange(n_points) %>%
    mutate(
        row_bottom = c(0, head(cumsum(row_height + row_gap), -1)),
        y0 = row_bottom - rel_min
    )

# Generate plotting data
dot_df <- dot_counts %>%
    left_join(row_layout %>% select(n_points, y0), by = "n_points") %>%
    uncount(count, .id = "stack_id") %>%
    mutate(
        y = y0 + lane_offset[design_chr] + (stack_id - 1) * stack_step
    )

y_limits <- range(dot_df$y) + c(-0.025, 0.025)

# Plot
a2 <- ggplot(dot_df, aes(x = level_plot, y = y, fill = design)) +
    geom_point(
        shape = 21,
        alpha = 0.75,
        size = point_size,
        colour = "white",
        stroke = 0.12
    ) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(
        limits = c(0, 25.5),
        breaks = c(0, 5, 10, 15, 20, 25),
        expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(
        breaks = row_layout$y0,
        labels = as.character(row_layout$n_points),
        limits = y_limits,
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        title = "b) Designs",
        x = expression(Nutrient~concentration~(mu*M)),
        y = "# of experimental units"
    ) +
    guides(fill = "none") +
    theme(
        plot.title = element_text(size = 15),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
    )
a2



#### 3. CRPS and ES ridgeline plots of differences

p_mumax <- ridgeline_diff(all_scores, 
                          CRPS_mumax, 
                          label_x = 0.16,
                          label_y_offset = 0.7,
                          labels = pct_opt_better,
                          label_col = "label_mumax",
                          title = expression("c)"~Estimation~score~`difference,`~mu[max])) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(day^{-1})))

p_k <- ridgeline_diff(all_scores, 
                      CRPS_k, 
                      label_x = 0.16,
                      label_y_offset = 0.7,
                      labels = pct_opt_better,
                      label_col = "label_k",
                      title = expression("d)"~Estimation~score~`difference,`~K)) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(mu*M)),
         y = "")

p_ES <- ridgeline_diff(all_scores, 
                       ES, 
                       label_x = 0.16,
                       label_y_offset = 0.7,
                       labels = pct_opt_better,
                       label_col = "label_ES",
                       title = expression("e)"~Estimation~score~`difference,`~all~parameters)) + 
    labs(x = expression(ES["opt"] - ES["unif"]~(dimensionless)))

agg_err <- mutate(agg_err, MAE = agg_MAE)
p_pred <- ridgeline_diff(agg_err, 
                         MAE, 
                         label_x = 0.16,
                         label_y_offset = 0.7,
                         labels = agg_opt_better,
                         label_col = "label_perc",
                         title = expression("f)"~Prediction~error~difference)) + 
    labs(x = expression(MAE["opt"] - MAE["unif"]~(day^{-1})),
         y = "")

a3 <- p_mumax
a4 <- p_k
a5 <- p_ES
a6 <- p_pred

###### Patching all main text figs together

section_header <- ggplot() +
    annotate("text", x = 0, y = 0.75,
             label = "Comparing designs' performance",
             hjust = 0.05, size = 7) +
    annotate("text", x = 0, y = 0.15,
             label = "Negative values indicate optimal design is better",
             hjust = 0.05, size = 4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(5, 0, 5, 0))


a1 + a2 + section_header + a3 + a4 + a5 + a6 + 
    plot_layout(
        design = "
AA
BB
CC
DE
FG
",
        heights = c(0.8, 1.2, 0.32, 0.9, 0.9)
    ) 


ggsave(
    here(
        "4_figures", 
        "nutrients_plot_grid.pdf"), 
    height = 24.62 * (4/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)



###### SUPPLEMENTARY FIGURES

#### S1. Priors 

s <- ggdistribution(dlnorm, 
                    seq(0, 5, 0.01), 
                    mean = 0.01, 
                    sd = 0.4, 
                    fill = "blue",
                    colour = NA,
                    ylab = "Probability density",
                    xlab = expression(italic(mu)[max]~(day^{-1}))) + 
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("a)"~italic(mu)[max])) 

t <- ggdistribution(dlnorm, 
                    seq(0, 5, 0.01), 
                    mean = 0.3, 
                    sd = 0.5, 
                    fill = "blue",
                    colour = NA,
                    ylab = "",
                    xlab = expression(K~(mu*M))) +
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("b)"~K)) 


u <- ggdistribution(dlnorm, 
                    seq(0, 0.2, 0.001), 
                    mean = -2.3, 
                    sd = 0.1, 
                    fill = "blue",
                    colour = NA,
                    ylab = "Probability density",
                    xlab = expression(sigma~(day^{-1}))) +
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("c)"~sigma)) 
u

s + t + u + 
    plot_layout(
        design = 
            "
AB
C#
        "
    )


ggsave(
    here(
        "4_figures", 
        "nutrients_SI_priors.pdf"), 
    height = 24.62 * (4/4.9) * 1, 
    width = 18 * 1,
    units = "cm"
)


#### S2. CRPS and ES on log scale (Violin plots)
cols_to_mean <- c("CRPS_mumax", "CRPS_k", "ES")

mean_score <- all_scores %>%
    group_by(n_points, design) %>%
    summarise(
        across(
            all_of(cols_to_mean),
            ~ mean(.x, na.rm = TRUE),
            .names = "mean_{.col}"
        ),
        .groups = "drop"
    )

mean_agg <- agg_err %>%
    group_by(n_points, design) %>%
    summarise(
        mean_MAE = mean(agg_MAE),
        .groups = "drop"
    )


# A. Max growth
s1 <- all_scores %>%
    ggplot(., aes(CRPS_mumax, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_mumax, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 0.5,
            y = factor(n_points),
            label = label_mumax
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab("# of experimental units") + 
    xlab(expression(CRPS~(day^{-1}))) +
    ggtitle(expression(paste("a) Estimation score, ", ~mu[max])),
            subtitle = "lower is better") + 
    guides(fill = "none") +
    theme(
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    ) 
s1

# B. Ks
s2 <- all_scores %>%
    ggplot(., aes(CRPS_k, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_k, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 7,
            y = factor(n_points),
            label = label_k
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab("") + 
    xlab(expression(CRPS~"("*mu*M*")")) + 
    ggtitle(expression(paste("b) Estimation score, ", ~K)),
            subtitle = "lower is better") + 
    guides(fill = "none") + 
    theme(
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    )
s2

# C. ES
s3 <- all_scores %>%
    ggplot(., aes(ES, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_ES, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 7,
            y = factor(n_points),
            label = label_ES
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() +
    ylab("# of experimental units") + 
    xlab(expression(Energy~Score~"(dimensionless)")) +
    ggtitle(expression("c) Estimation score, all parameters"),
            subtitle = "lower is better") + 
    guides(fill = "none") + 
    theme(
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    ) 
s3

# D. Aggregate prediction error
s4 <-
    agg_err %>%
    ggplot(., aes(agg_MAE, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_agg,
        aes(x = mean_MAE, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = agg_opt_better,
        aes(
            x = 0.3,
            y = factor(n_points),
            label = label_perc
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() +
    ylab("") + 
    xlab(expression(Mean~Absolute~Error~(day^{-1}))) +
    ggtitle("d) Prediction error",
            subtitle = "lower is better") + 
    guides(fill = "none") +
    theme( 
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    )
s4


s1 + s2 + s3 + s4 + 
    plot_layout(
        design = "
AB
CD
",
        heights = c(1, 1, 1, 1),
    )


ggsave(here(
    "4_figures", 
    "nutrients_SI_log_scale.pdf"), 
    height = 24.62 * (2/3) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)

#### S3. Pointwise prediction error

plot_mae_pointwise <- point_err %>%
    mutate(design = recode(design, opt = "Optimal", unif = "Uniform"),
           design = factor(design, levels = c("Optimal", "Uniform"))) %>%
    ggplot(aes(x = env, y = point_MAE, colour = design)) +
    geom_line(linewidth = 1.5) +
    facet_wrap(~forcats::fct_rev(as.factor(n_points)),
               ncol = 1, strip.position = "right")  +
    labs(
        x = expression(Nutrient~concentration~"("*mu*M*")"), 
        y = expression(Pointwise~Mean~Absolute~Error~(day^{-1})),
        colour = "Design"
    ) + 
    scale_colour_manual(values = c(Optimal = "purple", Uniform = "orange")) +
    guides(colour = guide_legend(
        keywidth       = unit(4, "lines"),
        keyheight      = unit(1, "lines")
    )) +
    theme(
        strip.text.y.right = element_text(angle = 0),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
    )


ggsave(
    here(
        "4_figures", 
        "nutrients_SI_pointwise_error.pdf"), 
    plot_mae_pointwise,
    height = 24.62 * (3/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)



##### S5. All simulated curves
driver_values <- seq(0, 25, 0.1)

# Plot curves from simulated parameter combinations
dat_plot <- params %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, k), monod_halfsat))

#Plot the function output for parameter combinations
ggplot(dat_plot, aes(x = driver, y = y, group = mu_max)) +
    geom_line(alpha = 0.05) +
    scale_x_continuous(
        limits = c(0, 25.5),
        breaks = c(0, 5, 10, 15, 20, 25),
        expand = expansion(mult = c(0.005, 0.005))
    ) +
    labs(
        title = "Nutrient curves simulated from priors",
        x = expression(Nutrient~concentration~(mu*M)),
        y = expression(Growth~rate~(day^{-1}))
    ) +
    theme(
        strip.text.y.right = element_text(angle = 0),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
    )


ggsave(
    here(
        "4_figures",
        "nutrients_SI_1000_simulated_curves.pdf"),
    height = 24.62 * (2/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)

