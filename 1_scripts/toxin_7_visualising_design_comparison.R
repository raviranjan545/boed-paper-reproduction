# Declare present directory
here::i_am("1_scripts/toxin_7_visualising_design_comparison.R")

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
        "1000_simulated_toxin_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(curve_id = row_number())

## Central curve parameter values
central_curve_params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "toxin_central_curve_mbd_parameters.csv")
) 

## Experimental designs

# Optimal
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

# Uniform
design_space <- c(0, 10^seq(log10(0.1), log10(1000), length.out = 101)) 

unif_tox5_levels <- design_space[round(seq(1, 102, length.out = 5), 0)]
unif_tox7_levels <- design_space[round(seq(1, 102, length.out = 7), 0)]
unif_tox10_levels <- design_space[round(seq(1, 102, length.out = 10), 0)]
unif_tox15_levels <- design_space[round(seq(1, 102, length.out = 15), 0)]

# Joining optimal and uniform designs into one data frame
designs <- data.frame(levels = c(opt_tox5_levels, opt_tox7_levels, 
                                 opt_tox10_levels, opt_tox15_levels,
                                 unif_tox5_levels, unif_tox7_levels, 
                                 unif_tox10_levels, unif_tox15_levels)) %>%
    mutate(design = c(rep("optimal", length(levels)/2), 
                      rep("uniform", length(levels)/2)),
           n_points = c(rep(5,5), rep(7,7), rep(10,10), rep(15,15),
                        rep(5,5), rep(7,7), rep(10,10), rep(15,15)))

# Removing objects no longer needed
rm(opt_tox5_levels, opt_tox7_levels, opt_tox10_levels, opt_tox15_levels,
   unif_tox5_levels, unif_tox7_levels, unif_tox10_levels, unif_tox15_levels,
   design_space)


## CRPS and ES scores
all_scores <- read.csv(here("3_simulation_result_summaries", 
                            "toxin_fit_opt_unif_scoring_rules_comparison.csv")) %>%
    mutate(design = factor(design, levels = c("unif", "opt")))

# Calculate percentage of simulations where optimal design is better for plotting
pct_opt_better <- all_scores %>%
    select(n_points, simulation, design, CRPS_mumax, CRPS_ec50, CRPS_slope, ES) %>%
    pivot_wider(
        names_from = design,
        values_from = c(CRPS_mumax, CRPS_ec50, CRPS_slope, ES)
    ) %>%
    group_by(n_points) %>%
    summarise(
        mumax_opt_better = 100 * mean(CRPS_mumax_opt < CRPS_mumax_unif, na.rm = TRUE),
        ec50_opt_better = 100 * mean(CRPS_ec50_opt < CRPS_ec50_unif, na.rm = TRUE),
        slope_opt_better = 100 * mean(CRPS_slope_opt < CRPS_slope_unif, na.rm = TRUE),
        ES_opt_better = 100 * mean(ES_opt < ES_unif, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        label_mumax = paste0(round(mumax_opt_better), "%"),
        label_ec50 = paste0(round(ec50_opt_better), "%"),
        label_slope = paste0(round(slope_opt_better), "%"),
        label_ES = paste0(round(ES_opt_better), "%")
    )

## Aggregate prediction error distributions
agg_err <- read.csv(here("3_simulation_result_summaries", 
                         "toxin_fit_opt_unif_prediction_error_comparison.csv")) %>%
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
                           "toxin_fit_opt_unif_pointwise_prediction_error_mbd.csv")) %>%
    group_by(design, n_points, env) %>%
    summarise(point_MAE = mean(abs_error),
              .groups = "drop_last") %>%
    mutate(design = factor(design, levels = c("unif", "opt")))




###### PLOTS 

#### 1. Curves

# Log-logistic toxin function from Ritz (2010)
loglogis <- function(toxin, mu_max, ec50, slope){
    gr <- mu_max / (1 + ((toxin/ec50)^slope))
    gr
}

# Set driver values 
driver_values <- c(0, 10^seq(log10(0.1), log10(1000), length.out = 101)) 

# Generate central curve predictions
dat_mbd_curve <- tibble(
    driver = driver_values,
    mu_max = central_curve_params$mu_max,
    ec50 = central_curve_params$ec50,
    slope = central_curve_params$slope
) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, ec50, slope), 
                               loglogis))

# Generate curve dataset for plotting
dat_plot <- params %>%
    filter(curve_id < 11) %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, ec50, slope), 
                               loglogis))

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
        limits = c(0, 1000),
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
        breaks = c(0, 1, 10, 100, 1000),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
        title = "a) Toxin curves",
        x = expression(Toxin~concentration~(mu*g~{L^-1})),
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
        limits = c(0, 1000),
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
        breaks = c(0, 1, 10, 100, 1000),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
        breaks = row_layout$y0,
        labels = as.character(row_layout$n_points),
        limits = y_limits,
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        title = "b) Designs",
        x = expression(Toxin~concentration~(mu*g~{L^-1})),
        y = "# of experimental units"
    ) +
    guides(fill = "none") +
    theme(
        plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 13.5),
        axis.title.y = element_text(size = 13.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )
a2


#### 3. CRPS and ES ridgeline plots of differences

p_mumax <- ridgeline_diff(all_scores, 
                          CRPS_mumax, 
                          label_x = 0.17,
                          label_y_offset = 0.7,
                          labels = pct_opt_better,
                          label_col = "label_mumax",
                          title = expression("c)"~Estimation~score~`difference,`~mu[max])) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(day^{-1})))

p_ec50 <- ridgeline_diff(all_scores, 
                         CRPS_ec50, 
                         label_x = 0.17,
                         label_y_offset = 0.7,
                         labels = pct_opt_better,
                         label_col = "label_ec50",
                         title = expression("d)"~Estimation~score~`difference,`~italic(e))) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(mu*g~{L^-1})),
         y = "")

p_slope <- ridgeline_diff(all_scores, 
                          CRPS_slope, 
                          label_x = 0.17,
                          label_y_offset = 0.7,
                          labels = pct_opt_better,
                          label_col = "label_slope",
                          title = expression("e)"~Estimation~score~`difference,`~italic(h))) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(dimensionless)))

p_ES <- ridgeline_diff(all_scores, 
                       ES, 
                       label_x = 0.17,
                       label_y_offset = 0.7,
                       labels = pct_opt_better,
                       label_col = "label_ES",
                       title = expression("f)"~Estimation~score~`difference,`~all~parameters)) + 
    labs(x = expression(ES["opt"] - ES["unif"]~(dimensionless)),
         y = "")

agg_err <- mutate(agg_err, MAE = agg_MAE)
p_pred <- ridgeline_diff(agg_err, 
                         MAE, 
                         label_x = 0.17,
                         label_y_offset = 0.7,
                         labels = agg_opt_better,
                         label_col = "label_perc",
                         title = expression("g)"~Prediction~error~difference)) + 
    labs(x = expression(MAE["opt"] - MAE["unif"]~(day^{-1})))

a3 <- p_mumax
a4 <- p_ec50
a5 <- p_slope
a6 <- p_ES
a7 <- p_pred

###### All main text figs together

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

a1 + a2 + section_header + a3 + a4 + a5 + a6 + a7 +
    plot_layout(
        design = "
AA
BB
CC
DE
FG
H#
",
        heights = c(0.8, 1.2, 0.32, 0.9, 0.9, 0.9)
    )


ggsave(
    here(
        "4_figures", 
        "toxin_plot_grid.pdf"), 
    height = 24.62 * (4.9/4.9) * 1.4, 
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
                    seq(0, 500, 0.1), 
                    mean = 4, 
                    sd = 1, 
                    fill = "blue",
                    colour = NA,
                    ylab = "",
                    xlab = expression(italic(e)~(mu*g~{L^-1}))) + 
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("b)"~italic(e)))

u <- ggdistribution(dlnorm, 
                    seq(0, 10, 0.1),
                    mean = 1, 
                    sd = 0.5, 
                    fill = "blue",
                    colour = NA,
                    ylab = "Probability density",
                    xlab = expression(italic(h)~(dimensionless))) + 
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("c)"~italic(h)))


v <- ggdistribution(dlnorm, 
                    seq(0, 0.2, 0.001), 
                    mean = -2.3, 
                    sd = 0.1, 
                    fill = "blue",
                    colour = NA,
                    ylab = "",
                    xlab = expression(sigma~(day^{-1}))) +
    theme(
        panel.grid.minor = element_blank()
    ) + 
    ggtitle(expression("d)"~sigma)) 

s + t + u + v + 
    plot_layout(
        design = "
AB
CD
")

ggsave(
    here(
        "4_figures", 
        "toxin_SI_priors.pdf"), 
    height = 24.62 * (4/4.9) * 1, 
    width = 18 * 1,
    units = "cm"
)

#### S2. CRPS and ES on log scale (Violin plots)
cols_to_mean <- c("CRPS_mumax", "CRPS_ec50", "CRPS_slope", "ES")

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
            x = 0.38,
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

# B. e
s2 <- all_scores %>%
    ggplot(., aes(CRPS_ec50, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_ec50, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 2500,
            y = factor(n_points),
            label = label_ec50
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab("") + 
    xlab(expression(CRPS~"("*mu*g~{L^-1}*")")) + 
    ggtitle(expression(paste("b) Estimation score, ", italic(e))),
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

# C. h
s3 <- all_scores %>%
    ggplot(., aes(CRPS_slope, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_slope, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 18,
            y = factor(n_points),
            label = label_slope
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab("# of experimental units") + 
    xlab(expression(CRPS~"(dimensionless)")) + 
    ggtitle(expression(paste("c) Estimation score, ", italic(h))),
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

# D. ES
s4 <- all_scores %>%
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
            x = 15,
            y = factor(n_points),
            label = label_ES
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() +
    ylab("") + 
    xlab(expression(Energy~Score~"(dimensionless)")) + 
    ggtitle(expression("d) Estimation score, all parameters"),
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

# E. Aggregate prediction error
s5 <- agg_err %>%
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
            x = 0.4,
            y = factor(n_points),
            label = label_perc
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() +
    ylab("# of experimental units") + 
    xlab(expression(Mean~Absolute~Error~(day^{-1}))) +
    ggtitle("e) Prediction error",
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
s5


s1 + s2 + s3 + s4 + s5 + 
    plot_layout(
        design = "
AB
CD
E#
",
        heights = c(1, 1, 1, 1, 1)
    )

ggsave(here(
    "4_figures", 
    "toxin_SI_log_scale.pdf"), 
    height = 24.62 * (3/3) * 1.4, 
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
        x = expression(Toxin~concentration~(mu*g~{L^-1})),
        y = expression(Pointwise~Mean~Absolute~Error~(day^{-1})),
        colour = "Design"
    ) + 
    scale_colour_manual(values = c(Optimal = "purple", Uniform = "orange")) +
    guides(colour = guide_legend(
        keywidth       = unit(4, "lines"),
        keyheight      = unit(1, "lines")
    )) +
    scale_x_continuous(
        limits = c(0, 1000),
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
        breaks = c(0, 1, 10, 100, 1000),
        expand = expansion(
            mult = c(0.01, 0.03))
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
        "toxin_SI_pointwise_error.pdf"), 
    plot_mae_pointwise,
    height = 24.62 * (3/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)

##### S4. All simulated curves

# Plot curves from simulated parameter combinations
dat_plot <- params %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, mu_max, ec50, slope), loglogis))

#Plot the function output for parameter combinations
ggplot(dat_plot, aes(x = driver, y = y, group = mu_max)) +
    geom_line(alpha = 0.05) +
    scale_x_continuous(
        limits = c(0, 1000),
        trans = scales::pseudo_log_trans(base = 10, sigma = 0.1),
        breaks = c(0, 1, 10, 100, 1000),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
        title = "Toxin curves simulated from priors",
        x = expression(Toxin~concentration~(mu*g~{L^-1})),
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
        "toxin_SI_1000_simulated_curves.pdf"),
    height = 24.62 * (2/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)

