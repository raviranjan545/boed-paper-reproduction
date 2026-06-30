# Declare present directory
here::i_am("1_scripts/temperature_7_visualising_design_comparison.R")

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
library(GGally)
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
        "1000_simulated_norberg_curve_parameter_combinations_original_scale.csv")
) %>%
    mutate(curve_id = row_number())

## Central curve parameter values
central_curve_params <- read.csv(
    here(
        "2_designs_and_other_simulation_inputs",
        "norberg_central_curve_mbd_parameters.csv")
) 

## Experimental designs

# Optimal
opt_temp5_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "temperature_SIG_optimal_design_5_points.rds")),
                         sort(final.d[[besti]]))
opt_temp7_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                      "temperature_SIG_optimal_design_7_points.rds")),
                         sort(final.d[[besti]]))
opt_temp10_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                       "temperature_SIG_optimal_design_10_points.rds")),
                          sort(final.d[[besti]]))
opt_temp15_levels <- with(readRDS(here("2_designs_and_other_simulation_inputs", 
                                       "temperature_SIG_optimal_design_15_points.rds")),
                          sort(final.d[[besti]]))

# Uniform
unif_temp5_levels <- round(seq(5, 35, length.out = 5) * 2) / 2
unif_temp7_levels <- round(seq(5, 35, length.out = 7) * 2) / 2
unif_temp10_levels <- round(seq(5, 35, length.out = 10) * 2) / 2
unif_temp15_levels <- round(seq(5, 35, length.out = 15) * 2) / 2

# Joining optimal and uniform designs into one data frame
designs <- data.frame(levels = c(opt_temp5_levels, opt_temp7_levels, 
                                 opt_temp10_levels, opt_temp15_levels,
                                 unif_temp5_levels, unif_temp7_levels, 
                                 unif_temp10_levels, unif_temp15_levels)) %>%
    mutate(design = c(rep('optimal', length(levels)/2), 
                      rep('uniform', length(levels)/2)),
           n_points = c(rep(5,5), rep(7,7), rep(10,10), rep(15,15),
                        rep(5,5), rep(7,7), rep(10,10), rep(15,15)))

# Removing objects no longer needed
rm(opt_temp5_levels, opt_temp7_levels, opt_temp10_levels, opt_temp15_levels,
   unif_temp5_levels, unif_temp7_levels, unif_temp10_levels, unif_temp15_levels)

## CRPS and ES scores 
all_scores <- read.csv(here("3_simulation_result_summaries", 
                            "norberg_fit_opt_unif_scoring_rules_comparison.csv")) %>%
    mutate(design = factor(design, levels = c('unif', 'opt')))

# Calculate percentage of simulations where optimal design is better for plotting
pct_opt_better <- all_scores %>%
    select(n_points, simulation, design, 
           CRPS_a, CRPS_b, CRPS_tmax, CRPS_tmin, ES_raw) %>%
    pivot_wider(
        names_from = design,
        values_from = c(CRPS_a, CRPS_b, CRPS_tmax, CRPS_tmin, ES_raw)
    ) %>%
    group_by(n_points) %>%
    summarise(
        a_opt_better = 100 * mean(CRPS_a_opt < CRPS_a_unif, na.rm = TRUE),
        b_opt_better = 100 * mean(CRPS_b_opt < CRPS_b_unif, na.rm = TRUE),
        tmax_opt_better = 100 * mean(CRPS_tmax_opt < CRPS_tmax_unif, na.rm = TRUE),
        tmin_opt_better = 100 * mean(CRPS_tmin_opt < CRPS_tmin_unif, na.rm = TRUE),
        ES_raw_opt_better = 100 * mean(ES_raw_opt < ES_raw_unif, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        label_a = paste0(round(a_opt_better), "%"),
        label_b = paste0(round(b_opt_better), "%"),
        label_tmax = paste0(round(tmax_opt_better), "%"),
        label_tmin = paste0(round(tmax_opt_better), "%"),
        label_ES_raw = paste0(round(ES_raw_opt_better), "%")
    )


## Aggregate prediction error distributions
agg_err <- read.csv(here("3_simulation_result_summaries", 
                         "norberg_fit_opt_unif_prediction_error_comparison.csv")) %>%
    select(design, n_points, simulation, mae_offdesign_mbd) %>%
    rename(agg_MAE = mae_offdesign_mbd) %>%
    mutate(design = factor(design, levels = c('unif', 'opt')))

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
                           "norberg_fit_opt_unif_pointwise_prediction_error_mbd.csv")) %>%
    group_by(design, n_points, env) %>%
    summarise(point_MAE = mean(abs_error),
              .groups = "drop_last") %>%
    mutate(design = factor(design, levels = c('unif', 'opt')))




###### PLOTS 

#### 1. Curves

# Temperature function (norberg)
norberg <- function(temp, a, b, tmax, tmin){
    gr <- (a*exp(b*temp))*(tmax - temp)*(temp - tmin)
    gr
}

# Set driver values 
driver_values <- seq(5, 35, 0.5)

# Generate central curve predictions
dat_mbd_curve <- tibble(
    driver = driver_values,
    a = central_curve_params$a,
    b = central_curve_params$b,
    tmax = central_curve_params$tmax,
    tmin = central_curve_params$tmin
) %>%
    mutate(y = purrr::pmap_dbl(list(driver,  a, b, tmax, tmin), 
                               norberg))

# Generate curve dataset for plotting
dat_plot <- params %>%
    filter(curve_id < 31) %>%
    filter(curve_id > 20) %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, a, b, tmax, tmin), 
                               norberg))

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
        limits = c(5, 35),
        breaks = seq(5, 35, 5),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
        title = "a) Temperature curves",
        x = expression(Temperature~(degree*C)),
        y = expression(Growth~rate~(day^{-1}))
    ) +
    scale_fill_manual(values = cols, name = NULL,
                      labels = c("Optimal", "Uniform")) +
    geom_hline(yintercept = 0, aes(colour = 'lightgrey')) + 
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
        limits = c(5, 35),
        breaks = seq(5, 35, 5),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
        breaks = row_layout$y0,
        labels = as.character(row_layout$n_points),
        limits = y_limits,
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        title = 'b) Designs',
        x = expression(Temperature~(degree*C)),
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

p_a <- ridgeline_diff(all_scores, 
                      CRPS_a, 
                      label_x = 0.16,
                      label_y_offset = 0.7,
                      labels = pct_opt_better,
                      label_col = 'label_a',
                      title = expression("c)"~Estimation~score~`difference,`~italic(a))) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(day^{-1}~degree*C^{-2})))

p_b <- ridgeline_diff(all_scores, 
                      CRPS_b, 
                      label_x = 0.16,
                      label_y_offset = 0.7,
                      labels = pct_opt_better,
                      label_col = 'label_b',
                      title = expression("d)"~Estimation~score~`difference,`~italic(b))) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(degree*C^{-1})),
         y = '')

p_tmin <- ridgeline_diff(all_scores, 
                         CRPS_tmin, 
                         label_x = 0.16,
                         label_y_offset = 0.7,
                         labels = pct_opt_better,
                         label_col = 'label_tmin',
                         title = expression("e)"~Estimation~score~`difference,`~italic(T)[min])) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(degree*C)))

p_tmax <- ridgeline_diff(all_scores, 
                         CRPS_tmax, 
                         label_x = 0.16,
                         label_y_offset = 0.7,
                         labels = pct_opt_better,
                         label_col = 'label_tmax',
                         title = expression("f)"~Estimation~score~`difference,`~italic(T)[max])) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(degree*C)),
         y = '')

p_ES <- ridgeline_diff(all_scores, 
                       ES_raw, 
                       label_x = 0.16,
                       label_y_offset = 0.7,
                       labels = pct_opt_better,
                       label_col = 'label_ES_raw',
                       title = expression("g)"~Estimation~score~`difference,`~all~parameters)) + 
    labs(x = expression(ES["opt"] - ES["unif"]~(dimensionless)))

agg_err <- mutate(agg_err, MAE = agg_MAE)
p_pred <- ridgeline_diff(agg_err, 
                         MAE, 
                         label_x = 0.16,
                         label_y_offset = 0.7,
                         labels = agg_opt_better,
                         label_col = 'label_perc',
                         title = expression("h)"~Prediction~error~difference)) + 
    labs(x = expression(MAE["opt"] - MAE["unif"]~(day^{-1})),
         y = '')

a3 <- p_a
a4 <- p_b
a5 <- p_tmin
a6 <- p_tmax
a7 <- p_ES
a8 <- p_pred


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


a1 + a2 + section_header + a3 + a4 + a5 + a6 + a7 + a8+
    plot_layout(
        design = "
AA
BB
CC
DE
FG
HI
",
        heights = c(0.8, 1.2, 0.32, 0.9, 0.9, 0.9, 0.9)
    )

ggsave(
    here(
        "4_figures", 
        "temperature_plot_grid.pdf"), 
    height = 24.62 * (4.9/4.9) * 1.4, 
    width = 18 * 1.4,
    units = 'cm'
)



###### SUPPLEMENTARY FIGURES

#### S1. Priors 

# Target marginals on ORIGINAL scale
mean_a <- 0.002; sd_a <- 0.0002
mean_b <- 0.025; sd_b <- 0.005
mean_tmax <- 35;  sd_tmax <- 2
mean_tmin <- 2;   sd_tmin <- 3

# Convert to lognormal params for a,b (sdlog first, then meanlog)
sdlog_a   <- sqrt(log(1 + (sd_a/mean_a)^2))
meanlog_a <- log(mean_a) - 0.5 * sdlog_a^2
sdlog_b   <- sqrt(log(1 + (sd_b/mean_b)^2))
meanlog_b <- log(mean_b) - 0.5 * sdlog_b^2

# Correlations on (log a, log b, tmax, tmin)
rho_ab     <- -0.6
rho_a_tmax <- -0.5
rho_a_tmin <-  0.3

C <- diag(4)
dimnames(C) <- list(c("loga","logb","tmax","tmin"), c("loga","logb","tmax","tmin"))
C["loga","logb"] <- C["logb","loga"] <- rho_ab
C["loga","tmax"] <- C["tmax","loga"] <- rho_a_tmax
C["loga","tmin"] <- C["tmin","loga"] <- rho_a_tmin
C <- as.matrix(Matrix::nearPD(C, corr = TRUE)$mat)  # ensure PSD

# MVN on (log a, log b, tmax, tmin)
mu_vec <- c(meanlog_a, meanlog_b, mean_tmax, mean_tmin)
sd_vec <- c(sdlog_a,   sdlog_b,   sd_tmax,  sd_tmin)
Sigma  <- diag(sd_vec) %*% C %*% diag(sd_vec)
L      <- t(chol(Sigma))

# Generate parameter values for supplementary plot
set.seed(93785)
draws_si <- MASS::mvrnorm(100000, mu = mu_vec, Sigma = Sigma)

params_si <- data.frame(
    log_a = draws_si[,1],
    log_b = draws_si[,2],
    Tmax = draws_si[,3],
    Tmin = draws_si[,4]
)

my_filled_density <- function(data, mapping, ..., bins = 6, alpha = 1) {
    ggplot(data = data, mapping = mapping) +
        stat_density_2d_filled(
            aes(fill = after_stat(level)),
            contour_var = "ndensity",
            bins = bins,
            alpha = alpha
        ) +
        guides(fill = "none")
}


my_cor <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    r <- cor(x, y, use = "pairwise.complete.obs")
    r <- round(r, 2)
    if (r == 0) r <- 0  # removes -0.0
    
    ggplot(data = data, mapping = mapping) +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("italic(r) == ", r),
                 parse = TRUE,
                 size = 4) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()
        )
}


ggpairs(
    params_si,
    columnLabels = c("log~a~(dimensionless)",
                     "log~b~(dimensionless)", 
                     "italic(T)[min]~(degree*C)", 
                     "italic(T)[max]~(degree*C)"),
    labeller = label_parsed,
    diag = list(
        continuous = wrap(
            "densityDiag",
            fill = "blue",
            colour = NA,
            alpha = 0.4,
            linewidth = 0.25
        )
    ),
    lower = list(
        continuous = wrap(
            my_filled_density,
            bins = 8,
            alpha = 1
        ),
        combo = "box_no_facet"
    ),
    upper = list(
        continuous = my_cor
    )
) + 
    theme(panel.grid.minor = element_blank())


ggsave(
    here(
        "4_figures", 
        "temperature_SI_priors.pdf"),
    height = 18 * 1, 
    width = 18 * 1,
    units = 'cm'
)

#### S2. CRPS and ES on log scale (Violin plots)
cols_to_mean <- c("CRPS_a", "CRPS_b", "CRPS_tmin", "CRPS_tmax", "ES_raw")

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

# A. a
s1 <- all_scores %>%
    ggplot(., aes(CRPS_a, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_a, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 0.00085,
            y = factor(n_points),
            label = label_a
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab("# of experimental units") + 
    xlab(expression(CRPS~(day^{-1}~degree*C^{-2}))) +
    ggtitle(expression(paste('a) Estimation score, ', italic(a))),
            subtitle = 'lower is better') + 
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

# B. b
s2 <- all_scores %>%
    ggplot(., aes(CRPS_b, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_b, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 0.025,
            y = factor(n_points),
            label = label_b
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab('') + 
    xlab(expression(CRPS~(degree*C^-1))) +
    ggtitle(expression(paste('b) Estimation score, ', italic(b))),
            subtitle = 'lower is better') + 
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

# C. Tmin
s3 <- all_scores %>%
    ggplot(., aes(CRPS_tmin, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_tmin, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 7,
            y = factor(n_points),
            label = label_tmin
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab('# of experimental units') +
    xlab(expression(CRPS~(degree*C))) + 
    ggtitle(expression(paste('c) Estimation score, ', italic(T)[min])),
            subtitle = 'lower is better') + 
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

# D. Tmax
s4 <- all_scores %>%
    ggplot(., aes(CRPS_tmax, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_CRPS_tmax, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 4,
            y = factor(n_points),
            label = label_tmax
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab('') +
    xlab(expression(CRPS~(degree*C))) + 
    ggtitle(expression(paste('d) Estimation score, ', italic(T)[max])),
            subtitle = 'lower is better') + 
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

# F. ES - a and b
s5 <-
    all_scores %>%
    ggplot(., aes(ES_raw, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score,
        aes(x = mean_ES_raw, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = pct_opt_better,
        aes(
            x = 5.5,
            y = factor(n_points),
            label = label_ES_raw
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() + 
    ylab('# of experimental units') +
    xlab(expression(Energy~Score~(dimensionless))) + 
    ggtitle('e) Estimation score, all parameters',
            subtitle = 'lower is better') + 
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



# F. Prediction error
s6 <-
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
    scale_fill_manual(values = c('orange', 'purple')) + 
    geom_text(
        data = agg_opt_better,
        aes(
            x = 0.24,
            y = factor(n_points),
            label = label_perc
        ),
        inherit.aes = FALSE,
        size = 4
    ) +
    scale_x_log10() +
    ylab('') + 
    xlab(expression(Mean~Absolute~Error~(day^{-1}))) +
    ggtitle('f) Prediction error',
            subtitle = 'lower is better') + 
    guides(fill = "none") +
    theme( 
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 12, hjust = 0.08, vjust = 1.8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    )
s6

s1 + s2 + s3 + s4 + s5 + s6 + 
    plot_layout(
        design = "
AB
CD
EF
",
        heights = c(1, 1, 1, 1, 1, 1)
    )

ggsave(
    here(
        "4_figures", 
        "temperature_SI_log_scale.pdf"), 
    height = 24.62 * (3/3) * 1.4, 
    width = 18 * 1.4,
    units = 'cm'
)

#### S3. Pointwise prediction error

plot_mae_pointwise <-  point_err %>%
    mutate(design = recode(design, opt = "Optimal", unif = "Uniform"),
           design = factor(design, levels = c("Optimal", "Uniform"))) %>%
    ggplot(aes(x = env, y = point_MAE, colour = design)) +
    geom_line(linewidth = 1.5) +
    facet_wrap(~forcats::fct_rev(as.factor(n_points)),
               ncol = 1, strip.position = 'right')  +
    labs(
        x = expression(Temperature~(degree*C)),
        y = expression(Pointwise~Mean~Absolute~Error~(day^{-1})),
        colour = "Design"
    ) + 
    scale_colour_manual(values = c(Optimal = "purple", Uniform = "orange")) +
    scale_x_continuous(
        limits = c(5, 35),
        breaks = seq(5, 35, 5),
        expand = expansion(
            mult = c(0.01, 0.01))
    ) +
    guides(colour = guide_legend(
        keywidth       = unit(4, "lines"),
        keyheight      = unit(1, "lines")
    )) +
    theme(
        strip.text.y.right = element_text(angle = 0),
        legend.position = 'top',
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
        "temperature_SI_pointwise_error.pdf"), 
    plot_mae_pointwise,
    height = 24.62 * (3/4.9) * 1.4, 
    width = 18 * 1.4,
    units = 'cm'
)



##### S4. All simulated curves

# Plot curves from simulated parameter combinations
dat_plot <- params %>%
    mutate(driver = list(driver_values)) %>%
    tidyr::unnest(driver) %>%
    mutate(y = purrr::pmap_dbl(list(driver, a, b, tmax, tmin), norberg))

#Plot the function output for parameter combinations
ggplot(dat_plot, aes(x = driver, y = y, group = a)) +
    geom_line(alpha = 0.05) +
    scale_x_continuous(
        limits = c(5, 35),
        breaks = seq(5, 35, 5),
        expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
        title = "Temperature curves simulated from priors",
        x = expression(Temperature~(degree*C)),
        y = expression(Growth~rate~(day^{-1}))
    ) +
    theme(
        strip.text.y.right = element_text(angle = 0),
        legend.position = 'top',
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
        '4_figures',
        'temperature_SI_1000_simulated_curves.pdf'),
    height = 24.62 * (2/4.9) * 1.4, 
    width = 18 * 1.4,
    units = 'cm'
)

