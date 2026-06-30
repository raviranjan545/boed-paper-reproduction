# Declare present directory
here::i_am("1_scripts/nutrients_11_bad_priors_visualising_design_comparison.R")

# Load libraries
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggridges)
theme_set(theme_minimal(base_size = 14))

###### Read in files

## Source plotting function
source(here(
    "1_scripts", 
    "00_support_file_ridgeline_plot_function.R")
)

## CRPS and ES scores - bad priors
badprior_scores <- read.csv(here("3_simulation_result_summaries",
                                 "monod_badprior_fit_opt_unif_scoring_rules_comparison.csv"), 
) %>%
    mutate(design = factor(design, levels = c("unif", "opt")))

badprior_pct_opt_better <- badprior_scores %>%
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

###### PLOTS 

#### Bad prior CRPS and ES ridgeline plots

b1 <- ridgeline_diff(badprior_scores, 
                     CRPS_mumax, 
                     label_x = 0.2,
                     label_y_offset = 0.7,
                     labels = badprior_pct_opt_better,
                     label_col = "label_mumax",
                     title = expression("a)"~Estimation~score~`difference,`~mu[max])) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(day^{-1})))

b2 <- ridgeline_diff(badprior_scores, 
                     CRPS_k, 
                     label_x = 0.2,
                     label_y_offset = 0.7,
                     labels = badprior_pct_opt_better,
                     label_col = "label_k",
                     title = expression("b)"~Estimation~score~`difference,`~K)) + 
    labs(x = expression(CRPS["opt"] - CRPS["unif"]~(mu*M)),
         y = "")

b3 <- ridgeline_diff(badprior_scores, 
                     ES, 
                     label_x = 0.2,
                     label_y_offset = 0.7,
                     labels = badprior_pct_opt_better,
                     label_col = "label_ES",
                     title = expression("c)"~Estimation~score~`difference,`~all~parameters)) + 
    labs(x = expression(ES["opt"] - ES["unif"]~(dimensionless)))


b1 + b2 + b3 +
    plot_layout(
        design = "
AB
C#
",
        heights = c(0.9, 0.9)
    )


ggsave(
    here(
        "4_figures", 
        "nutrients_bad_priors_plot_grid.pdf"), 
    height = 24.62 * (1.8/4.9) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)



###### SUPPLEMENTARY FIGURES

#### S4. Bad priors violin plots

cols_to_mean <- c("CRPS_mumax", "CRPS_k", "ES")

mean_score_badprior <- badprior_scores %>%
    group_by(n_points, design) %>%
    summarise(
        across(
            all_of(cols_to_mean),
            ~ mean(.x, na.rm = TRUE),
            .names = "mean_{.col}"
        ),
        .groups = "drop"
    )

# A. Max growth
sb1 <- badprior_scores %>%
    ggplot(., aes(CRPS_mumax, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score_badprior,
        aes(x = mean_CRPS_mumax, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = badprior_pct_opt_better,
        aes(
            x = 0.4,
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
sb1

# B. Ks
sb2 <- badprior_scores %>%
    ggplot(., aes(CRPS_k, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score_badprior,
        aes(x = mean_CRPS_k, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +  
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = badprior_pct_opt_better,
        aes(
            x = 5.5,
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
sb2

# C. ES
sb3 <- badprior_scores %>%
    ggplot(., aes(ES, factor(n_points), group = interaction(design, n_points), 
                  fill = design)) + 
    geom_violin(quantile.linetype = 1,
                quantiles = 0.5) + 
    geom_point(
        data = mean_score_badprior,
        aes(x = mean_ES, y = factor(n_points)),
        size = 4, colour = "red",
        position = position_dodge2(width = 0.9)
    ) +  
    scale_fill_manual(values = c("orange", "purple")) + 
    geom_text(
        data = badprior_pct_opt_better,
        aes(
            x = 6,
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
sb3


sb1 + sb2 + sb3 + 
    plot_layout(
        design = "
AB
CD
",
        heights = c(1, 1, 1, 1),
    )


ggsave(here(
    "4_figures", 
    "nutrients_SI_bad_priors_log_scale.pdf"), 
    height = 24.62 * (2/3) * 1.4, 
    width = 18 * 1.4,
    units = "cm"
)

