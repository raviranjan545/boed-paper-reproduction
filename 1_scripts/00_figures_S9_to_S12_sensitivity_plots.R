# Declare present directory and create directory structure if not done
here::i_am("1_scripts/00_figures_S9_to_S12_sensitivity_plots.R")
source(here::here("1_scripts", "00_support_file_create_directories_and_install_packages.R"))

# Load libraries
library(here)
library(ggplot2)
library(patchwork)

# Sensitivity analysis figures

run_sensitivity_figures <- function(project_dir = NULL) {
    # If you use here, this will usually resolve to the project root.
    # Otherwise, edit project_dir manually when calling this function.
    if (is.null(project_dir)) {
        if (requireNamespace("here", quietly = TRUE)) {
            project_dir <- here::here()
        } else {
            project_dir <- getwd()
        }
    }
    
    base_path  <- file.path(project_dir, "4_figures")
    params_dir <- file.path(project_dir, "2_designs_and_other_simulation_inputs")
    dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
    
    read_params <- function(filename, n_params) {
        path <- file.path(params_dir, filename)
        dat <- utils::read.csv(path, check.names = FALSE)
        vals <- suppressWarnings(as.numeric(dat[1, seq_len(n_params)]))
        
        if (anyNA(vals)) {
            stop("Could not read numeric parameter values from: ", path, call. = FALSE)
        }
        
        vals
    }
    
    make_curve <- function(fun, x_min, x_max, n = 1001) {
        x <- seq(x_min, x_max, length.out = n)
        data.frame(x = x, sensitivity = fun(x))
    }
    
    sensitivity_theme <- function() {
        theme_bw(base_size = 16, base_family = "Helvetica") +
            theme(
                panel.grid = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
                axis.title = element_text(colour = "black"),
                axis.text = element_text(colour = "black"),
                plot.title = element_text(hjust = 0.5, colour = "black", size = 16),
                plot.background = element_rect(fill = "white", colour = NA),
                panel.background = element_rect(fill = "white", colour = NA),
                aspect.ratio = 1 / 1.61803398875,
                plot.margin = margin(8, 8, 8, 8)
            )
    }
    
    make_sensitivity_plot <- function(curve, title, x_label, x_limits,
                                      x_breaks = waiver(), x_labels = waiver(),
                                      x_trans = "identity") {
        ggplot(curve, aes(x = x, y = sensitivity)) +
            geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
            geom_line(linewidth = 0.7, colour = "black") +
            scale_x_continuous(
                limits = x_limits,
                breaks = x_breaks,
                labels = x_labels,
                trans = x_trans,
                expand = expansion(mult = 0.02)
            ) +
            scale_y_continuous(expand = expansion(mult = 0.08)) +
            labs(x = x_label, y = "Sensitivity", title = title) +
            sensitivity_theme()
    }
    
    save_pdf <- function(plot, filename, width, height) {
        ggplot2::ggsave(
            filename = file.path(base_path, filename),
            plot = plot,
            width = width,
            height = height,
            units = "cm",
            device = grDevices::pdf,
            bg = "white"
        )
    }
    
    # ---------------------------------------------------------------------------
    # 1. MONOD FUNCTION (NUTRIENTS)
    # ---------------------------------------------------------------------------
    monod_params <- read_params("monod_central_curve_mbd_parameters.csv", 2)
    mumax_m <- monod_params[1]
    k_m     <- monod_params[2]
    
    d_mumax_monod <- function(x) x / (x + k_m)
    d_k_monod     <- function(x) -(mumax_m * x) / (x + k_m)^2
    
    x_label_monod <- expression(paste("Nutrient concentration (", mu, "M)"))
    
    p_m1 <- make_sensitivity_plot(
        make_curve(d_mumax_monod, 0, 25),
        title = expression(mu[max]),
        x_label = x_label_monod,
        x_limits = c(0, 25)
    )
    
    p_m2 <- make_sensitivity_plot(
        make_curve(d_k_monod, 0, 25),
        title = "K",
        x_label = x_label_monod,
        x_limits = c(0, 25)
    )
    
    plot_monod <- (p_m1 + p_m2) +
        plot_annotation(
            title = "Sensitivity of Monod parameters",
            theme = theme(
                plot.title = element_text(size = 24, hjust = 0.5, family = "Helvetica"),
                plot.background = element_rect(fill = "white", colour = NA)
            )
        )
    
    save_pdf(plot_monod, 
             "nutrients_SI_sensitivity.pdf", 
             width = 11.2*2.54, 
             height = 4.2*2.54
    )
    
    # ---------------------------------------------------------------------------
    # 2. EILERS-PEETERS FUNCTION (LIGHT)
    # ---------------------------------------------------------------------------
    light_params <- read_params("eilerspeeters_central_curve_mbd_parameters.csv", 3)
    mumax_l <- light_params[1]
    alpha_l <- light_params[2]
    iopt_l  <- light_params[3]
    
    ep_denominator <- function(light) {
        (mumax_l / (alpha_l * iopt_l^2)) * light^2 +
            (1 - (2 * mumax_l) / (alpha_l * iopt_l)) * light +
            (mumax_l / alpha_l)
    }
    
    d_mumax_ep <- function(light) {
        q <- ep_denominator(light)
        dq_dmumax <- light^2 / (alpha_l * iopt_l^2) -
            2 * light / (alpha_l * iopt_l) +
            1 / alpha_l
        
        (light * q - mumax_l * light * dq_dmumax) / q^2
    }
    
    d_alpha_ep <- function(light) {
        q <- ep_denominator(light)
        dq_dalpha <- -mumax_l * light^2 / (alpha_l^2 * iopt_l^2) +
            2 * mumax_l * light / (alpha_l^2 * iopt_l) -
            mumax_l / alpha_l^2
        
        -(mumax_l * light * dq_dalpha) / q^2
    }
    
    d_iopt_ep <- function(light) {
        q <- ep_denominator(light)
        dq_diopt <- -2 * mumax_l * light^2 / (alpha_l * iopt_l^3) +
            2 * mumax_l * light / (alpha_l * iopt_l^2)
        
        -(mumax_l * light * dq_diopt) / q^2
    }
    
    x_label_light <- expression(paste("Light intensity (", mu, "E ", m^{-2}, " ", s^{-1}, ")"))
    
    p_l1 <- make_sensitivity_plot(
        make_curve(d_mumax_ep, 0, 1000),
        title = expression(mu[max]),
        x_label = x_label_light,
        x_limits = c(0, 1000)
    )
    
    p_l2 <- make_sensitivity_plot(
        make_curve(d_alpha_ep, 0, 1000),
        title = expression(alpha),
        x_label = x_label_light,
        x_limits = c(0, 1000)
    )
    
    p_l3 <- make_sensitivity_plot(
        make_curve(d_iopt_ep, 0, 1000),
        title = expression(I[opt]),
        x_label = x_label_light,
        x_limits = c(0, 1000)
    )
    
    plot_light <- (p_l1 + p_l2 + p_l3) +
        plot_layout(design = "AB\nC#") +
        plot_annotation(
            title = "Sensitivity of Eilers-Peeters parameters",
            theme = theme(
                plot.title = element_text(size = 24, hjust = 0.5, family = "Helvetica"),
                plot.background = element_rect(fill = "white", colour = NA)
            )
        )
    
    save_pdf(plot_light, 
             "light_SI_sensitivity.pdf", 
             width = 11.2*2.54, 
             height = 7.6*2.54
    )
    
    # ---------------------------------------------------------------------------
    # 3. NORBERG FUNCTION (TEMPERATURE)
    # ---------------------------------------------------------------------------
    temp_params <- read_params("norberg_central_curve_mbd_parameters.csv", 4)
    a_m    <- temp_params[1]
    b_m    <- temp_params[2]
    tmax_m <- temp_params[3]
    tmin_m <- temp_params[4]
    
    d_a <- function(temp) {
        exp(b_m * temp) * (tmax_m - temp) * (temp - tmin_m)
    }
    
    d_b <- function(temp) {
        a_m * temp * exp(b_m * temp) * (tmax_m - temp) * (temp - tmin_m)
    }
    
    d_tmax <- function(temp) {
        a_m * exp(b_m * temp) * (temp - tmin_m)
    }
    
    d_tmin <- function(temp) {
        -a_m * exp(b_m * temp) * (tmax_m - temp)
    }
    
    x_label_temp <- expression(paste("Temperature (", degree, "C)"))
    
    p_t1 <- make_sensitivity_plot(
        make_curve(d_a, 5, 34.9),
        title = "a",
        x_label = x_label_temp,
        x_limits = c(5, 34.9)
    )
    
    p_t2 <- make_sensitivity_plot(
        make_curve(d_b, 5, 34.9),
        title = "b",
        x_label = x_label_temp,
        x_limits = c(5, 34.9)
    )
    
    p_t3 <- make_sensitivity_plot(
        make_curve(d_tmax, 5, 34.9),
        title = expression(T[max]),
        x_label = x_label_temp,
        x_limits = c(5, 34.9)
    )
    
    p_t4 <- make_sensitivity_plot(
        make_curve(d_tmin, 5, 34.9),
        title = expression(T[min]),
        x_label = x_label_temp,
        x_limits = c(5, 34.9)
    )
    
    plot_temp <- (p_t1 + p_t2 + p_t4 + p_t3) +
        plot_layout(design = "AB\nCD") +
        plot_annotation(
            title = "Sensitivity of Norberg parameters",
            theme = theme(
                plot.title = element_text(size = 24, hjust = 0.5, family = "Helvetica"),
                plot.background = element_rect(fill = "white", colour = NA)
            )
        )
    
    save_pdf(plot_temp, 
             "temperature_SI_sensitivity.pdf", 
             width = 11.2*2.54, 
             height = 7.6*2.54
    )
    
    # ---------------------------------------------------------------------------
    # 4. LOG-LOGISTIC FUNCTION (TOXINS), WITH PSEUDO-LOG X SCALE
    # ---------------------------------------------------------------------------
    toxin_params <- read_params("toxin_central_curve_mbd_parameters.csv", 3)
    d_m <- toxin_params[1]
    e_m <- toxin_params[2]
    h_m <- toxin_params[3]
    
    dose_z <- function(x) (x / e_m)^h_m
    
    d_d <- function(x) {
        z <- dose_z(x)
        1 / (1 + z)
    }
    
    d_h <- function(x) {
        z <- dose_z(x)
        out <- -d_m * z * log(x / e_m) / (1 + z)^2
        out[x == 0] <- 0
        out
    }
    
    d_e <- function(x) {
        z <- dose_z(x)
        d_m * h_m * z / (e_m * (1 + z)^2)
    }
    
    x_label_toxin <- expression(paste("Toxin concentration (", mu, "g ", L^{-1}, ")"))
    pseudo_log_x <- scales::pseudo_log_trans(sigma = 0.1, base = 10)
    toxin_breaks <- c(0, 1, 10, 100, 1000)
    
    p_dr1 <- make_sensitivity_plot(
        make_curve(d_d, 0, 1000),
        title = expression(mu[max]),
        x_label = x_label_toxin,
        x_limits = c(0, 1000),
        x_breaks = toxin_breaks,
        x_labels = as.character(toxin_breaks),
        x_trans = pseudo_log_x
    )
    
    p_dr2 <- make_sensitivity_plot(
        make_curve(d_h, 0, 1000),
        title = "h",
        x_label = x_label_toxin,
        x_limits = c(0, 1000),
        x_breaks = toxin_breaks,
        x_labels = as.character(toxin_breaks),
        x_trans = pseudo_log_x
    )
    
    p_dr3 <- make_sensitivity_plot(
        make_curve(d_e, 0, 1000),
        title = "e",
        x_label = x_label_toxin,
        x_limits = c(0, 1000),
        x_breaks = toxin_breaks,
        x_labels = as.character(toxin_breaks),
        x_trans = pseudo_log_x
    )
    
    plot_dr <- (p_dr1 + p_dr3 + p_dr2) +
        plot_layout(design = "AB\nC#") +
        plot_annotation(
            title = "Sensitivity of dose-response parameters",
            theme = theme(
                plot.title = element_text(size = 24, hjust = 0.5, family = "Helvetica"),
                plot.background = element_rect(fill = "white", colour = NA)
            )
        )
    
    save_pdf(plot_dr, 
             "toxin_SI_sensitivity.pdf", 
             width = 11.2*2.54, 
             height = 7.6*2.54
    )
    
    invisible(list(
        monod = plot_monod,
        light = plot_light,
        temperature = plot_temp,
        dose_response = plot_dr
    ))
}

# Run to generate and save figures
run_sensitivity_figures(project_dir = NULL)
