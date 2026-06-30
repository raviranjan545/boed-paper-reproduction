# ============================================================================
# ridgeline_diff(): paired (opt - unif) difference ridgeline for one parameter
# ----------------------------------------------------------------------------
#   data           : long data frame with columns n_points, design, simulation,
#                    and one column per parameter (e.g. CRPS_mumax, ES, ...)
#   param          : unquoted or quoted column name to plot the difference of
#   clip_probs     : lower/upper quantiles for the clip range (default 2/98%)
#   fill_cols      : c(low, mid, high) for the diverging gradient
#                    (low = optimal-better side, high = uniform-better side)
#   scale_ridge    : ggridges `scale` (1 = zero overlap)
#   title          : plot title (expression or string); NULL = no title
#   labels         : data frame with n_points + one label column per parameter
#   label_col      : REQUIRED if labels given: which label column to use
#   label_x        : horizontal label position, fraction of x-range (0=L, 1=R)
#   label_y_offset : vertical label offset above each ridge baseline
#   label_size     : label text size
#
# The x-axis title is intentionally left to the default; set it after the
# call with + labs(x = ...).
# ============================================================================

# If running separately, load library
# library(ggridges)

ridgeline_diff <- function(data,
                           param,
                           clip_probs     = c(.02, .98),
                           fill_cols      = c("purple", "grey92", "orange"),
                           scale_ridge    = 1,
                           title          = NULL,
                           # --- right-side text labels (one per n_points) ---
                           labels         = NULL,        # data frame with n_points + label column(s)
                           label_col      = NULL,        # REQUIRED if labels given: name of the label-text column to use
                           label_x        = 0.9,         # horizontal position, fraction of visible x-range (0=left, 1=right)
                           label_y_offset = 0.5,         # vertical offset above each ridge baseline, in ridge-height units
                           label_size     = 4) {         # text size
    
    param <- rlang::ensym(param)        # accept mumax or "CRPS_mumax" style
    pname <- rlang::as_string(param)
    
    # --- paired differences: opt - unif, per (n_points, simulation) ---
    diff_df <- data %>%
        select(n_points, design, simulation, value = !!param) %>%
        pivot_wider(names_from = design, values_from = value) %>%
        mutate(diff = opt - unif,
               n_points = factor(n_points, levels = sort(unique(n_points))))
    
    n_na <- sum(is.na(diff_df$diff))
    if (n_na > 0)
        warning(sprintf("%s: %d NA differences (check for missing/duplicate design cells)",
                        pname, n_na))
    diff_df <- diff_df %>% filter(!is.na(diff))
    
    # --- clip range derived from THIS parameter's own distribution ---
    qlim <- quantile(diff_df$diff, clip_probs, names = FALSE)
    
    # --- per-level mean (pre-clip), clamped to the visible range for plotting ---
    mean_df <- diff_df %>%
        group_by(n_points) %>%
        summarise(mean_diff = mean(diff), .groups = "drop") %>%
        mutate(mean_plot = pmin(pmax(mean_diff, qlim[1]), qlim[2]),
               y_num = as.integer(n_points))
    
    # Build the ridges layer FIRST so we can extract the actual rendered
    # heights that ggridges drew (it uses a joint bandwidth and a global-max
    # normalisation that aren't trivial to replicate). Then we add the mean
    # segments using those exact heights, guaranteeing the line tops snap to
    # the curve regardless of bandwidth choices, grid resolution, or scale.
    p_base <- ggplot(diff_df, aes(diff, n_points)) +
        geom_density_ridges_gradient(
            aes(fill = after_stat(x)),
            scale = scale_ridge, rel_min_height = 0.004,
            colour = "white", linewidth = 0.25) +
        scale_x_continuous(limits = qlim, oob = scales::oob_squish,
                           expand = expansion(mult = c(.01, .01)))
    
    ridges_built <- suppressMessages(
        suppressWarnings(ggplot2::ggplot_build(p_base))
    )$data[[1]]
    mean_df <- mean_df %>%
        rowwise() %>%
        mutate(
            y_top = {
                rows <- ridges_built[ridges_built$ymin == y_num, ]
                rows <- rows[!duplicated(rows$x), ]
                rows <- rows[order(rows$x), ]
                if (nrow(rows) == 0 || mean_plot < min(rows$x) || mean_plot > max(rows$x)) {
                    y_num
                } else {
                    stats::approx(rows$x, rows$ymax, xout = mean_plot)$y
                }
            }
        ) %>%
        ungroup()
    
    p <- p_base +
        geom_vline(xintercept = 0, linetype = "solid",
                   colour = "black", linewidth = 0.8) +
        geom_segment(
            data = mean_df,
            aes(x = mean_plot, xend = mean_plot,
                y = y_num, yend = y_top),
            colour = "red", linewidth = 1, inherit.aes = FALSE
        ) +
        scale_y_discrete(expand = expansion(mult = c(0.04, .18))) +
        scale_fill_gradient2(low = fill_cols[1], mid = fill_cols[2],
                             high = fill_cols[3], midpoint = 0, guide = "none") +
        labs(title = title, y = "# of experimental units") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(size = 15),
              plot.subtitle = element_text(size = 10, colour = "grey35"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_text(size = 13.5),
              axis.title.y = element_text(size = 13.5))
    
    # --- optional right-side text labels, one per n_points ---
    if (!is.null(labels)) {
        if (!"n_points" %in% names(labels))
            stop("`labels` must contain an `n_points` column to join on.")
        if (is.null(label_col))
            stop("`label_col` must be specified when `labels` is given (name the label-text column for this parameter, e.g. label_col = \"label_mumax\").")
        if (!label_col %in% names(labels))
            stop(sprintf("`labels` has no column '%s'. Available columns: %s",
                         label_col, paste(setdiff(names(labels), "n_points"), collapse = ", ")))
        
        # horizontal position: a fraction of the visible x-range, so it adapts
        # to each parameter's own clip range.
        x_at <- qlim[1] + label_x * (qlim[2] - qlim[1])
        
        # JOIN BY VALUE, not position. Match the label's n_points to the factor
        # levels actually plotted; warn about any that don't line up.
        lvl <- levels(diff_df$n_points)
        lab_df <- labels %>%
            mutate(n_points_chr = as.character(n_points)) %>%
            filter(n_points_chr %in% lvl) %>%
            mutate(y_num = match(n_points_chr, lvl),
                   x = x_at,
                   y = y_num + label_y_offset,
                   label = .data[[label_col]])
        
        missing_lvls <- setdiff(lvl, as.character(labels$n_points))
        if (length(missing_lvls) > 0)
            warning(sprintf("No label provided for n_points: %s",
                            paste(missing_lvls, collapse = ", ")))
        extra_lbls <- setdiff(as.character(labels$n_points), lvl)
        if (length(extra_lbls) > 0)
            warning(sprintf("Labels given for n_points not in the data (ignored): %s",
                            paste(extra_lbls, collapse = ", ")))
        
        p <- p +
            geom_text(data = lab_df,
                      aes(x = x, y = y, label = label),
                      inherit.aes = FALSE, size = label_size)
    }
    
    p
}