# Install needed packages if not present already
local({
    packages <- c(
        "acebayes",
        "brms",
        "conflicted",
        "dplyr",
        "here",
        "GGally",
        "ggplot2",
        "ggdist",
        "ggfortify",
        "ggridges",
        "MASS",
        "Matrix",
        "parallel",
        "patchwork",
        "posterior",
        "purrr",
        "Rcpp",
        "rlang",
        "roahd",
        "rstan",
        "scales",
        "scoringRules",
        "this.path",
        "tidyr"
    )
    
    missing <- packages[
        !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
    ]
    
    if (length(missing)) {
        install.packages(
            missing,
            repos = "https://cloud.r-project.org",
            quiet = TRUE
        )
    }
})

setwd(this.path::this.dir())   # move into 1_scripts, where this file lives
here::i_am("1_scripts/00_support_file_create_directories_and_install_packages.R")

# Directories to be created
dirs <- c(
    "2_designs_and_other_simulation_inputs",
    "3_simulation_result_summaries",
    "4_figures",
    "nutrient_curves_posteriors",
    "nutrient_curves_bad_priors_posteriors",
    "light_curves_posteriors",
    "temperature_curves_posteriors",
    "toxin_curves_posteriors",
    "nutrient_curves_simulations_and_fits",
    "nutrient_curves_bad_priors_simulations_and_fits",
    "light_curves_simulations_and_fits",
    "temperature_curves_simulations_and_fits",
    "toxin_curves_simulations_and_fits"
)

invisible(lapply(
    dirs,
    \(x) dir.create(here::here(x), recursive = TRUE, showWarnings = FALSE)
))

rm(dirs)
