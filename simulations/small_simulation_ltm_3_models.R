rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)

fit_all_three_models <- function(data, model, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, debug = TRUE, force = FALSE,
                                 path_prefix = "") {

  stan_data_orig <- data_2_stan(data, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  stan_data_skew <- data_2_stan(data, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  path_orig <- file.path("fitted_objects", sprintf("%s3_models_orig_thresholds.rds", path_prefix))
  path_skew <- file.path("fitted_objects", sprintf("%s3_models_skew_thresholds.rds", path_prefix))
  path_free <- file.path("fitted_objects", sprintf("%s3_models_free_thresholds.rds", path_prefix))

  cat("Fitting original threshold model\n")
  fit_orig <- save_or_run_model(
    model$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples),
    path_orig, force
  )
  cat("Fitting skew threshold model\n")
  fit_skew <- save_or_run_model(
    model$variational(data = stan_data_skew, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples),
    path_skew, force
  )
  cat("Fitting free threshold model\n")
  fit_free <- save_or_run_model(
    model$variational(data = stan_data_free, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples),
    path_free, force
  )

  return(list(
    fit_orig = fit_orig,
    fit_skew = fit_skew,
    fit_free = fit_free
  ))
}

get_performance_tib <- function(fit, dat, threshold_name, nc, nr) {
  get_means_tib(fit) |>
    add_implied_thresholds(nc, nr) |>
    add_true_values_to_tib(dat) |>
    mutate(parameter = if_else(parameter == "free_thresholds", threshold_name, parameter))
}

get_scatterplot_retrieval <- function(fit, tib) {

  idx <- grep("thresholds", tib$parameter, fixed = TRUE)[1]
  threshold_nm <- tib$parameter[idx]
  tib |>
    scatterplot_retrieval() +
    ggtitle(sprintf("Variational inference for %s took %.3f seconds", threshold_nm, fit$time())) +
    labs(x = "True value", y = "Variational Estimate")
}

np <- 50             # no patients
ni <- 5              # no items
nr <- 10             # no raters
nc <- 5              # no response categories
no_rater_groups <- 1 # no groups of rater

mod_ltm <- compile_stan_model("stanmodels/LTM_3_models.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")

set.seed(1234)
dat_orig <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_skew_thresholds = FALSE, use_free_thresholds = FALSE)
dat_skew <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_skew_thresholds = TRUE,  use_free_thresholds = FALSE)
dat_free <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_skew_thresholds = FALSE, use_free_thresholds = TRUE)

fits_orig <- fit_all_three_models(dat_orig, mod_ltm, path_prefix = "data_orig")
system("beep_finished.sh 0")

perf_orig_dat_orig_tib <- get_performance_tib(fits_orig$fit_orig, dat_orig, "2-parameter_thresholds", nc, nr)
perf_skew_dat_orig_tib <- get_performance_tib(fits_orig$fit_skew, dat_orig, "skew_thresholds", nc, nr)
perf_free_dat_orig_tib <- get_performance_tib(fits_orig$fit_free, dat_orig, "free_thresholds", nc, nr)

g1 <- get_scatterplot_retrieval(fits_orig$fit_orig, perf_orig_dat_orig_tib)
g2 <- get_scatterplot_retrieval(fits_orig$fit_skew, perf_skew_dat_orig_tib)
g3 <- get_scatterplot_retrieval(fits_orig$fit_free, perf_free_dat_orig_tib)

ff <- gridExtra::grid.arrange(g1, g2, g3, ncol = 3)
gg <- cowplot::plot_grid(g1, g2, g3, nrow = 1)

ggplot2::benchplot(gg)
system.time(grid::grid.draw(ff))



# skew model
fits_skew <- fit_all_three_models(dat_skew, mod_ltm, path_prefix = "data_skew")
system("beep_finished.sh 0")

perf_orig_dat_skew_tib <- get_performance_tib(fits_skew$fit_orig, dat_skew, "2-parameter_thresholds", nc, nr)
perf_skew_dat_skew_tib <- get_performance_tib(fits_skew$fit_skew, dat_skew, "skew_thresholds", nc, nr)
perf_free_dat_skew_tib <- get_performance_tib(fits_skew$fit_free, dat_skew, "free_thresholds", nc, nr)

get_scatterplot_retrieval(fits_skew$fit_orig, perf_orig_dat_skew_tib)
get_scatterplot_retrieval(fits_skew$fit_skew, perf_skew_dat_skew_tib)
get_scatterplot_retrieval(fits_skew$fit_free, perf_free_dat_skew_tib)

# free model
fits_free <- fit_all_three_models(dat_free, mod_ltm, path_prefix = "data_free")
system("beep_finished.sh 0")

perf_orig_dat_free_tib <- get_performance_tib(fits_free$fit_orig, dat_free, "2-parameter_thresholds", nc, nr)
perf_skew_dat_free_tib <- get_performance_tib(fits_free$fit_skew, dat_free, "skew_thresholds", nc, nr)
perf_free_dat_free_tib <- get_performance_tib(fits_free$fit_free, dat_free, "free_thresholds", nc, nr)

get_scatterplot_retrieval(fits_free$fit_orig, perf_orig_dat_free_tib)
get_scatterplot_retrieval(fits_free$fit_skew, perf_skew_dat_free_tib)
get_scatterplot_retrieval(fits_free$fit_free, perf_free_dat_free_tib)
