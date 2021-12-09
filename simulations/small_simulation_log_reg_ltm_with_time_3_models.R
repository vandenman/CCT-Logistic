rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(gridExtra)
library(cmdstanr)
library(dplyr)
library(purrr)

# functions ----
fit_all_three_models <- function(data, model, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, debug = TRUE, force = FALSE,
                                 path_prefix = "", store_predictions = TRUE) {

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  stan_data_orig <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  stan_data_skew <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

  path_orig <- file.path("fitted_objects", sprintf("%slog_reg_3_models_orig_thresholds.rds", path_prefix))
  path_skew <- file.path("fitted_objects", sprintf("%slog_reg_3_models_skew_thresholds.rds", path_prefix))
  path_free <- file.path("fitted_objects", sprintf("%slog_reg_3_models_free_thresholds.rds", path_prefix))

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

split_slopes_on_covariate_vs_lt <- function(tib, no_covariates) {

  idx_intercept <- which(tib$parameter == "log_reg_intercept")
  idx_slopes    <- which(tib$parameter == "log_reg_slopes")

  tib$parameter[idx_intercept]       <- "log_reg_slopes (covariates)"
  tib$parameter[idx_slopes[1:no_covariates]]    <- "log_reg_slopes (covariates)"
  tib$parameter[idx_slopes[-(1:no_covariates)]] <- "log_reg_slopes (latent truth)"
  tib
}


get_performance_tib <- function(fit, dat, threshold_name, nc, nr, no_covariates) {
  nc            <- dat$ltm_data$nc
  nr            <- dat$ltm_data$nr
  no_covariates <- ncol(dat$log_reg$design_matrix) - dat$ltm_data$ni
  get_means_tib(fit) |>
    add_implied_thresholds(nc, nr) |>
    add_true_values_to_tib(dat) |>
    split_slopes_on_covariate_vs_lt(no_covariates) |>
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

get_predictions_plot <- function(fit, name, data) {

  probs <- fit$draws("log_reg_predictions", format = "draws_matrix")
  mean_probs <- colMeans(unname(probs))

  tib <- tibble(
    x1 = factor(dat_log_reg_orig$log_reg$y),
    x2 = c(dat_log_reg_orig$log_reg$true_probs),
    y  = mean_probs
  )

  g_left <- ggplot(tib, aes(x = x1, y = y)) +
    geom_point() +
    labs(x = "Observed outcome", y = "Posterior mean") +
    ggtitle(name) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()

  g_right <- ggplot(tib, aes(x = x2, y = y)) +
   jaspGraphs::geom_abline2(slope = 1, intercept = 0, color = "gray") +
    geom_point() +
    labs(x = "True probability", y = "Posterior mean") +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()

  list(g_left, g_right)

}

fit_single_dataset <- function(dataset, model, data_name) {

  path_prefix <- paste0("data_", data_name)
  fits <- fit_all_three_models(dataset, model, path_prefix = path_prefix)
  perf_tib <- purrr::map2(fits, c("2-parameter_thresholds", "skew_thresholds", "free_thresholds"), get_performance_tib, dat = dataset)

  graphs_retrieval <- map2(fits, perf_tib, get_scatterplot_retrieval)
  graphs_predictions <- purrr::map2(fits, names(fits), get_predictions_plot, data = dataset)

  # TODO: these should not plot directly
  graphs_retrieval_joined   <- grid.arrange(grobs = graphs_retrieval, ncol = 3)
  graphs_predictions_joined <- grid.arrange(grobs = matrix(unlist(graphs_predictions, recursive = FALSE), 3, 2, TRUE), nrow = 2, ncol = 3)

  return(list(
    fits                      = fits,
    perf_tib                  = perf_tib,
    graphs_retrieval          = graphs_retrieval,
    graphs_predictions        = graphs_predictions,
    graphs_retrieval_joined   = graphs_retrieval_joined,
    graphs_predictions_joined = graphs_predictions_joined
  ))
}

# simulate all data sets ----
np <- 50             # no patients
ni <- 5              # no items
nr <- 10             # no raters
nc <- 5              # no response categories
no_rater_groups <- 1 # no groups of rater
no_covariates <- 5   # no additional covariates
no_time_points <- 2

set.seed(1234)

dat_orig <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points, use_skew_thresholds = FALSE, use_free_thresholds = TRUE)
dat_skew <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points, use_skew_thresholds = TRUE,  use_free_thresholds = FALSE)
dat_free <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points, use_skew_thresholds = FALSE, use_free_thresholds = TRUE)

dat_log_reg_orig <- simulate_logistic_regression(dat_orig, no_covariates = no_covariates, slopes = seq_len(no_covariates))
dat_log_reg_skew <- simulate_logistic_regression(dat_skew, no_covariates = no_covariates, slopes = seq_len(no_covariates))
dat_log_reg_free <- simulate_logistic_regression(dat_free, no_covariates = no_covariates, slopes = seq_len(no_covariates))

dataset_list <- list(
  orig = dat_log_reg_orig,
  skew = dat_log_reg_skew,
  free = dat_log_reg_free
)

# compile the stan model
mod_log_reg_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")

# fit all 9 models
results <- pmap(list(dataset = dataset_list, data_name = names(dataset_list)), fit_single_dataset, model = mod_log_reg_ltm)

system("beep_finished.sh 0")

grid.arrange(results$orig$graphs_retrieval_joined)
grid.arrange(results$orig$graphs_retrieval_joined)
grid.arrange(results$skew$graphs_retrieval_joined)
grid.arrange(results$skew$graphs_retrieval_joined)
grid.arrange(results$free$graphs_retrieval_joined)
grid.arrange(results$free$graphs_retrieval_joined)


# fit original model ----

fits_orig <- fit_all_three_models(dat_log_reg_orig, mod_log_reg_ltm, path_prefix = "data_orig")

perf_dat_orig_tib <- purrr::map2(fits_orig, c("2-parameter_thresholds", "skew_thresholds", "free_thresholds"), \(x, name) get_performance_tib(x, dat_log_reg_orig, name, nc, nr, no_covariates))
# individual graphs
graphs_dat_orig <- map2(fits_orig, perf_dat_orig_tib, get_scatterplot_retrieval)
# joined graph of parameter retrieval
g_joined <- grid.arrange(grobs = graphs_dat_orig, ncol = 3)

# predictions
g_predictions <- purrr::map2(fits_orig, names(fits_orig), get_predictions_plot, data = dat_log_reg_orig)

g_predictions_joined <- grid.arrange(grobs = matrix(unlist(g_predictions, recursive = FALSE), 3, 2, TRUE), nrow = 2, ncol = 3)
system("beep_finished.sh 0")

# fit skew model ----
fits_skew <- fit_all_three_models(dat_log_reg_skew, mod_log_reg_ltm, path_prefix = "data_skew")

perf_dat_orig_tib <- purrr::map2(fits_skew, c("2-parameter_thresholds", "skew_thresholds", "free_thresholds"), \(x, name) get_performance_tib(x, dat_log_reg_skew, name, nc, nr, no_covariates))
# individual graphs
graphs_dat_orig <- map2(fits_orig, perf_dat_orig_tib, get_scatterplot_retrieval)
# joined graph of parameter retrieval
g_joined <- grid.arrange(grobs = graphs_dat_orig, ncol = 3)

# predictions
g_predictions <- purrr::map2(fits_orig, names(fits_orig), get_predictions_plot, data = dat_log_reg_orig)

g_predictions_joined <- grid.arrange(grobs = matrix(unlist(g_predictions, recursive = FALSE), 3, 2, TRUE), nrow = 2, ncol = 3)
system("beep_finished.sh 0")


perf_orig_dat_orig_tib <- get_performance_tib(fits_orig$fit_orig, dat_log_reg_orig, "2-parameter_thresholds", nc, nr, no_covariates)
perf_skew_dat_orig_tib <- get_performance_tib(fits_orig$fit_skew, dat_log_reg_orig, "skew_thresholds",        nc, nr, no_covariates)
perf_free_dat_orig_tib <- get_performance_tib(fits_orig$fit_free, dat_log_reg_orig, "free_thresholds",        nc, nr, no_covariates)

g1 <- get_scatterplot_retrieval(fits_orig$fit_orig, perf_orig_dat_orig_tib)
g2 <- get_scatterplot_retrieval(fits_orig$fit_skew, perf_skew_dat_orig_tib)
g3 <- get_scatterplot_retrieval(fits_orig$fit_free, perf_free_dat_orig_tib)

g_joined <- gridExtra::grid.arrange(g1, g2, g3, ncol = 3)


iter           <- 3e4
adapt_iter     <- 500
output_samples <- 2e3
grad_samples   <- 5
elbo_samples   <- 5
debug          <- TRUE
force          <- TRUE
path_prefix    <- ""
path_orig <- file.path("fitted_objects", sprintf("%s3_log_reg_models_orig_thresholds.rds", path_prefix))

purrr::map2


fit_orig <- save_or_run_model(
  mod_log_reg_ltm$variational(data = stan_data, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples),
  path_orig, force
)
system("beep_finished.sh 0")
# fit_orig <- save_or_run_model(
#   mod_log_reg_ltm$sample(data = stan_data, iter_sampling = output_samples, chains = 3),
#   path_orig, force
# )

# TODO: split log_reg_slopes based on lt vs covariate
get_means_tib(fit_orig) |>
  add_true_values_to_tib(dat_log_reg) |>
  split_slopes_on_covariate_vs_lt(no_covariates) |>
  scatterplot_retrieval() +
  ggtitle(sprintf("Variational inference for %s took %.3f seconds", "free thresholds", fit_orig$time())) +
  labs(x = "True value", y = "Variational Estimate")

probs <- fit_orig$draws("log_reg_predictions", format = "draws_matrix")

qplot(
  factor(dat_log_reg$log_reg$y),
  colMeans(unname(probs)),
  xlab = "Observed outcome",
  ylab = "Posterior mean",
)
qplot(
  c(dat_log_reg$log_reg$true_probs),
  colMeans(unname(probs)),
  xlab = "True probability",
  ylab = "Posterior mean",
) + geom_abline(slope = 1, intercept = 0)



# assuming we "know" the true latent truths
np <- 90             # no patients
ni <- 20             # no items
nr <- 10             # no raters
nc <- 5              # no response categories
no_rater_groups <- 1 # no groups of rater
no_covariates <- 3

set.seed(1234)

dat_orig <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_skew_thresholds = FALSE, use_free_thresholds = TRUE)
# debugonce(simulate_logistic_regression)
dat_log_reg <- simulate_logistic_regression(dat_orig, no_covariates = no_covariates, slopes = seq_len(no_covariates))

mod_simple_log_reg <- compile_stan_model("stanmodels/simple_logistic_regression.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")
design_mat <- dat_log_reg$log_reg$design_matrix[, -1L]
stan_data_log_reg <- list(
  np               = nrow(design_mat),
  no_covariates    = ncol(design_mat),
  design_matrix    = design_mat,
  log_reg_outcomes = dat_log_reg$log_reg$y
)

log_reg_fit <- mod_simple_log_reg$variational(stan_data_log_reg, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples)
log_reg_samps <- log_reg_fit$draws(format = "draws_matrix")[, -(1:2)]

log_reg_fit <- mod_simple_log_reg$sample(stan_data_log_reg, chains = 4)
log_reg_samps <- log_reg_fit$draws(format = "draws_matrix")[, -1]

qplot(
  c(dat_log_reg$log_reg$log_reg_intercept, dat_log_reg$log_reg$log_reg_slopes),
  colMeans(log_reg_samps),
  xlab = "Observed outcome",
  ylab = "Posterior mean",
) + geom_abline(intercept = 0, slope = 1)


x <- list(1, 1, 1)
y <- list(10, 20, 30)
z <- list(100, 200, 300)

