rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)

np <- 40             # no patients
ni <- 20             # no items
nr <- 30             # no raters
nc <- 5              # no response categories
no_rater_groups <- 3 # no groups of rater

set.seed(1234)
dat <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_free_thresholds = FALSE)

mod_sparse <- compile_stan_model("stanmodels/LTM.stan",                 pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")
mod_free   <- compile_stan_model("stanmodels/LTM_free_thresholds.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")

stan_data <- ltm_data_2_stan(dat, debug = TRUE)

path_fit_vb1 <- file.path("fitted_objects", "simulation_ltm_sparse_fixed_thresholds.rds")
path_fit_vb2 <- file.path("fitted_objects", "simulation_ltm_sparse_free_thresholds.rds")

if (!file.exists(path_fit_vb1)) {
  fit_sparse <- mod_sparse$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5)
  fit_sparse$save_object(path_fit_vb1)
} else {
  fit_sparse <- readRDS(path_fit_vb1)
}

if (!file.exists(path_fit_vb2)) {
  fit_free <- mod_free$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3)
  fit_free$save_object(path_fit_vb2)
} else {
  fit_free <- readRDS(path_fit_vb2)
}
system("beep_finished.sh 0")

get_means_tib(fit_sparse) |>
  add_implied_thresholds(nc, nr) |>
  add_true_values_to_tib(dat) |>
  mutate(parameter = if_else(parameter == "free_thresholds", "sparse_thresholds", parameter)) |>
  scatterplot_retrieval() +
  ggtitle(sprintf("variational inference for 2-parameter thresholds took %.3f seconds", fit_sparse$time())) +
  labs(x = "Variational Estimate", y = "True value")

get_means_tib(fit_free) |>
  add_true_values_to_tib(dat) |>
  scatterplot_retrieval() +
  ggtitle(sprintf("variational inference for free thresholds took %.3f seconds", fit_free$time())) +
  labs(x = "Variational Estimate", y = "True value")


log_probs_sparse           <- fit_sparse$draws("log_probs", format = "draws_matrix")
mean_probs_sparse          <- matrix(apply(log_probs_sparse, 2L, function(x) mean(exp(x))), nc, nrow(dat$df))
model_probabilities_sparse <- rowMeans(mean_probs_sparse)
rm(log_probs_sparse); gc()

log_probs_free           <- fit_free$draws("log_probs", format = "draws_matrix")
mean_probs_free          <- matrix(apply(log_probs_free, 2L, function(x) mean(exp(x))), nc, nrow(dat$df))
model_probabilities_free <- rowMeans(mean_probs_free)
rm(log_probs_free); gc()

observed_proportions <- obs_proportions(dat$df$score, nc = nc)

inspect_probability_plot(c(observed_proportions), c(model_probabilities_sparse)) + ggtitle("sparse thresholds")
inspect_probability_plot(c(observed_proportions), c(model_probabilities_free))   + ggtitle("free thresholds")


# now with free thresholds as true model
set.seed(1234)
dat_free <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, use_free_thresholds = TRUE)

stan_data <- ltm_data_2_stan(dat_free, debug = TRUE)

path_fit_vb1 <- file.path("fitted_objects", "simulation_ltm_free_fixed_thresholds.rds")
path_fit_vb2 <- file.path("fitted_objects", "simulation_ltm_free_free_thresholds.rds")

if (!file.exists(path_fit_vb1)) {
  fit_sparse <- mod_sparse$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5)
  fit_sparse$save_object(path_fit_vb1)
} else {
  fit_sparse <- readRDS(path_fit_vb1)
}

if (!file.exists(path_fit_vb2)) {
  fit_free <- mod_free$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3)
  fit_free$save_object(path_fit_vb2)
} else {
  fit_free <- readRDS(path_fit_vb2)
}

get_means_tib(fit_sparse) |>
  add_implied_thresholds(nc, nr) |>
  add_true_values_to_tib(dat) |>
  mutate(parameter = if_else(parameter == "free_thresholds", "sparse_thresholds", parameter)) |>
  scatterplot_retrieval() +
  ggtitle(sprintf("variational inference for 2-parameter thresholds took %.3f seconds", fit_sparse$time())) +
  labs(x = "Variational Estimate", y = "True value")

get_means_tib(fit_free) |>
  add_true_values_to_tib(dat) |>
  scatterplot_retrieval() +
  ggtitle(sprintf("variational inference for free thresholds took %.3f seconds", fit_free$time())) +
  labs(x = "Variational Estimate", y = "True value")


log_probs_sparse           <- fit_sparse$draws("log_probs", format = "draws_matrix")
mean_probs_sparse          <- matrix(apply(log_probs_sparse, 2L, function(x) mean(exp(x))), nc, nrow(dat$df))
model_probabilities_sparse <- rowMeans(mean_probs_sparse)
rm(log_probs_sparse); gc()

log_probs_free           <- fit_free$draws("log_probs", format = "draws_matrix")
mean_probs_free          <- matrix(apply(log_probs_free, 2L, function(x) mean(exp(x))), nc, nrow(dat$df))
model_probabilities_free <- rowMeans(mean_probs_free)
rm(log_probs_free); gc()

observed_proportions <- obs_proportions(dat$df$score, nc = nc)

inspect_probability_plot(c(observed_proportions), c(model_probabilities_sparse)) + ggtitle("sparse thresholds")
inspect_probability_plot(c(observed_proportions), c(model_probabilities_free))   + ggtitle("free thresholds")
