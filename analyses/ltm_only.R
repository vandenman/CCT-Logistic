rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)

normalize_factor <- function(x) {
  as.integer(droplevels(x))
}

all_data <- read_long_data()
data_2_analyze <- all_data |>
  as_tibble() |>
  filter(as.integer(time) == 1 & !is.na(score)) |>
  select(-c(time, patient_age_group, violent_before, violent_between, violent_after,
                 treatement_duration_group, diagnosis_group, crime_group)) |>
  mutate(
    patient     = normalize_factor(patient),
    item        = normalize_factor(item),
    rater       = normalize_factor(rater),
    rater_group = normalize_factor(rater_group),
    score       = score + 1 # transform score from 0 -17 to 1 - 18
  )

np <- max(data_2_analyze$patient)
ni <- max(data_2_analyze$item)
nr <- max(data_2_analyze$rater)
nc <- 18L

# TODO: rewrite this with fit_all_three_models

stan_data <- ltm_data_2_stan(data_2_analyze, nc = 18, debug = TRUE, vary_lambda_across_patients = FALSE)

# mod <- compile_stan_model("stanmodels/LTM.stan", pedantic = TRUE, quiet = FALSE)
mod <- compile_stan_model("stanmodels/LTM_free_thresholds.stan", pedantic = TRUE, quiet = FALSE)

path_for_fitted_results <- file.path("fitted_objects", "fitted_ltm_only_fit_lambda_fixed_across_patients_vb.rds")
# path_for_fitted_results <- file.path("fitted_objects", "fitted_ltm_only_fit_vb.rds")
if (!file.exists(path_for_fitted_results)) {

  fit_vb <- mod$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 2, elbo_samples = 2)
  failure <- any(fit_vb$return_codes() != 0)

  if (!failure)
    fit_vb$save_object(path_for_fitted_results)

  system(paste("beep_finished.sh", as.integer(failure)), wait = FALSE)

} else {
  fit_vb <- readRDS(path_for_fitted_results)
  # samples <- readRDS(file.path("fitted_objects", "fitted_ltm_only_samples.rds"))
}

# tib_vb <- fit_vb |> get_means_tib()

log_probs <- fit_vb$draws("log_probs", format = "draws_matrix")
# sum(exp(log_probs[5, 1:18 + 36]))
dim(log_probs)
head(log_probs[, 1:20])

mean_probs <- matrix(apply(log_probs, 2L, function(x) mean(exp(x))), nc, nrow(data_2_analyze))
max(abs(colSums(mean_probs) - 1))

observed_proportions <- obs_proportions(data_2_analyze$score)
model_probabilities  <- rowMeans(mean_probs)
pplot(c(observed_proportions), c(model_probabilities))
inspect_probability_plot(c(observed_proportions), c(model_probabilities))

# aggregated over raters and items
observed_proportions <- obs_proportions(data_2_analyze$score, data_2_analyze$patient)
model_probabilities  <- sapply(seq_len(np), \(p) rowMeans(mean_probs[, data_2_analyze$patient == p, drop = FALSE]) / np)
# inspect_probability_plot(c(observed_proportions), c(model_probabilities))
pplot(c(observed_proportions), c(model_probabilities))

observed_proportions <- obs_proportions(data_2_analyze$score, data_2_analyze$rater)
model_probabilities  <- sapply(seq_len(nr), \(r) rowMeans(mean_probs[, data_2_analyze$rater == r, drop = FALSE]) / nr)
pplot(c(observed_proportions), c(model_probabilities))

# aggregated over patients and raters
observed_proportions <- obs_proportions(data_2_analyze$score, data_2_analyze$item)
model_probabilities  <- sapply(seq_len(ni), \(i) rowMeans(mean_probs[, data_2_analyze$item == i, drop = FALSE]) / ni)
pplot(c(observed_proportions), c(model_probabilities))

# aggregated over patients and raters
observed_proportions <- obs_proportions(data_2_analyze$score, data_2_analyze$item)
model_probabilities  <- sapply(seq_len(ni), \(i) rowMeans(mean_probs[, data_2_analyze$item == i, drop = FALSE]) / ni)
pplot(c(observed_proportions), c(model_probabilities))

data_one_item <- data_2_analyze |> filter(item == 1)
which.max(tapply(data_one_item$score, data_one_item$rater, length))

data_rater_54 <- data_2_analyze |> filter(item == 1, rater == 54)
observed_proportions <- obs_proportions(data_rater_54$score)
model_probabilities  <- rowMeans(sapply(seq_len(ni), \(i) rowMeans(mean_probs[, data_2_analyze$item == 1 & data_2_analyze$rater == 54, drop = FALSE])))
pplot(c(observed_proportions), c(model_probabilities))


eee <- runif(100)
mean(eee)
plogis(mean(qlogis(eee)))

eee <- colSums(mean_log_probs)
max(abs(eee - 1))

prob_incorrect       <- numeric(ncol(mean_log_probs))
mode_obs             <- integer(ncol(mean_log_probs))
prob_close_incorrect <- numeric(ncol(mean_log_probs))
dist <- 4

score_to_interval <- function(score, distance, lower = 1, upper = 18) {
  interval <- seq(score - distance, score + distance)
  if (interval[1L] < lower)
    interval <- interval + (lower - interval[1L])

  n <- length(interval)
  if (interval[n] > upper)
    interval <- interval + (upper - interval[n])

  interval
}


for (i in seq_along(prob_incorrect)) {
  prob_incorrect[i] <- 1 - mean_log_probs[data_2_analyze$score[i], i]
  interval <- score_to_interval(data_2_analyze$score[i], dist)

  prob_close_incorrect[i] <- 1 - sum(mean_log_probs[interval, i])
  mode_obs[i] <- which.max(mean_log_probs[, i]) == data_2_analyze$score[i]
}

mean(mode_obs)
mean(prob_incorrect)
hist(prob_incorrect)

mean(prob_close_incorrect)
hist(prob_close_incorrect)

