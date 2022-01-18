rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(purrr)

cpp_options <- list(
  "CXXFLAGS+= -oFast -march=native -mtune=native"
)
cmdstan_make_local(cpp_options = cpp_options)
rebuild_cmdstan(cores = 8)

fit_all_three_models <- function(data, model, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, threads = 8, debug = FALSE, force = FALSE,
                                 path_prefix = "", store_predictions = TRUE) {

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  stan_data_orig <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  stan_data_skew <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

  path_orig <- file.path("fitted_objects", sprintf("%stime_3_models_orig_thresholds.rds", path_prefix))
  path_skew <- file.path("fitted_objects", sprintf("%stime_3_models_skew_thresholds.rds", path_prefix))
  path_free <- file.path("fitted_objects", sprintf("%stime_3_models_free_thresholds.rds", path_prefix))

  fit_model <- function(data) {
    model$variational(data = data, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples,
                      threads = threads)
  }

  cat("Fitting original threshold model\n")
  fit_orig <- save_or_run_model(fit_model(stan_data_orig), path_orig, force)

  cat("Fitting skew threshold model\n")
  fit_skew <- save_or_run_model(fit_model(stan_data_skew), path_skew, force)

  cat("Fitting free threshold model\n")
  fit_free <- save_or_run_model(fit_model(stan_data_free), path_free, force)

  return(list(
    fit_orig = fit_orig,
    fit_skew = fit_skew,
    fit_free = fit_free
  ))
}

normalize_factor <- function(x) {
  as.integer(droplevels(x))
}

# TODO: this data-cleaning step should happen elsewhere, so that the data can be easily switched out for something else
all_data <- read_long_data()
data_2_analyze <- all_data |>
  as_tibble() |>
  filter(!is.na(score)) |>
  select(-c(patient_age_group, violent_before, violent_between, violent_after,
                 treatement_duration_group, diagnosis_group, crime_group)) |>
  mutate(
    patient     = normalize_factor(patient),
    item        = normalize_factor(item),
    rater       = normalize_factor(rater),
    time        = normalize_factor(time),
    rater_group = normalize_factor(rater_group),
    score       = score + 1 # transform score from 0 -17 to 1 - 18
  ) |>
  arrange(rater_group, patient, item, rater, time)

np <- max(data_2_analyze$patient)
ni <- max(data_2_analyze$item)
nr <- max(data_2_analyze$rater)
nt <- max(data_2_analyze$time)
nc <- 18L

mod_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels",
                              cpp_options = list(stan_threads=TRUE))

# debugonce(CCTLogistic:::save_or_run_model)
fits <- fit_all_three_models(data_2_analyze, mod_ltm, path_prefix = "mesdag_ltm")
# system("beep_finished.sh 0")

observed_proportions <- obs_proportions(data_2_analyze$score)
rmse <- function(x, y) sqrt(mean((x - y)^2))

mean_probs <- map(fits, function(fit) {
  log_probs <- fit$draws("log_probs", format = "draws_matrix")
  mean_probs <- matrix(apply(log_probs, 2L, function(x) mean(exp(x))), nc, nrow(data_2_analyze))
  mean_probs
})

plts <- map(mean_probs, function(m_prob) {
  model_probabilities  <- rowMeans(m_prob)
  inspect_probability_plot(c(observed_proportions), c(model_probabilities))
})

rmses <- map(mean_probs, function(m_prob) {
  model_probabilities  <- rowMeans(m_prob)
  c(
    "rmsea" = rmse(observed_proportions, model_probabilities),
    "wasserstein" = transport::wasserstein1d(a = 1:nc, b = 1:nc, wa = observed_proportions, wb = model_probabilities)
  )
})


plts <- map(setNames(names(mean_probs), names(mean_probs)), function(nm) {
  model_probabilities  <- rowMeans(mean_probs[[nm]])
  inspect_probability_plot2(c(observed_proportions), c(model_probabilities), yLimLeft = c(-0.06, 0.08), title = nm)
})
# system("beep_finished.sh 0")

pdf("simulation_figures/ltm_only_joined.pdf", width = 20, height = 10)
gridExtra::grid.arrange(plts$fit_orig)
gridExtra::grid.arrange(plts$fit_skew)
gridExtra::grid.arrange(plts$fit_free)
dev.off()

map(names(plts), \(nm) {
  filename <- file.path("simulation_figures", sprintf("ltm_only_%s.pdf", nm))
  pdf(filename, width = 20, height = 10)
  on.exit(dev.off())
  gridExtra::grid.arrange(plts[[nm]])
  # print(plts[[nm]] + ggtitle(nm))
  invisible(NULL)
})


transport::wasserstein1d(1:nc, 1:nc, wa = observed_proportions, wb = observed_proportions)
transport::wasserstein1d(sample(nc, 1e6, TRUE, prob = observed_proportions), sample(nc, 1e6, TRUE, prob = observed_proportions))


inspect_probability_plot2 <- function(observed_proportions, model_probabilities, yLimLeft = NULL, yLimRight = NULL, title = NULL) {

  nc <- length(observed_proportions)
  diff <- observed_proportions - model_probabilities
  tib <- tibble(
    category    = factor(rep(1:nc, 3)),
    probability = c(observed_proportions, model_probabilities, diff),
    group       = factor(rep(c("observed", "model", "diff"), each = nc)),
    panel       = rep(c("raw", "difference"), c(2*nc, nc))
  )

  tibAbline <- tibble(
    intercept = 0,
    slope     = 0,
    panel     = "difference"
  )

  yBreaksLeft <- jaspGraphs::getPrettyAxisBreaks(c(yLimLeft, 0, tib |> filter(panel == "difference") |> select(probability) |> unlist()))
  pLeft <- ggplot(data = tib |> filter(panel == "difference"), aes(x = category, y = probability, group = group, fill = group)) +
    geom_bar(position="dodge", stat="identity", width = .4) +
    scale_y_continuous(breaks = yBreaksLeft, limits = range(yBreaksLeft)) +
    labs(fill = NULL, x = "Category", y = "Probability") +
    facet_wrap(~panel, scales = "free_y") +
    ggtitle(title) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right")

  yBreaksRight <- jaspGraphs::getPrettyAxisBreaks(c(yLimRight, 0, tib |> filter(panel == "raw") |> select(probability) |> unlist()))
  pRight <- ggplot(data = tib |> filter(panel == "raw"), aes(x = category, y = probability, group = group, fill = group)) +
    geom_bar(position="dodge", stat="identity", width = .4) +
    scale_y_continuous(breaks = yBreaksRight, limits = range(yBreaksRight)) +
    labs(fill = NULL, x = "Category", y = "Probability") +
    facet_wrap(~panel, scales = "free_y") +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right")

  return(gridExtra::arrangeGrob(grobs = list(pLeft, pRight), nrow = 1L))
}







log_probs <- fit_vb$draws("log_probs", format = "draws_matrix")
# sum(exp(log_probs[5, 1:18 + 36]))
dim(log_probs)
head(log_probs[, 1:20])

mean_probs <- matrix(apply(log_probs, 2L, function(x) mean(exp(x))), nc, nrow(data_2_analyze))
max(abs(colSums(mean_probs) - 1))



stan_data <- data_2_stan(data_2_analyze, nc = 18, debug = TRUE, vary_lambda_across_patients = FALSE)

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

