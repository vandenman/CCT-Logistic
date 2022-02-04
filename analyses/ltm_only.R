rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(purrr)

# cpp_options <- list(
#   "CXXFLAGS+= -oFast -march=native -mtune=native"
# )
# cmdstan_make_local(cpp_options = cpp_options)
# rebuild_cmdstan(cores = 8)

fit_all_three_models <- function(data, model, logistic_dat = NULL, logistic_target = NULL, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, threads = 8, debug = FALSE, force = FALSE,
                                 path_prefix = "", store_predictions = TRUE) {

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  stan_data_orig <- data_2_stan(data, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  stan_data_skew <- data_2_stan(data, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

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
    fit = list(
      orig = fit_orig,
      skew = fit_skew,
      free = fit_free
    ),
    stan_data = list(
      orig = stan_data_orig,
      skew = stan_data_skew,
      free = stan_data_free
    )
  ))
}

normalize_factor <- function(x) {
  if (is.factor(x))
    as.integer(droplevels(x))
  else if (is.integer(x))
    match(x, unique(x))
}

compute_mean_probs <- function(fits, cache_filename, force = FALSE) {
  nms <- names(fits$fit)
  save_or_run_model(
    setNames(pmap(list(
      nms, fits$fit, fits$stan_data
    ), \(name, fit, stan_data) {
      cat(sprintf("computing probs for %s\n", name))
      variable_sizes <- fit$metadata()$stan_variable_sizes
      draws          <- fit$draws(format = "draws_matrix")
      compute_probs(draws, stan_data, variable_sizes, return_mean = TRUE)
    }), nms),
    path = cache_filename,
    force = force
  )
}

# TODO: this data-cleaning step should happen elsewhere, so that the data can be easily switched out for something else
all_data <- read_long_data()
data_2_analyze <- all_data |>
  as_tibble() |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
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

data_violence <- all_data |>
  as_tibble() |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
  select(c(patient, patient_age_group, violent_before, violent_between, violent_after, treatement_duration_group, diagnosis_group, crime_group)) |>
  mutate(
    patient = normalize_factor(patient),
    violent_after = as.integer(violent_after) - 1L
  ) |>
  filter(!duplicated(patient))

# # omit 2 patients with missing values
# bad_patients <- unique(which(is.na(data_violence),arr.ind = TRUE)[, 1L])
#
# data_2_analyze <- data_2_analyze |>
#   filter(!(patient %in% bad_patients)) |>
#   mutate(patient = normalize_factor(patient))
#
# data_violence  <- data_violence |>
#   filter(!(patient %in% bad_patients)) |>
#   mutate(patient = normalize_factor(patient))

# ensure that both datasets contain the same patients
all(unique(data_2_analyze$patient) == unique(data_violence$patient))

colnames(data_2_analyze)
# [1] "time"        "patient"     "rater"       "item"        "score"       "rater_group"
colnames(data_violence)

np <- max(data_2_analyze$patient)
ni <- max(data_2_analyze$item)
nr <- max(data_2_analyze$rater)
nt <- max(data_2_analyze$time)
nc <- 18L

mod_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels",
                              cpp_options = list(stan_threads=TRUE))

fits_without_lr <- fit_all_three_models(data_2_analyze, mod_ltm, NULL,          NULL,            path_prefix = "mesdag_ltm_without_logistic_regression")
fits_with_lr    <- fit_all_three_models(data_2_analyze, mod_ltm, data_violence, "violent_after", path_prefix = "mesdag_ltm_with_logistic_regression")
# system("beep_finished.sh 0")

mean_probs_without_lr <- compute_mean_probs(fits_without_lr, "fitted_objects/mesdag_ltm_probs_without_logistic_regression.rds")
mean_probs_with_lr    <- compute_mean_probs(fits_with_lr,    "fitted_objects/mesdag_ltm_probs_with_logistic_regression.rds")

# TODO: the code below should be done for both sets of fits!!
observed_proportions <- get_obs_proportions(data_2_analyze$score)
rmse <- function(x, y) sqrt(mean((x - y)^2))

raw_probabilities_without_lr <- unlist(map(mean_probs_without_lr, rowMeans), use.names = FALSE)
raw_probabilities_with_lr    <- unlist(map(mean_probs_with_lr,    rowMeans), use.names = FALSE)

plot_ltm_fit <- function(raw_probabilities, observed_proportions, ylim = NULL) {
  nc <- length(observed_proportions)
  betterNames <- c("Original", "Skew", "Free")
  tib <- tibble(
    panel = rep(c("Raw", "Diff"), c(4L * nc, 3L * nc)),
    fill  = c(rep(c(betterNames, "Observed"), each = nc), rep(betterNames, each = nc)),
    y     = c(raw_probabilities, observed_proportions, raw_probabilities - rep(observed_proportions, 3L)),
    x     = rep(1:nc, 7L)
  )
  cols <- hcl.colors(3L, "Set 2")
  colors <- c("Observed" = "black", setNames(cols, betterNames))

  if (is.null(ylim)) {
    breaks_fun <- jaspGraphs::getPrettyAxisBreaks
    limits_fun <- \(x) range(jaspGraphs::getPrettyAxisBreaks(x))
  } else {
    stopifnot(is.list(ylim) && length(ylim) == 2L)
    counter_breaks <<- 0L
    counter_limits <<- 0L
    limits_fun <- \(x) {res <- ylim[[counter_limits+1]]; counter_limits <<- (counter_limits+1L) %% 2; res}
    breaks_fun <- \(x) {res <- jaspGraphs::getPrettyAxisBreaks(ylim[[counter_breaks+1L]]); counter_breaks <<- (counter_breaks + 1L) %% 2; res}
  }

  fit_plot <- ggplot(data = tib, mapping = aes(x = x, y = y, group = fill, fill = fill)) +
    geom_bar(position="dodge", stat="identity", width = .5) +
    facet_wrap(~panel, scales = "free_y") +
    scale_x_continuous(name = "Score", breaks = 1:18, limits = c(0, 19)) +
    scale_y_continuous(name = "Diff (1) / Probability (2)", breaks = breaks_fun, limits = limits_fun) +
    scale_fill_manual(name = NULL, values = colors) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right") +
    theme(strip.text = element_text(size = 28))
  return(fit_plot)
}

pdf("simulation_figures/ltm_only_joined_ggplot22.pdf", width = 16, height = 8)
ylim <- list(c(-.2, .4), c(0, .4))
plot_ltm_fit(raw_probabilities_without_lr, observed_proportions, ylim = ylim) + labs(title = "Without logistic regression")
plot_ltm_fit(raw_probabilities_with_lr, observed_proportions,    ylim = ylim) + labs(title = "With logistic regression")
dev.off()

# one 4 panel plot
betterNames <- c("Original", "Skew", "Free")
tib <- tibble(
  cols  = rep(rep(c("Raw", "Diff"),     c(4L * nc, 3L * nc)), 2),
  rows  = rep(c("without", "with"), each = 4L * nc + 3L * nc),
  fill  = rep(c(rep(c(betterNames, "Observed"), each = nc), rep(betterNames, each = nc)), 2),
  y     = c(c(raw_probabilities_without_lr, observed_proportions, raw_probabilities_without_lr - rep(observed_proportions, 3L)),
            c(raw_probabilities_with_lr, observed_proportions, raw_probabilities_with_lr - rep(observed_proportions, 3L))),
  x     = rep(1:nc, 14L)
)
cols <- hcl.colors(3L, "Set 2")
colors <- c("Observed" = "black", setNames(cols, betterNames))

# ylim <- list(c(-.2, .4), c(0, .4))
# counter_breaks <<- 0L
# counter_limits <<- 0L
# limits_fun <- \(x) {res <- ylim[[counter_limits+1]]; counter_limits <<- (counter_limits+1L) %% 2; res}
# breaks_fun <- \(x) {res <- jaspGraphs::getPrettyAxisBreaks(ylim[[counter_breaks+1L]]); counter_breaks <<- (counter_breaks + 1L) %% 2; res}

fit_plot <- ggplot(data = tib, mapping = aes(x = x, y = y, group = fill, fill = fill)) +
    geom_bar(position="dodge", stat="identity", width = .5) +
    facet_wrap(rows~cols, scales = "free_y") +
    # facet_grid(rows = vars(rows), cols = vars(cols), scales = "free") +
    scale_x_continuous(name = "Score", breaks = 1:18, limits = c(0, 19)) +
    scale_y_continuous(name = "Diff (1) / Probability (2)", breaks = breaks_fun, limits = limits_fun) +
    scale_fill_manual(name = NULL, values = colors) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right") +
    theme(strip.text = element_text(size = 28))

pdf("simulation_figures/ltm_only_joined_ggplot2_4plts.pdf", width = 16, height = 8)
fit_plot
dev.off()

plotly::ggplotly(fit_plot)

pp <- plotly::ggplotly(plot_ltm_fit(raw_probabilities_without_lr, observed_proportions, ylim = ylim) + labs(title = "Without logistic regression"))
pp


observed_violence <- fits_with_lr$stan_data$orig$log_reg_outcomes
prediction_results <- map(fits_with_lr$fit, \(fit) {
  log_reg_predictions <- unname(colMeans(fit$draws("log_reg_predictions", format = "draws_matrix")))
  discrete_prediction <- 1 * (log_reg_predictions > .5)
  confusion_matrix <- table(observed_violence, discrete_prediction)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  return(list(
    log_reg_predictions = log_reg_predictions,
    discrete_prediction = discrete_prediction,
    confusion_matrix    = confusion_matrix,
    accuracy            = accuracy
  ))
})


prediction_results$orig
prediction_results$skew

idx_incorrect <- which(prediction_results$free$discrete_prediction != observed_violence)
prediction_results$free$log_reg_predictions[idx_incorrect] # reasonably close to 0.5

# mean_probs <- map(fits, function(fit) {
#   log_probs <- fit$draws("log_probs", format = "draws_matrix")
#   mean_probs <- matrix(apply(log_probs, 2L, function(x) mean(exp(x))), nc, nrow(data_2_analyze))
#   mean_probs
# })

plts00 <- map(mean_probs, function(m_prob) {
  model_probabilities  <- rowMeans(m_prob)
  inspect_probability_plot(c(observed_proportions), c(model_probabilities))
})

rawDiffs <- do.call(rbind, map(mean_probs, function(m_prob) {
  model_probabilities  <- rowMeans(m_prob)
  observed_proportions - model_probabilities
}))


rmses <- map(mean_probs, function(m_prob) {
  model_probabilities  <- rowMeans(m_prob)
  c(
    "rmsea" = rmse(observed_proportions, model_probabilities)#,
    # "wasserstein" = transport::wasserstein1d(a = 1:nc, b = 1:nc, wa = observed_proportions, wb = model_probabilities)
  )
})

yLimLeft <- range(rawDiffs)
plts <- map(setNames(names(mean_probs), names(mean_probs)), function(nm) {
  model_probabilities  <- rowMeans(mean_probs[[nm]])
  inspect_probability_plot2(c(observed_proportions), c(model_probabilities), yLimLeft = c(-0.06, 0.08), title = nm)
})
# system("beep_finished.sh 0")

pdf("simulation_figures/ltm_only_joined.pdf", width = 20, height = 10)
gridExtra::grid.arrange(plts$orig)
gridExtra::grid.arrange(plts$skew)
gridExtra::grid.arrange(plts$free)
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

