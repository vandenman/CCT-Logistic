rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(purrr)

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

extract_numeric <- function(character_vector) {
  res <- regmatches(character_vector, gregexpr("[[:digit:]]+", character_vector))
  if (all(lengths(res) == 1L))
    return(unlist(res))
  return(res)
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
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(-c(age, violent_before, violent_between, violent_after,
                 treatment_duration, diagnosis, crime)) |>
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
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(c(patient, age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
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

mean_probs_without_lr <- compute_mean_probs(fits_without_lr, "fitted_objects/mesdag_ltm_probs_without_logistic_regression.rds")
mean_probs_with_lr    <- compute_mean_probs(fits_with_lr,    "fitted_objects/mesdag_ltm_probs_with_logistic_regression.rds")
# system("beep_finished.sh 0")

# TODO: the code below should be done for both sets of fits!!
# inspect LTM fit ----
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

ylim <- list(c(-.2, .4), c(0, .4))
counter_breaks <<- 0L
counter_limits <<- 0L
limits_fun <- \(x) {res <- ylim[[counter_limits+1]]; counter_limits <<- (counter_limits+1L) %% 2; res}
breaks_fun <- \(x) {res <- jaspGraphs::getPrettyAxisBreaks(ylim[[counter_breaks+1L]]); counter_breaks <<- (counter_breaks + 1L) %% 2; res}

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

plot_for_presentation <- tib |>
  filter(cols == "Raw" & rows == "without") |>
  ggplot(mapping = aes(x = x, y = y, group = fill, fill = fill)) +
    geom_bar(position="dodge", stat="identity", width = .5) +
    scale_x_continuous(name = "Score", breaks = 1:18, limits = c(0, 19)) +
    scale_y_continuous(name = "Probability", breaks = seq(0, .4, .1), limits = c(0, .4)) +
    scale_fill_manual(name = NULL, values = colors) +
    # jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right") +
    theme(strip.text = element_text(size = 28)) +
    geom_segment(x = 1, xend = 18, y = -Inf, yend = -Inf) +
    geom_segment(x = -Inf, xend = -Inf, y = 0.0, yend = 0.4)
saveRDS(plot_for_presentation, file = file.path("figure_r_objs", "ltm_fit.rds"))
pdf("simulation_figures/ltm_only_joined_ggplot2_4plts.pdf", width = 8, height = 8)
fit_plot
dev.off()

plotly::ggplotly(plot_for_presentation)

# inspect fit logisitc regression ----
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


map_dbl(prediction_results, `[[`, "accuracy")

idx_incorrect <- which(prediction_results$free$discrete_prediction != observed_violence)
prediction_results$free$log_reg_predictions[idx_incorrect] # reasonably close to 0.5

modelNames <- names(fits_with_lr$fit)
modelNames <-  setNames(as.list(modelNames), modelNames)

log_reg_results <- map(modelNames, \(nm) {
  log_reg_slopes_samps <- fits_with_lr$fit[[nm]]$draws("log_reg_slopes", format = "draws_matrix")
  log_reg_slopes_means <- colMeans(log_reg_slopes_samps)
  log_reg_slopes_cris  <- matrixStats::colQuantiles(log_reg_slopes_samps, probs = c(0.025, 0.5, 0.975))
  colnames(log_reg_slopes_cris) <- c("lower", "median", "upper")
  design_mat <- fits_with_lr$stan_data[[nm]]$design_matrix_covariates
  return(list(log_reg_slopes_means = log_reg_slopes_means, log_reg_slopes_cris = log_reg_slopes_cris, design_mat = design_mat))
})

no_slopes <- length(log_reg_results$orig$log_reg_slopes_means)

covariate_groups <- gsub("[[:digit:]]+", "", colnames(log_reg_results$orig$design_mat))
covariate_groups <- c(covariate_groups, paste("item", rep(1:ni, 2), "time", rep(1:2, each = ni)))
length(covariate_groups)

log_reg_tib <- tibble(
  fit    = rep(names(modelNames), each = no_slopes),
  item   = rep(seq_len(no_slopes), 3),
  group  = rep(covariate_groups, 3),
  type   = rep(ifelse(seq_len(no_slopes) <= 18, "covariate", "item"), 3),
  mean   = unlist(map(log_reg_results, `[[`, "log_reg_slopes_means"), use.names = FALSE),
  lower  = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "lower"]), use.names = FALSE),
  upper  = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "upper"]), use.names = FALSE),
  median = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "median"]), use.names = FALSE)
)

# This figure can become A LOT more informative by labelling the covariates more properly!
# for the covariates, group by categorical level (e.g., crime).
# for the items, group by (1) the questionnaire they come from, (2) time point.
yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(log_reg_tib$lower, log_reg_tib$upper))
ggplot(data = log_reg_tib |> filter(item < 18),
       aes(y = mean, group = interaction(fit, item, group), color = fit, x = group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(.5), width = .1) +
  geom_point(position=position_dodge(.5)) +
  scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw() +
  theme(axis.text.x = element_text(angle = 90))

# Preliminary conclusions:
#  1. It looks like the uncertainty kills any chance of meaningful conclusions.
#  2. Some covariates really stand out.
#  3. Maybe repeat this with only the last time point?

largest_betas <- do.call(rbind, map(log_reg_results, \(x) {
  o <- order(abs(x$log_reg_slopes_means), decreasing = TRUE)
  x$log_reg_slopes_means[o]
}))

largest_betas_inds <- do.call(rbind, map(log_reg_results, \(x) {
  o <- order(abs(x$log_reg_slopes_means), decreasing = TRUE)
  as.numeric(extract_numeric(names(x$log_reg_slopes_means[o])))
}))

item_index <- ifelse(largest_betas_inds > 18, largest_betas_inds - 18L, NA_integer_)

tib <- tibble(
  value  = c(largest_betas),
  fit    = rep(names(modelNames), ncol(largest_betas)),
  type   = c(ifelse(largest_betas_inds <= 18L, "covariate", "item")),
  shape  = factor(item_index)
)


log_reg_results$orig$log_reg_slopes_means
