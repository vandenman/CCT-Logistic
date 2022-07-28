rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(gridExtra)
library(cmdstanr)
library(dplyr)
library(purrr)

# TODO:

# - [ ] try vectorizing over raters and do
#  target += ordered_logistic_lpdf(scores | locations[idx_rater[r]] / scales[idx_rater[r]], thresholds[r])
#
# the above does not work because it's thresholds[r] / scales and there is no vectorized cutpoints implementation
# target += ordered_logistic_lpdf(scores | locations[idx_rater[r]] / scales[idx_rater[r]], thresholds[r])
#
# - [ ] try to reorder and avoid the if-else branches in ordered_logistic_lpdf?
# perhaps precompute three vectors, isOne, isNc, and isOther
# then loop over isOne with log1m_inv_logit, isNc with log_inv_logit, and the rest with log_inv_logit_diff.

# - [ ] try to rewrite this to decorrelate log_E and log_a
# log_inv_logit((location - threshold_shift - threshold_scale * default_thresholds[x - 1]) / scale)
# where scale = exp(log_lambda[i] - log_E[r]),
#       threshold_scale = exp(log_a[r]),
#       threshold_shift = b[r]
#       location = lt[p, i] +/- lt_offset[p, i]
#
# perhaps
# (lt[p, i] +/- lt_offset[p, i] - b[r] - a[r] * default_thresholds[x - 1]) / (lambda[i] / E[r])
# ((lt[p, i] +/- lt_offset[p, i]) / E[r] - b[r] - a[r] * default_thresholds[x - 1]) / lambda[i]
# ?



# functions ----
fit_all_three_models <- function(data, model, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, debug = TRUE, force = FALSE,
                                 path_prefix = "", store_predictions = TRUE, threads = 8L) {

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  stan_data_orig <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  stan_data_skew <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

  path_orig <- file.path("fitted_objects", sprintf("%soptimized_orig_thresholds.rds", path_prefix))
  path_skew <- file.path("fitted_objects", sprintf("%soptimized_skew_thresholds.rds", path_prefix))
  path_free <- file.path("fitted_objects", sprintf("%soptimized_free_thresholds.rds", path_prefix))

  cat("Fitting original threshold model\n")
  fit_orig <- save_or_run_model(
    model$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads),
    path_orig, force
  )
  cat("Fitting skew threshold model\n")
  fit_skew <- save_or_run_model(
    model$variational(data = stan_data_skew, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads),
    path_skew, force
  )
  cat("Fitting free threshold model\n")
  fit_free <- save_or_run_model(
    model$variational(data = stan_data_free, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads),
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

get_performance_tib <- function(dat, fit, threshold_name = "threshold", ...) {
  UseMethod("get_performance_tib", dat)
}

get_performance_tib.logistic_regression_ltm_data <- function(dat, fit, threshold_name) {
  nc <- get_nc(dat)
  nr <- get_nr(dat)
  no_covariates <- get_no_covariates(dat)
  get_means_tib(fit) |>
    add_implied_thresholds(nc, nr) |>
    add_true_values_to_tib(dat) |>
    split_slopes_on_covariate_vs_lt(no_covariates) |>
    mutate(parameter = if_else(parameter == "free_thresholds", threshold_name, parameter))
}

get_performance_tib.ltm_data <- function(dat, fit, threshold_name) {
  nc <- get_nc(dat)
  nr <- get_nr(dat)
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

get_fit_ltm_plot <- function(fit, name, data) {
  nc <- get_nc(data)

  obs_proportions   <- get_obs_proportions(get_score(data), nc = nc)
  log_probs <- fit$draws("log_probs", format = "draws_matrix")
  mean_probs <- matrix(apply(log_probs, 2L, function(x) mean(exp(x))), nrow = nc)
  model_probabilities  <- rowMeans(mean_probs)
  inspect_probability_plot2(obs_proportions, model_probabilities, title = name)

}

fit_single_dataset <- function(dataset, model, data_name, debug = TRUE, store_predictions = TRUE) {

  path_prefix <- paste0("data_", data_name)
  fits <- fit_all_three_models(dataset, model, path_prefix = path_prefix, debug = debug, store_predictions = store_predictions)

  perf_tib <- purrr::map2(fits, c("2-parameter_thresholds", "skew_thresholds", "free_thresholds"), get_performance_tib, dat = dataset)

  graphs_retrieval <- map2(fits, perf_tib, get_scatterplot_retrieval)
  graphs_retrieval_joined   <- arrangeGrob(grobs = graphs_retrieval, ncol = 3)

  graphs_fit_ltm <- graphs_fit_ltm_joined <- NULL
  if (debug) {
    graphs_fit_ltm        <- purrr::map2(fits, names(fits), get_fit_ltm_plot, data = dataset)
    graphs_fit_ltm_joined <- arrangeGrob(grobs = matrix(unlist(lapply(graphs_fit_ltm, `[[`, "grobs"), recursive = FALSE), 3, 2, TRUE), nrow = 2, ncol = 3)
  }

  graphs_predictions_log_reg <- graphs_predictions_log_reg_joined <- NULL
  if (store_predictions && is.logistic_regression_ltm_data(dataset)) {
    graphs_predictions_log_reg <- purrr::map2(fits, names(fits), get_predictions_plot, data = dataset)
    graphs_predictions_log_reg_joined <- arrangeGrob(grobs = matrix(unlist(graphs_predictions_log_reg, recursive = FALSE), 3, 2, TRUE), nrow = 2, ncol = 3)
  }

  return(list(
    fits                              = fits,
    perf_tib                          = perf_tib,

    graphs_retrieval                  = graphs_retrieval,
    graphs_fit_ltm                    = graphs_fit_ltm,
    graphs_predictions_log_reg        = graphs_predictions_log_reg,

    graphs_retrieval_joined           = graphs_retrieval_joined,
    graphs_fit_ltm_joined             = graphs_fit_ltm_joined,
    graphs_predictions_log_reg_joined = graphs_predictions_log_reg_joined
  ))
}

# simulate all data sets ----
np <- 70             # no patients
ni <- 8              # no items
nr <- 20             # no raters
nc <- 10              # no response categories
no_rater_groups <- 1 # no groups of rater
no_covariates   <- 5   # no additional covariates
no_time_points  <- 2
# for debugging
# np <- 3              # no patients
# ni <- 4              # no items
# nr <- 3              # no raters
# nc <- 5              # no response categories
# no_rater_groups <- 1 # no groups of rater
# no_covariates   <- 5   # no additional covariates
# no_time_points  <- 2
# compile the stan model
mod_log_reg_ltm1 <- compile_stan_model("stanmodels/LTM_optimized.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE),
                                      compile = TRUE)

mod_log_reg_ltm2 <- compile_stan_model("stanmodels/LTM_optimized_noLogE.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE),
                                      compile = TRUE)

mod_log_reg_ltm3 <- compile_stan_model("stanmodels/LTM_loop.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE),
                                       compile = TRUE)
mod_log_reg_ltm4 <- compile_stan_model("stanmodels/LTM_loop_testfile.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE),
                                       compile = TRUE)


set.seed(1234)
threshold_types <- c("logistic", "skew_logistic", "free")
dataset_list0 <- map(c("logistic", "skew_logistic", "free"),
                    \(threshold_type) simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points, threshold_type = threshold_type)
)
names(dataset_list0) <- threshold_types

n_total <- np * ni * nr * no_time_points
n_missing <- round(0.170753 * n_total) # same percentage of missing as in the actual dataset
idx_missing <- sample(n_total, n_missing)
dataset_list <- dataset_list0
dataset_list <- map(dataset_list, \(x) {
  x$missing_data <- x$df[idx_missing, ]
  x$df <- x$df[-idx_missing, ]
  x
})

iter = 3e4; adapt_iter = 500; output_samples = 2e3; grad_samples = 5; elbo_samples = 5; debug = TRUE; force = FALSE;
path_prefix = ""; store_predictions = TRUE; threads = 8L

stan_data_orig <- data_2_stan(dataset_list$logistic, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)

stan_data_orig2 <- stan_data_orig

idx_1     <- which( stan_data_orig2$x == 1L)
idx_nc    <- which( stan_data_orig2$x == nc)
idx_other <- which(!stan_data_orig2$x %in% c(1L, nc))

stan_data_orig2$n_observed_1     <- length(idx_1)
stan_data_orig2$idx_rater_1      <- stan_data_orig2$idx_rater[idx_1]
stan_data_orig2$idx_item_1       <- stan_data_orig2$idx_item[idx_1]
stan_data_orig2$idx_patient_1    <- stan_data_orig2$idx_patient[idx_1]
stan_data_orig2$idx_time_point_1 <- stan_data_orig2$idx_time_point[idx_1]

stan_data_orig2$n_observed_nc     <- length(idx_nc)
stan_data_orig2$idx_rater_nc      <- stan_data_orig2$idx_rater[idx_nc]
stan_data_orig2$idx_item_nc       <- stan_data_orig2$idx_item[idx_nc]
stan_data_orig2$idx_patient_nc    <- stan_data_orig2$idx_patient[idx_nc]
stan_data_orig2$idx_time_point_nc <- stan_data_orig2$idx_time_point[idx_nc]

stan_data_orig2$n_observed     <- length(idx_other)
stan_data_orig2$x              <- stan_data_orig2$x[idx_other]
stan_data_orig2$idx_rater      <- stan_data_orig2$idx_rater[idx_other]
stan_data_orig2$idx_item       <- stan_data_orig2$idx_item[idx_other]
stan_data_orig2$idx_patient    <- stan_data_orig2$idx_patient[idx_other]
stan_data_orig2$idx_time_point <- stan_data_orig2$idx_time_point[idx_other]



# fit1 <- mod_log_reg_ltm1$variational(data = stan_data_orig,  iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
fit3 <- mod_log_reg_ltm3$variational(data = stan_data_orig,  iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
fit4 <- mod_log_reg_ltm4$variational(data = stan_data_orig2,  iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
# fit1$time()
fit3$time()
fit4$time()
system("sound.sh 0")

samps3 <- fit3$draws(format = "draws_matrix")
samps4 <- fit4$draws(format = "draws_matrix")

layout(1)
plot(colMeans(samps3[, -1]), colMeans(samps4[, -1]))

h <- .1
probs <- seq(h, 1 - h, h)
nn <- 6
start <- 2
idx <- start:(start+nn-1)
layout(matrix(1:(nn*nn), nn, nn))
for (i in seq(2, by = 1, length.out = nn*nn)) {
  plot(quantile(samps3[, i], probs = probs), quantile(samps4[, i], probs = probs))
  abline(0, 1)
}

plot_cor_log_E_log_a <- function(fit) {
  samps <- exp(fit$draws(variables = c("log_E[1]", "log_a[1]"), format = "draws_matrix"))
  plot(samps[, 1], samps[, 2], main = sprintf("cor(log_E[1], log_a[1]) = %.3f", cor(samps[, 1], samps[, 2])))
}
plot_cor_log_E_log_a(fit3)
plot_cor_log_E_log_a(fit4)
colnames(sss)
plot(sss[, "log_E[1]"], sss[, "log_a[1]"])
cor(sss[, "log_E[1]"], sss[, "log_a[1]"])



starting_values_logistic <- function(data, threshold_type = c("logistic", "free")) {

  sscale <- function(x) {
    scale(x)
    attr(x, "scaled:center") <- attr(x, "scaled:scale") <- NULL
    if (NCOL(x) == 1L || NROW(x) == 1L)
      x <- c(x)
    x
  }

  nc <- vctrs::vec_unique_count(data$score)
  nr <- vctrs::vec_unique_count(data$rater)
  if (!is.factor(data$score))
    data$score <- factor(data$score, levels = seq_len(nc))

  start_lt <- tapply(as.integer(data$score), list(data$patient, data$item), mean)
  start_mu_lt <- mean(start_lt)
  start_sd_lt <- sd(start_lt)
  start_lt <- sscale(start_lt)

  if ("time_point" %in% colnames(data) && vctrs::vec_unique_count(data$time_point) == 2L) {
    temp0 <- tapply(as.integer(data$score), list(data$patient, data$item, data$time_point), mean)
    start_offset_lt <- temp0[, , 2] - start_lt
    # start_lt + start_offset_lt == temp0[, , 2]

  } else {
    start_offset_lt <- NULL
  }
  start_log_lambda <- sscale(log(tapply(as.integer(data$score), list(data$item),  sd)))

  probs_obs <- do.call(cbind, tapply(data$score, data$rater, \(x) { tb <- table(x); tb / sum(tb) }))

  gamma <- 1:(nc - 1)
  gamma <- log(gamma / (nc - gamma))


  start_log_E <- numeric(nr)
  start_log_a <- numeric(nr)
  start_b     <- numeric(nr)

  data_split <- split(data, data$rater)
  mean_location <- sapply(data_split, \(x) mean(start_lt[cbind(x$patient, x$item)]))
  mean_scale    <- exp(sapply(data_split, \(x) mean(start_log_lambda[x$item])))
  # r <- 1
  # x <- c(starting_values$log_E[r], starting_values$log_a[r], starting_values$b[r])
  # mean_location_r <- mean_location[r]
  # mean_scale_r <- mean_scale[r]
  for (r in seq_len(nr)) {
    opt <- optim(c(0, 0, 1), \(x, probs_obs, mean_location_r, mean_scale_r, gamma) {
      E_r <- exp(x[1])
      a_r <- exp(x[2])
      b_r <- x[3]

      delta <- a_r * gamma + b_r

      probs_est <- CCTLogistic:::get_probs_ordered(delta, mean_location_r, mean_scale_r / E_r)

      sum((probs_obs - probs_est)^2)
    }, probs_obs = probs_obs[, r], mean_location = mean_location[r], mean_scale = mean_scale[r], gamma = gamma)

    start_log_E[r] <- opt$par[1L]
    start_log_a[r] <- opt$par[2L]
    start_b[r]     <- opt$par[3L]

  }

  start_sd_log_lambda <- sd(start_log_lambda)

  data_rater_group   <- split(data, data$rater_group)
  start_mu_b         <- sapply(data_rater_group, \(d) mean(start_b[unique(d$rater)]))
  start_sd_b         <- sapply(data_rater_group, \(d) sd(start_b[unique(d$rater)]))
  start_mu_log_E     <- sapply(data_rater_group, \(d) mean(start_log_E[unique(d$rater)]))
  start_log_sd_log_E <- sapply(data_rater_group, \(d) sd(start_log_E[unique(d$rater)]))
  start_log_sd_log_a <- sapply(data_rater_group, \(d) mean(start_log_a[unique(d$rater)]))

  return(list(
    lt             = start_lt,
    log_lambda     = start_log_lambda,
    offset_lt      = start_offset_lt,
    log_E          = start_log_E,
    log_a          = start_log_a,
    b              = start_b,
    # hyperparameters
    mu_lt          = start_mu_lt,
    sd_lt          = start_sd_lt,
    sd_log_lambda  = start_sd_log_lambda,
    mu_b           = start_mu_b,
    sd_b           = start_sd_b,
    mu_log_E       = start_mu_log_E,
    log_sd_log_E   = start_log_sd_log_E,
    log_sd_log_a   = start_log_sd_log_a,
    free_thresholds = numeric(),
    threshold_shape = numeric(),

    log_reg_intercept = numeric(),
    log_reg_slopes    = numeric()
  ))
}

starting_values_2_raw <- function(starting_values) {

  np <- nrow(starting_values$lt)
  ni <- ncol(starting_values$lt)
  nr <- length(starting_values$log_E)

  QR_np <- CCTLogistic:::Q_sum_to_zero_QR_R(np)
  QR_ni <- CCTLogistic:::Q_sum_to_zero_QR_R(ni)
  QR_nr <- CCTLogistic:::Q_sum_to_zero_QR_R(nr)

  starting_values$lt_raw <- matrix(NA, np - 1, ni)
  for (i in seq_len(ni))
    starting_values$lt_raw[, i]  <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$lt[, i],    QR_np)
  starting_values$log_lambda_raw <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$log_lambda, QR_ni)
  starting_values$log_E_raw      <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$log_E,      QR_nr)
  starting_values$log_a_raw      <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$log_a,      QR_nr)

  if (length(starting_values$mu_log_E) > 1L) {
    starting_values$mu_log_E_raw     <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$mu_log_E)
    starting_values$log_sd_log_E_raw <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$log_sd_log_E)
    starting_values$log_sd_log_a_raw <- CCTLogistic:::sum_to_zero_QR_inv_R(starting_values$log_sd_log_a)
  } else {
    starting_values$mu_log_E_raw     <- numeric()
    starting_values$log_sd_log_E_raw <- numeric()
    starting_values$log_sd_log_a_raw <- numeric()
  }

  starting_values

}



debugonce(starting_values_logistic)

starting_values <- starting_values_logistic(dataset_list$logistic$df)
diag(cor(starting_values$lt, dataset_list$logistic$parameters$lt))
diag(cor(starting_values$offset_lt, dataset_list$logistic$parameters$offset_lt))
cor(starting_values$log_lambda, dataset_list$logistic$parameters$log_lambda)
cor(starting_values$log_E, dataset_list$logistic$parameters$log_E)
cor(starting_values$log_a, dataset_list$logistic$parameters$log_a)
cor(starting_values$b, dataset_list$logistic$parameters$b)


ddd <- dataset_list$logistic
CCTLogistic::log_likelihood_ltm(
  ddd,
  log_a           = ddd$parameters$log_a,
  b               = ddd$parameters$b,
  lt              = ddd$parameters$lt,
  log_E           = ddd$parameters$log_E,
  log_lambda      = ddd$parameters$log_lambda,
  offset_lt       = ddd$parameters$offset_lt,
  threshold_shape = ddd$parameters$threshold_shape,
  free_thresholds = ddd$parameters$free_thresholds,
  threshold_type  = ddd$threshold_type
)
CCTLogistic::log_likelihood_ltm(
  ddd,
  log_a           = starting_values$log_a,
  b               = starting_values$b,
  lt              = starting_values$lt,
  log_E           = starting_values$log_E,
  log_lambda      = starting_values$log_lambda,
  offset_lt       = starting_values$offset_lt,
  threshold_shape = NULL,
  free_thresholds = NULL,
  threshold_type  = ddd$threshold_type
)

#
# ddd <- dataset_list$logistic
# CCTLogistic::log_likelihood_ltm(
#   ddd,
#   log_a           = ddd$parameters$log_a,
#   b               = ddd$parameters$b,
#   lt              = ddd$parameters$lt,
#   log_E           = ddd$parameters$log_E,
#   log_lambda      = ddd$parameters$log_lambda,
#   offset_lt       = ddd$parameters$offset_lt,
#   threshold_shape = ddd$parameters$threshold_shape,
#   free_thresholds = ddd$parameters$free_thresholds,
#   threshold_type  = ddd$threshold_type
# )
# ddd2 <- ddd
# a_r
# CCTLogistic::log_likelihood_ltm(
#   ddd2,
#   log_a           = log(exp(2 * ddd2$parameters$log_a)),
#   b               = ddd2$parameters$b,
#   lt              = ddd2$parameters$lt,
#   log_E           = ddd2$parameters$log_E,
#   log_lambda      = ddd2$parameters$log_lambda,
#   offset_lt       = ddd2$parameters$offset_lt,
#   threshold_shape = ddd2$parameters$threshold_shape,
#   free_thresholds = ddd2$parameters$free_thresholds,
#   threshold_type  = ddd2$threshold_type
# )


starting_values_raw <- starting_values_2_raw(starting_values)

fit0 <- mod_log_reg_ltm$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
fit1 <- mod_log_reg_ltm$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads,
                            init = list(starting_values_raw))
fit0$time()
fit1$time()

fit00 <- mod_log_reg_ltm$sample(data = stan_data_orig, iter_sampling = output_samples, chains = 1, threads_per_chain = 16)
fit00$time()

sss <- fit00$draws(format = "draws_matrix")
colnames(sss)
plot(sss[, "log_E[1]"], sss[, "log_a[1]"])
cor(sss[, "log_E[1]"], sss[, "log_a[1]"])

fit2 <- mod_log_reg_ltm2$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
fit2$time()

fit3 <- mod_log_reg_ltm3$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads)
fit4 <- mod_log_reg_ltm3$variational(data = stan_data_orig, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples, threads = threads,
                                     init = list(starting_values_raw[lengths(starting_values_raw) > 0]))
fit3$time()
fit4$time()
fit4$output()
# debugonce(fit_all_three_models)
fit_test <- fit_all_three_models(dataset_list$logistic, mod_log_reg_ltm, path_prefix = "test", force = TRUE)
# new (does not use threads)
# 47.3, 98.4, 36.5
# old (with threads)
# 72.9, 78.5, 36.7

mod_log_reg_ltm_profile <- compile_stan_model("stanmodels/LTM_optimized_with_profile.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE),
                                      compile = TRUE)
fit_test <- fit_all_three_models(dataset_list$logistic, mod_log_reg_ltm_profile, path_prefix = "test_profile", force = TRUE)

profile_tib <- rbind(
  fit_test$fit_orig$profiles()[[1]] |> as_tibble() |> mutate(fit = "fit_orig") |> relocate(fit) |> select(-thread_id),
  fit_test$fit_skew$profiles()[[1]] |> as_tibble() |> mutate(fit = "fit_skew") |> relocate(fit) |> select(-thread_id),
  fit_test$fit_free$profiles()[[1]] |> as_tibble() |> mutate(fit = "fit_free") |> relocate(fit) |> select(-thread_id)
)
print(profile_tib, n = 30)

compute_imputation_performance <- function(fit, data) {
  draws_orig <- fit$draws()
  cnms <- colnames(draws_orig)
  idx <- which(startsWith(cnms, "predicted_missings["))

  orig_missing_posterior_modes <- numeric(length(idx))
  orig_missing_posterior_probs <- matrix(NA_real_, nc, length(idx))
  for (i in seq_along(idx)) {
    tb <- table(factor(draws_orig[, idx[i]], levels = 1:nc))
    orig_missing_posterior_modes[i]   <- which.max(tb)
    orig_missing_posterior_probs[, i] <- tb / sum(tb)
  }
  tb_est <- table(data$missing_data$score, orig_missing_posterior_modes)
  brier_score <- sum(orig_missing_posterior_probs[cbind(data$missing_data$score, 1:ncol(orig_missing_posterior_probs))]) / length(orig_missing_posterior_probs)
  c("classification" = sum(diag(tb_est)) / sum(tb_est),
    "Brier score" = brier_score)
}
map2(fit_test, rep(dataset_list["logistic"], 3), compute_imputation_performance)
# $fit_orig
# classification    Brier score
#     0.50936768     0.08416276
#
# $fit_skew
# classification    Brier score
#     0.50468384     0.08451967
#
# $fit_free
# classification    Brier score
#      0.5035129      0.0838363

# classification of ~0.5 is good (a lot better than  guessing, 0.2)

# fit all 9 models
results <- pmap(list(dataset = dataset_list, data_name = names(dataset_list)), fit_single_dataset, model = mod_log_reg_ltm)

grid.arrange(results$logistic$graphs_retrieval_joined)
grid.arrange(results$logistic$graphs_fit_ltm_joined)
grid.arrange(results$skew_logistic$graphs_retrieval_joined)
grid.arrange(results$skew_logistic$graphs_fit_ltm_joined)
grid.arrange(results$free$graphs_retrieval_joined)
grid.arrange(results$free$graphs_fit_ltm_joined)

# add logistic regression
set.seed(5678)
dataset_log_reg_list <- map(dataset_list, \(dat) simulate_logistic_regression(dat, no_covariates = no_covariates, slopes = seq_len(no_covariates)))
names(dataset_log_reg_list) <- paste0("log_reg_", names(dataset_log_reg_list))

# fit all 9 models
# debugonce(fit_single_dataset)
results <- pmap(list(dataset = dataset_log_reg_list, data_name = names(dataset_log_reg_list)), fit_single_dataset, model = mod_log_reg_ltm)


# first fit the 9 models without logistic regression
dataset_list <- list(
  orig = dat_log_reg_orig,
  skew = dat_log_reg_skew,
  free = dat_log_reg_free
)

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



# fit model with missings
mod_log_reg_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels")

set.seed(9123)
missing_idx <- sort(sample(np, 10))
non_missing_idx <- (1:np)[-missing_idx]

dat_log_reg_logistic_true <- dataset_log_reg_list$log_reg_logistic
dat_log_reg_logistic <- dat_log_reg_logistic_true

debugonce(data_2_stan)
stan_dat_with_missing <- data_2_stan(dat_log_reg_logistic, missing_idx = missing_idx, store_predictions = TRUE)
str(stan_dat_with_missing[c("log_reg_outcomes", "log_reg_np_idx", "np_log_reg", "design_matrix_covariates")])

fit_with_missing <- mod_log_reg_ltm$variational(data = stan_dat_with_missing,
  iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5
)

log_reg_predictions <- fit_with_missing$draws("log_reg_predictions", format = "draws_matrix")
mean_pred_probs <- unname(colMeans(log_reg_predictions))
post_preds <- as.integer(mean_pred_probs >= .5)
table(dat_log_reg_logistic_true$log_reg$y, post_preds)
table(dat_log_reg_logistic_true$log_reg$y[missing_idx], post_preds[missing_idx])
table(dat_log_reg_logistic_true$log_reg$y[non_missing_idx], post_preds[non_missing_idx])
system("beep_finished.sh")


# same experiment, equal patients and items as in actual data set
np <- 103            # no patients
ni <- 23             # no items
nr <- 10             # no raters
nc <- 5              # no response categories
no_rater_groups <- 1 # no groups of rater
no_covariates   <- 5 # no additional covariates
no_time_points  <- 2

set.seed(4567)
dat_missing_new <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points, threshold_type = "free")
dat_missing_new_log_reg <- simulate_logistic_regression(dat_missing_new, no_covariates = no_covariates, slopes = seq_len(no_covariates))
missing_idx <- sort(sample(np, floor(0.2*np)))
non_missing_idx <- (1:np)[-missing_idx]

stan_dat_with_missing <- data_2_stan(dat_missing_new_log_reg, missing_idx = missing_idx, store_predictions = TRUE)

fit_with_missing <- mod_log_reg_ltm$variational(data = stan_dat_with_missing,
  iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5
)

log_reg_predictions <- fit_with_missing$draws("log_reg_predictions", format = "draws_matrix")
mean_pred_probs <- unname(colMeans(log_reg_predictions))
post_preds <- as.integer(mean_pred_probs >= .5)

compute_accuracy <- function(x, y) {
  tb <- table(x, y)
  sum(diag(tb)) / sum(tb)
}

table(dat_missing_new_log_reg$log_reg$y, post_preds)
compute_accuracy(dat_missing_new_log_reg$log_reg$y, post_preds)
table(dat_missing_new_log_reg$log_reg$y[missing_idx], post_preds[missing_idx])
compute_accuracy(dat_missing_new_log_reg$log_reg$y[missing_idx], post_preds[missing_idx])
table(dat_missing_new_log_reg$log_reg$y[non_missing_idx], post_preds[non_missing_idx])
compute_accuracy(dat_missing_new_log_reg$log_reg$y[non_missing_idx], post_preds[non_missing_idx])
system("beep_finished.sh")

