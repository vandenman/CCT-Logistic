rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)
library(dplyr)
library(purrr)

options("cores" = as.integer(parallel::detectCores() / 2))
# options("cores" = parallel::detectCores())
fit_all_three_models <- function(data, model, logistic_dat = NULL, logistic_target = NULL, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 25, elbo_samples = 25, threads = getOption("cores", 1L), debug = FALSE, force = FALSE,
                                 path_prefix = "", store_predictions = TRUE, max_retries = 5L) {

  if (path_prefix != "" && !endsWith(path_prefix, "_"))
    path_prefix <- paste0(path_prefix, "_")

  data_nonmissing <- data[!is.na(data$score), ]
  data_missing    <- data[ is.na(data$score), ]
  data_missing$score <- NULL

  stan_data_orig <- data_2_stan(data_nonmissing, missing_data = data_missing, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)
  # stan_data_skew <- data_2_stan(data_nonmissing, missing_data = data_missing, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = TRUE,  use_free_logistic_thresholds = FALSE)
  stan_data_free <- data_2_stan(data_nonmissing, missing_data = data_missing, logistic_dat = logistic_dat, logistic_target = logistic_target, store_predictions = store_predictions, debug = debug, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE)

  path_orig <- file.path("fitted_objects", sprintf("%stime_3_models_orig_thresholds.rds", path_prefix))
  # path_skew <- file.path("fitted_objects", sprintf("%stime_3_models_skew_thresholds.rds", path_prefix))
  path_free <- file.path("fitted_objects", sprintf("%stime_3_models_free_thresholds.rds", path_prefix))

  fit_model <- function(data) {
    model$variational(data = data, iter = iter, adapt_iter = adapt_iter, output_samples = output_samples, grad_samples = grad_samples, elbo_samples = elbo_samples,
                      threads = threads)
  }

  cat(sprintf("Using %s threads\n", threads))
  cat("Fitting original threshold model\n")
  fit_orig <- save_or_run_model(fit_model(stan_data_orig), path_orig, force, max_retries = max_retries)

  # cat("Fitting skew threshold model\n")
  # fit_skew <- save_or_run_model(fit_model(stan_data_skew), path_skew, force, max_retries = max_retries)

  cat("Fitting free threshold model\n")
  fit_free <- save_or_run_model(fit_model(stan_data_free), path_free, force, max_retries = max_retries)

  return(list(
    fit = list(
      orig = fit_orig,
      # skew = fit_skew,
      free = fit_free
    ),
    stan_data = list(
      orig = stan_data_orig,
      # skew = stan_data_skew,
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

data_2_analyze <- read_ifte_data()
data_violence <- read_violence_data()

data_nonmissing <- data_2_analyze[!is.na(data_2_analyze$score), ]
data_missing    <- data_2_analyze[ is.na(data_2_analyze$score), ]
data_missing$score <- NULL

# ensure that both datasets contain the same patients
all(unique(data_2_analyze$patient) == unique(data_violence$patient))

colnames(data_2_analyze)
# [1] "time"        "patient"     "rater"       "item"        "score"       "rater_group"
colnames(data_violence)

np <- nlevels(droplevels(data_2_analyze$patient))
ni <- nlevels(droplevels(data_2_analyze$item))
nr <- nlevels(droplevels(data_2_analyze$rater))
nt <- nlevels(droplevels(data_2_analyze$time))
nc <- 17L

mod_ltm <- compile_stan_model("stanmodels/LTM_loop.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels",
                              cpp_options = list(stan_threads=TRUE))
# mod_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time_and_missing.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels",
#                               cpp_options = list(stan_threads=TRUE))

fits_without_lr <- fit_all_three_models(data_2_analyze, mod_ltm, NULL,          NULL,            path_prefix = file.path("ltm_only", "mesdag_ltm_without_logistic_regression"))
# debugonce(fit_all_three_models)
fits_with_lr    <- fit_all_three_models(data_2_analyze, mod_ltm, data_violence, "violent_after", path_prefix = file.path("ltm_only", "mesdag_ltm_with_logistic_regression"))

mean_probs_without_lr <- compute_mean_probs(fits_without_lr, file.path("fitted_objects", "ltm_only", "mesdag_ltm_probs_without_logistic_regression.rds"))
mean_probs_with_lr    <- compute_mean_probs(fits_with_lr,    file.path("fitted_objects", "ltm_only", "mesdag_ltm_probs_with_logistic_regression.rds"))
# system("sound.sh 0")

# inspect LTM fit ----
observed_proportions <- get_obs_proportions(data_2_analyze$score)
rmse <- function(x, y) sqrt(mean((x - y)^2))

raw_probabilities_without_lr <- unlist(map(mean_probs_without_lr, rowMeans), use.names = FALSE)
raw_probabilities_with_lr    <- unlist(map(mean_probs_with_lr,    rowMeans), use.names = FALSE)

plot_ltm_fit <- function(raw_probabilities, observed_proportions, ylim = NULL) {
  nc <- length(observed_proportions)
  fitIndex <- which(c("orig", "skew", "free") %in% names(mean_probs_without_lr))
  nFits <- length(fitIndex)
  betterNames <- c("Original", "Skew", "Free")[fitIndex]
  tib <- tibble(
    panel = rep(c("Raw", "Diff"), c((nFits + 1L) * nc, nFits * nc)),
    fill  = c(rep(c(betterNames, "Observed"), each = nc), rep(betterNames, each = nc)),
    y     = c(raw_probabilities, observed_proportions, raw_probabilities - rep(observed_proportions, nFits)),
    x     = rep(1:nc, 2L * nFits + 1L)
  )
  cols <- hcl.colors(3L, "Set 2")[fitIndex] # keeps the same colors even if we drop one fit
  colors <- c("Observed" = "black", setNames(cols, betterNames))

  if (is.null(ylim)) {
    breaks_fun <- jaspGraphs::getPrettyAxisBreaks
    limits_fun <- \(x) range(jaspGraphs::getPrettyAxisBreaks(x))
  } else {
    stopifnot(is.list(ylim) && length(ylim) == 2L)
    counter_breaks <<- 0L
    counter_limits <<- 0L
    limits_fun <- \(x) {res <- ylim[[counter_limits + 1]]; counter_limits <<- (counter_limits + 1L) %% 2; res}
    breaks_fun <- \(x) {res <- jaspGraphs::getPrettyAxisBreaks(ylim[[counter_breaks+1L]]); counter_breaks <<- (counter_breaks + 1L) %% 2; res}
  }

  fit_plot <- ggplot(data = tib, mapping = aes(x = x, y = y, group = fill, fill = fill)) +
    geom_bar(position="dodge", stat = "identity", width = .5) +
    facet_wrap(~panel, scales = "free_y") +
    scale_x_continuous(name = "Score", breaks = 1:nc, limits = c(0, 19)) +
    scale_y_continuous(name = "Diff (1) / Probability (2)", breaks = breaks_fun, limits = limits_fun) +
    scale_fill_manual(name = NULL, values = colors) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right") +
    theme(strip.text = element_text(size = 28))
  return(fit_plot)
}

# NOTE: the fit to the IFTE items is approximately identical for the models with and without logistic regression
# pdf("simulation_figures/ltm_only_joined_ggplot22.pdf", width = 16, height = 8)
ylim <- list(c(-.2, .4), c(0, .4))
plot_ltm_fit(raw_probabilities_without_lr, observed_proportions, ylim = ylim) + labs(title = "Without logistic regression")
plot_ltm_fit(raw_probabilities_with_lr, observed_proportions,    ylim = ylim) + labs(title = "With logistic regression")
# dev.off()

# inspect fit logistic regression ----
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

# key_diagnosis <- setNames(
#   c("As_1_Other", "Autism spectrum disorder", "Personality disorder cluster B", "Personality disorder other", "Schizophrenia and other psychotic disorders"),
#   levels(data_violence$diagnosis)
# )
# key_crime <- setNames(
#   c("Arson", "Manslaughter", "Murder", "Power with violence", "Power/medium violence", "Sex offense", "Heavy violence"),
#   levels(data_violence$crime)
# )
#
# key_treatment_duration <- setNames(
#   c("0-2 years", "2-4 years", "4-6 years", "6+ years"),
#   levels(data_violence$treatment_duration)
# )
key_diagnosis <- setNames(
  c("Axis 1", "Autism spectrum disorder", "Personality disorder cluster B", "Other personality disorders", "Schizophrenia and other psychotic disorders"),
  levels(data_violence$diagnosis)
)
key_crime <- setNames(
  c("Arson", "Manslaughter", "Murder", "Violent property crime", "Moderate violence / property crime", "Sex offense", "Aggrevated assault"),
  levels(data_violence$crime)
)
key_treatment_duration <- setNames(
  c("0-2 years", "2-4 years", "4-6 years", "6+ years"),
  levels(data_violence$treatment_duration)
)


covariate_groups <-  c("age", "violent_before", "violent_between", "treatment_duration",  "diagnosis", "crime")
# covariate_groups <- gsub("[[:digit:]]+", "", colnames(log_reg_results$orig$design_mat))
covariate_groups <- c(covariate_groups, "item T1", "item T2")
covariate_levels <- c(colnames(log_reg_results$orig$design_mat), paste("item", rep(1:ni, 2), "time", rep(1:2, each = ni)))
covariate_groups_with_repeats <- rep(
  covariate_groups,
  c(sapply(data_violence[covariate_groups[1:6]], nlevels)-1L, ni, ni)
)
covariate_levels_stripped <- sapply(seq_along(covariate_levels), \(i) {
  gsub(covariate_groups_with_repeats[i], "", covariate_levels[i], fixed = TRUE)
  })
covariate_levels_stripped[19:64] <- 1:ni

covariate_levels_stripped <- covariate_levels_stripped |>
  recode(!!!key_diagnosis) |>
  recode(!!!key_crime) |>
  recode(!!!key_treatment_duration)

nFits <- length(fits_with_lr)
log_reg_tib <- tibble(
  fit    = rep(names(modelNames), each = no_slopes),
  item   = rep(seq_len(no_slopes), nFits),
  group  = factor(rep(covariate_groups_with_repeats, nFits), levels = covariate_groups),
  level  = rep(covariate_levels_stripped, nFits),
  type   = rep(ifelse(seq_len(no_slopes) <= nc, "covariate", "item"), nFits),
  mean   = unlist(map(log_reg_results, `[[`, "log_reg_slopes_means"), use.names = FALSE),
  lower  = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "lower"]), use.names = FALSE),
  upper  = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "upper"]), use.names = FALSE),
  median = unlist(map(log_reg_results, \(x) x[["log_reg_slopes_cris"]][, "median"]), use.names = FALSE)
)
# all(as.character(log_reg_tib$group) == rep(covariate_groups_with_repeats, nFits))

# This figure can become A LOT more informative by labelling the covariates more properly!
# for the covariates, group by categorical level (e.g., crime).
# for the items, group by (1) the questionnaire they come from, (2) time point.
yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(log_reg_tib$lower, log_reg_tib$upper))
ggplot(data = log_reg_tib |> filter(fit == "free"),# |> filter(item < nc),
       aes(y = mean, group = interaction(fit, item, group), x = group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(.5), width = .1) +
  geom_point(position=position_dodge(.5)) +
  scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw() +
  theme(axis.text.x = element_text(angle = 90))

idx_violence_before <- which(rep(covariate_levels  == "violent_before1", nFits))
idx_violence_between <- which(rep(covariate_levels == "violent_between1", nFits))
log_reg_tib$group

log_reg_tib2 <- log_reg_tib |>
  mutate(
    group = droplevels(recode_factor(group, violent_before = "Violent", violent_between = "Violent",
                                     treatment_duration = "Treatment\n\nduration",
                                     crime = "Offense", age = "Age", diagnosis = "Diagnosis",
                                     `item T1` = "IFTE T1", `item T2` = "IFTE T2"
                                     )
                       )
  )
log_reg_tib2$level[idx_violence_before]  <- "before"
log_reg_tib2$level[idx_violence_between] <- "between"
log_reg_tib2$level <- factor(log_reg_tib2$level, unique(log_reg_tib2$level))
log_reg_tib2$level <- recode_factor(log_reg_tib2$level,
  `Autism spectrum disorder`                    = "autism",
  `Personality disorder cluster B`              = "personality B",
  `Personality disorder other`                  = "personality Other",
  `Schizophrenia and other psychotic disorders` = "schizophrenia",

                                                  #Manslaughter
                                                  #Power/medium violence
  `Power with violence`                         = "VPC",
  `Power/medium violence`                       = "MPC",
  `Heavy violence`                              = "AA",
  `Manslaughter`                                = "manslaughter",
  `Murder`                                      = "murder",
  `Sex offense`                                 = "sex offense"
)
tib_segment <- tibble(x = -Inf, xend = -Inf, y = yBreaks[1], yend = yBreaks[length(yBreaks)], group="age")


posterior_mean_95CRI <- ggplot(data = log_reg_tib2 |> filter(fit == "free"),# |> mutate(level = abbreviate(level, minlength = 14)),# |> filter(item < nc),
       aes(y = mean, group = interaction(fit, item, group), x = level)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(.5), width = .3) +
  geom_point(position=position_dodge(.5), size = 3) +
  scale_y_continuous(name = "Posterior mean", breaks = yBreaks, limits = range(yBreaks)) +
  xlab(NULL) +
  facet_grid(.~group,space="free_x", scales="free_x", switch="x") +
  jaspGraphs::geom_rangeframe() +
  # jaspGraphs::geom_rangeframe(sides = "b") +
  # geom_segment(data = tib_segment, mapping = aes(x=x,xend=xend,y=y,yend=yend), inherit.aes = FALSE) +
  # geom_segment(x = -Inf, xend = -Inf, y = yBreaks[1], yend = yBreaks[length(yBreaks)]) +
  jaspGraphs::themeJaspRaw(fontsize = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "grey80"),
    # panel.grid.minor.y = element_line(colour = "grey40"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
    panel.spacing.x = unit(0, "cm"),
    strip.placement = "outside"
  )

# save_figure(figure = posterior_mean_95CRI, file = "posterior_mean_95CRI.svg", width = 15, height = 7.5)

gt = ggplot_gtable(ggplot_build(posterior_mean_95CRI))
idx_null <- which(vapply(gt$widths, \(x) grid::unitType(x) == "null", FUN.VALUE = logical(1)))
gt$widths[idx_null] <- gt$widths[idx_null] * c(3 / 2.2, 4 / 3.2, 7 / 6.2, 4 / 3.2, 5/ 4.2,  7 / 23.2, 7 / 23.2)

# required because Latex renders the text in the svg
gt$heights[8] <- unit(8, "cm")

gt$grobs[[23]]$children$axis$grobs$`1`$children[[1]]$rot <- 0
# for (i in 16:20)
#   gt$grobs[[i]]$children$axis$grobs$`1`$children[[1]]$rot <- 45

gt$grobs[[21]]$children$axis$grobs$`1`$children[[1]]$rot <- 0
gt$grobs[[22]]$children$axis$grobs$`1`$children[[1]]$rot <- 0
gt$grobs[[21]]$children$axis$grobs$`1`$children[[1]]$hjust <- 0.5
gt$grobs[[22]]$children$axis$grobs$`1`$children[[1]]$hjust <- 0.5
gt$grobs[[21]]$children$axis$grobs$`1`$children[[1]]$vjust <- 1.0
gt$grobs[[22]]$children$axis$grobs$`1`$children[[1]]$vjust <- 1.0
gt$grobs[[21]]$children$axis$grobs$`1`$children[[1]]$label[(1:23)[-c(1, 5, 10, 15, 20)]] <- ""
gt$grobs[[22]]$children$axis$grobs$`1`$children[[1]]$label[(1:23)[-c(1, 5, 10, 15, 20)]] <- ""

prefix_newline <- function(x) paste0("\n", x)
for (i in c(21, 22))
  gt$grobs[[i]]$children$axis$grobs$`1`$children[[1]]$label[c(1, 5, 10, 15, 20)] <- prefix_newline(gt$grobs[[i]]$children$axis$grobs$`1`$children[[1]]$label[c(1, 5, 10, 15, 20)])

grid::grid.newpage()
grid::grid.draw(gt)
# save_figure(figure = gt, file = "posterior_mean_95CRI2.svg", width = 15, height = 15)

ggplot(data = log_reg_tib2 |> filter(fit == "free"),# |> filter(item < nc),
       aes(y = mean, group = interaction(fit, item, group), x = group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(.5), width = .1) +
  geom_point(position=position_dodge(.5)) +
  scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw() +
  theme(axis.text.x = element_text(angle = 90))





# multiply slopes of with value of lt
raw_samples <- fits_with_lr$fit$free$draws(format = "draws_matrix")
no_covariates <- fits_with_lr$stan_data$free$no_covariates
colnames(raw_samples)
derived <- array(NA, c(np, ni, nt))
derived2 <- matrix(0.0, nrow(raw_samples), np)
for (p in seq_len(np)) for (i in seq_len(ni)) for (t in seq_len(nt)) {
  lt        <- raw_samples[, sprintf("lt[%d,%d]", p, i)]
  offset_lt <- raw_samples[, sprintf("offset_lt[%d,%d]", p, t)]
  slope     <- raw_samples[, sprintf("log_reg_slopes[%d]", no_covariates + i)]
  tt <- 2*t - 1
  derived[p, i, t] <- mean((lt + tt*offset_lt) * slope)
  derived2[, p] <- derived2[, p] + (lt + tt*offset_lt) * slope
}

all.equal(rowSums(derived), colMeans(derived2)) # consistency check

derived2_cris <- apply(derived2, 2L, quantile, probs = c(0.025, 0.975))
ifte_aggregate_tib <- tibble(
  x        = seq_len(nrow(data_violence)),
  y        = colMeans(derived2),
  ci_lower = derived2_cris["2.5%", ],
  ci_upper = derived2_cris["97.5%", ],
  violent  = ifelse(as.integer(data_violence$violent_after) == 2L, "Violent", "Nonviolent")
)
range(derived2_cris)
rowMeans(apply(derived2_cris, 2L, range))

# with CRIs -- unreadable
yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(ifte_aggregate_tib$y, ifte_aggregate_tib$ci_lower, ifte_aggregate_tib$ci_upper))
ggplot(data = ifte_aggregate_tib, aes(x = x, y = y, ymin = ci_lower, ymax = ci_upper, shape = violent, fill = violent)) +
  geom_errorbar() +
  geom_point(size = 5) +
  scale_y_continuous("Posterior mean", breaks = yBreaks, limits = range(yBreaks)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Patient") +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw(legend.position = "right") +
  scale_y_continuous("Posterior mean", breaks = seq(-4, 4, 2)) +
  ggplot2::coord_cartesian(ylim = c(-4, 4))

# without CRIs (as in manuscript) --
# yBreaks <- jaspGraphs::getPrettyAxisBreaks(ifte_aggregate_tib$y)
yBreaks <- jaspGraphs::getPrettyAxisBreaks(ifte_aggregate_tib$y)
posterior_mean_IFTE_aggregate <- ggplot(data = ifte_aggregate_tib, aes(x = x, y = y, shape = violent, fill = violent)) +
  geom_point(size = 5) +
  scale_y_continuous("Posterior mean", breaks = yBreaks, limits = range(yBreaks)+c(0, 0.1)) +
  scale_shape_manual(name = NULL, values = c(21, 22)) +
  # scale_shape_manual(name = NULL, values = c("N", "V")) +
  scale_fill_manual(name = NULL, values = c("grey", "white")) +
  scale_x_continuous("Patient", breaks = c(1, 25, 50, 75, 104)) +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw(legend.position = c(0.15, 0.995))
posterior_mean_IFTE_aggregate
save_figure(figure = posterior_mean_IFTE_aggregate, file = "posterior_mean_IFTE_aggregate.svg", width = 7, height = 7)


no_betas <- length(log_reg_results$orig$log_reg_slopes_means)
no_beta_covariates <- ncol(log_reg_results$orig$design_mat)
no_beta_items <- no_betas - no_beta_covariates

tib <- tibble(
  value      = unlist(lapply(log_reg_results, \(x) x$log_reg_slopes_means)),
  fit        = rep(names(modelNames), each = no_betas),
  beta_type  = rep(rep(c("covariate", "item"), c(no_beta_covariates, no_beta_items)), 2),
  item_index = rep(c(rep(NA_integer_, no_beta_covariates), 1:23, 1:23), 2),
  time_index = rep(rep(c(NA_integer_, 1L, 2L), c(no_beta_covariates, no_beta_items / 2, no_beta_items / 2)), 2)
)

# order items by effect
tib_item <- tib |> filter(beta_type == "item" & fit == "free")
print(tib_item[order(tib_item$value), ], n = 50)

# TODO:
# lt * beta == (-lt) * (-beta) so the sign might not be identified!
# the same possibly goes for lt and the thresholds, check if
# ordered_logistic_lpmf(x[o] | location / scale, free_thresholds[r] ./ scale)
# is equal to
# ordered_logistic_lpmf(x[o] | -location / scale, reverse(free_thresholds[r]) ./ scale);
# but free_thresholds[r] are ordered and scale is positive so this cannot happen?

# 1  Does the patient show problem insight?                                                                  &   Protective behaviors    &   HKT-R                   \\
# 2  Does the patient cooperate with your treatment?                                                         &   Protective behaviors    &   HKT-R                   \\
# 3  Does the patient admit and take responsibility for the crime(s)?                                        &   Protective behaviors    &   HKT-R                   \\
# 4  Does the patient show adequate coping skills?                                                           &   Protective behaviors    &   HKT-R                   \\
# 5  Does the patient have balanced daytime activities?                                                      &   Resocialization Skills  &   HKT-R                   \\
# 6  Does the patient show sufficient labor skills?                                                          &   Resocialization Skills  &   HKT-R                   \\
# 7  Does the patient show sufficient common social skills?                                                  &   Resocialization Skills  &   HKT-R                   \\
# 8  Does the patient show sufficient skills to take care of oneself?                                        &   Resocialization Skills  &   HKT-R                   \\
# 9  Does the patient show sufficient financial skills?                                                      &   Resocialization Skills  &   Proposed by clinicians  \\
# 10 Does the patient show impulsive behavior?                                                               &   Problematic behavior    &   HKT-R                   \\
# 11 Does the patient show antisocial behavior?                                                              &   Problematic behavior    &   HKT-R                   \\
# 12 Does the patient show hostile behavior?                                                                 &   Problematic behavior    &   HKT-R                   \\
# 13 Does the patient show sexual deviant behavior?                                                          &   Problematic behavior    &   Proposed by clinicians  \\
# 14 Does the patient show manipulative behavior?                                                            &   Problematic behavior    &   Proposed by clinicians  \\
# 15 Does the patient comply with the rules\\ and conditions of the center and/or the treatment?             &   Problematic behavior    &   HKT-R                   \\
# 16 Does the patient have antisocial associates?                                                            &   Problematic behavior    &   HKT-R                   \\
# 17 Does the patient use his medication in a consistent and adequate manner?                                &   Protective behaviors    &   Proposed by clinicians  \\
# 18 Does the patient have psychotic symptoms?                                                               &   Problematic behavior    &   HKT-R                   \\
# 19 Does the patient show skills to prevent drug and alcohol use?                                           &   Protective behaviors    &   ASP                     \\
# 20 Does the patient use any drug or alcohol?                                                               &   Problematic behavior    &   HKT-R                   \\
# 21 Does the patient show skills to prevent physical aggressive behavior?                                   &   Protective behaviors    &   ASP                     \\
# 22 Does the patient show skills to prevent sexual deviant behavior?                                        &   Protective behaviors    &   ASP                     \\

# plot posterior mean of lt vs sample mean of items ----

samples_lt <- fits_with_lr$fit$free$draws(variables = "lt", format = "draws_matrix")
tail(colnames(samples_lt)) # 2nd dim is ni
post_means_lt <- matrix(colMeans(samples_lt), np, ni)

observed_data_means_tib <- data_2_analyze |>
  group_by(patient, item) |>
  summarise(mean = mean(score, na.rm = TRUE)) |>
  ungroup(item) |>
  mutate(mean2 = scale(mean)) |> # - mean(mean, na.rm = TRUE)) |>
  ungroup()

observed_data_means_mat <- matrix(observed_data_means_tib$mean2, np, ni, byrow = TRUE)

pplot <- function(x, y) {
  plot(x, y, main = sprintf("cor = %.3f", cor(x, y, use = "pairwise.complete.obs")))
  abline(0, 1)
}
pplot(post_means_lt[1, ], observed_data_means_mat[1, ])
pplot(post_means_lt[2, ], observed_data_means_mat[2, ])

pplot(c(post_means_lt), c(observed_data_means_mat))

observed_data_means_tib$post_mean_lt <- c(t(post_means_lt))
tib_text <- tibble(
  # x = -2,
  # y = 1.5,
  x = 2,
  y = 2.8,
  label = sprintf("$\\rho = %.3f$", cor(observed_data_means_tib$mean2, observed_data_means_tib$post_mean_lt, use = "pairwise.complete.obs"))
)
observed_data_means_tib$colors <- cut(as.integer(observed_data_means_tib$item), breaks = 5)

prot_idx <- c(1, 2, 3, 4, 17, 19, 21, 22)
reso_idx <- c(5, 6, 7, 8, 9)
prob_idx <- setdiff(1:23, c(prot_idx, reso_idx))
observed_data_means_tib$group <- ifelse(
         observed_data_means_tib$item %in% prot_idx, "Protective behaviors",
  ifelse(observed_data_means_tib$item %in% reso_idx, "Resocialization skills",
                                                     "Problematic behavior")
)

breaks <- -3:3
limits <- range(breaks)
post_mean_lt_vs_item_sample_means <- ggplot(data = observed_data_means_tib, aes(x = post_mean_lt, y = mean2, shape = group, fill = group, color = group)) +
  jaspGraphs::geom_point(alpha = .85, size = 2.5) +
  jaspGraphs::geom_abline2(color = "grey80") +
  geom_text(data = tib_text, aes(x = x, y = y, label = label), size = .35*jaspGraphs::graphOptions("fontsize"), inherit.aes = FALSE) +
  scale_shape_manual(values = 21:23) +
  labs(x = "Posterior mean of $\\theta$", y = "Standardized item mean", shape = NULL, fill = NULL, color = NULL) +
  scale_x_continuous(breaks = breaks, limits = limits) +
  scale_y_continuous(breaks = breaks, limits = limits) +
  jaspGraphs::scale_JASPcolor_discrete() +
  jaspGraphs::scale_JASPfill_discrete() +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw(legend.position = c(.125, .99))
post_mean_lt_vs_item_sample_means
save_figure(figure = post_mean_lt_vs_item_sample_means, file = "post_mean_lt_vs_item_sample_means.svg", width = 7, height = 7)
