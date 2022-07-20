rm(list = ls())
library(CCTLogistic)
library(dplyr)
library(purrr)
library(ranger)
library(gbm)
library(brms)
library(tidyr)
library(pROC)
library(ggplot2)
library(mice)

# functions ----
get_prediction_accuracy <- function(table) sum(diag(table)) / sum(table)

get_confusion_table <- function(predictions, observed, maxCategories = 2, round = FALSE) {

  if (round)
    predictions <- round(predictions)
  r <- factor(levels = 0:(maxCategories - 1))
  # we add e.g., 1:5, 1:5 to both observations and predictions so that all outcomes appear in the table
  # next we subtract 1 from the diagonal to remove these observations.
  join <- function(a, b) unlist(list(a, b))
  return(table(join(r, predictions), join(r, observed)))# - diag(1, maxCategories, maxCategories))

}

fit_rf <- function(data, target = "violent_after") {
  rf_obj <- ranger::ranger(
    dependent.variable.name = target,
    data = data, oob.error = TRUE, verbose = TRUE, num.trees = 1e4, importance = "impurity", probability = TRUE)
  return(rf_obj)
}

fit_gbm <- function(data, target = "violent_after") {

  data[[target]] <- as.character(data[[target]])
  f <- as.formula(sprintf("%s ~ .", target))
  return(gbm::gbm(formula = f, data = data, distribution = "bernoulli",
                  cv.folds = 5, n.cores = 6, shrinkage = 0.001, interaction.depth = 3,
                  n.trees = 1000))

}

fit_logistic_regression_bay <- function(data, target = "violent_after", i, ...) {

  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  # data[[target]] <- as.character(data[[target]])
  f <- as.formula(sprintf("%s ~ .", target))
  return(brms::brm(formula = f, data = data, family = brms::bernoulli(link = "logit"),
                   iter       = 5000,
                   warmup     = 2000,
                   chains     = 6,
                   cores      = 8,
                   silent     = 2,
                   refresh    = 1000,
                   save_model = "saved_stan_models/logistic_fit.stan",
                   backend    = "cmdstanr",
                   ...
  ))

}

fit_logistic_regression_freq <- function(data, target = "violent_after") {

  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  f <- as.formula(sprintf("%s ~ .", target))
  return(glm(f, data = data, family = binomial()))

}

fit_baseline_model <- function(data, target = "violent_after") {
  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  f <- as.formula(sprintf("%s ~ 1", target)) # intercept only
  return(glm(f, data = data, family = binomial()))
}

fit_baseline_violence_only_model <- function(data, target = "violent_after") {
  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  f <- as.formula(sprintf("%s ~ 1 + violent_before + violent_between", target)) # intercept only
  return(glm(f, data = data, family = binomial()))
}

fit_baseline_no_item_model <- function(data, target = "violent_after") {
  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  predictors <- colnames(data)
  predictors <- predictors[!startsWith(predictors, "IFBE_")]
  predictors <- setdiff(predictors, target)
  f <- as.formula(sprintf("%s ~ 1 + %s", target, paste(predictors, collapse = " + "))) # intercept only
  return(glm(f, data = data, family = binomial()))
}

fit_baseline_no_violence <- function(data, target = "violent_after") {
  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  predictors <- colnames(data)
  predictors <- setdiff(predictors, c("violent_before", "violent_between"))
  predictors <- setdiff(predictors, target)
  f <- as.formula(sprintf("%s ~ 1 + %s", target, paste(predictors, collapse = " + "))) # intercept only
  return(glm(f, data = data, family = binomial()))
}

do_ml <- function(data_train, data_test,
                  method = c("random_forest", "gbm", "bayesian_logistic", "frequentist_logistic", "baseline_model", "baseline_violence_only", "baseline_no_item", "baseline_no_violence"),
                  target = "violent_after") {

  method <- match.arg(method)

  fit <- switch(method,
    "random_forest"          = fit_rf(data_train, target),
    "gbm"                    = fit_gbm(data_train, target),
    "bayesian_logistic"      = fit_logistic_regression_bay(data_train, target, prior = brms::prior(normal(0, 3), class = "b")),
    "frequentist_logistic"   = fit_logistic_regression_freq(data_train, target),
    "baseline_model"         = fit_baseline_model(data_train, target),
    "baseline_violence_only" = fit_baseline_violence_only_model(data_train, target),
    "baseline_no_item"       = fit_baseline_no_item_model(data_train, target),
    "baseline_no_violence"   = fit_baseline_no_violence(data_train, target),
  )
  perf <- assess_performance(fit, data_train, data_test)
  list(fit = fit, perf = perf)
}

summarize_predictive_performance <- function(predictions, observed) {
  if (is.list(predictions))
    predictions <- predictions[["predictions"]]
  confusion_table     <- get_confusion_table(predictions, observed)
  prediction_accuracy <- get_prediction_accuracy(confusion_table)
  return(list(confusion_table = confusion_table, prediction_accuracy = prediction_accuracy))
}

compute_mse <- function(predicted_probs, outcomes_factor) {
  outcomes <- as.integer(outcomes_factor) - 1L
  return(mean((predicted_probs - outcomes)^2))
}

compute_log_loss <- function(predicted_probs, outcomes_factor) {
  outcomes <- as.integer(outcomes_factor) - 1L
  return(sum(dbinom(outcomes, 1, predicted_probs, log = TRUE)))
}
# TODO: refactor into a prediction S3 method that returns the probabilities and use one function for assess_performance
assess_performance <- function(obj, data_train, data_test, target = "violent_after", ...) {
  UseMethod("assess_performance", obj)
}

assess_performance2 <- function(predicted_probs_train, predicted_probs_test, data_train, data_test, target) {

  predicted_outcome_train <- factor(as.integer(predicted_probs_train >= 0.5), levels = levels(data_train[[target]]))
  predicted_outcome_test  <- factor(as.integer(predicted_probs_test  >= 0.5), levels = levels(data_test[[target]]))

  performance_train <- summarize_predictive_performance(predicted_outcome_train, data_train[[target]])
  performance_test  <- summarize_predictive_performance(predicted_outcome_test,  data_test[[target]])

  mse_train <- compute_mse(predicted_probs_train, data_train[[target]])
  mse_test  <- compute_mse(predicted_probs_test,   data_test[[target]])

  log_loss_train <- compute_log_loss(predicted_probs_train, data_train[[target]])
  log_loss_test  <- compute_log_loss(predicted_probs_test,   data_test[[target]])

  return(list(
    train = c(performance_train, list(
      predicted_probs = predicted_probs_train,
      mse             = mse_train,
      log_loss        = log_loss_train
    )),
    test = c(performance_test, list(
      predicted_probs = predicted_probs_test,
      mse             = mse_test,
      log_loss        = log_loss_test
    ))
  ))
}

assess_performance.gbm <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_train_temp <- data_train
  data_test_temp  <- data_test

  data_train_temp[[target]] <- as.character(data_train[[target]])
   data_test_temp[[target]] <- as.character( data_test[[target]])

  predicted_probs_train <- predict(obj, data_train_temp, type = "response")
  predicted_probs_test  <- predict(obj, data_test_temp,  type = "response")

  return(assess_performance2(predicted_probs_train, predicted_probs_test, data_train, data_test, target))

}

assess_performance.ranger <- function(obj, data_train, data_test, target = "violent_after", ...) {

  raw_predictions_rf_train <- predict(obj, data_train)
  raw_predictions_rf_test  <- predict(obj, data_test)

  predicted_probs_train <- raw_predictions_rf_train[["predictions"]][, 2L]
  predicted_probs_test  <- raw_predictions_rf_test[["predictions"]][, 2L]

  return(assess_performance2(predicted_probs_train, predicted_probs_test, data_train, data_test, target))

}

assess_performance.brmsfit <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_test  <- data_test[,  which(!names(data_train) %in% c("rater", "patient", "item"))]
  data_train <- data_train[, which(!names(data_train) %in% c("rater", "patient", "item"))]

  predicted_probs_train <- predict(obj, data_train)[, "Estimate"]
  predicted_probs_test  <- predict(obj, data_test)[, "Estimate"]

  return(assess_performance2(predicted_probs_train, predicted_probs_test, data_train, data_test, target))

}

assess_performance.glm <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_test  <- data_test[,  which(!names(data_train) %in% c("rater", "patient", "item"))]
  data_train <- data_train[, which(!names(data_train) %in% c("rater", "patient", "item"))]

  predicted_probs_train <- predict(obj, data_train, type = "response")
  predicted_probs_test  <- predict(obj, data_test, type = "response")

  return(assess_performance2(predicted_probs_train, predicted_probs_test, data_train, data_test, target))

}

fitRandomforest <- function(dat_train, dat_test) {

  tStart <- Sys.time()
  dat_train <- as.data.frame(lapply(dat_train, factor))
  dat_test  <- as.data.frame(lapply(dat_test, factor))

  rfObj <- ranger::ranger(dependent.variable.name = "Outcome", data = dat_train,
                          oob.error = TRUE, verbose = TRUE, num.trees = 1e4)

  preds <- predict(rfObj, data = dat_test[, colnames(dat_test) != "Outcome"])
  confusionTable <- get_confusion_table(preds$predictions, dat_test[, colnames(dat_test) == "Outcome"])
  predAccuracy <- get_prediction_accuracy(confusionTable)
  tStop <- Sys.time()
  return(list(
    # trainAccuracy   = sum(diag(rfObj$confusion.matrix)) / sum(rfObj$confusion.matrix),
    predAccuracy    = predAccuracy,
    confusionTable  = confusionTable,
    rfObj           = rfObj,
    time            = tStop - tStart
  ))
}

fitBoosting <- function(datTrain, datTest, n.trees = 500, cv.folds = 5, shrinkage = 1e-3,
                        train = FALSE) {

  tStart <- Sys.time()

  # don't attempt to optimize hyperparameters
  gbmObj <- gbm(
    formula = Outcome ~ .,
    data    = datTrain,
    distribution = "multinomial"
  )

  predsProbs <- predict(gbmObj, newdata = datTest[, colnames(datTest) != "Outcome"], n.trees = gbmObj$n.trees,
                    type = "response")
  preds <- apply(predsProbs, 1L, which.max)
  confusionTable <- makeConfusionTable(preds, datTest[, colnames(datTest) == "Outcome"])
  predAccuracy <- getPredictionAccuracy(confusionTable)
  tStop <- Sys.time()
  return(list(
    predAccuracy    = predAccuracy,
    confusionTable  = confusionTable,
    gbmObj          = gbmObj,
    time            = tStop - tStart
  ))
}

# load data ----
all_data <- read_long_data()
data_2_analyze <- all_data |>
  filter(
    # !is.na(score) &
    !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(-c(age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  arrange(rater_group, patient, item, rater, time)

data_violence <- all_data |>
  filter(
    # !is.na(score) &
    !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(c(patient, age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  filter(!duplicated(patient))

measure_vars <- paste0("IFBE_", 0:22)
data_wide <- read_wide_data()
data_wider <- data_wide |>
  filter(!is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  pivot_wider(names_from = time, values_from = all_of(measure_vars), names_sep = "_time_") |>
  group_by(patient) |>
  mutate(
    across(
      all_of(paste0(measure_vars, "_time_", rep(1:2, length(measure_vars)))),
      ~ mean(.x, na.rm = TRUE)
    )
  ) |>
  ungroup() |>
  distinct(patient, .keep_all = TRUE) |>
  select(-c(rater, rater_group))

set.seed(31415)
mice_results <- mice::mice(data_wider, m = 1, blocks = colnames(data_wider)[-1])
data_wider_imputed <- complete(mice_results) |> as_tibble()
anyNA(complete(data_wider_imputed)) # FALSE

# usedToBeNa <- which(is.na(data_wider), arr.ind = TRUE)
# data_wider[usedToBeNa[1, 1], usedToBeNa[1, 2]]
# data_wider_imputed[usedToBeNa[1, 1], usedToBeNa[1, 2]]
# data_wider_imputed[usedToBeNa[2, 1], usedToBeNa[2, 2]]


# assert that data_wider and data_violence are identical, except that data_wider contains the IFBE items in wide format
validate_datasets <- function(data_wider_imputed, data_violence, measure_vars) {
  shared_colnames <- intersect(colnames(data_wider_imputed), colnames(data_violence))
  assertthat::assert_that(
    nrow(data_wider_imputed) == nrow(data_violence),
    ncol(data_wider_imputed) == ncol(data_violence) + 2 * length(measure_vars),
    identical(data_wider_imputed[shared_colnames], data_violence[shared_colnames]), # data are identical except for the IFBE items in data_wider
    !anyNA(data_wider_imputed[setdiff(colnames(data_wider_imputed), shared_colnames)]) # no NAs in the IFBE items
  )
}
validate_datasets(data_wider_imputed, data_violence, measure_vars)

data_nonmissing <- data_2_analyze[!is.na(data_2_analyze$score), ]
data_missing    <- data_2_analyze[ is.na(data_2_analyze$score), ]
data_missing$score <- NULL
# factor model doesn't fit very well
# data_for_pca <- data_wider[, seq(9, 8+23*2, 2)]
# ff <- psych::fa.parallel(data_for_pca)
# data_for_cfa <- data_for_pca[, -1]
# cfa_fit <- lavaan::cfa(model = "
#   prot =~ IFBE_1_time_1 + IFBE_2_time_1 + IFBE_3_time_1 + IFBE_4_time_1 + IFBE_17_time_1 + IFBE_19_time_1 + IFBE_21_time_1 + IFBE_22_time_1
#   reso =~ IFBE_5_time_1 + IFBE_6_time_1 + IFBE_7_time_1 + IFBE_8_time_1 + IFBE_9_time_1
#   prob =~ IFBE_10_time_1 + IFBE_11_time_1 + IFBE_12_time_1 + IFBE_13_time_1 + IFBE_14_time_1 + IFBE_15_time_1 + IFBE_16_time_1 + IFBE_18_time_1 + IFBE_20_time_1
#                        ", data = data_for_cfa)
#
# lavaan::fitmeasures(cfa_fit) # rmsea 0.164, cfi 0.533, tli 0.476

# fit ML methods ----

# meta information
n_obs      <- nrow(data_wider_imputed)
all_obs    <- seq_len(n_obs)
unique_obs <- n_obs
n_holdout  <- floor(0.2 * unique_obs)

cols <- c("time", "patient", "rater", "item", "score", "age",
          "treatment_duration",
          "violent_before", "violent_between",
          "violent_after",
          "diagnosis", "crime", "rater_group"
)
cols <- setdiff(colnames(data_wider_imputed), c("time", "Aantal_Patienten", "patient"))

mod_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time_and_missing.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads=TRUE))

# set.seed(123)
set.seed(42)
no_cross_validations <- 10

idx_violent    <- which(data_wider_imputed$violent_after == 1)
idx_nonviolent <- which(data_wider_imputed$violent_after == 0)
idx_violent_split    <- split(idx_violent,    cut(seq_along(idx_violent),    no_cross_validations, labels = FALSE))
idx_nonviolent_split <- split(idx_nonviolent, cut(seq_along(idx_nonviolent), no_cross_validations, labels = FALSE))
cv_indices <- unname(Map(\(x, y) sort(c(x, y)), idx_violent_split, idx_nonviolent_split))


ml_methods <- c("random_forest", "gbm", "frequentist_logistic", "baseline_model", "baseline_violence_only", "baseline_no_item", "baseline_no_violence")
objs <- matrix(list(), length(ml_methods) + 1L, no_cross_validations, dimnames = list(c(ml_methods, "CCT"), NULL))
data_test_objs <- vector("list", length = no_cross_validations)

seeds <- matrix(sample(100000, no_cross_validations * nrow(objs)), nrow(objs), no_cross_validations, dimnames = dimnames(objs))
force <- FALSE

# prediction_fitting_imputed fits all methods to the imputed data
fitted_objects_dir <- file.path("fitted_objects", "prediction_fitting_imputed")
if (!dir.exists(fitted_objects_dir)) dir.create(fitted_objects_dir)

for (i in seq_len(no_cross_validations)) {

  idx_holdout <- cv_indices[[i]]

  data_train <- data_wider_imputed[-idx_holdout, cols]
  data_test  <- data_wider_imputed[ idx_holdout, cols]
  data_test_objs[[i]] <- data_test
  cat(sprintf("nrow(data_train) = %d\nnrow(data_test)  = %d\nfold = %s", nrow(data_train), nrow(data_test), i), sep = "\n")

  for (method in ml_methods) {
    cat("method:", method, "\n")
    set.seed(seeds[method, i])
    objs[[method, i]] <- CCTLogistic::save_or_run_model(
      do_ml(data_train, data_test, method),
      path = file.path(fitted_objects_dir, sprintf("%s-%d.rds", method, i)),
      force = force
    )
  }
  cat("method: ltm\n")

  ltm_data <- data_2_stan(
    data_nonmissing,
    logistic_dat = data_violence, logistic_target = "violent_after",
    missing_data = data_missing,
    store_predictions = TRUE, debug = FALSE, nc = 17,
    use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE,
    missing_idx = sort(idx_holdout)
  )

  ltm_retries <- 5
  for (r in 0:ltm_retries) {
    set.seed(seeds["CCT", i] + r)
    ltm_fit <- try(CCTLogistic::save_or_run_model(
      mod_ltm$variational(data = ltm_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, threads = 8),
      path = file.path(fitted_objects_dir, sprintf("cct-%d.rds", i)), force = force
    ))

    # browser()
    # TODO: this check doesn't work!
    if (!inherits(ltm_fit, "try-error") && (identical(ltm_fit$return_codes(), 0) || identical(ltm_fit$return_codes(), NA_integer_))) {
      system("beep_finished.sh 0")
      break
    } else if (r == ltm_retries) {
      system("beep_finished.sh 1")
      stop("maximum ltm_retries reached!")
    } else {
      system("beep_finished.sh 1")
    }
  }

  log_reg_predictions <- ltm_fit$draws("log_reg_predictions", format = "draws_matrix")
  mean_pred_probs <- unname(colMeans(log_reg_predictions))

  ltm_perf <- assess_performance2(
    predicted_probs_train = mean_pred_probs[-idx_holdout],
    predicted_probs_test  = mean_pred_probs[ idx_holdout],
    data_train            = data_train,
    data_test             = data_test,
    target                = "violent_after"
  )

  ltm_obj <- list(fit = ltm_fit, perf = ltm_perf)

  objs[["CCT", i]] <- ltm_obj

}

system("beep_finished.sh")

# names used in manuscript ----
name_map <- c(baseline_model         = "LR-Intercept",
              baseline_no_item       = "LR-No IFTE",
              baseline_no_violence   = "LR-No violence",
              baseline_violence_only = "LR-Violence",
              frequentist_logistic   = "LR",
              gbm                    = "GBM",
              random_forest          = "Random forest",
              CCT                    = "LR-LTM")

method_order <- c(
  "LR-LTM",
  "LR-Violence",
  "LR-No IFTE",
  "LR-No violence",
  "Random forest",
  "GBM",
  "LR",
  "LR-Intercept"
)

# plot ROC curves ----
get_tpr_tnr <- function(observed, predicted) {
  tb <- table(observed, predicted)
  tnr <- tb[1L, 1L] / (tb[1L, 1L] + tb[1L, 2L])
  tpr <- tb[2L, 2L] / (tb[2L, 1L] + tb[2L, 2L])
  return(c(tpr, tnr))
}

compute_tpr_tnr_mat <- function(observed, predicted_probs, thresholds) {
  mm <- matrix(NA_real_, length(thresholds), 2)
  for (i in seq_along(thresholds)) {
    predicted <- factor(as.integer(thresholds[i] <= predicted_probs), levels = levels(observed))
    mm[i, ] <- get_tpr_tnr(observed, predicted)
  }

  se_mm <- mm[, 1]#sort(mm[, 1], decreasing = TRUE)
  sp_mm <- mm[, 2]#sort(mm[, 2], decreasing = FALSE)
  tibble(tpr = se_mm, tfr = 1 - sp_mm, thresholds = thresholds)
}

compute_tpr_tnr_mat_from_obj <- function(objs, data_test_objs, thresholds, method, cv_fold) {
  probs <- objs[[method, cv_fold]]$perf$test$predicted_probs
  labels   <- data_test_objs[[cv_fold]]$violent_after
  compute_tpr_tnr_mat(labels, probs, thresholds)
}

# validate against pROC::roc
rf_probs1 <- objs[["random_forest", 1]]$perf$test$predicted_probs
labels1   <- data_test_objs[[1]]$violent_after
roc_obj_rf1 <- roc(as.integer(labels1), rf_probs1, plot = TRUE, legacy.axes=TRUE, xlab = "False positive rate", ylab = "True positive rate")
roc_obj_rf1_remade <- compute_tpr_tnr_mat_from_obj(objs, data_test_objs, roc_obj_rf1$thresholds, "random_forest", 1)
roc_obj_rf1$sensitivities - roc_obj_rf1_remade$tpr
roc_obj_rf1$specificities - (1 - roc_obj_rf1_remade$tfr)
plot(roc_obj_rf1)
auc(roc_obj_rf1)

plot_roc <- function(objs, data_test_objs, method, cv_fold) {
  probs <- objs[[method, cv_fold]]$perf$test$predicted_probs
  labels   <- data_test_objs[[cv_fold]]$violent_after
  roc_obj <- roc(as.integer(labels), probs, plot = TRUE, legacy.axes=TRUE, xlab = "False positive rate", ylab = "True positive rate", direction = "<")
  invisible(roc_obj)
}
oo <- plot_roc(objs, data_test_objs, "random_forest", 1)

no_reps <- map_int(data_test_objs, \(x) length(x$violent_after))


objs_tib <- tibble(
  # method     = rep(rownames(objs), each = no_reps, ncol(objs)),
  # iter       = rep(seq_len(ncol(objs)), each = nrow(objs), no_reps),
  method     = unlist(lapply(no_reps, \(r) rep(rownames(objs), each = r))),
  iter       = rep(seq_len(ncol(objs)), nrow(objs) * no_reps),
  test_probs = unlist(lapply(objs, \(x) x$perf$test$predicted_probs)),
  test_label = unlist(lapply(seq_along(data_test_objs), \(i) rep(data_test_objs[[i]]$violent_after, nrow(objs))))
)

# TODO: compute roc once
roc_thresholds_method <- objs_tib |>
  group_by(method) |>
  summarize(thresholds = roc(test_label, test_probs, plot = FALSE)$thresholds)

roc_auc_method_iter <- objs_tib |>
  group_by(method, iter) |>
  summarize(auc = c(roc(test_label, test_probs, plot = FALSE)$auc))

roc_auc_method <-roc_auc_method_iter |>
  group_by(method) |>
  summarize(mean_auc = mean(auc))

# roc_thresholds <- roc_obj$thresholds
# roc_thresholds <- c(-Inf, seq(0.01, 0.99, length.out = 40), Inf)

# roc_tib2 <- objs_tib |>
#   group_by(method, iter) |>
#   summarise(
#     temp = roc(
#       as.integer(data_test_objs[[iter[[1L]]]]$violent_after),
#       objs[[method[[1L]], iter[[1L]]]]$perf$test$predicted_probs,
#       plot = FALSE,
#       direction = "<"
#     )[c("sensitivities", "specificities", "auc")],
#     tpr = temp$sensitivities,
#     tfr = 1 - temp$specificities,
#     auc = temp$auc
#   ) |>
#   select(-temp)
#   ungroup() |>
#   unnest_wider(temp)

roc_tib <- objs_tib |>
  group_by(method, iter) |>
  summarise(
    temp = compute_tpr_tnr_mat(test_label, test_probs, roc_thresholds_method$thresholds[roc_thresholds_method$method == method[1L]]),
  ) |>
  ungroup() |>
  unnest_wider(temp) |>
  group_by(method, iter) |>
  arrange(tfr, tpr)


roc_mean_tib <- roc_tib |>
  group_by(method, thresholds) |>
  summarise(mean_tpr = mean(tpr), mean_tfr = mean(tfr))

# filt <- function(data) data |> filter(method %in% c("CCT", "baseline_violence_only"))
# filt <- function(data) data |> filter(!(method %in% c("CCT", "baseline_violence_only")))
filt <- identity

roc_tib_to_plt_tib <- function(tib, method_order, name_map2, name_map) {
  skip <- c("LR-No violence", "LR-Intercept")
  idx1 <- match(skip, names(name_map2))
  skip <- name_map2[idx1]

  method_order_idx <- match(method_order, names(name_map2))
  method_order <- name_map2[method_order_idx]

  legend_order <- setdiff(method_order, skip)
  names(name_map2) <- names(name_map)
  tib |>
    mutate(method = recode(method, !!!name_map2)) |>
    filter(!(method %in% skip)) |>
    mutate(method = factor(method, levels = legend_order))
}


auc_idx <- match(names(name_map), roc_auc_method$method)
name_map2 <- setNames(sprintf("%s (%.3f)", name_map, roc_auc_method$mean_auc[auc_idx]), name_map)

roc_tib_plt <- roc_tib_to_plt_tib(roc_tib,           method_order, name_map2, name_map)
roc_mean_tib_plt <- roc_tib_to_plt_tib(roc_mean_tib, method_order, name_map2, name_map)

roc_plot <- ggplot(
  data = roc_tib_plt |> filt(),
  # aes(x = tfr, y = tpr, group = interaction(method, iter), color = factor(iter))
  aes(x = tfr, y = tpr, group = interaction(method, iter), color = method)
  # aes(x = tfr, y = tpr, group = interaction(method, iter), color = factor(iter), linetype = method)
) +
  jaspGraphs::geom_abline2(intercept = 0, slope = 1, color = "grey") +
  # geom_line(alpha = .3) + # <- showing the individual curves does not help
  geom_line(data = roc_mean_tib_plt |> filt(), aes(x = mean_tfr, y = mean_tpr, color = method, group = method), size = 0.9) +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::scale_JASPcolor_discrete(name = "Method (AUC)") +
  labs(x = "False positive rate", y = "True positive rate") +
  jaspGraphs::themeJaspRaw(legend.position = c(0.9, 0.50))
roc_plot

# a rough double check for the AUC (values do not exactly match roc_auc_method, nor do they need to)
# roc_mean_tib |>
#   group_by(method) |>
#   summarise(
#     temp = {
#       yy <- mean_tpr
#       xx <- mean_tfr
#       sy <- split(yy, xx)
#       xx <- as.numeric(names(sy))
#       yy <- sapply(sy, max)
#       list(list(xx = xx, yy = yy))
#     },
#     af  = list(approxfun(temp[[1]]$xx, temp[[1]]$yy, yleft = 0, yright = 1)),
#     auc = integrate(af[[1]], 0, 1)$value
#   )
roc_auc_method

save_figure(roc_plot, "roc_curve_cv.svg")

# roc_tibble <- tibble()
# for (i in seq_along(data_train_objs)) {
#   test_labels <- as.integer(data_train_objs[[i]]$violent_after)
#   for (j in seq_len(nrow(objs))) {
#     temp_roc_obj <- pROC::roc(test_labels, objs[[j, i]]$perf$test$predicted_probs)
#     temp_tib <- tibble(
#       method    = rownames(objs)[j],
#       iteration = i,
#       tpr       =     sort(temp_roc_obj$sensitivities, decreasing = TRUE),
#       fpr       = 1 - sort(temp_roc_obj$specificities, decreasing = FALSE)
#     )
#     roc_tibble <- rbind(roc_tibble, temp_tib)
#   }
# }
#
# i <- j <- 1
# test_labels <- as.integer(data_train_objs[[i]]$violent_after)
# temp_roc_obj <- pROC::roc(test_labels, objs[[j, i]]$perf$test$predicted_probs)
# tt1 <- tibble(
#   method    = rownames(objs)[j],
#   iteration = i,
#   tpr       =     sort(temp_roc_obj$sensitivities, decreasing = TRUE),
#   fpr       = 1 - sort(temp_roc_obj$specificities, decreasing = FALSE)
# )
# temp_roc_obj2 <- cutpointr(x = objs[[j, i]]$perf$test$predicted_probs, class = test_labels, direction = ">=", pos_class = 2)
# tt2 <- temp_roc_obj2 |> select(c("roc_curve")) |>
#     tidyr::unnest(.data$roc_curve) |>
#     select(c("tpr", "tnr")) |>
#     mutate(fnr = 1 - tnr) |>
#     select(-tnr)
#
# tt1$tpr - rev(tt2$tpr)
# tt1$fpr - rev(tt2$fnr)
#
# roc_tibble$tpr - rev(individual_roc_tib$tpr)
#
# ggplot(data = roc_tibble, aes(x = fpr, y = tpr, group = interaction(method, iteration), color = method)) +
#   jaspGraphs::geom_abline2(intercept = 0, slope = 1, color = "grey") +
#   geom_line(alpha = .5) +
#   jaspGraphs::geom_rangeframe() +
#   jaspGraphs::scale_JASPcolor_discrete() +
#   jaspGraphs::themeJaspRaw(legend.position = "right")
#
# library(ggplot2)
# library(cutpointr)
#
# predictions_100_samples <- data.frame(
#     Sample = rep(c(1:100), times = 195),
#     PredictionValues = c(rnorm(n = 9750), rnorm(n = 9750, mean = 1)),
#     RealClass = c(rep("benign", times = 9750), rep("pathogenic", times = 9750))
# )
#
# roc_tib <- tibble()
# individual_roc_tib <- tibble()
# mean_roc_tib <- tibble()
# cutpoints <- seq(0, 1, .05)
# for (j in seq_len(nrow(objs))) {
#
#   method <- rownames(objs)[j]
#   temp_tib <- tibble(
#     probs  = unlist(lapply(seq_len(nrow(objs)), \(i) objs[[i, j]]$perf$test$predicted_probs)),
#     class  = unlist(lapply(seq_len(nrow(objs)), \(i) as.integer(data_train_objs[[i]]$violent_after))),
#     cv     = unlist(lapply(seq_len(nrow(objs)), \(i) rep(i, length(data_train_objs[[i]]$violent_after)))),
#     method = method
#   )
#   roc_tib <- rbind(roc_tib, temp_tib)
#
#   temp_individual_roc_tib <- cutpointr(data = temp_tib, x = probs, class = class, subgroup = cv, direction = ">=", pos_class = 2)
#   plot_roc(temp_individual_roc_tib)
#
#   temp_individual_roc_tib$method <- method
#   temp_individual_roc_tib <- temp_individual_roc_tib |>
#     select(c("roc_curve", "subgroup", "method")) |>
#     tidyr::unnest(.data$roc_curve) |>
#     select(c("tpr", "tnr", "subgroup", "method")) |>
#     mutate(
#       tpr =     sort(tpr, decreasing = TRUE),
#       fpr = 1 - sort(tnr, decreasing = FALSE)
#     ) |>
#     select(-tnr)
#   individual_roc_tib <- rbind(individual_roc_tib, temp_individual_roc_tib)
#
#   temp_mean_roc_tib <- map_df(cutpoints, function(cp) {
#     out <- cutpointr(data = temp_tib, x = probs, class = class, method = oc_manual, cutpoint = cp, direction = ">=", subgroup = cv, pos_class = 2)
#     tibble(cutoff = cp,
#            sensitivity = mean(out$sensitivity),
#            specificity = mean(out$specificity))
#   })
#   temp_mean_roc_tib$method <- method
#   mean_roc_tib <- rbind(mean_roc_tib, temp_mean_roc_tib)
#
# }
#
# ggplot(data = individual_roc_tib, aes(x = fpr, y = tpr, group = interaction(method, subgroup), color = method)) +
#   jaspGraphs::geom_abline2() +
#   geom_line(alpha = .5) +
#   geom_line(data = mean_roc_tib, mapping = aes(x = 1 - specificity, y = sensitivity, group = method, color = method), inherit.aes = FALSE) +
#   jaspGraphs::geom_rangeframe() +
#   jaspGraphs::themeJaspRaw(legend.position = "right")
#

# performance tables ----
rowSds <- matrixStats::rowSds
train_accuracy <- matrix(unlist(map(objs, \(x) x$perf$train$prediction_accuracy)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
test_accuracy  <- matrix(unlist(map(objs, \(x) x$perf$test $prediction_accuracy)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
rowMeans(train_accuracy)
rowMeans(test_accuracy)
rowSds(train_accuracy)
rowSds(test_accuracy)

train_mse <- matrix(unlist(map(objs, \(x) x$perf$train$mse)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
test_mse  <- matrix(unlist(map(objs, \(x) x$perf$test $mse)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
rowMeans(train_mse)
rowMeans(test_mse)
matrixStats::rowSds(train_mse)
matrixStats::rowSds(test_mse)


rownames <- rownames(test_mse)# c("random forest", "gbm", "Bayesian LR", "Frequentist LR", "Intercept only LR", "Bayesian LR + CCT")
order <- order(rowMeans(test_mse))# c(5, 4, 3, 6, 1, 2)

method_names <- name_map[match(rownames, names(name_map))]
order <- match(method_order, method_names)
method_names[order]
skip <- c("LR-No violence")
order <- order[!(method_names[order] %in% skip)]

tb_accuracy <- tibble(
  # temp                = rownames,
  method              = method_names,
  mean_train_accuracy = rowMeans(train_accuracy),
    sd_train_accuracy = rowSds  (train_accuracy),
  mean_test_accuracy  = rowMeans(test_accuracy ),
    sd_test_accuracy  = rowSds  (test_accuracy ),
)[order, ]
tb_mse <- tibble(
  # temp                = rownames,
  method              = method_names,
  mean_train_accuracy = rowMeans(train_mse),
    sd_train_accuracy = rowSds  (train_mse),
  mean_test_accuracy  = rowMeans(test_mse ),
    sd_test_accuracy  = rowSds  (test_mse ),
)[order, ]

tb_joined <- cbind(tb_accuracy, tb_mse[-1])

format_tb <- function(tb) {
  tibble(
    method = tb[[1]],
    train_acc  = sprintf("%.3f (%.3f)", tb[[2]], tb[[3]]),
    test_acc   = sprintf("%.3f (%.3f)", tb[[4]], tb[[5]]),
    train_mse  = sprintf("%.3f (%.3f)", tb[[6]], tb[[7]]),
    test_mse   = sprintf("%.3f (%.3f)", tb[[8]], tb[[9]])
  )
}

do_xtable <- function(tb) {
  print(xtable::xtable(tb), include.rownames = FALSE, floating = FALSE, sanitize.text.function = identity)
}

tb_joined |>
  format_tb() |>
  do_xtable()

# tb_mse |>
#   format_tb() |>
#   do_xtable()

print(xtable::xtable(format_tb(tb_mse)),      include.rownames = FALSE, floating = FALSE, sanitize.text.function = identity)
