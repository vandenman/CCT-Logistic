rm(list = ls())
library(CCTLogistic)
library(data.table) # do we need this one?
library(dplyr)
library(purrr)
library(ranger)
library(gbm)
library(brms) # <- do this manually

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

do_ml <- function(data_train, data_test,
                  method = c("random_forest", "gbm", "bayesian_logistic", "frequentist_logistic", "baseline_model"),
                  target = "violent_after") {

  method <- match.arg(method)

  fit <- switch(method,
    "random_forest"        = fit_rf(data_train, target),
    "gbm"                  = fit_gbm(data_train, target),
    "bayesian_logistic"    = fit_logistic_regression_bay(data_train, target, prior = brms::prior(normal(0, 3), class = "b")),
    "frequentist_logistic" = fit_logistic_regression_freq(data_train, target),
    "baseline_model"       = fit_baseline_model(data_train, target)
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
    patient_orig = patient,
    patient = normalize_factor(patient),
    violent_after = as.integer(violent_after) - 1L
  ) |>
  filter(!duplicated(patient))

wider_imputed_data <- read_wider_imputed_data()
data_wider <- mice::complete(wider_imputed_data, 1) |>
  filter(patient %in% data_violence$patient_orig)

# fit ML methods ----

# meta information
n_obs      <- nrow(data_wider)
all_obs    <- seq_len(n_obs)
unique_obs <- n_obs
n_holdout  <- floor(0.2 * unique_obs)

cols <- c("time", "patient", "rater", "item", "score", "patient_age_group",
          "treatement_duration_group",
          "violent_before", "violent_between",
          "violent_after",
          "diagnosis_group", "crime_group", "rater_group"
)
cols <- setdiff(colnames(data_wider), c("time", "Aantal_Patienten", "patient"))

mod_ltm <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels", cpp_options = list(stan_threads=TRUE))

set.seed(123)
no_cross_validations <- 10
seeds <- sample(100000, no_cross_validations)
force <- FALSE

ml_methods <- c("random_forest", "gbm", "bayesian_logistic", "frequentist_logistic", "baseline_model")
objs <- matrix(list(), length(ml_methods) + 1L, no_cross_validations, dimnames = list(c(ml_methods, "CCT"), NULL))

for (i in seq_len(no_cross_validations)) {

  set.seed(seeds[i])
  idx_unique_obs <- sample(unique_obs, n_holdout)
  idx_holdout    <- which(all_obs %in% idx_unique_obs)

  data_train <- data_wider[-idx_holdout, cols]
  data_test  <- data_wider[ idx_holdout, cols]
  cat(sprintf("nrow(data_train) = %d\nnrow(data_test)  = %d\nfold = %s", nrow(data_train), nrow(data_test), i), sep = "\n")

  for (method in ml_methods) {
    cat("method:", method, "\n")
    objs[[method, i]] <- CCTLogistic::save_or_run_model(
      do_ml(data_train, data_test, method),
      path = sprintf("fitted_objects/prediction_fitting/%s-%d.rds", method, i),
      force = force
    )
  }
  cat("method: ltm\n")

  ltm_data <- data_2_stan(
    data_2_analyze,
    logistic_dat = data_violence, logistic_target = "violent_after",
    store_predictions = TRUE, debug = FALSE, nc = 18,
    use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE,
    missing_idx = sort(idx_holdout)
  )

  ltm_fit <- CCTLogistic::save_or_run_model(
    mod_ltm$variational(data = ltm_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, threads = 8),
    path = sprintf("fitted_objects/prediction_fitting/cct-%d.rds", i), force = force
  )

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

train_accuracy <- matrix(unlist(map(objs, \(x) x$perf$train$prediction_accuracy)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
test_accuracy  <- matrix(unlist(map(objs, \(x) x$perf$test $prediction_accuracy)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
rowMeans(train_accuracy)
rowMeans(test_accuracy)

train_mse <- matrix(unlist(map(objs, \(x) x$perf$train$mse)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
test_mse  <- matrix(unlist(map(objs, \(x) x$perf$test $mse)), nrow(objs), ncol(objs), FALSE, dimnames(objs))
rowMeans(train_mse)
rowMeans(test_mse)


debugonce(brms::brm)
ff <- brms::brm(formula = violent_after ~ ., data = data_train,
          family = brms::bernoulli(link = "logit"),
          iter       = 5000,
          warmup     = 2000,
          chains     = 6,
          backend    = "cmdstanr"
          # file       = "saved_stan_models/logistic_fit.RData",
          # save_model = "saved_stan_models/logistic_fit.stan"
)
