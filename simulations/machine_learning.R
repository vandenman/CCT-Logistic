rm(list = ls())
library(data.table)
library(ranger)
library(gbm)
library(brms)

source(file.path("simulations", "utils.R"))
# source(file.path("simulations", "mlhelpers.R")) # <- TODO make this file and put all common ml stuff in there

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
    data = data, oob.error = TRUE, verbose = TRUE, num.trees = 1e4, importance = "impurity")
  return(rf_obj)
}

fit_gbm <- function(data, target = "violent_after") {

  data[[target]] <- as.character(data[[target]])
  f <- as.formula(sprintf("%s ~ .", target))
  return(gbm::gbm(formula = f, data = data, distribution = "bernoulli",
                  cv.folds = 5, n.cores = 6, shrinkage = 0.001, interaction.depth = 3,
                  n.trees = 1000))

}

fit_logistic_regression_bay <- function(data, target = "violent_after", ...) {

  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  # data[[target]] <- as.character(data[[target]])
  f <- as.formula(sprintf("%s ~ .", target))
  return(brms::brm(formula = f, data = data, family = brms::bernoulli(link = "logit"),
                   iter       = 5000,
                   warmup     = 2000,
                   chains     = 6,
                   file       = "saved_stan_models/logistic_fit.RData",
                   save_model = "saved_stan_models/logistic_fit.stan",
                   ...
  ))

    # BAS::bas.glm(f, data = data, family = stats::binomial()))

}

fit_logistic_regression_freq <- function(data, target = "violent_after") {

  data <- data[, which(!names(data) %in% c("rater", "patient", "item"))]
  f <- as.formula(sprintf("%s ~ .", target))
  return(glm(f, data = data, family = binomial()))

}

summarize_predictive_performance <- function(predictions, observed) {
  if (is.list(predictions))
    predictions <- predictions[["predictions"]]
  confusion_table     <- get_confusion_table(predictions, observed)
  prediction_accuracy <- get_prediction_accuracy(confusion_table)
  return(list(confusion_table = confusion_table, prediction_accuracy = prediction_accuracy))
}

assess_performance <- function(obj, data_train, data_test, target = "violent_after", ...) {
  UseMethod("assess_performance", obj)
}

assess_performance.gbm <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_train[[target]] <- as.character(data_train[[target]])
  data_test[[target]] <- as.character(data_test[[target]])

  predictions_gbm_train <- as.factor(round(predict(obj, data_train, type = "response")))
  predictions_gbm_test  <- as.factor(round(predict(obj, data_test, type = "response")))

  performance_gbm_train <- summarize_predictive_performance(predictions_gbm_train, data_train[[target]])
  performance_gbm_test  <- summarize_predictive_performance(predictions_gbm_test,  data_test[[target]])

  return(list(performance_train = performance_gbm_train,
              performance_test  = performance_gbm_test))

}

assess_performance.ranger <- function(obj, data_train, data_test, target = "violent_after", ...) {
  predictions_rf_train <- predict(obj, data_train)
  predictions_rf_test  <- predict(obj, data_test)

  performance_rf_train <- summarize_predictive_performance(predictions_rf_train, data_train[[target]])
  performance_rf_test  <- summarize_predictive_performance(predictions_rf_test,  data_test[[target]])

  return(list(performance_train = performance_rf_train,
              performance_test  = performance_rf_test))

}

assess_performance.brmsfit <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_test  <- data_test[,  which(!names(data_train) %in% c("rater", "patient", "item"))]
  data_train <- data_train[, which(!names(data_train) %in% c("rater", "patient", "item"))]

  predictions_brmsfit_train <- predict(obj, data_train)
  predictions_brmsfit_test  <- predict(obj, data_test)

  predictions_brmsfit_train <- factor(1 * (predictions_brmsfit_train[, "Estimate"] >= 0.5), levels = levels(data_train[[target]]))
  predictions_brmsfit_test  <- factor(1 * (predictions_brmsfit_test[, "Estimate"] >= 0.5),  levels = levels(data_train[[target]]))

  performance_brmsfit_train <- summarize_predictive_performance(predictions_brmsfit_train, data_train[[target]])
  performance_brmsfit_test  <- summarize_predictive_performance(predictions_brmsfit_test,  data_test[[target]])

  return(list(performance_train = performance_brmsfit_train,
              performance_test  = performance_brmsfit_test))

}

assess_performance.glm <- function(obj, data_train, data_test, target = "violent_after", ...) {

  data_test  <- data_test[,  which(!names(data_train) %in% c("rater", "patient", "item"))]
  data_train <- data_train[, which(!names(data_train) %in% c("rater", "patient", "item"))]

  predictions_glm_train <- predict(obj, data_train, type = "response")
  predictions_glm_test  <- predict(obj, data_test,  type = "response")

  performance_glm_train <- summarize_predictive_performance(1 * (predictions_glm_train >= 0.5), data_train[[target]])
  performance_glm_test  <- summarize_predictive_performance(1 * (predictions_glm_test  >= 0.5),  data_test[[target]])

  return(list(performance_train = performance_glm_train,
              performance_test  = performance_glm_test))

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



# reshape so that there is a single outcome per patient! This doesn't make any sense!
# data_long <- read_long_data()
# data_long <- data_long[complete.cases(data_long), ]

# data_wide <- read_wide_data()
wider_imputed_data <- read_wider_imputed_data()
data_wider <- mice::complete(wider_imputed_data, 1)
nrow(data_wider)
# length(unique(data_wide$patient))

mice::md.pattern(data_wider, plot = FALSE)

# data_wider <- data_wider[complete.cases(data_wider), ]
# data_wider_imputed <- mice::mice(data_wider, m = 5, maxit = 50, method = 'pmm', seed = 500)


# create a training data set and a test data set
set.seed(123)
n_obs <- nrow(data_wider)

# all_obs <- paste(data_long$time, data_long$patient)
all_obs <- seq_len(n_obs)
# unique_obs <- unique(all_obs)
unique_obs <- n_obs
n_holdout <- floor(0.2 * unique_obs)

idx_unique_obs <- sample(unique_obs, n_holdout)
idx_holdout    <- which(all_obs %in% idx_unique_obs)

cols <- c("time", "patient", "rater", "item", "score", "patient_age_group",
          "treatement_duration_group",
          "violent_before", "violent_between",
          "violent_after",
          "diagnosis_group", "crime_group", "rater_group"
)
cols <- setdiff(colnames(data_wider), c("time", "Aantal_Patienten", "patient"))

data_train <- data_wider[-idx_holdout, cols]
data_test  <- data_wider[ idx_holdout, cols]
cat(sprintf("nrow(data_train) = %d\nnrow(data_test)  = %d", nrow(data_train), nrow(data_test)), sep = '\n')

rf_obj <- fit_rf(data_train)
sort(importance(rf_obj), TRUE)

rf_perf <- assess_performance(rf_obj, data_train, data_test)

gbm_obj <- fit_gbm(data_train)
sort(gbm::relative.influence(gbm_obj), TRUE)
gbm_perf <- assess_performance(gbm_obj, data_train, data_test)

options(mc.cores = 6)
brms_obj <- fit_logistic_regression_bay(data_train, prior = brms::prior(normal(0, 3), class = "b"))
summary(brms_obj)
brms_perf <- assess_performance(brms_obj, data_train, data_test)


glm_obj <- fit_logistic_regression_freq(data_train)
summary(glm_obj)
glm_perf <- assess_performance(glm_obj, data_train, data_test)
beepr::beep(5)

perf_list <- list(rf = rf_perf, gbm = gbm_perf, brms = brms_perf, glm = glm_perf)

do.call(cbind, lapply(perf_list, \(x) x[["performance_train"]][["confusion_table"]]))
do.call(cbind, lapply(perf_list, \(x) x[["performance_test"]][["confusion_table"]]))
sapply(perf_list, \(x) x[["performance_train"]][["prediction_accuracy"]])
sapply(perf_list, \(x) x[["performance_test"]][["prediction_accuracy"]])

# what happened below?
# cont_vars <- sapply(data_wider, is.numeric)
# diffs <- aggregate(data_wider[, cont_vars], list(data_wider$violent_after), mean)[, -2]
# o <- order(abs(apply(diffs[, -1L], 2L, diff)), decreasing = TRUE)
# diffs[, c(1L, o + 1L)]
# apply(diffs[, c(1L, o + 1L)][, -1L], 2L, diff)
#
#
# cc <- coef(glm_obj)
# names(cc)[1] <- "Intercept"
# oo <- order(abs(cc), decreasing = TRUE)
# cc <- cc[oo]
#
# ff <- fixef(brms_obj)
# oo <- order(abs(ff[, "Estimate"]), decreasing = TRUE)
# ff <- ff[oo, "Estimate"]
#
# length(intersect(names(cc), names(ff))) == length(cc)
# length(intersect(names(cc), names(ff))) == length(ff)
#
# options(scipen = 100)
# print(cbind(cc, ff[names(cc)]), digits = 3)
# print(cbind(ff, cc[names(ff)]), digits = 3)
#
#
# cat_vars <- which(sapply(data_wider, \(x) !is.numeric(x) && (nlevels(x) <= 20)))
# cat_vars <- cat_vars[setdiff(names(cat_vars), c("violent_after", "time"))]
# sapply(data_wider[, ..cat_vars], nlevels)
#
# result <- setNames(vector("list", length(cat_vars)), names(cat_vars))
# for (i in seq_along(cat_vars)) {
#   tb <- table(data_wider[["violent_after"]], data_wider[[cat_vars[i]]])
#   result[[i]] <- tb / sum(tb)
# }
# result
#
# diffs <- aggregate(data_wider[, ..cat_vars], list(data_wider[, ..cat_vars]), \(x) {tb <- table(x); tb / sum(tb)}, na.rm = TRUE)[, -2]
# o <- order(abs(apply(diffs[, -1L], 2L, diff)), decreasing = TRUE)
# diffs[, c(1L, o + 1L)]
# apply(diffs[, c(1L, o + 1L)][, -1L], 2L, diff)
#
# (\(x) {
#   tb <- table(x)
#   tb / sum(tb)
# })(c(1, 2, 3))
#
# (function(x) {
#   tb <- table(x)
#   tb / sum(tb)
# })(c(1, 2, 3))
