#' @export
simulate_logistic_regression <- function(np, no_covariates, slopes = NULL, design_matrix = NULL, ...) {
  UseMethod("simulate_logistic_regression", np)
}

#' @export
simulate_logistic_regression.default <- function(np, no_covariates, slopes = NULL, design_matrix = NULL, ...) {

  assert_that(is.count(np), is.count(no_covariates))

  design_matrix <- design_matrix %||% cbind(1, center(matrix(rnorm(np * (no_covariates - 1L)), np, no_covariates - 1L)))
  slopes        <- slopes        %||% rnorm(no_covariates)

  assert_that(
    is.matrix(design_matrix),
    ncol(design_matrix) == no_covariates,
    nrow(design_matrix) == np,
    length(slopes) == no_covariates
  )

  true_probs <- plogis(design_matrix %*% slopes)
  y <- as.integer(runif(np) <= true_probs)

  res <- list(y = y, design_matrix = design_matrix, log_reg_intercept = slopes[1L], log_reg_slopes = slopes[-1L], true_probs = true_probs)
  class(res) <- "logistic_regression_data"

  return(res)

}

#' @export
is.logistic_regression_data <- function(x) {
  inherits(x, "logistic_regression_data")
}

#' @export
simulate_logistic_regression.ltm_data <- function(np, no_covariates, slopes = NULL, design_matrix = NULL, slopes_items = NULL, ...) {

  # ensure names make sense
  ltm_data <- np
  np <- ltm_data$np
  ni <- ltm_data$ni
  no_time_points <- ltm_data$no_time_points

  assert_that(is.count(no_covariates))

  design_matrix <- design_matrix %||% cbind(1, center(matrix(rnorm(np * (no_covariates - 1L)), np, no_covariates - 1L)))
  slopes        <- slopes        %||% rnorm(no_covariates)
  slopes_items  <- slopes_items  %||% rnorm(ni * no_time_points)

  assert_that(
    is.matrix(design_matrix),
    ncol(design_matrix)  == no_covariates,
    nrow(design_matrix)  == np,
    length(slopes)       == no_covariates,
    length(slopes_items) == ni * no_time_points
  )

  if (no_time_points == 1L)
    design_matrix <- cbind(design_matrix, ltm_data$parameters$lt)
  else
    design_matrix <- cbind(design_matrix,
                           ltm_data$parameters$lt - ltm_data$parameters$offset_lt,
                           ltm_data$parameters$lt + ltm_data$parameters$offset_lt)

  slopes <- c(slopes, slopes_items)

  log_reg <- simulate_logistic_regression(np, ncol(design_matrix), slopes, design_matrix)

  res <- list(
    log_reg     = log_reg,
    ltm_data    = ltm_data
  )
  class(res) <- "logistic_regression_ltm_data"
  return(res)

}

#' @export
is.logistic_regression_ltm_data <- function(x) {
  inherits(x, "logistic_regression_ltm_data")
}

#' @export
data_2_stan.logistic_regression_ltm_data <- function(dat, store_predictions = FALSE, missing_idx = integer(), ...) {

  assert_that(is.flag(store_predictions), inherits(dat, "logistic_regression_ltm_data"))

  ni             <- dat$ltm_data$ni
  no_time_points <- dat$ltm_data$no_time_points
  no_covariates  <- ncol(dat$log_reg$design_matrix) - 1L - ni * no_time_points

  data <- data_2_stan(dat$ltm_data, ...)

  data$store_predictions        <- store_predictions
  data$log_reg_outcomes         <- dat$log_reg$y
  # we start from 2 because we drop the intercept since Stan's bernoulli_logit_glm has a separate argument for the intercept
  # the final columns contain the true values for lt, which we also drop
  data$design_matrix_covariates <- dat$log_reg$design_matrix[, seq(2L, 1L + no_covariates), drop = FALSE]
  data$no_covariates            <- no_covariates
  data$fit_logistic             <- 1L

  data <- data_2_stan_incorporate_missing_log_reg(data, missing_idx)

  return(data)

}

data_2_stan_incorporate_missing_log_reg <- function(data, missing_idx) {

  if (length(missing_idx) == 0L) {
    data$np_log_reg     <-   data$np
    data$log_reg_np_idx <- 1:data$np
    return(data)
  }

  np <- data$np
  assertthat::assert_that(
    all(missing_idx > 0L & missing_idx <= np)
  )
  non_missing_idx <- setdiff(1:np, missing_idx)

  if (length(data$log_reg_outcomes) == np)
    data$log_reg_outcomes <- data$log_reg_outcomes[non_missing_idx]

  data$log_reg_np_idx <- non_missing_idx
  data$np_log_reg     <- length(non_missing_idx)
  return(data)
}
