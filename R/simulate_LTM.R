#' @export
simulate_data_ltm <- function(np, ni, nr, nc, no_rater_groups = 4L, no_time_points = 1L,
                          lt = NULL, log_lambda = NULL, log_a = NULL, b = NULL, log_E = NULL, threshold_shape = NULL,
                          mu_log_a = NULL, sd_log_a = NULL, mu_b = NULL, sd_b = NULL, mu_log_E = NULL, sd_log_E = NULL,
                          return_mus_sds = TRUE, vary_lambda_across_patients = NULL, store_probabilities = FALSE,
                          threshold_type = c("logistic", "skew_logistic", "free"), thresholds = NULL) {

  threshold_type <- match.arg(threshold_type)
  use_skew_thresholds <- threshold_type == "skew"
  use_free_thresholds <- threshold_type == "free"

  assert_counts(np, ni, nr, no_rater_groups)
  assert_that(
    is.flag(store_probabilities),
    is.flag(use_skew_thresholds),
    is.flag(use_free_thresholds),
    !(use_skew_thresholds && use_free_thresholds),
    is.count(no_time_points) && no_time_points >= 1L && no_time_points <= 2L,
    is.null(thresholds) || (is.matrix(thresholds) && is.numeric(thresholds) && nrow(thresholds) == nr && ncol(thresholds) == nc - 1)
  )


  vary_lambda_across_patients <- vary_lambda_across_patients %||% (if (is.null(log_lambda)) FALSE else NCOL(log_lambda) == np)

  if (no_rater_groups == 1L) {
    mu_log_a <- 0
    mu_b     <- 0
    mu_log_E <- 0
    sd_log_a <- 1
    sd_b     <- 1
    sd_log_E <- 1
    rater_group_assignment <- rep(1L, nr)
  } else {
    mu_log_a <- mu_log_a %||% center(rnorm(no_rater_groups))
    mu_b     <- mu_b     %||% center(rnorm(no_rater_groups))
    mu_log_E <- mu_log_E %||% center(rnorm(no_rater_groups))
    sd_log_a <- sd_log_a %||% rgamma(no_rater_groups, 5, 20)
    sd_b     <- sd_b     %||% rgamma(no_rater_groups, 5, 10)
    sd_log_E <- sd_log_E %||% rgamma(no_rater_groups, 5, 20)
    rater_group_assignment <- sort(rep(seq_len(no_rater_groups), length.out = nr))
  }

  assert_that(
    length(mu_log_a) == no_rater_groups,
    length(mu_b)     == no_rater_groups,
    length(mu_log_E) == no_rater_groups,
    length(sd_log_a) == no_rater_groups && all(sd_log_a > 0),
    length(sd_b)     == no_rater_groups && all(sd_b     > 0),
    length(sd_log_E) == no_rater_groups && all(sd_log_E > 0)
  )

  np_log_lambda <- if (vary_lambda_across_patients) np else 1L
  lt         <- lt         %||% center(matrix(rnorm(np * ni), np, ni))
  log_lambda <- log_lambda %||% center(matrix(log(rlnorm(ni*np_log_lambda, 0, .8)), ni, np_log_lambda))
  log_a      <- log_a      %||% center(log(rlnorm(nr, mu_log_a[rater_group_assignment], sd_log_a[rater_group_assignment])))
  b          <- b          %||%             rnorm(nr, mu_b    [rater_group_assignment], sd_b    [rater_group_assignment])
  log_E      <- log_E      %||% center(log(rlnorm(nr, mu_log_E[rater_group_assignment], sd_log_E[rater_group_assignment])))

  threshold_shape <- if (!use_skew_thresholds) numeric(nr) else {
    threshold_shape %||% rnorm(nr, 1)
  }

  offset_lt <- if (no_time_points == 2L) matrix(rnorm(np*ni), np, ni) else matrix(NA_real_, 0, 0)

  validate_parameter_dimensions(np, ni, nr, no_time_points, log_a, b, log_E, lt, log_lambda, offset_lt, vary_lambda_across_patients)

  df <- tibble::tibble(
    patient     = rep(seq_len(np),                times = ni * nr),
    item        = rep(seq_len(ni), each = np,     times = nr     ),
    rater       = rep(seq_len(nr), each = np * ni                ),
    rater_group = rater_group_assignment[rater]
  )

  if (no_time_points == 2L) {
    df <- rbind(df, df)
    df[["time_point"]] <- rep(1:2, each = np * ni * nr)
  }

  assertthat::assert_that(nrow(df) == np * ni * nr * no_time_points)

  threshold_probs <- seq_len(nc-1L) / nc
  gamma <- (1:(nc - 1))
  gamma <- log(gamma / (nc - gamma))

  score <- integer(nrow(df))
  probs <- if (store_probabilities) matrix(NA_real_, nrow(df), nc) else NULL

  # free_thresholds <- NULL
  if (!is.null(thresholds)) {
    free_thresholds <- thresholds
  } else {
    free_thresholds <- matrix(NA, nr, nc - 1)
    if (use_free_thresholds) {
      for (r in seq_len(nr))
        free_thresholds[r, ] <- sort(rlogis(nc - 1))
    } else {
      for (r in seq_len(nr))
        free_thresholds[r, ] <- qELGW(threshold_probs, location = b[r], scale = exp(log_a[r]), shape = threshold_shape[r])
    }
  }

  for (o in seq_len(nrow(df))) {

    p <- df$patient[o]
    i <- df$item[o]
    r <- df$rater[o]

    delta <- if (use_free_thresholds) {
      free_thresholds[r, ]
    } else {
      qELGW(threshold_probs, location = b[r], scale = exp(log_a[r]), shape = threshold_shape[r])
    }

    location <- lt[p, i]
    scale    <- exp(log_lambda[i] - log_E[r])

    if (no_time_points == 2L) {
      t <- df$time_point[o]
      # (2 * t - 3) maps 1:2 to c(-1, 1)
      location <- location + (2 * t - 3) * offset_lt[p, i]
    }

    prob <- get_probs_ordered(delta, location, scale)
    score[o] <- sample(nc, 1L, prob = prob)

    if (store_probabilities)
      probs[o, ] <- prob

  }

  df$score <- score

  parameters_used <- if (use_free_thresholds) {
    c("lt", "log_lambda", "log_E", "free_thresholds")
  } else if (use_skew_thresholds) {
    c("lt", "log_lambda", "log_E", "log_a", "b", "threshold_shape")
  } else {
    c("lt", "log_lambda", "log_E", "log_a", "b")
  }

  if (no_time_points == 2L)
    parameters_used <- c(parameters_used, "offset_lt")

  res <- list(
    df = df,
    np = np, ni = ni, nr = nr, nc = nc, no_rater_groups = no_rater_groups, no_time_points = no_time_points,
    rater_group_assignment = rater_group_assignment, threshold_type = threshold_type,
    parameters = list(
      lt = lt, log_lambda = log_lambda,
      log_a = log_a, b = b, log_E = log_E,
      threshold_shape = threshold_shape,
      free_thresholds = free_thresholds,
      probabilities = probs,
      offset_lt = offset_lt
    ),
    parameters_used = parameters_used
  )
  class(res) <- "ltm_data"
  return(res)
}

#' @export
is.ltm_data <- function(x) inherits(x, "ltm_data")

guess_nc <- function(dat, score = "score") {
  nc <- if (is.factor(dat[[score]])) {
    nlevels(dat[[score]])
  } else {
    vctrs::vec_unique_count(dat[[score]])
  }
  message("Guessing nc at ", nc)
  nc
}

#' @export
add_implied_thresholds <- function(tib, nc, nr) {

  if ("free_thresholds" %in% tib$parameter)
    return(tib)

  tib_to_add <- tibble(
    parameter = rep("free_thresholds", nr * (nc - 1)),
    estimate  = numeric(nr * (nc - 1))
  )
  default_thresholds <- log(seq_len(nc-1L) / (nc - seq_len(nc-1L)))
  shifts <- tib |> dplyr::filter(parameter == "b")     |> dplyr::select(estimate) |> unlist(use.names = FALSE)
  scales <- tib |> dplyr::filter(parameter == "log_a") |> dplyr::select(estimate) |> unlist(use.names = FALSE)

  for (r in seq_len(nr)) {
    thresholds <- shifts[r] + scales[r] * default_thresholds
    # tib_to_add[["estimate"]][seq(1 + (nc - 1) * (r - 1), (nc - 1) * r)] <- thresholds
    tib_to_add[["estimate"]][seq(r, (nc-1) * nr, nr)] <- thresholds

  }

  return(rbind(tib, tib_to_add))
}


#' @export
log_likelihood_ltm <- function(x, nc = NULL, log_a, b, log_E, lt, log_lambda,...) {
  UseMethod("log_likelihood_ltm", x)
}

#' @export
log_likelihood_ltm.ltm_data <- function(x, nc = NULL, log_a, b, log_E, lt, log_lambda) {

  if (!is.null(nc))
    message("nc is ignored when passing an object of class 'ltm_data'")

  log_likelihood_ltm.data.frame(x$df, x$nc, log_a, b, log_E, lt, log_lambda)
}

#' @export
log_likelihood_ltm.data.frame <- function(x, nc = NULL, log_a, b, log_E, lt, log_lambda) {

  gamma <- 1:(nc - 1)
  gamma <- log(gamma / (nc - gamma))

  nc <- nc %||% guess_nc(x)

  log_likelihood <- 0.0

  for (o in seq_len(nrow(x))) {

    p <- x$patient[o]
    i <- x$item[o]
    r <- x$rater[o]

    delta <- exp(log_a[r]) * gamma + b[r]

    location <- lt[p, i]
    scale    <- exp(log_lambda[i] - log_E[r])

    log_prob <- get_probs_ordered(delta, location, scale)[x$score[o]]
    log_likelihood <- log_likelihood + log_prob

  }
  log_likelihood
}

#' @export
data_2_stan <- function(dat, ...) {
  UseMethod("data_2_stan", dat)
}

#' @export
data_2_stan.ltm_data <- function(dat, nc = NULL, prior_only = FALSE, debug = TRUE, vary_lambda_across_patients = FALSE,
                           use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE,
                           mu_log_lambda = 0, mu_log_a = 0, mu_b = 0,
                           a_sd_lt         = 1,   b_sd_lt         = 1,
                           a_sd_log_lambda = 1.1, b_sd_log_lambda = 1.1,
                           a_sd_log_E      = 1,   b_sd_log_E      = 1,
                           a_sd_log_a      = 1,   b_sd_log_a      = 1,
                           a_sd_b          = 1,   b_sd_b          = 1,
                           ...
                           ) {

  x_df <- dat$df
  np   <- dat$np
  ni   <- dat$ni
  nr   <- dat$nr
  nc   <- dat$nc
  no_rater_groups <- dat$no_rater_groups
  rater_group_assignment <- dat$rater_group_assignment
  no_time_points <- dat$no_time_points
  idx_time_point <- as.integer(x_df$time_point)

  data <- stan_data_ltm_inner(
    np, ni, nr, no_rater_groups, nc,
    a_sd_lt, b_sd_lt, a_sd_log_lambda, b_sd_log_lambda, a_sd_log_E, b_sd_log_E, a_sd_log_a, b_sd_log_a,
    a_sd_b, b_sd_b, mu_log_lambda, mu_log_a, mu_b,
    prior_only, vary_lambda_across_patients,  use_skew_logistic_thresholds, use_free_logistic_thresholds,
    rater_group_assignment, x_df, no_time_points, idx_time_point, debug
  )

  return(data)
}

stan_data_ltm_inner <- function(np, ni, nr, no_rater_groups, nc, a_sd_lt, b_sd_lt, a_sd_log_lambda, b_sd_log_lambda, a_sd_log_E, b_sd_log_E, a_sd_log_a, b_sd_log_a, a_sd_b, b_sd_b, mu_log_lambda, mu_log_a, mu_b, prior_only, vary_lambda_across_patients, use_skew_logistic_thresholds, use_free_logistic_thresholds, rater_group_assignment, x_df, no_time_points, idx_time_point, debug) {
  assert_counts(np, ni, nr, no_rater_groups)
  assert_that(is.positive_int(nc, larger_than = 2L))

  positive_reals <- c(
    a_sd_lt             = a_sd_lt,
    b_sd_lt             = b_sd_lt,
    a_sd_log_lambda     = a_sd_log_lambda,
    b_sd_log_lambda     = b_sd_log_lambda,
    a_sd_log_E          = a_sd_log_E,
    b_sd_log_E          = b_sd_log_E,
    a_sd_log_a          = a_sd_log_a,
    b_sd_log_a          = b_sd_log_a,
    a_sd_b              = a_sd_b,
    b_sd_b              = b_sd_b
  )

  assert_that(
    is.number(mu_log_lambda),
    is.number(mu_log_a),
    is.number(mu_b),
    all(is.numeric(positive_reals) & positive_reals > 0),
    is.number(b_sd_log_lambda) && b_sd_log_lambda > 0,
    is.flag(prior_only),
    is.flag(vary_lambda_across_patients),
    is.flag(use_skew_logistic_thresholds),
    is.flag(use_free_logistic_thresholds),
    !is.unsorted(rater_group_assignment),
    !(use_skew_logistic_thresholds && use_free_logistic_thresholds)
  )

  return(list(
    x                   = x_df$score,
    np                  = np,
    ni                  = ni,
    nr                  = nr,
    nc                  = nc,
    no_time_points      = no_time_points,
    no_rater_groups     = no_rater_groups,

    n_observed          = length(x_df$score),
    idx_patient         = as.integer(x_df$patient),
    idx_item            = as.integer(x_df$item),
    idx_rater           = as.integer(x_df$rater),
    idx_rater_group     = as.integer(rater_group_assignment),
    idx_time_point      = idx_time_point,
    rater_group_counts  = as.integer(unname(c(table(rater_group_assignment)))),

    n_missing           = 0L,
    idx_patient_missing = integer(),
    idx_item_missing    = integer(),
    idx_rater_missing   = integer(),
    idx_longer_lt       = rep(1:(np * ni), nr),

    mu_log_lambda       = mu_log_lambda,
    mu_log_a            = rep(0, no_rater_groups),#mu_log_a,
    mu_b                = mu_b,

    a_sd_lt             = a_sd_lt,
    b_sd_lt             = b_sd_lt,
    a_sd_log_lambda     = a_sd_log_lambda,
    b_sd_log_lambda     = b_sd_log_lambda,
    a_sd_log_E          = a_sd_log_E,
    b_sd_log_E          = b_sd_log_E,
    a_sd_log_a          = a_sd_log_a,
    b_sd_log_a          = b_sd_log_a,
    a_sd_b              = a_sd_b,
    b_sd_b              = b_sd_b,


    debug                        = as.integer(debug),
    prior_only                   = as.integer(prior_only),
    vary_lambda_across_patients  = as.integer(vary_lambda_across_patients),
    # vectorized                   = as.integer(vectorized),
    use_skew_logistic_thresholds = as.integer(use_skew_logistic_thresholds),
    use_free_logistic_thresholds = as.integer(use_free_logistic_thresholds),

    # these may be overwritten by data_2_stan.logistic_regression_ltm_data
    fit_logistic             = 0L,
    store_predictions        = 0L,
    log_reg_outcomes         = integer(),
    no_covariates            = 0L,
    design_matrix_covariates = matrix(NA_real_, 0L, 0L),
    np_log_reg               = 0L,
    log_reg_np_idx           = integer()

  ))
}


#TODO: this function need individual documentation! maybe it shouldn't even be an S3 method?
#' @export
data_2_stan.data.frame <- function(dat, score = "score", patient = "patient", item = "item",
                               rater = "rater", time = "time", rater_group = "rater_group",
                               nc = NULL,
                               logistic_dat = NULL, logistic_target = NULL, logistic_patient = "patient",
                               missing_idx = integer(),
                               store_predictions = FALSE,
                               prior_only = FALSE, debug = TRUE, vary_lambda_across_patients = FALSE,
                               use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = TRUE,
                               mu_log_lambda = 0, mu_log_a = 0, mu_b = 0,
                               a_sd_lt         = 1,   b_sd_lt         = 1,
                               a_sd_log_lambda = 1.1, b_sd_log_lambda = 1.1,
                               a_sd_log_E      = 1,   b_sd_log_E      = 1,
                               a_sd_log_a      = 1,   b_sd_log_a      = 1,
                               a_sd_b          = 1,   b_sd_b          = 1,
                               ...) {

  # TODO: add safety checks!

  # this is not very clean but it gets the job done
  dat <- dat |>
    mutate(
      across(all_of(c(patient, item, rater, rater_group, time)),
             normalize_factor
      ),
    )

  nc <- nc %||% guess_nc(dat)

  np <- max(as.integer(dat[[patient]]))
  ni <- max(as.integer(dat[[item]]))
  nr <- max(as.integer(dat[[rater]]))
  no_rater_groups <- max_or_1(dat[[rater_group]])
  rater_group_assignment <- dat[[rater_group]][!duplicated(dat[[rater]])]

  idx_time_point <- as.integer(dat[[time]])
  no_time_points <- max_or_1(idx_time_point)

  data <- stan_data_ltm_inner(
    np, ni, nr, no_rater_groups, nc,
    a_sd_lt, b_sd_lt, a_sd_log_lambda, b_sd_log_lambda, a_sd_log_E, b_sd_log_E, a_sd_log_a, b_sd_log_a,
    a_sd_b, b_sd_b, mu_log_lambda, mu_log_a, mu_b,
    prior_only, vary_lambda_across_patients,  use_skew_logistic_thresholds, use_free_logistic_thresholds,
    rater_group_assignment, dat, no_time_points, idx_time_point, debug
  )

  if (!is.null(logistic_dat)) {

    # TODO: should this be a separate function so it can be done without ltm? the model doesnt support that though
    f <- paste0(logistic_target, "~ 1")
    if (ncol(logistic_dat) > 2L)
      f <- paste(f, "+", paste(setdiff(colnames(logistic_dat), c(logistic_target, logistic_patient)), collapse = " + "))
    f <- as.formula(f)
    design_matrix <- model.matrix(f, logistic_dat)

    no_covariates <- ncol(design_matrix) - 1L

    data$store_predictions        <- store_predictions
    data$log_reg_outcomes         <- normalize_factor(logistic_dat[[logistic_target]]) - 1L
    # we start from 2 because we drop the intercept since Stan's bernoulli_logit_glm has a separate argument for the intercept
    data$design_matrix_covariates <- design_matrix[, -1L, drop = FALSE]
    data$no_covariates            <- no_covariates
    data$fit_logistic             <- 1L

    data <- data_2_stan_incorporate_missing_log_reg(data, missing_idx)

  }

  return(data)

}

max_or_1 <- function(x) {
  if (is.null(x))
    1L
  else
    vctrs::vec_unique_count(x)
  # if (is.factor(x))
  #   nlevels(x)
  # else
  #   max(x)
}
