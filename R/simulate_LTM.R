#' @export
simulate_data_ltm <- function(np, ni, nr, nc, no_rater_groups = 4L,
                          lt = NULL, log_lambda = NULL, log_a = NULL, b = NULL, log_E = NULL,
                          mu_log_a = NULL, sd_log_a = NULL, mu_b = NULL, sd_b = NULL, mu_log_E = NULL, sd_log_E = NULL,
                          return_mus_sds = TRUE, vary_lambda_across_patients = NULL, store_probabilities = FALSE) {

  assert_counts(np, ni, nr, no_rater_groups)
  assert_that(is.count(nc), is.flag(store_probabilities))


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
  log_lambda <- log_lambda %||% center(matrix(log(rlnorm(ni*np_log_lambda, 0, .2)), ni, np_log_lambda))
  log_a      <- log_a      %||% center(log(rlnorm(nr, mu_log_a[rater_group_assignment], sd_log_a[rater_group_assignment])))
  b          <- b          %||%             rnorm(nr, mu_b    [rater_group_assignment], sd_b    [rater_group_assignment])
  log_E      <- log_E      %||% center(log(rlnorm(nr, mu_log_E[rater_group_assignment], sd_log_E[rater_group_assignment])))

  validate_parameter_dimensions(np, ni, nr, log_a, b, log_E, lt, log_lambda)

  df <- tibble::tibble(
    patient     = rep(seq_len(np),                times = ni * nr),
    item        = rep(seq_len(ni), each = np,     times = nr     ),
    rater       = rep(seq_len(nr), each = np * ni                ),
    rater_group = rater_group_assignment[rater]
  )

  assertthat::assert_that(nrow(df) == np * ni * nr)

  gamma <- (1:(nc - 1))
  gamma <- log(gamma / (nc - gamma))

  score <- integer(nrow(df))
  probs <- if (store_probabilities) matrix(NA_real_, nrow(df), nc) else NULL

  for (o in seq_len(nrow(df))) {

    p <- df$patient[o]
    i <- df$item[o]
    r <- df$rater[o]

    delta <- exp(log_a[r]) * gamma + b[r]

    location <- lt[p, i]
    scale    <- exp(log_lambda[i] - log_E[r])

    prob <- get_probs_ordered(delta, location, scale)
    score[o] <- sample(nc, 1L, prob = prob)

    if (store_probabilities)
      probs[o, ] <- prob

  }

  df$score <- score

  res <- list(
    df = df,
    np = np, ni = ni, nr = nr, nc = nc, no_rater_groups = no_rater_groups,
    rater_group_assignment = rater_group_assignment,
    parameters = list(
      lt = lt, log_lambda = log_lambda,
      log_a = log_a, b = b, log_E = log_E,
      probabilities = probs
    )
  )
  class(res) <- "ltm_data"
  return(res)
}

guess_nc <- function(df) {
  nc <- if (is.factor(df$score)) {

  } else {
    length(unique())
  }
  message("Guessing nc at ", nc)
  nc
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
ltm_data_2_stan <- function(dat, nc, prior_only = FALSE, debug = TRUE, vary_lambda_across_patients = FALSE, vectorized = FALSE,
                       mu_log_lambda = 0, mu_log_a = 0, mu_b = 0,
                       a_sd_lt         = 1,   b_sd_lt         = 1,
                       a_sd_log_lambda = 1.1, b_sd_log_lambda = 1.1,
                       a_sd_log_E      = 1,   b_sd_log_E      = 1,
                       a_sd_log_a      = 1,   b_sd_log_a      = 1,
                       a_sd_b          = 1,   b_sd_b          = 1
                       ) {


  if (!is.data.frame(dat)) {

    x_df <- dat$df
    np   <- dat$np
    ni   <- dat$ni
    nr   <- dat$nr
    nc   <- dat$nc
    no_rater_groups <- dat$no_rater_groups
    rater_group_assignment <- dat$rater_group_assignment

  } else {
    x_df <- dat
    np <- max(as.integer(x_df$patient))
    ni <- max(as.integer(x_df$item))
    nr <- max(as.integer(x_df$rater))
    no_rater_groups <- if (!is.null(x_df$rater_group)) max(x_df$rater_group) else 1L
    stop("define rater_group_assignment!")
  }

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

  assertthat::assert_that(
    assertthat::is.number(mu_log_lambda),
    assertthat::is.number(mu_log_a),
    assertthat::is.number(mu_b),
    all(is.numeric(positive_reals) & positive_reals > 0),
    assertthat::is.number(b_sd_log_lambda) && b_sd_log_lambda > 0,
    assertthat::is.flag(prior_only),
    assertthat::is.flag(vary_lambda_across_patients),
    assertthat::is.flag(vectorized)
  )

  data <- list(
    x                   = x_df$score,
    np                  = np,
    ni                  = ni,
    nr                  = nr,
    nc                  = nc,
    no_rater_groups     = no_rater_groups,

    n_observed          = length(x_df$score),
    idx_patient         = as.integer(x_df$patient),
    idx_item            = as.integer(x_df$item),
    idx_rater           = as.integer(x_df$rater),
    idx_rater_group     = as.integer(rater_group_assignment),

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


    debug                       = as.integer(debug),
    prior_only                  = as.integer(prior_only),
    vary_lambda_across_patients = as.integer(vary_lambda_across_patients),
    vectorized                  = as.integer(vectorized)
  )
  return(data)
}

get_probs_ordered <- function(delta, location, scale, cdf = plogis) {

  # replica of what Stan does
  nc <- length(delta) + 1L
  prob <- numeric(nc)

  vals <- cdf(location - delta, 0, scale)
  prob[1] <- 1 - vals[1]
  for (c in 2:(nc - 1)) {
    prob[c] <- vals[c - 1L] - vals[c]
  }
  prob[nc] <- vals[nc - 1]

  return(prob)
}

#' @export
rordered <- function(n, delta, location = 0, scale = 1, cdf = stats::plogis) {
  assert_that(is.integer(n) && n >= 0L, length(delta) >= 1L, scale > 0L, is.function(cdf))
  if (n == 0L) return(integer())

  prob <- get_probs_ordered(delta, location, scale, cdf = cdf)
  sample(length(delta) + 1L, size = n, replace = TRUE, prob = prob)
}

#' @export
dordered <- function(x, delta, location = 0, scale = 1, log = TRUE, cdf = stats::plogis) {
  prob <- get_probs_ordered(delta, location, scale, cdf = cdf)
  if (log)
    return(log(prob[x]))
  else return(prob[x])
}
