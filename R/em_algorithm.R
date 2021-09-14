em_lt_log_lambda <- function(np, ni, nr, df, log_a, b, log_E, log_lambda, lt) {

  prior_shape <- 1
  prior_rate  <- 1

  shape <- np * nr / 2 + 1
  for (i in seq_len(ni)) {
    idx_i <- df$item == i


    for (p in seq_len(np)) {

      idx_pi <- which(idx_i & df$patient == p)

      mus <- c(0, (df$score[idx_pi] - b) / exp(log_a))
      sds <- c(1, exp(log_E + log_lambda[i] - log_a))
      prop <- product_n_univariate_gaussian_pdfs(mus, sds)
      lt[p, i] <- prop[1]
    }

    rate <- sum(((df$score[idx_i] - exp(log_a[df$rater[idx_i]]) * lt[df$patient[idx_i], i] - b[df$rater[idx_i]]) / exp(log_E[df$rater[idx_i]]))^2) / 2
    log_lambda[i] <- -log(shape / rate) / 2

  }
  return(list(lt = lt, log_lambda = log_lambda))
}

em_log_a_log_E_b <- function(nr, df, log_a, b, log_E, log_lambda, lt) {

  prior_shape <- 1
  prior_rate  <- 1

  shape <- np * ni / 2 + 1 + prior_shape
  for (r in seq_len(nr)) {
    idx_r <- which(df$rater == r)

    # b
    mus <- c(0, (df$score[idx_r] - exp(log_a[r]) * lt[cbind(df$patient[idx_r], df$item[idx_r])]))
    sds <- c(1, exp(log_E[r] + log_lambda[df$item[idx_r]]))

    prop <- product_n_univariate_gaussian_pdfs(mus, sds)
    b[r] <- prop[1]

    # log_a
    mus <- c(0, (b[r] - df$score[idx_r]) / lt[cbind(df$patient[idx_r], df$item[idx_r])])
    sds <- sds / c(1, abs(lt[cbind(df$patient[idx_r], df$item[idx_r])]))

    prop <- product_n_univariate_gaussian_pdfs(mus, sds)
    log_a[r] <- log(abs(prop[1]))

    rate <- sum(((df$score[idx_r] - exp(log_a[r]) * lt[cbind(df$patient[idx_r], df$item[idx_r])] - b[r]) / exp(log_lambda[df$item[idx_r]]))^2) / 2
    log_E[r] <- -log((shape + prior_shape) / (rate + prior_rate)) / 2
    if (abs(log_E[r]) > 20)
      browser()

  }

  return(list(log_a = log_a, b = b, log_E = log_E))

}

#'@export
results_to_tib <- function(..., iteration = NULL) {

  assertthat::assert_that(
    !any(is.na(...names())),
    is.null(iteration) ||assertthat::is.count(iteration)
  )

  dots <- list(...)
  tibble::tibble(
    parameter = rep(...names(), lengths(dots)),
    estimate  = unlist(dots, use.names = FALSE),
    iteration = if (!is.null(iteration)) iteration
  )
}

#'@export
fit_em <- function(df, n_iter = 5, record_all_iterations = 1L, progress = TRUE, compute_locations_scales = FALSE, init_method = c("data", "random"),
                   log_a = NULL, b = NULL, log_E = NULL, log_lambda = NULL, lt = NULL) {

  init_method <- match.arg(init_method)
  assertthat::assert_that(
    assertthat::is.count(n_iter),
    record_all_iterations == 0 || assertthat::is.count(record_all_iterations) || assertthat::is.flag(record_all_iterations),
    assertthat::is.flag(compute_locations_scales)
  )

  tmp <- get_np_ni_nr_from_df(df)
  np <- tmp[1L]; ni <- tmp[2L]; nr <- tmp[3L]

  init <- data_2_init(df, return_tibble = FALSE)
  log_a      <- log_a      %||% if (init_method == "data") init$log_a      else rnorm(nr)
  b          <- b          %||% if (init_method == "data") init$b          else rnorm(nr)
  log_E      <- log_E      %||% if (init_method == "data") init$log_E      else rnorm(nr)
  lt         <- lt         %||% if (init_method == "data") init$lt         else matrix(rnorm(np*ni), np, ni)
  log_lambda <- log_lambda %||% if (init_method == "data") init$log_lambda else rnorm(ni)

  validate_parameter_dimensions(np, ni, nr, log_a, b, log_E, lt, log_lambda)

  if (record_all_iterations)
    result <- tibble::tibble()

  prog <- if (progress)
    progress::progress_bar$new(total = n_iter, format = "[:bar] :percent eta: :eta")

  # browser()
  for (it in seq_len(n_iter)) {

    tmp_pi  <- em_lt_log_lambda(np, ni, nr, df, log_a, b, log_E, log_lambda, lt)

    lt         <- tmp_pi$lt
    log_lambda <- tmp_pi$log_lambda

    tmp_r <- em_log_a_log_E_b(nr, df, log_a, b, log_E, log_lambda, lt)

    # normalize the scale of log_a and lt, the mean of lt and b, and the sum of log_E and log_lambda
    b     <- tmp_r$b
    log_a <- tmp_r$log_a
    log_E <- tmp_r$log_E

    mu_lt <- mean(lt)
    b     <- b + mu_lt
    lt    <- lt - mean(lt)

    sd_lt <- sd(lt)
    log_a <- log_a + log(sd_lt)
    lt    <- lt / sd_lt

    mu_log_lambda <- mean(log_lambda)
    log_lambda    <- log_lambda - mean(log_lambda)
    log_E         <- log_E + mu_log_lambda

    # results_to_tib(log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda, lt = lt, iteration = it) |>
    #   add_true_values_to_em_tib(dat) |>
    #   scatterplot_retrieval(facets = iteration~parameter)

    if (it %% record_all_iterations == 0L) {
      toAdd <- if (compute_locations_scales) {
        locations <- exp(log_a[df$rater]) * c(lt[cbind(df$patient, df$item)]) + b[df$rater]
        scales    <- exp(log_E[df$rater] + log_lambda[df$item])
        results_to_tib(log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda, lt = lt, locations = locations, scales = scales, iteration = it)
      } else {
        results_to_tib(log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda, lt = lt, iteration = it)
      }
      result <- rbind(result,  toAdd)
    }

    if (progress)
      prog$tick()

  }

  if (record_all_iterations)
    return(result)

  if (compute_locations_scales) {
    locations <- exp(log_a[df$rater]) * c(lt[cbind(df$patient, df$item)]) + b[df$rater]
    scales    <- exp(log_E[df$rater] + log_lambda[df$item])
    return(results_to_tib(log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda, lt = lt, locations = locations, scales = scales))
  } else {
    return(results_to_tib(log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda, lt = lt))
  }

}

#'@export
add_true_values_to_em_tib <- function(em_tib, dat) {
  param_nms <- unique(em_tib$parameter)
  nrep <- if ("iteration" %in% names(em_tib)) vctrs::vec_unique_count(em_tib$iteration) else 1L
  em_tib$true_value <- rep(unlist(dat$parameters[param_nms], use.names = FALSE), nrep)

  if (dat$no_rater_groups != 1L)
    em_tib$rater_group <- rep(factor(unlist(Map(function(nm, len) {
    if (nm %in% c("log_a", "b", "log_E")) {
      dat$rater_group_assignment
    } else {
      rep(NA_integer_, len)
    }
  }, param_nms, lengths(dat$parameters[param_nms])), use.names = FALSE)), nrep)

  em_tib
}


get_np_ni_nr_from_df <- function(df) {
  np <- vctrs::vec_unique_count(df$patient)
  ni <- vctrs::vec_unique_count(df$item)
  nr <- vctrs::vec_unique_count(df$rater)
  return(c(np, ni, nr))
}

#'@export
data_2_init <- function(df, return_tibble = TRUE) {

  tmp <- get_np_ni_nr_from_df(df)
  np <- tmp[1]; ni <- tmp[2]; nr <- tmp[3]

  np <- vctrs::vec_unique_count(df$patient)
  ni <- vctrs::vec_unique_count(df$item)
  nr <- vctrs::vec_unique_count(df$rater)

  b          <-        unlist(tapply(df$score, df$rater,                  mean, simplify = FALSE), use.names = FALSE)
  log_a      <-    log(unlist(tapply(df$score, df$rater,                  sd,   simplify = FALSE), use.names = FALSE))
  log_E      <-    log(unlist(tapply(df$score, df$rater,                  sd,   simplify = FALSE), use.names = FALSE))
  log_lambda <-    log(unlist(tapply(df$score, df$item,                   sd,   simplify = FALSE), use.names = FALSE))
  lt         <- matrix(unlist(tapply(df$score, list(df$patient, df$item), mean, simplify = FALSE), use.names = FALSE), np, ni)

  mu_lt <- mean(lt)
  b     <- b# + mu_lt
  lt    <- lt - mean(lt)

  sd_lt <- sd(lt)
  log_a <- log_a + log(sd_lt)
  lt    <- lt / sd_lt

  mu_log_lambda <- mean(log_lambda)
  log_lambda    <- log_lambda - mean(log_lambda)
  log_E         <- log_E + mu_log_lambda

  log_a      <- center(log_a)
  log_E      <- center(log_E)
  log_lambda <- center(log_lambda)

  if (return_tibble)
    results_to_tib(lt = lt, log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda)
  else
    list(lt = lt, log_a = log_a, b = b, log_E = log_E, log_lambda = log_lambda)
}

compute_hyperparams <- function(x, dat) {
  x$mu_lt         <- mean(x$lt)
  x$sd_lt         <- sd(x$lt)

  x$mu_b          <- tapply(x$b, dat$rater_group_assignment, mean)
  x$sd_b          <- tapply(x$b, dat$rater_group_assignment, sd)

  x$mu_log_E      <- tapply(x$log_E, dat$rater_group_assignment, mean)
  x$sd_log_E      <- tapply(x$log_E, dat$rater_group_assignment, sd)

  x$sd_log_lambda <- sd(x$log_lambda)
  x$sd_log_a      <- tapply(x$log_a, dat$rater_group_assignment, sd)

  return(x)
}

compute_raw_params <- function(x) {
  x$log_lambda_raw <- CRMLogistic:::sum_to_zero_QR_inv_R(x$log_lambda)
  x$log_lambda     <- CRMLogistic:::center(x$log_lambda)

  x$log_a_raw <- CRMLogistic:::sum_to_zero_QR_inv_R(x$log_a)
  x$log_a <- CRMLogistic:::center(x$log_a)
  x
}

#' @export
em_tib_to_init <- function(tib, dat) {
  em_tib |>
    subset(parameter %in% c("log_a", "b", "log_E", "log_lambda", "lt") & iteration == max(iteration)) |>
    split(~parameter) |>
    lapply(`[[`, "estimate") |>
    compute_hyperparams(dat) |>
    compute_raw_params()
}
