get_params <- function(fit, filter = "lp__|lp_approx__|_raw|sd_|mu_") {
  params <- fit$metadata()$stan_variables
  params[!grepl(filter, params)]
}

#'@export
get_means_tib <- function(fit, dat) {
  params <- get_params(fit)
  draws <- fit$draws(params, format = "draws_matrix")

  tib <- tibble(
    variable = colnames(draws),
    value    = colMeans(draws)
  )

  true_parameters <- dat$parameters[params]
  tib <- tibble(
    estimate       = tib$value,
    true_value     = unlist(true_parameters, use.names = FALSE),
    parameter      = rep(names(true_parameters), lengths(true_parameters)),
    bias           = true_value - estimate,
    param_index    = unlist(lapply(lengths(true_parameters), seq_len), use.names = FALSE),
  )

  tib$rater_group <- factor(unlist(Map(function(nm, len) {
    if (nm %in% c("log_a", "b", "log_E")) {
      dat$rater_group_assignment
    } else {
      rep(NA_integer_, len)
    }
  }, names(true_parameters), lengths(true_parameters)), use.names = FALSE))
  return(tib)
}

#'@export
benchmark_table <- function(...) {
  nms <- ...names()
  dots <- list(...)
  tib <- tibble(
    method   = ...names(),
    time     = bench::as_bench_time(unlist(lapply(dots, \(x) x$time()$total))),
    relative = c(time / min(time))
  )
  return(tib)
}
