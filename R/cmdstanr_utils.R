#' @title compile a stan model
#' @inheritParams cmdstanr::cmdstan_model
#' @param dir directory for compiled executables
#' @details This is a helper for this project so you don't have to adjust any directories.
#' It's not recommended for use outside of this project, just use cmdstanr::cmdstan_model directly.
#' @export
compile_stan_model <- function(stan_file, compile = TRUE, dir = file.path("stanmodels", "compiled_stan_models"), ...) {
  cmdstanr::cmdstan_model(stan_file = stan_file, compile = compile, dir = dir, ...)
}

#' @title run a stan model
#'
#' @param expr some R expression
#' @param path path to saved results
#' @param force if TRUE, ignore any saved results the
#' @details if path exists and force is FALSE it loads the results with `readRDS`. Otherwise, evaluate `expr`.
#' @export
save_or_run_model <- function(expr, path, force) {

  if (file.exists(path) && !force) {
    return(readRDS(path))
  } else {
    fit <- force(expr)
    if (inherits(fit, "CmdStanFit")) {
      fit$save_object(path)
    } else {
      saveRDS(fit, file = path)
    }
    return(fit)
  }
}

get_params <- function(fit, filter = "lp__|lp_approx__|_raw|sd_|mu_|log_probs") {
  params <- fit$metadata()$stan_variables
  params[!grepl(filter, params)]
}

#'@export
get_means_tib <- function(fit) {
  params <- get_params(fit)
  draws <- fit$draws(params, format = "draws_matrix")

  meta <- fit$metadata()
  tibble(
    parameter = unlist(lapply(params, \(x) rep(x, prod(meta$stan_variable_sizes[[x]]))), use.names = FALSE),
    estimate  = unname(colMeans(draws))
  )
}

#'@export
add_true_values_to_tib <- function(tib, dat, filter_unused = TRUE, keep_free_thresholds = TRUE) {

  if (filter_unused)
    tib <- tib |> dplyr::filter(parameter %in% c(get_params_used(dat), if (keep_free_thresholds) "free_thresholds"))

  param_nms <- unique(tib$parameter)

  tib$true_value <- get_true_values(dat, param_nms)
  tib

}

get_params_used <- function(dat) {
  UseMethod("get_params_used", dat)
}

get_params_used.default <- function(dat) {
  dat$parameters_used
}

get_params_used.logistic_regression_ltm_data <- function(dat) {
  return(c(get_params_used(dat$log_reg), get_params_used(dat$ltm_data)))
}

get_params_used.logistic_regression_data <- function(dat) {
  nms <- c("log_reg_intercept", "log_reg_slopes")
  assert_that(all(nms %in% names(dat)))
  return(nms)
}

get_true_values <- function(dat, param_nms) {
  UseMethod("get_true_values", dat)
}
get_true_values.default <- function(dat, param_nms) {
  unlist(dat$parameters[param_nms], use.names = FALSE)
}

get_true_values.logistic_regression_ltm_data <- function(dat, param_nms) {
  unlist(c(dat$log_reg, dat$ltm_data$parameters)[param_nms], use.names = FALSE)
}

#
#   true_parameters <- dat$parameters[params]
#   tib <- tibble(
#     estimate       = tib$value,
#     true_value     = unlist(true_parameters, use.names = FALSE),
#     parameter      = rep(names(true_parameters), lengths(true_parameters)),
#     bias           = true_value - estimate,
#     param_index    = unlist(lapply(lengths(true_parameters), seq_len), use.names = FALSE),
#   )
#
#   tib$rater_group <- factor(unlist(Map(function(nm, len) {
#     if (nm %in% c("log_a", "b", "log_E")) {
#       dat$rater_group_assignment
#     } else {
#       rep(NA_integer_, len)
#     }
#   }, names(true_parameters), lengths(true_parameters)), use.names = FALSE))
#   return(tib)
# }

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
