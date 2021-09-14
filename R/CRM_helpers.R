assert_counts <- function (np, ni, nr, no_rater_groups) {
  assertthat::assert_that(
    assertthat::is.count(np),
    assertthat::is.count(ni),
    assertthat::is.count(nr),
    assertthat::is.count(no_rater_groups),
    nr > no_rater_groups
  )
}


get_means <- function(samples, paramName) {
  paramNames <- dimnames(samples)$parameters
  idx <- grep(sprintf("^%s\\[", paramName), paramNames)
  if (length(idx) == 0L) {
    idx <- which(paramNames == paramName)
    if (length(idx) != 1L)
      stop("did not understand paramName ", paramName)
  }
  means <- apply(samples[, , idx, drop = FALSE], 3, mean)
  return(means)
}

trimBrackets <- function(x) {
  gsub("(\\[.*)", "",x)
}


stack_chain <- function(samps) {
  # TODO: optimize this ugly code
  do.call(rbind, lapply(seq_len(ncol(samps)), \(x) samps[, x, ]))
}

center <- function(x) {
  y <- scale(x, scale = FALSE)
  attr(y, "scaled:center") <- NULL
  if (is.vector(x))
    setNames(c(y), names(x))
  else
    y
}

prod_to_1 <- function(x) exp(center(log(x)))#x / prod(x)

cor_2_tib <- function(m) {
  ut <- upper.tri(m)
  tib <- tibble::tibble(row = rownames(m)[row(m)[ut]], col = colnames(m)[col(m)[ut]], corr = m[ut])
  tib <- tib[order(abs(tib$corr), decreasing = TRUE), ]
  tib
}

log_likelihood_crm <- function(stan_data, params, nr, idx_obs) {

  locations <- exp(params$log_a[stan_data$idx_rater]) * c(t(params$lt)) + params$b[stan_data$idx_rater]
  scales    <- exp(params$log_E[stan_data$idx_rater] + params$log_lambda[stan_data$idx_item])
  sum(dnorm(stan_data$x, locations, scales, log = TRUE))

}

get_samples_tib <- function(stan_fit, dat) {

  samples <- rstan::extract(stan_fit, permuted = FALSE)
  obs_param_names <- unique(trimBrackets(dimnames(samples)[[3L]]))
  all_param_names <- names(dat$parameters)
  all_param_names <- intersect(all_param_names, obs_param_names)

  post_mean_lst <- setNames(lapply(all_param_names, get_means, samples = samples), all_param_names)

  true_parameters <- dat$parameters[all_param_names]
  return(tibble::tibble(
    posterior_mean = unlist(lapply(all_param_names, get_means, samples = samples), use.names = FALSE),
    true_value     = unlist(true_parameters, use.names = FALSE),
    parameter      = rep(names(true_parameters), lengths(true_parameters)),
    # color          = color,
    bias           = posterior_mean - true_value,
    param_index    = unlist(lapply(lengths(post_mean_lst), seq_len), use.names = FALSE)
  ))
}

