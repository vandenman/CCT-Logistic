# #' @details see Asgharzadeh, A., Esmaeili, L., Nadarajah, S., & Shih, S. H. (2013). A generalized skew logistic distribution. REVSTAT–Statistical Journal, 11(3), 317-338.
# #' @rdname skew_logistic_distribution
# #' @export
# dslogis <- function(x, location, scale, skew = 0, log = FALSE) {
#   # definition 1.5 of Asgharzadeh, Esmaily, and Nadarajah (2015).
#   value <- log(2) - log(scale) + dlogis(x, location, scale, log = TRUE) + plogis(skew * (x - location) / scale, log.p = TRUE)
#
#   if (log)
#     return(value)
#   return(exp(value))
#
# }
#
# #' @rdname skew_logistic_distribution
# #' @export
# rslogis <- function(x, location, scale, skew = 0) {
#
# }
#
# #' @rdname skew_logistic_distribution
# #' @export
# pslogis <- function(x, location, scale, skew = 0) {
#
# }
#
#
# #' @rdname skew_logistic_distribution
# #' @export
# qslogis <- function(x, location, scale, skew = 0) {
#
# }


# mathematica ----

# SKL[y_, l_, s_, k_] :=(2/s) *PDF[LogisticDistribution[0, 1], (y - l) / s]*CDF[LogisticDistribution[0, 1], k*(y - l) / s]
# Plot[{SKL[y, 0, 1, 1], SKL[y, 0, 1, 5], SKL[y, 0, 1, 10]},{y, -5, 5}]

# Exponential Logistic Generalized Weibull ----
log1pexp <- function(x) {
  # maybe swap this out for an improved implementation
  log1p(exp(x))
}

#' @title The exponential logistic generalized Weibull distribution
#' @rdname exponential-logistic-generalized-Weibull
#' @details see Aljarrah, M.A., Famoye, F. & Lee, C. Generalized logistic distribution and its regression model. J Stat Distrib App 7, 7 (2020). https://doi.org/10.1186/s40488-020-00107-8
#' @export
dELGW <- function(x, location, scale, shape, log = FALSE) {

  # corollary 1.a
  if (shape == 0)
    return(dlogis(x, location, scale, log))

  z <- (x - location) / scale * sign(shape)
  shape <- abs(shape)
  res <- -log(scale) + z + (shape - 1) * log1p(exp(z)) + (1 - (1 + exp(z))^shape) / shape

  if (log)
    return(res)
  return(exp(res))
}

#' @rdname exponential-logistic-generalized-Weibull
#' @export
pELGW <- function(q, location, scale, shape, log = FALSE) {

  if (shape == 0)
    return(plogis(q, location, scale, TRUE, log))

  z <- (q - location) / scale * sign(shape)
  # tricky to compute this on a log scale by default
  direction <- sign(shape) # sign(scale)
  shape <- abs(shape)
  res <- 0.5 + direction * (0.5 - exp((1 - (1 + exp(z))^shape) / shape))

  if (log)
    return(log(res))
  return(res)

}

#' @rdname exponential-logistic-generalized-Weibull
#' @export
qELGW <- function(p, location, scale, shape, log = FALSE) {

  if (shape == 0)
    return(qlogis(p, location, scale, TRUE, log))

  # tricky to compute this on a log scale by default
  direction <- sign(shape) # sign(scale)
  shape <- abs(shape)
  res <- location + direction * scale * log((1 - shape * log(0.5 - direction * (p - 0.5))) ^ (1 / shape) - 1)

  if (log)
    return(log(res))
  return(res)

}

#' @rdname exponential-logistic-generalized-Weibull
#' @export
rELGW <- function(n, location, scale, shape) {
  if (shape == 0)
    return(rlogis(n, location, scale))
  return(qELGW(runif(n), location, scale, shape))
}
