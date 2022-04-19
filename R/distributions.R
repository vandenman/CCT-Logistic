get_probs_ordered <- function(thresholds, location, scale, cdf = plogis) {

  # replica of what Stan does
  nc <- length(thresholds) + 1L
  prob <- numeric(nc)

  vals <- cdf(location - thresholds, 0, scale)
  prob[1L] <- 1 - vals[1L]
  for (c in 2:(nc - 1)) {
    prob[c] <- vals[c - 1L] - vals[c]
  }
  prob[nc] <- vals[nc - 1]

  return(prob)
}

get_probs_ordered_logistic <- function(thresholds, location, scale) {
  # optimized implementation of get_probs_ordered for logistic cdf
  vals <- 1 / (1 + exp((thresholds - location) / scale))
  c(
    1 - vals[1L],
    vals[-length(vals)] - vals[-1L],
    vals[length(vals)]
  )
}


#' @export
cutpoints_to_probs <- function(cutpoints) {

  # replica of what Stan does
  nc <- length(cutpoints) + 1L
  prob <- numeric(nc)

  prob[1L] <- 1 - cutpoints[1L]
  for (c in 2:(nc - 1)) {
    prob[c] <- cutpoints[c - 1L] - cutpoints[c]
  }
  prob[nc] <- cutpoints[nc - 1]

  return(prob)
}


cutpoints_to_probs_fast <- function(cutpoints) {

  c(
    1 - cutpoints[1L],
    cutpoints[-length(cutpoints)] - cutpoints[-1L],
    cutpoints[length(cutpoints)]
  )
}

#' @export
dordered_logistic <- function(x, thresholds, location = 0, scale = 1, log = TRUE) {
  dordered(x, log, stats::plogis(location - thresholds, 0, scale))
}

#' @export
dordered_skew_logistic <- function(x, thresholds, location = 0, scale = 1, shape = 0, log = TRUE) {
  dordered(x, log, pELGW(location - thresholds, 0, scale, shape))
}

#' @export
rordered_logistic <- function(n, thresholds, location = 0, scale = 1) {
  rordered(n, stats::plogis(location - thresholds, 0, scale))
}

#' @export
rordered_skew_logistic <- function(x, thresholds, location = 0, scale = 1, shape = 0) {
  rordered(x, pELGW(location - thresholds, 0, scale, shape))
}

#' @export
rordered <- function(n, cutpoints) {
  n == 0L && return(integer())
  prob <- cutpoints_to_probs_fast(cutpoints)
  sample(length(cutpoints) + 1L, size = n, replace = TRUE, prob = prob)
}

#' @export
dordered <- function(x, cutpoints, log = TRUE) {
  probs <- cutpoints_to_probs_fast(cutpoints)
  if (log)
    return(log(prob[x]))
  else return(prob[x])
}

