#' @export
compute_probs <- function(draws, stan_data, stan_variable_sizes, return_mean = TRUE) {

  # this is slower than having stan do the work, but it's MUCH less memory intense to do it this way

  nc             <- stan_data$nc
  n_observed     <- stan_data$n_observed
  no_time_points <- stan_data$no_time_points

  idx_patient    <- stan_data$idx_patient
  idx_item       <- stan_data$idx_item
  idx_rater      <- stan_data$idx_rater
  idx_time_point <- stan_data$idx_time_point

  use_free_logistic_thresholds <- stan_data$use_free_logistic_thresholds
  use_skew_logistic_thresholds <- stan_data$use_skew_logistic_thresholds
  vary_lambda_across_patients  <- stan_data$vary_lambda_across_patients

  cnms <- colnames(draws)

  idx_lt              <- cnms[startsWith(cnms, "lt[")]
  idx_b               <- cnms[startsWith(cnms, "b[")]
  idx_log_a           <- cnms[startsWith(cnms, "log_a[")]
  idx_log_E           <- cnms[startsWith(cnms, "log_E[")]
  idx_log_lambda      <- cnms[startsWith(cnms, "log_lambda[")]
  idx_free_thresholds <- cnms[startsWith(cnms, "free_thresholds[")]
  idx_threshold_shape <- cnms[startsWith(cnms, "threshold_shape[")]
  idx_offset_lt       <- cnms[startsWith(cnms, "offset_lt[")]

  default_thresholds      <- (1:(nc - 1))
  default_thresholds      <- log(default_thresholds / (nc - default_thresholds))
  default_threshold_probs <- (1:(nc - 1)) / nc

  no_draws <- nrow(draws)
  if (return_mean) {
    mean_log_probs <- matrix(0, nc, n_observed)
  } else {
    individual_log_probs <- matrix(NA, no_draws, nc * n_observed,
                                   dimnames = list(NULL, paste0("log_probs[", rep(1:nc, n_observed), ", ", rep(1:n_observed, each = nc), "]")))
  }

  for (d in seq_len(no_draws)) {

    lt              <- array(draws[d, idx_lt             ], stan_variable_sizes[["lt"]]              %||% 1)
    b               <- array(draws[d, idx_b              ], stan_variable_sizes[["b"]]               %||% 1)
    log_a           <- array(draws[d, idx_log_a          ], stan_variable_sizes[["log_a"]]           %||% 1)
    log_E           <- array(draws[d, idx_log_E          ], stan_variable_sizes[["log_E"]]           %||% 1)
    log_lambda      <- array(draws[d, idx_log_lambda     ], stan_variable_sizes[["log_lambda"]]      %||% 1)
    free_thresholds <- array(draws[d, idx_free_thresholds], stan_variable_sizes[["free_thresholds"]] %||% 1)
    threshold_shape <- array(draws[d, idx_threshold_shape], stan_variable_sizes[["threshold_shape"]] %||% 1)
    offset_lt       <- array(draws[d, idx_offset_lt      ], stan_variable_sizes[["offset_lt"]]       %||% 1)

    log_probs <- matrix(NA, nc, n_observed)

    for (o in seq_len(n_observed)) {

      p <- idx_patient[o]
      i <- idx_item[o]
      r <- idx_rater[o]

      location <- lt[p, i]
      scale <- exp(log_lambda[i, if(vary_lambda_across_patients) p else 1] - log_E[r])

      if (no_time_points == 2) {

        t <- idx_time_point[o]
        # (2 * t - 3) maps {1, 2} to {-1, 1}
        location <- location + (2 * t - 3) * offset_lt[p, i]

      }

      if (use_free_logistic_thresholds) {

        log_probs[, o] <- get_probs_ordered_logistic(free_thresholds[r, ], location, scale)

      } else {

        threshold_scale = exp(log_a[r])
        threshold_shift = b[r]

         if (use_skew_logistic_thresholds) {
           log_probs[, o] <- get_probs_ordered_logistic(
             qELGW(default_threshold_probs, threshold_shift, threshold_scale, threshold_shape[r]),
             location,
             scale
           )
         } else {
           log_probs[, o] <- get_probs_ordered_logistic(
             qlogis(default_threshold_probs, threshold_shift, threshold_scale),
             location,
             scale
           )
         }

        # log_probs[, o] <- get_probs_ordered(
        #   qlogis(default_threshold_probs, threshold_shift, threshold_scale),
        #   location, scale
        # )

        # log_probs[, o] <- get_probs_ordered_logistic(
        #   qELGW(default_threshold_probs, threshold_shift, threshold_scale, if (use_skew_logistic_thresholds) threshold_shape[r] else 0),
        #   location,
        #   scale
        # )

      }
    }

    if (return_mean) {
      mean_log_probs <- mean_log_probs + log_probs / no_draws
    } else {
      individual_log_probs[d, ] <- c(log_probs)
    }

  }

  if (return_mean) {
    return(mean_log_probs)
  } else {
    return(individual_log_probs)
  }
}
