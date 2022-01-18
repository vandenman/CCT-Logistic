saveRDS(object = cnms, "~/test/cnms.rds")

grep(pattern = "(.*)\\[{0,1}", x = cnms[1:20])

cnms2 <- gsub(pattern = "(\\[.*\\])", replacement = "", x = cnms)
table(cnms2)

log_probs_idx <- which(startsWith(cnms, "log_probs"))
not_log_probs_idx <- which(!startsWith(cnms, "log_probs"))

aa <- fits$fit_orig$.__enclos_env__$private$draws_[, log_probs_idx]
bb <- fits$fit_orig$.__enclos_env__$private$draws_[, not_log_probs_idx]
format(object.size(aa), units = "auto")
format(object.size(bb), units = "auto")

dim(bb)
colnames(bb)

stan_data <- data_2_stan(data_2_analyze, store_predictions = FALSE, debug = FALSE, use_skew_logistic_thresholds = FALSE, use_free_logistic_thresholds = FALSE)


n <- 5e4
asd <- rnorm(n)
print(mean(asd), 20)
print(Reduce(`+`, asd / n, 0), 20)

compute_probs <- function(draws, stan_data) {

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
  mean_log_probs <- matrix(0, nc, n_observed)

  for (d in seq_len(no_draws)) {

    lt              <- array(draws[d, idx_lt             ], meta$stan_variable_sizes[["lt"]]              %||% 1)
    b               <- array(draws[d, idx_b              ], meta$stan_variable_sizes[["b"]]               %||% 1)
    log_a           <- array(draws[d, idx_log_a          ], meta$stan_variable_sizes[["log_a"]]           %||% 1)
    log_E           <- array(draws[d, idx_log_E          ], meta$stan_variable_sizes[["log_E"]]           %||% 1)
    log_lambda      <- array(draws[d, idx_log_lambda     ], meta$stan_variable_sizes[["log_lambda"]]      %||% 1)
    free_thresholds <- array(draws[d, idx_free_thresholds], meta$stan_variable_sizes[["free_thresholds"]] %||% 1)
    threshold_shape <- array(draws[d, idx_threshold_shift], meta$stan_variable_sizes[["threshold_shape"]] %||% 1)
    offset_lt       <- array(draws[d, idx_offset_lt      ], meta$stan_variable_sizes[["offset_lt"]]       %||% 1)

    log_probs <- matrix(NA, nc, n_observed)

    for (o in seq_len(n_observed)) {

      p <- idx_patient[o];
      i <- idx_item[o];
      r <- idx_rater[o];

      location <- lt[p, i];
      scale <- exp(log_lambda[i, if(vary_lambda_across_patients) p else 1] - log_E[r]);

      if (no_time_points == 2) {

        t <- idx_time_point[o];
        # (2 * t - 3) maps {1, 2} to {-1, 1}
        location <- location + (2 * t - 3) * offset_lt[p, i]

      }

      if (use_free_logistic_thresholds) {

        log_probs[, o] <- cutpoints_to_probs((location - free_thresholds[r]) / scale)

      } else {

        threshold_scale = exp(log_a[r]);
        threshold_shift = b[r];

        if (use_skew_logistic_thresholds) {

          log_probs[, o] <- cutpoints_to_probs(location - qELGW(default_threshold_probs, threshold_shift, threshold_scale, threshold_shape)) / scale

        } else {

          log_probs[, o] <- cutpoints_to_probs((location - threshold_shift - threshold_scale * default_thresholds) / scale)

        }
      }
    }

    mean_log_probs <- mean_log_probs + log_probs / no_draws

  }
}
