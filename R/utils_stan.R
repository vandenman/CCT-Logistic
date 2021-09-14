#' @title compile a stan model
#' @inheritParams cmdstanr::cmdstan_model
#' @param dir
#' @details This is a helper for this project so you don't have to adjust any directories.
#' It's not recommended for use outside of this project, just use cmdstanr::cmdstan_model directly.
#' @export
compile_stan_model <- function(stan_file, compile = TRUE, dir = file.path("stanmodels", "compiled_stan_models"), ...) {
  cmdstanr::cmdstan_model(stan_file, compiler, dir = dir, ...)
}


Q_sum_to_zero_QR_R <- function(N) {
  Q_r <- numeric(2*N)

  for (i in seq_len(N)) {
    Q_r[i]     = -sqrt((N-i)/(N-i+1.0));
    Q_r[i + N] = 1 / sqrt((N-i) * (N-i+1));
  }
  Q_r
}

sum_to_zero_QR_R <- function(x_raw, Q_r) {
  N <- length(x_raw) + 1
  x <- numeric(N)
  x_aux = 0;

  for(i in 1:(N-1)) {
    x[i]  = x_aux + x_raw[i] * Q_r[i];
    x_aux = x_aux + x_raw[i] * Q_r[i+N];
  }
  x[N] = x_aux;
  x
}

sum_to_zero_QR_inv_R <- function(x, Q_r = Q_sum_to_zero_QR_R(length(x))) {
  N <- length(x)
  x_raw <- numeric(N - 1)
  x_aux <- 0;

  for (i in 1:(N-1)) {
    x_raw[i] = (x[i] - x_aux) / Q_r[i]
    x_aux    = x_aux + x_raw[i] * Q_r[i + N]
  }
  x_raw
}

