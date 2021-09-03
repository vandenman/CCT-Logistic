rm(list = ls())
library(rstan)
library(tibble)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = 8)

# functions ----
pplot <- function(x, y, ...) {
  plot(x, y, ...)
  abline(0, 1)
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

scatterplot_retrieval <- function(tib) {
  ggplot(data = tib, mapping = aes(x = true_value, y = posterior_mean)) +
    geom_abline() +
    geom_point() +
    facet_wrap(facets = ~parameter, scales = "free") +
    theme_bw()
}

simulate_data <- function(np, ni, nr, lt = NULL, log_lambda = NULL, log_a = NULL, b = NULL, log_E = NULL, return_mus_sds = TRUE) {

  lt     <- matrix(rnorm(ni*np), ni, np)
  lambda <- matrix(rgamma(ni*np, 5, 2), ni, np)
  # lambda <- matrix(1, ni, np)
  log_lambda <- log(lambda)

  log_a <- rnorm(nr, 1)
  b     <- rnorm(nr)
  log_E <- rnorm(nr, 1)

  x <- array(NA, c(np, ni, nr), dimnames = list(patient = seq_len(np), item = seq_len(ni), rater = seq_len(nr)))

  if (return_mus_sds) {
    mus <- array(NA, c(np, ni, nr), dimnames = dimnames(x))
    sds <- array(NA, c(np, ni, nr), dimnames = dimnames(x))
  } else {
    mus <- sds <- NULL
  }

  # TODO: this can probably be vectorized over some dimension
  for (p in seq_len(np)) {
    for (i in seq_len(ni)) {
      for (r in seq_len(nr)) {

        mu <- exp(log_a[r]) * lt[i, p] + b[r]
        sd <- exp(log_E[r] + log_lambda[i, p])
        if (return_mus_sds) {
          mus[p, i, r] <- mu
          sds[p, i, r] <- sd
        }
        x[p, i, r] <- rnorm(1, mu, sd)
      }
    }
  }
  return(list(
    x = x,
    np = np, ni = ni, nr = nr,
    parameters = list(
      lt = lt, log_lambda = log_lambda,
      log_a = log_a, b = b, log_E = log_E,
      mus = mus, sds = sds
    )
  ))
}

dat_2_stan <- function(dat, prior = FALSE, debug = TRUE) {

  x_df <- as.data.frame.table(dat$x, responseName = "score")
  all(x_df$score == dat$x[cbind(x_df$patient, x_df$item, x_df$rater)])

  data <- list(
    x                   = x_df$score,
    np                  = dat$np,
    ni                  = dat$ni,
    nr                  = dat$nr,

    n_observed          = length(x_df$score),
    idx_patient         = as.integer(x_df$patient),
    idx_item            = as.integer(x_df$item),
    idx_rater           = as.integer(x_df$rater),

    n_missing           = 0L,
    idx_patient_missing = integer(),
    idx_item_missing    = integer(),
    idx_rater_missing   = integer(),

    debug      = as.integer(debug),
    prior_only = as.integer(prior)
  )
  return(data)
}

# inspect prior ----

np <- 3
ni <- 20
nr <- 5

dat <- simulate_data(np, ni, nr)

mod <- stan_model(file = "stanmodels/extendedCRM_long.stan", verbose = TRUE)
fit <- sampling(mod, data = dat_2_stan(dat, prior = TRUE), iter = 5e3, chains = 6)

samples <- extract(fit, permuted = FALSE)
dim(samples)

mu_log_E = get_means(samples, "mu_log_E")
sd_log_E = get_means(samples, "tau_log_E")
log_E_prior_samples = samples[, , "log_E[1]"]

plot(density(log_E_prior_samples))
plot(\(x) dnorm(x, mu_log_E, sd_log_E), from = -5, to = 5, col = 2, add = TRUE)

mu_log_a = 0.0# fixed get_means(samples, "mu_log_a")
sd_log_a = get_means(samples, "tau_log_a")
log_a_prior_samples = samples[, , "log_a[2]"]

plot(density(log_a_prior_samples))
plot(\(x) dnorm(x, mu_log_a, sd_log_a), from = -5, to = 5, col = 2, add = TRUE)

# simulations with data ----

np <- 3
ni <- 20
nr <- 5

dat <- simulate_data(np, ni, nr)

x <- dat$x
tibble(
  true_value = c(dat$parameters$lt),
  sum_score  = c(t(apply(x, 1:2, mean))),
  patient    = factor(rep(seq_len(np), each = ni))
) |> ggplot(mapping = aes(x = true_value, y = sum_score)) +
  geom_abline() +
  geom_point() +
  facet_wrap(~patient) +
  theme_bw()


mod <- stan_model(file = "stanmodels/extendedCRM_long.stan", verbose = TRUE)
fit <- sampling(mod, data = dat_2_stan(dat), iter = 5e3, chains = 6)
# fit <- vb(mod, data = data, iter = 5e3 * 6)
# system("beep_finished.sh")

samples <- extract(fit, permuted = FALSE)
dim(samples)

get_means(samples, "mu_lt")
get_means(samples, "sd_lt")
get_means(samples, "sd_log_lambda")


obsParamNames <- unique(trimBrackets(dimnames(samples)[[3]]))
allParamNames <- names(dat$parameters)
allParamNames <- intersect(allParamNames, obsParamNames)

post_means <- tibble(
  posterior_mean = unlist(lapply(allParamNames, get_means, samples = samples), use.names = FALSE),
  true_value     = unlist(dat$parameters, use.names = FALSE),
  parameter      = rep(names(dat$parameters), lengths(dat$parameters))
)

post_means |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()
post_means |> subset(parameter %in% c("mus", "sds")) |> scatterplot_retrieval()

system("beep_finished.sh")



samples_mat <- do.call(rbind, lapply(1:ncol(samples), \(x) samples[, x, ]))
dim(samples_mat)
cor_samps <- cor(samples_mat)
cor_2_tib <- function(m) tibble(row = rownames(m)[row(m)[upper.tri(m)]], col = colnames(m)[col(m)[upper.tri(m)]], corr = m[upper.tri(m)])

tib_cor <- cor_2_tib(cor_samps)
tib_cor[order(abs(tib_cor$corr), decreasing = TRUE), ]

plot(samples_mat[, "log_E[8]"], samples_mat[, "log_a[8]"])
plot(samples_mat[, "lt[1,1]"], samples_mat[, "lt[2,1]"])

