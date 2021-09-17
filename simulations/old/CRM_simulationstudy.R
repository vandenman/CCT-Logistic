rm(list = ls())
# library(rstan)
library(tibble)
library(ggplot2)
library(cmdstanr)
source(file.path("simulations", "utils.R"))
source(file.path("simulations", "CRM_helpers.R"))
# rstan_options(auto_write = TRUE)
# options(mc.cores = 8)

# inspect prior ----

np <- 3
ni <- 10
nr <- 5

dat <- simulate_data(np, ni, nr)
sum(lengths(dat$parameters)[1:5])
mod <- stan_model(file = "stanmodels/extendedCRM_long.stan", verbose = TRUE)
fit <- sampling(mod, data = dat_2_stan(dat, prior_only = TRUE), iter = 5e3, chains = 6)

samples <- extract(fit, permuted = FALSE)
dim(samples)

n_points <- 2^10
mu_log_E = get_means(samples, "mu_log_E")
sd_log_E = get_means(samples, "sd_log_E")
mu_log_a = 0.0# fixed get_means(samples, "mu_log_a")
sd_log_a = get_means(samples, "sd_log_a")
dens_log_E <- density(samples[, , "log_E[1]"], n = n_points)
dens_log_a <- density(samples[, , "log_a[1]"], n = n_points)
xrange <- range(dens_log_E$x, dens_log_a$x)
xx <- seq(xrange[1], xrange[2], length.out = n_points)
tibble(
  x = c(dens_log_E$x, dens_log_a$x, xx, xx),
  y = c(dens_log_E$y, dens_log_a$y, dnorm(xx), dnorm(xx)),
  g = rep(c("Sample", "Density"), each = 2*n_points),
  parameter = rep(rep(c("log_E", "log_a"), each = n_points), 2)
) |> ggplot(mapping = aes(x = x, y = y, group = g, color = g)) +
  geom_line() + facet_wrap(~parameter, scales = "free_y", strip.position = "bottom") + labs(y = "Density", x = NULL, color = NULL) +
  ggplot2::scale_x_continuous(breaks = jaspGraphs::getPrettyAxisBreaks(xrange), limits = range(jaspGraphs::getPrettyAxisBreaks(xrange))) +
  jaspGraphs::geom_rangeframe() + jaspGraphs::themeJaspRaw(legend.position = c(.9, .99)) + theme(strip.placement = "outside")

pars_2_plot <- c(
  paste0("lt[", 1:2, ",", rep(1:1, each = 2), "]"),
  paste0("log_lambda[", 1:2, ",", rep(1:1, each = 2), "]"),
  paste0("b[", 1:2, "]"),
  paste0("log_E[", 1:2, "]"),
  paste0("log_a[", 1:2, "]")
)
all(pars_2_plot %in% names(fit))
w <- 7*length(pars_2_plot)
plt <- bayesplot::mcmc_pairs(samples, pars = pars_2_plot, diag_fun = "dens", off_diag_fun = "hex")
class(plt)

pdf("figures/prior_plot.pdf", width = w, height = w)
plot(plt)
dev.off()

paramNames <- dimnames(samples)[[3]]
idxSkipMusAndSds <- paramNames[!grepl("(mus|sds|lp__|lt\\[|log_a\\[|log_lambda\\[)", paramNames)]
samples_mat <- do.call(rbind, lapply(1:ncol(samples), \(x) samples[, x, idxSkipMusAndSds]))
dim(samples_mat)
cor_samps <- cor(samples_mat)
tib_cor <- cor_2_tib(cor_samps)

# plot(samples_mat[, log])


# simulations with data ----

np <- 10
ni <- 15
nr <- 20

# debugonce(simulate_data)
dat <- simulate_data(np, ni, nr)#, log_lambda = matrix(rnorm(np*ni, 0, .25), ni, np))
dat$parameters$log_lambda
mean(dat$parameters$log_a)
mean(dat$parameters$log_E)

tibble(
  true_value = c(dat$parameters$lt),
  sum_score  = c(t(apply(dat$x, 1:2, mean))),
  patient    = factor(rep(seq_len(np), each = ni))
) |> ggplot(mapping = aes(x = true_value, y = sum_score)) +
  geom_abline() +
  geom_point() +
  facet_wrap(~patient, labeller = \(x) data.frame(patient = paste("patient:", x$patient))) +
  theme_bw()


# mcmc_iter   <- 5e3
# mcmc_chains <- 6
# stan_data <- dat_2_stan(dat, debug = FALSE)

# benchmark_model <- function(file, stan_data, mcmc_iter, mcmc_chains) {
#   mod <- stan_model(file = file, verbose = TRUE)
#   fit <- sampling(mod, data = stan_data, iter = mcmc_iter, chains = mcmc_chains)
#   timings <- get_elapsed_time(fit)
#   return(list(mod = mod, fit = fit, timings = timings))
# }
#
# rr1 <- benchmark_model("stanmodels/extendedCRM_long_nonvectorized.stan", stan_data, mcmc_iter, mcmc_chains)
# rr2 <- benchmark_model("stanmodels/extendedCRM_long_vectorized.stan",    stan_data, mcmc_iter, mcmc_chains)

# benches <- list(
#   nonvectorized = rr1,
#   vectorized    = rr2
# )
# rstan::check_divergences(rr1$fit)
# rstan::check_divergences(rr2$fit)
# lapply(benches, `[[`, "timings")
# sapply(benches, \(x) colMeans(x$timings))
#
# mod2 <- stan_model(file = "stanmodels/extendedCRM_long_nonvectorized.stan", verbose = TRUE)
# fit2 <- sampling(mod2, data = stan_data, iter = mcmc_iter, chains = mcmc_chains)
# timings2 <- get_elapsed_time(fit2)
#
# mod3 <- stan_model(file = "stanmodels/extendedCRM_long_vectorized_2.stan", verbose = TRUE)
# stan_data3 <- stan_data
# fit3 <- sampling(mod3, data = stan_data3, iter = mcmc_iter, chains = mcmc_chains)
# timings3 <- get_elapsed_time(fit3)
#
# colMeans(timings)
# colMeans(timings2)
# colMeans(timings3)
#
# tib1 <- get_samples_tib(rr1$fit, dat)
# tib1 |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()
#
# tib2 <- get_samples_tib(rr2$fit, dat)
# tib2 |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()

mod <- stan_model(file = "stanmodels/extendedCRM_long_new.stan", verbose = TRUE)

mcmc_iter   <- 5e3
mcmc_chains <- 6
stan_data <- dat_2_stan(dat, debug = FALSE)
# fit <- sampling(mod, data = stan_data, iter = mcmc_iter, chains = mcmc_chains)

point <- optimizing(mod, data = stan_data)

fit <- vb(mod, data = stan_data, iter = 1e4, output_samples = mcmc_iter, grad_samples = 10)

samples_tib <- get_samples_tib(fit, dat)
samples_tib |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()
samples_tib |> subset( parameter %in% c("mus", "sds")) |> scatterplot_retrieval()


samples_tib |> subset( parameter %in% c("mus", "sds")) |> scatterplot_retrieval(aes(x = param_index, y = bias))

cpp_options = list(stan_threads = TRUE)
mod <- cmdstan_model("stanmodels/extendedCRM_long_new.stan", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_threads = TRUE))

# system("beep_finished.sh")
#
# fit <- rr2$fit
# samples <- extract(fit, permuted = FALSE)
# dim(samples)
#
# get_means(samples, "mu_log_E")
#
# idx_log_E <- which(startsWith(dimnames(samples)[[3]], "log_E["))
# samps_log_E <- stack_chain(samples[, , idx_log_E, drop = FALSE])
# rowMeans(samps_log_E)
# apply(samples[, , idx_log_E], 3, mean) == get_means(samples, "log_E")
#
# means_log_E <- center(get_means(samples, "log_E"))
# pplot(dat$parameter$log_E, means_log_E)
#
# idx_log_lambda <- which(startsWith(dimnames(samples)[[3]], "log_lambda["))
# samps_log_lambda <- stack_chain(samples[, , idx_log_lambda, drop = FALSE])
# rowMeans(samps_log_lambda)
# colMeans(samples[, , idx_log_E]) == get_means(samples, "mu_log_lambda")
#
#
# get_means(samples, "mu_lt")
# get_means(samples, "sd_lt")
# get_means(samples, "sd_log_lambda")
#
# # dat$parameters$lt <- t(dat$parameters$lt)
# obs_param_names <- unique(trimBrackets(dimnames(samples)[[3]]))
# all_param_names <- names(dat$parameters)
# all_param_names <- intersect(all_param_names, obs_param_names)
#
# true_parameters <- dat$parameters[all_param_names]
# post_means <- tibble(
#   posterior_mean = unlist(lapply(all_param_names, get_means, samples = samples), use.names = FALSE),
#   true_value     = unlist(true_parameters, use.names = FALSE),
#   parameter      = rep(names(true_parameters), lengths(true_parameters))
# )
#
# post_means |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()
# post_means |> subset( parameter %in% c("mus", "sds")) |> scatterplot_retrieval()
#
# system("beep_finished.sh")
#
# paramNames <- dimnames(samples)[[3]]
# idxSkipMusAndSds <- paramNames[!grepl("(mus|sds|lp__|_raw\\[)", paramNames)]
# samples_mat <- do.call(rbind, lapply(1:ncol(samples), \(x) samples[, x, idxSkipMusAndSds]))
# dim(samples_mat)
# cor_samps <- cor(samples_mat)
# cor_2_tib <- function(m) {
#   ut <- upper.tri(m)
#   tib <- tibble(row = rownames(m)[row(m)[ut]], col = colnames(m)[col(m)[ut]], corr = m[ut])
#   tib <- tib[order(abs(tib$corr), decreasing = TRUE), ]
#   tib
# }
#
# tib_cor <- cor_2_tib(cor_samps)
#
# print(tib_cor[500:800,], n = 150)
#
# idx_no_lt <- idxSkipMusAndSds[!startsWith(idxSkipMusAndSds, "lt[")]
# tib_cor_no_lt <- cor_2_tib(cor_samps[idx_no_lt, idx_no_lt])
# print(tib_cor_no_lt, n = 150)
# get_means(samples, "log_lambda[17,3]")
# plot(samples_mat[, "log_lambda[1,1]"])
#
# layout(matrix(1:16, 4, 4))
# for (i in 1:16) {
#   xx <- samples_mat[, tib_cor$row[i]]
#   yy <- samples_mat[, tib_cor$col[i]]
#   plot(yy, xx, xlab = tib_cor$row[i], ylab = tib_cor$col[i])
#   # plot(lm(yy ~ xx), add = TRUE)
# }
#
# debugonce(get_means)
#
#
# max(tib_cor$corr)
#
# plot(samples_mat[, "lt[1,1]"], samples_mat[, "lt[2,1]"])


# simulation study with similar properties to data set ----

# ddd <- read_long_data()
# length(unique(ddd$patient)) # 105
# length(unique(ddd$item))    # 23
# length(unique(ddd$rater))   # 188
# # patient rater combinations
# length(unique(paste(ddd$patient, ddd$rater, sep = "-"))) # 611
# mean(table(ddd$rater, ddd$patient) != 0L) # 0.03095238, so on average a rater rates about 3.1% of patients
#
# tt <- table(ddd$rater) / 23 # no patients rated by the raters
# obs_props <- table(tt) / sum(table(tt)) # normalized to probabilities
# obs_props["40"] <- 0.0 # nobody rated more than 30 patients so we set the max at 40
# pdf <- approx(as.numeric(names(obs_props)), obs_props, xout = 1:40)$y # linearly interpolate missing probabilities.
# pdf <- pdf / sum(pdf) # ensure sum to 0
#
# plot(as.numeric(names(obs_props)), obs_props)
# lines(pdf)

# dput(pdf)

pdf <- c(0.264058679706601, 0.176039119804401, 0.0929095354523227, 0.0929095354523227,
         0.0440097799511002, 0.107579462102689, 0.0293398533007335, 0.0244498777506112,
         0.0244498777506112, 0.00488997555012225, 0.00488997555012225,
         0.00488997555012225, 0.00488997555012225, 0.00488997555012225,
         0.00488997555012225, 0.00733496332518337, 0.0097799511002445,
         0.0097799511002445, 0.0097799511002445, 0.00733496332518337,
         0.00488997555012225, 0.00488997555012225, 0.00488997555012225,
         0.00488997555012225, 0.00488997555012225, 0.00488997555012225,
         0.00488997555012225, 0.00488997555012225, 0.00488997555012225,
         0.00488997555012225, 0.00440097799511002, 0.0039119804400978,
         0.00342298288508557, 0.00293398533007335, 0.00244498777506112,
         0.0019559902200489, 0.00146699266503667, 0.00097799511002445,
         0.000488997555012225, 0)

np <- 105L
ni <- 23
nr <- 188
dat <- simulate_data(np, ni, nr)
x_df <- as.data.frame.table(dat$x, responseName = "score", stringsAsFactors = TRUE)
for (v in c("patient", "item", "rater"))
  x_df[[v]] <- as.integer(x_df[[v]])
mean(table(x_df$rater, x_df$patient) != 0L)

# TODO this is also not optimal!
# for each rater, sample the number of patients they saw and then sample those patients
all_patients <- unique(x_df$patient)
patients_not_sampled <- all_patients
combs_keep <- c()
for (r in seq_len(nr)) {
  nkeep <- sample(1:40, nr, prob = pdf, replace = TRUE)
  if (length(patients_not_sampled) > nkeep)
    patients_sampled <- sample(patients_not_sampled, nkeep, FALSE)
  else if (length(patients_not_sampled) == 0L) {
    patients_sampled <- sample(all_patients, nkeep, FALSE)
  } else {
    patients_sampled <- c(patients_not_sampled, sample(all_patients, nkeep - length(patients_not_sampled)))
  }
  rater_combs_keep <- paste(patients_sampled, r, sep = " - ")
  combs_keep <- c(combs_keep, rater_combs_keep)
}

obs_combinations <- paste(x_df$patient, x_df$rater, sep = " - ")
idx_keep <- obs_combinations %in% combs_keep

x_df_remaining <- x_df[idx_keep, ]
mean(table(x_df_remaining$rater, x_df_remaining$patient) != 0L)


# # TODO: this is closer to the empirical data set but the while loop sucks! improve this!
# obs_combinations <- paste(x_df$patient, x_df$rater, sep = " - ")
# unique_combinations <- apply(expand.grid(unique(x_df$patient), unique(x_df$rater)), 1L, paste, collapse = " - ")
# n_idx_keep <- floor(3.095238 / 100 * length(unique_combinations))
# x_df_remaining <- NULL
# while (is.null(x_df_remaining) || !(length(unique(x_df_remaining$rater)) == nr && length(unique(x_df_remaining$patient)) == np)) {
#   idx_remaining <- sample(length(unique_combinations), n_idx_keep, replace = FALSE)
#   combs_remaining  <- unique_combinations[idx_remaining]
#   idx_df_remaining <- which(obs_combinations %in% combs_remaining)
#   x_df_remaining <- x_df[idx_df_remaining, ]
# }

assertthat::assert_that(
  all(apply(x_df_remaining[, c("patient", "rater", "item")], 2L, \(x) all(table(x) > 0L))),
  all(apply(x_df_remaining[, c("patient", "rater", "item")], 2L, \(x) length(unique(x)) == max(x))),
  assertthat::are_equal(mean(table(x_df_remaining$rater, x_df_remaining$patient) != 0L), 0.03095238, tolerance = .1) # ~ 0.03 as desired
)

# mod <- stan_model(file = "stanmodels/extendedCRM_long.stan", verbose = TRUE)

stan_data <- dat_2_stan(x_df_remaining, debug = FALSE)
mcmc_iter   <- 5e3
mcmc_chains <- 6
# fit_mcmc <- sampling(mod, data = stan_data, iter = mcmc_iter, chains = mcmc_chains)

# point <- optimizing(mod, data = stan_data)

mod2   <- stan_model(file = "stanmodels/extendedCRM_long_notransformedparameters.stan", verbose = TRUE)
fit_vb <- vb(mod2, data = stan_data, output_samples = mcmc_iter, iter = 2e4)

samples <- extract(fit_vb, permuted = FALSE)
dim(samples)

# count the number of observations per item and per rater
patient_count_by_rater <- rowSums(table(x_df_remaining$rater, x_df_remaining$patient) / ni)
table(patient_count_by_rater)

rater_count_by_patient <- colSums(table(x_df_remaining$rater, x_df_remaining$patient) / ni)
table(rater_count_by_patient)

patient_params <- c("lt")
rater_params   <- c("log_a", "b", "log_E")

obs_param_names <- unique(trimBrackets(dimnames(samples)[[3]]))
all_param_names <- names(dat$parameters)
all_param_names <- intersect(all_param_names, obs_param_names)

post_mean_lst <- setNames(lapply(all_param_names, get_means, samples = samples), all_param_names)
color_lst <- lapply(names(post_mean_lst), \(x) {
  if (x %in% rater_params) {
    patient_count_by_rater
  } else if (x %in% patient_params) {
    rep(rater_count_by_patient, length.out = length(post_mean_lst[[x]]))
  } else {
    rep(NA_integer_, length(post_mean_lst[[x]]))
  }
})
assertthat::assert_that(all(lengths(color_lst) == lengths(post_mean_lst)))

lengths(dat$parameters)
lengths(post_mean_lst)

true_parameters <- dat$parameters[all_param_names]
post_means <- tibble(
  posterior_mean = unlist(lapply(all_param_names, get_means, samples = samples), use.names = FALSE),
  true_value     = unlist(true_parameters, use.names = FALSE),
  parameter      = rep(names(true_parameters), lengths(true_parameters)),
  # color          = color,
  bias           = posterior_mean - true_value,
  param_index    = unlist(lapply(lengths(post_mean_lst), seq_len), use.names = FALSE)
)


post_means |> subset(!parameter %in% c("mus", "sds")) |> scatterplot_retrieval()#mapping = aes(x = true_value, y = posterior_mean, color = color))
post_means |> subset( parameter %in% c("mus", "sds")) |> scatterplot_retrieval()

system("beep_finished.sh")

post_means |>
  subset(!parameter %in% c("mus", "sds")) |>
  ggplot(mapping = aes(x = param_index, y = bias, color = color)) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point() +
  facet_wrap(facets = ~parameter, scales = "free") +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw()

post_means |>
  subset(!parameter %in% c("mus", "sds")) |>
  ggplot(mapping = aes(x = color, y = abs(bias))) +
  geom_smooth(se = FALSE) +
  geom_point() +
  facet_wrap(facets = ~parameter, scales = "free") +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw()

by(post_means[, c("bias", "color")], post_means$parameter, \(x) cor(abs(x$bias), x$color))

tapply(post_means$bias, post_means$parameter, )

system("beep_finished.sh")


# simulation study comparing MCMC to vb ----

ddd <- read_long_data()
dim(ddd)
length(unique(ddd$patient)) # 105
length(unique(ddd$item))    # 23
length(unique(ddd$rater))   # 188
# patient rater combinations
length(unique(paste(ddd$patient, ddd$rater, sep = "-"))) # 611
tb <- table(ddd$rater, ddd$patient)
100 * mean(tb != 0L) # 3.095238, so on average a rater rates about 3.1% of patients


opts <- tibble(expand.grid(nps = c(50L, 100L, 150L), nis = c(10L, 20L, 30L), nrs = c(5L, 20L, 40L), repeats = 1:5))
opts
