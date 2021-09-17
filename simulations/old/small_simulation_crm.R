rm(list = ls())
library(tibble)
library(ggplot2)
library(cmdstanr)
library(CRMLogistic)

# analysis ----
mod <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new.stan",
                     dir = "stanmodels", force_recompile = TRUE,
                     pedantic = TRUE, quiet = FALSE)#, cpp_options = list(stan_opencl = TRUE))

np <- 10
ni <- 15
nr <- 40
no_rater_groups <- 2
debugonce(simulate_data)
dat <- simulate_data(np, ni, nr, no_rater_groups = no_rater_groups)


tmp_mu  <- 1
tmp_var <- 0.1

# debugonce(dat_2_stan)
stan_data <- dat_2_stan(
  dat#,
  # required to constrain the mle, but not needed for vb or mcmc
  # a_sd_log_a      = (tmp_mu / tmp_var)^2, b_sd_log_a      = tmp_mu / tmp_var^2,
  # a_sd_log_lambda = (tmp_mu / tmp_var)^2, b_sd_log_lambda = tmp_mu / tmp_var^2
)
stan_data

hyperparams <- c("mu_lt", "mu_log_a", "mu_b", "mu_log_E", "sd_lt", "sd_log_lambda", "sd_log_E", "sd_log_a", "sd_b")

# fit_mle <- mod$optimize(data = stan_data, threads = 4, iter = 2e4)
# tib_mle <- get_means_tib(fit_mle, dat)
# tib_mle |> scatterplot_retrieval() + ggtitle(sprintf("mod$optimize: %.3f seconds", fit_mle$time())) + labs(x = "Estimate", y = "True value")


stan_data <- dat_2_stan(dat, vectorized = FALSE)
fit_vb <- mod$variational(data = stan_data, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3#,
                          # threads = 6,
                          # opencl_ids = c(1, 0), refresh = 0
)
tib_vb <- get_means_tib(fit_vb, dat)
tib_vb |> scatterplot_retrieval() + ggtitle(sprintf("mod$variational: %.3f seconds", fit_vb$time())) + labs(x = "Variational Estimate", y = "True value")
system("beep_finished.sh")


tib_log_a_vb <- tib_vb |> subset(parameter == "log_a")
mean(tib_log_a_vb$estimate)

sum_vb <- fit_vb$summary()
sum_hyperparams_vb <- sum_vb |> subset(grepl(paste0("^(", paste(hyperparams, collapse = "|"), ")"), variable))
sum_hyperparams_vb

fit_mcmc <- mod$sample(data = stan_data, parallel_chains = 6, iter_sampling = 4e3, chains = 6, threads_per_chain = 1)
tib_mcmc <- get_means_tib(fit_mcmc, dat)
tib_mcmc |> scatterplot_retrieval() + ggtitle(sprintf("mod$sample: %.3f seconds", fit_mcmc$time()$total)) + labs(x = "MCMC Estimate", y = "True value")
system("beep_finished.sh")

# inform vb with results from the mle --  this does not need to improve things since the mle may diverge
# params <- get_params(fit_mle, filter = "lp__")
# init <- list(setNames(lapply(params, \(x) unlist(fit_mle$draws(x, format = "draws_list"), use.names = FALSE)), params))
# fit_vb2 <- mod$variational(data = stan_data, threads = 4, init = init, iter = 2e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
# tib_vb2 <- get_means_tib(fit_vb2, dat)
# tib_vb2 |> scatterplot_retrieval() + ggtitle(sprintf("mod$variational: %.3f seconds", fit_vb$time())) + labs(x = "Variational Estimate", y = "True value")
# system("beep_finished.sh")


mod_cpu <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new.stan", dir = "stanmodels", pedantic = TRUE, quiet = FALSE)
mod_gpu <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new.stan", dir = "stanmodels", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_opencl = TRUE))

stan_data_vec     <- dat_2_stan(dat, vectorized = FALSE)
stan_data_non_vec <- dat_2_stan(dat, vectorized = FALSE)
fit_vb_non_vec <- mod_cpu$variational(data = stan_data_non_vec, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb_vec_cpu <- mod_cpu$variational(data = stan_data_vec,     iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb_vec_gpu <- mod_gpu$variational(data = stan_data_vec,     iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3, opencl_ids = c(1, 0), refresh = 0)

timings <- tibble(
  method = c("cpu_loop", "cpu_vectorized", "gpu_vectorized"),
  time   = bench::as_bench_time(c(fit_vb_non_vec$time()$total, fit_vb_vec_cpu$time()$total, fit_vb_vec_gpu$time()$total)),
  rel    = c(time / min(time))
)
timings


benchmark_table(non_vec = fit_vb_non_vec, vec_cpu = fit_vb_vec_cpu, vec_gpu = fit_vb_vec_gpu)

mod_cpu0 <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new.stan",         dir = "stanmodels", pedantic = TRUE, quiet = FALSE)
mod_cpu  <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new_no_genq.stan", dir = "stanmodels", pedantic = TRUE, quiet = FALSE)
mod_gpu  <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new_no_genq.stan", dir = "stanmodels", pedantic = TRUE, quiet = FALSE, cpp_options = list(stan_opencl = TRUE))

np <- 20
ni <- 25
nr <- 50
no_rater_groups <- 1
dat <- simulate_data(np, ni, nr, no_rater_groups = no_rater_groups)
dim(dat$df)

stan_data_vec     <- dat_2_stan(dat, vectorized = FALSE)
stan_data_non_vec <- dat_2_stan(dat, vectorized = FALSE)

fit_vb_non_vec0 <- mod_cpu0$variational(data = stan_data_non_vec, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb_non_vec  <- mod_cpu$variational (data = stan_data_non_vec, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb_vec_cpu  <- mod_cpu$variational (data = stan_data_vec,     iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb_vec_gpu  <- mod_gpu$variational (data = stan_data_vec,     iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3, opencl_ids = c(1, 0), refresh = 0)

timings <- benchmark_table(cpu_loop0 = fit_vb_non_vec0, cpu_loop = fit_vb_non_vec, cpu_vectorized = fit_vb_vec_cpu, gpu_vectorized = fit_vb_vec_gpu)
timings
system("beep_finished.sh")



np <- 20
ni <- 25
nr <- 10
no_rater_groups <- 2
dat <- simulate_data(np, ni, nr, no_rater_groups = no_rater_groups)

data_2_init(dat$df) |> add_true_values_to_em_tib(dat) |> scatterplot_retrieval()

em_tib <- fit_em(dat$df, record_all_iterations = 2, n_iter = 10, compute_locations_scales = TRUE)
em_tib |> add_true_values_to_em_tib(dat) |> scatterplot_retrieval(facets = iteration~parameter)




qq <- Q_sum_to_zero_QR_R(5)
ee <- rnorm(4)
ff <- sum_to_zero_QR_R(ee, qq)
sum(ff)
gg <- sum_to_zero_QR_inv_R(ff, qq)
ee - gg


em_tib |>
  subset(parameter %in% c("log_a", "b", "log_E", "log_lambda", "lt") & iteration == max(iteration)) |>
  split(~parameter) |>
  lapply(\(x) x$estimate) |>
  compute_raw_params() |>
  compute_hyperparams()

zapsmall(tapply(em_tib$estimate, list(em_tib$parameter, em_tib$iteration), mean))
tapply(em_tib$estimate, list(em_tib$parameter, em_tib$iteration), sd)
