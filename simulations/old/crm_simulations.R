rm(list = ls())
library(cmdstanr)
library(CRMLogistic)

np <- 50
ni <- 20
nr <- 30
no_rater_groups <- 4
dat <- simulate_data_crm(np, ni, nr, no_rater_groups = no_rater_groups)

data_2_init(dat$df) |> add_true_values_to_em_tib(dat) |> scatterplot_retrieval() +
  ggtitle("parameter estimates based on data heuristics")


em_tib <- fit_em(dat$df, record_all_iterations = 5, n_iter = 10, compute_locations_scales = TRUE)
em_tib |> add_true_values_to_em_tib(dat) |> scatterplot_retrieval(facets = iteration~parameter)


mod <- cmdstanr::cmdstan_model("stanmodels/extendedCRM_long_new_no_genq.stan",
                               dir = "stanmodels",
                               pedantic = TRUE, quiet = FALSE)

stan_data <- crm_data_2_stan(dat, vectorized = FALSE)

vb_init_fit <- mod$variational(data = stan_data, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3,
                          init = list(em_tib_to_init(em_tib, dat)))
vb_init_tib <- get_means_tib(vb_init_fit, dat)
vb_init_tib |> scatterplot_retrieval() + ggtitle(sprintf("mod$variational informed by EM: %.3f seconds", vb_init_fit$time())) + labs(x = "Variational Estimate", y = "True value")


vb_init_no_adapt_fit <- mod$variational(data = stan_data, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3,
                          init = list(init_from_em), eta = 1, adapt_engaged = FALSE)
vb_init_no_adapt_tib <- get_means_tib(vb_init_no_adapt_fit, dat)
vb_init_no_adapt_tib |> scatterplot_retrieval() + ggtitle(sprintf("mod$variational informed by EM without adaptation: %.3f seconds", vb_init_no_adapt_fit$time())) + labs(x = "Variational Estimate", y = "True value")
