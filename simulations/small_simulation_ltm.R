rm(list = ls())
library(CCTLogistic)
library(tibble)
library(ggplot2)
library(cmdstanr)

np <- 50             # no patients
ni <- 20             # no items
nr <- 30             # no raters
nc <- 5              # no response categories
no_rater_groups <- 1 # no groups of rater

set.seed(1234)
# debugonce(simulate_data_ltm)
dat <- simulate_data_ltm(np, ni, nr, nc, no_rater_groups = no_rater_groups)

# debugonce(log_likelihood_ltm)
log_likelihood_ltm(dat,
                   log_a      = dat$parameters$log_a,
                   b          = dat$parameters$b,
                   log_E      = dat$parameters$log_E,
                   lt         = dat$parameters$lt,
                   log_lambda = dat$parameters$log_lambda)

mod <- compile_stan_model("stanmodels/extendedLTM_long_new_no_genq.stan", pedantic = TRUE, quiet = FALSE)


stan_data <- ltm_data_2_stan(dat)
fit_vb <- mod$variational(data = stan_data, iter = 3e4, grad_samples = 3, elbo_samples = 2, adapt_iter = 500, output_samples = 5e3)
fit_vb <- mod$variational(data = stan_data, iter = 3e4, adapt_iter = 500, output_samples = 5e3)
tib_vb <- get_means_tib(fit_vb) |> add_true_values_to_tib(dat)

tib_vb |>
  scatterplot_retrieval() +
  ggtitle(sprintf("mod$variational: %.3f seconds", fit_vb$time())) +
  labs(x = "Variational Estimate", y = "True value")

system("beep_finished.sh")

lst <- fit_vb$draws(format = "draws_list")
plot(lst$`1`$mu_lt)

plot(lst$`1`$)

mean(dat$parameters$b) - mean(dat$parameters$lt)
diff(sum_vb[sum_vb$variable %in% c("mu_lt", "mu_b[1]"), ][["mean"]])

lt_lst <- fit_vb$draws(variables = "lt", format = "draws_matrix")
b_lst  <- fit_vb$draws(variables = "b",  format = "draws_matrix")
cbind(rowMeans(lt_lst),    lst$`1`$mu_lt)
cbind(rowMeans(b_lst),     lst$`1`$`mu_b[1]`)
cbind(apply(b_lst, 1, sd), lst$`1`$`sd_b[1]`)

sum_vb <- fit_vb$summary()

sum_vb[sum_vb$variable %in% c("mu_lt", "mu_b"), ]

str(sum_vb)
sum_vb$variable, 200)
