# alternative rstan:::rstudio_stanc()
rstudio_stanc_cmdstanr <- function(filename) {
  model <- cmdstanr::cmdstan_model(stan_file = filename, compile = FALSE)
  model$check_syntax(stanc_options = list("allow-undefined"), pedantic = FALSE)
}
