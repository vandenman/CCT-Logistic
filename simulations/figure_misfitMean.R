rm(list = ls())
library(tibble)
library(ggplot2)
source(file.path("simulations", "utils.R"))
source(file.path("simulations", "CRM_helpers.R"))

np <- 1
ni <- 23
nr <- 10

set.seed(42)
dat <- simulate_data(np, ni, nr)

figure <- tibble(
  true_value = c(dat$parameters$lt),
  sum_score  = c(t(apply(dat$x, 1:2, mean))),
  patient    = factor(rep(seq_len(np), each = ni))
) |> ggplot(mapping = aes(x = true_value, y = sum_score)) +
  jaspGraphs::geom_abline2(color = "grey20") +
  jaspGraphs::geom_point() +
  jaspGraphs::scale_x_continuous(name = "True item score") +
  jaspGraphs::scale_y_continuous(name = "Sample mean") +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw()

save_figure(figure, "misfitMean.svg")
