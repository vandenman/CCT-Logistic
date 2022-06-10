rm(list = ls())
library(tibble)
library(dplyr)
library(ggplot2)
library(CCTLogistic)

np <- 100
ni <- 10
nr <- 10
nc <- 7

patient_idx <- 1:3

set.seed(42)
# thresholds <- t(apply((matrix(rnorm(n = (nc - 1) * nr, 2, 2), nr, nc - 1)), 1, sort))
data_obj <- simulate_data_ltm(np, ni, nr, nc, threshold_type = "free")#, thresholds = thresholds)
observed_data <- data_obj$df |>
  filter(patient %in% patient_idx) |>
  group_by(patient, item) |>
  summarise(mean = mean(score)) |>
  ungroup(item) |>
  summarise(mean2 = mean - mean(mean))
latent_truth  <- data_obj$parameters$lt[patient_idx, ]
patient_grp   <-  factor(rep(patient_idx, each = ni))

figure <- tibble(
  true_value = c(t(latent_truth)),
  mean_score  = observed_data$mean2,
  patient    = patient_grp
) |> ggplot(mapping = aes(x = true_value, y = mean_score, color = patient, fill = patient, shape = patient)) +
  jaspGraphs::geom_abline2(color = "grey20") +
  jaspGraphs::geom_point(size = 5, alpha = 0.8) +
  jaspGraphs::scale_x_continuous(name = "True item score") +
  jaspGraphs::scale_y_continuous(name = "Sample mean") +
  jaspGraphs::scale_JASPcolor_discrete() +
  labs(color = "Patient", fill = "Patient", shape = "Patient") +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw(legend.position = c(0.07, .95), )

figure

save_figure(figure, "misfitMean.svg")

ltm_fit <- CCTLogistic::save_or_run_model({
  ltm_model <- compile_stan_model("stanmodels/LTM_3_models_with_logistic_regression_with_time.stan", pedantic = TRUE, quiet = FALSE, include_paths = "stanmodels", cpp_options = list(stan_threads=TRUE))
  ltm_data  <- data_2_stan(data_obj)
  ltm_model$variational(data = ltm_data, iter = 3e4, adapt_iter = 500, output_samples = 2e3, grad_samples = 5, elbo_samples = 5, threads = 8)
}, force = FALSE, path = "fitted_objects/figure_misfitMean2.rds")

lt_samples <- ltm_fit$draws("lt", format = "draws_matrix")
lt_post_means <- matrix(colMeans(lt_samples), nrow = np, ncol = ni)

additional_figure_data <- tibble(
  true_value = c(t(latent_truth), t(latent_truth)),
  mean_score = c(observed_data$mean2, t(lt_post_means[patient_idx, ])),
  patient    = c(patient_grp, patient_grp),
  what       = rep(c("means", "model"), each = length(patient_idx) * ni)
)

additional_figure <- ggplot(additional_figure_data, mapping = aes(x = true_value, y = mean_score, color = patient, fill = patient, shape = what)) +
  jaspGraphs::geom_abline2(color = "grey20") +
  jaspGraphs::geom_point(size = 5, alpha = 0.8) +
  jaspGraphs::scale_x_continuous(name = "True item score") +
  jaspGraphs::scale_y_continuous(name = "Sample mean") +
  jaspGraphs::scale_JASPcolor_discrete() +
  jaspGraphs::scale_JASPfill_discrete() +
  labs(color = "Patient", fill = "Patient", shape = "Method") +
  jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw(legend.position = "right")
additional_figure

rmse <- function(x, y) sqrt(mean((x - y)^2))
additional_figure_data |>
  group_by(what) |>
  summarise(rmse = rmse(true_value, mean_score))
# # A tibble: 2 × 2
#   what   rmse
#   <chr> <dbl>
# 1 means 0.672
# 2 model 0.333


