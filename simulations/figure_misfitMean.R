rm(list = ls())
library(tibble)
library(dplyr)
library(ggplot2)
library(CCTLogistic)
library(patchwork)

create_figure_tib <- function(data_obj) {
  observed_data <- data_obj$df |>
    filter(patient %in% patient_idx) |>
    group_by(patient, item) |>
    summarise(mean = mean(score)) |>
    ungroup(item) |>
    summarise(mean2 = mean - mean(mean))
  latent_truth  <- data_obj$parameters$lt[patient_idx, ]
  patient_grp   <-  factor(rep(patient_idx, each = ni))

  tibble(
    true_value = c(t(latent_truth)),
    mean_score  = observed_data$mean2,
    patient    = patient_grp
  )
}

create_figure2 <- function(tib, x = c(0, 0), y = c(0, 0)) {

  cors <- tib_joined |> group_by(type) |> summarise(cor = cor(true_value, mean_score))
  tib_text <- tibble(
    x = x,
    y = y,
    label = sprintf("ρ = %.3f", cors$cor),
    type = cors$type
  )

  ggplot(tib, mapping = aes(x = true_value, y = mean_score, color = patient, fill = patient, shape = patient)) +
    jaspGraphs::geom_abline2(color = "grey20") +
    jaspGraphs::geom_point(size = 5, alpha = 0.8) +
    jaspGraphs::scale_x_continuous(name = "True item score") +
    jaspGraphs::scale_y_continuous(name = "Standardised item mean") +
    jaspGraphs::scale_JASPcolor_discrete() +
    geom_text(data = tib_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, show.legend = FALSE, size = .35*jaspGraphs::graphOptions("fontsize")) +
    labs(color = "Patient", fill = "Patient", shape = "Patient") +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = c(0.07, .95)) +
    facet_wrap(~type)
}


create_figure <- function(data_obj, x = 0, y = 0) {
  observed_data <- data_obj$df |>
    filter(patient %in% patient_idx) |>
    group_by(patient, item) |>
    summarise(mean = mean(score)) |>
    ungroup(item) |>
    summarise(mean2 = mean - mean(mean))
  latent_truth  <- data_obj$parameters$lt[patient_idx, ]
  patient_grp   <-  factor(rep(patient_idx, each = ni))

  tib <- tibble(
    true_value = c(t(latent_truth)),
    mean_score  = observed_data$mean2,
    patient    = patient_grp
  )

  # x_range <- range(jaspGraphs::getPrettyAxisBreaks(tib$true_value))
  # y_range <- range(jaspGraphs::getPrettyAxisBreaks(tib$mean_score))
  tib_text <- tibble(
    x = x,
    y = y,
    label = sprintf("Cor = %.3f", cor(tib$true_value, tib$mean_score))
  )

  ggplot(tib, mapping = aes(x = true_value, y = mean_score, color = patient, fill = patient, shape = patient)) +
    jaspGraphs::geom_abline2(color = "grey20") +
    jaspGraphs::geom_point(size = 5, alpha = 0.8) +
    jaspGraphs::scale_x_continuous(name = "True item score") +
    jaspGraphs::scale_y_continuous(name = "Standardised item mean") +
    jaspGraphs::scale_JASPcolor_discrete() +
    geom_text(data = tib_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, show.legend = FALSE, size = .35*jaspGraphs::graphOptions("fontsize")) +
    labs(color = "Patient", fill = "Patient", shape = "Patient") +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = c(0.07, .95))# +
    # ggtitle(sprintf("Cor = %.3f", cor(tib$true_value, tib$mean_score))) +
    # theme(plot.title = element_text(hjust = 0.1))
}

np <- 100
ni <- 10
nr <- 10
nc <- 7

patient_idx <- 1:3

# TODO: side by side of possibly under model vs ideal scenario for mean

set.seed(42)
thresholds <- t(apply((matrix(abs(rnorm(n = (nc - 1) * nr, 2, 2)), nr, nc - 1)), 1, sort))
log_E <- log(abs(rnorm(nr)))
data_obj_bad  <- simulate_data_ltm(np, ni, nr, nc, threshold_type = "free", thresholds = thresholds, log_E = log_E)
set.seed(42)
data_obj_good  <- simulate_data_ltm(np, ni, nr, nc, threshold_type = "free")

tib_joined <- rbind(
  create_figure_tib(data_obj_bad)  |> mutate(type = "heterogenous"),
  create_figure_tib(data_obj_good) |> mutate(type = "homogenous")
)
figure <- create_figure2(tib_joined, c(1.0, 1.0), c(2.5, 2.5)) + theme(strip.text = element_blank())

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


