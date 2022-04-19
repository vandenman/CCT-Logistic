rm(list = ls())
library(CCTLogistic)
library(dplyr)
library(ggplot2)


all_data <- read_long_data()
data_2_analyze <- all_data |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(-c(age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  arrange(rater_group, patient, item, rater, time)

data_violence <- all_data |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(c(patient, age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  filter(!duplicated(patient))

# data_2_analyze <- all_data |>
#   filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
#   select(-c(patient_age_group, violent_before, violent_between, violent_after,
#                  treatement_duration_group, diagnosis, crime)) |>
  # mutate(
  #   across(c(patient, item, rater, time, rater_group), \(x) normalize_factor(factor(x)))
  #   score       = score + 1 # transform score from 0 -17 to 1 - 18
  # ) |>
  # arrange(rater_group, patient, item, rater, time)

# data_violence <- all_data |>
#   filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
#   select(c(patient, patient_age_group, violent_before, violent_between, violent_after, treatement_duration_group, diagnosis_group, crime_group)) |>
  # mutate(
  #   patient = normalize_factor(factor(patient)),
  #   violent_after = as.integer(violent_after) - 1L
  # ) |>
  # filter(!duplicated(patient))

h <- hist(data_2_analyze$score, breaks = 0:18, plot = FALSE)
all(h$counts == table(data_2_analyze$score))

tib_scores <- data_2_analyze |>
  summarise(
    counts        = hist(x = score, plot = FALSE, breaks = 0:18)$counts,
    unique_scores = factor(1:18),
  )


yBreaks <- jaspGraphs::getPrettyAxisBreaks(tib_scores$counts)
plot_scores <-
  ggplot(data = tib_scores, aes(x = unique_scores, y = counts)) +
  # jaspGraphs::geom_abline2(data = data.frame(i = seq(500, 3500, 500), s = rep(0, 7)), mapping = aes(intercept = i, slope = s), color = "grey80") +
  geom_col(fill = "grey", col = "black") +
  scale_y_continuous(name = "Count", breaks = yBreaks, limits = range(yBreaks)) +
  jaspGraphs::geom_rangeframe(sides = "bl") +
  labs(x = "Score") +
  jaspGraphs::themeJaspRaw() +
  theme(panel.grid.major.y = element_line(colour = "grey80"))
plot_scores

tib_scores_by_items <- data_2_analyze |>
  group_by(item, time) |>
  summarise(
    counts        = hist(x = score, plot = FALSE, breaks = 0:18)$counts,
    unique_scores = factor(1:18),
    .groups = "keep"
  ) |>
  ungroup() |>
  mutate(
    item = factor(item, levels = 23:1),
    time = factor(time)
  )


yLimits <- range(jaspGraphs::getPrettyAxisBreaks(tib_scores_by_items$counts))
plot_scores_by_items <-
  ggplot(data = tib_scores_by_items, aes(x = unique_scores, y = counts, group = time, fill = time)) +
  geom_col(position = position_dodge()) +
  jaspGraphs::geom_rangeframe(sides = "bl") +
  # scale_y_continuous(name = "Item", breaks = NULL, limits = yLimits) +
  # scale_y_continuous(name = NULL, breaks = yLimits, limits = yLimits, position = "right",
  #                    sec.axis = dup_axis(name = "Item", breaks = NULL, labels = NULL)) +
  scale_y_continuous(name = "Items", breaks = yLimits, limits = yLimits) +
  jaspGraphs::scale_JASPfill_discrete() +
  facet_grid(rows = vars(item), switch = "y") +
  labs(fill = "Time", x = "Score") +
  jaspGraphs::themeJaspRaw() +
  theme(
    panel.spacing = 2*theme_grey()$panel.spacing,
    strip.text.y.left = element_text(angle = 0),
    legend.position = "right",
    # axis.text.y.right = element_blank(),
    # axis.ticks.y.right = element_blank()
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank()
  )
plot_scores_by_items

tib_raters_by_patients <- data_2_analyze |>
  group_by(patient, rater) |>
  summarise(
    count = factor(table(item)[1]),
    .groups = "keep"
  )# |>
  # ungroup() |>
  # mutate(
  #   patient = factor(patient),
  #   rater   = factor(rater)
  # )

range_fact <- function(x) range(as.numeric(x))
colors <- c("0" = "#FFFFFF", setNames(jaspGraphs::JASPcolors(palette = "colorblind", asFunction = TRUE)(2), 1:2))
grey <- "grey95"
plot_raters_by_patients <- ggplot(tib_raters_by_patients, aes(x = as.integer(rater), y = as.integer(patient), color = count, fill = factor(count))) +
  geom_raster() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  # jaspGraphs::scale_JASPfill_discrete() +
  # scale_y_continuous(breaks = range(tib_raters_by_patients$patient), labels = range(tib_raters_by_patients$patient), minor_breaks = seq_along(tib_raters_by_patients$patient)-0.5, expand = expansion(), limits = range(tib_raters_by_patients$patient) + c(-0.9, 0.9)) +
  # scale_x_continuous(breaks = range(tib_raters_by_patients$rater),   labels = range(tib_raters_by_patients$rater),   minor_breaks = seq_along(tib_raters_by_patients$rater)  -0.5, expand = expansion(), limits = range(tib_raters_by_patients$rater) + c(-0.9, 0.9)) +
  scale_y_continuous(breaks = range_fact(tib_raters_by_patients$patient), labels = range_fact(tib_raters_by_patients$patient), minor_breaks = seq_along(tib_raters_by_patients$patient)-0.5, expand = expansion(), limits = range_fact(tib_raters_by_patients$patient) + c(-0.5, 0.5)) +
  scale_x_continuous(breaks = range_fact(tib_raters_by_patients$rater),   labels = range_fact(tib_raters_by_patients$rater),   minor_breaks = seq_along(tib_raters_by_patients$rater)  -0.5, expand = expansion(), limits = range_fact(tib_raters_by_patients$rater) + c(-0.5, 0.5)) +
  labs(x = "Rater", y = "Patient", fill = "Frequency") +
  coord_fixed(expand = FALSE) +
  # jaspGraphs::geom_rangeframe() +
  jaspGraphs::themeJaspRaw() +
  theme(
    legend.position = "right",
    panel.border = element_blank(),
    # panel.border = element_rect(colour = "grey80", fill=NA, size=0.5),
    # panel.border = element_rect(colour = "grey80", fill=NA, size=1),
    panel.grid.minor = element_line(colour = grey, size = 0.5),
    # axis.title.x = element_text(margin = margin(t = -20), debug = FALSE),
    axis.title.y = element_text(margin = margin(r = -20), debug = FALSE),
    plot.margin = margin(l = 10),
    legend.key = element_rect(colour = grey, size = 2),
    axis.ticks.length = unit(1, "mm")
  )
plot_raters_by_patients
# TODO: add histograms on top and right?
temp <- table(all_data$patient, all_data$rater) / 23
colSums(temp)
rowSums(temp)


combined_plt <- patchwork::wrap_plots(
  plot_scores             + ggtitle("A") + theme(plot.title = element_text(hjust = 0)),
  plot_raters_by_patients + ggtitle("B") + theme(plot.title = element_text(hjust = 0))
)
combined_plt
save_figure_obj(combined_plt, file = "descriptives_combined_plt.rds")
save_figure(figure = combined_plt, file = "descriptives_combined_plt.svg", width = 15, height = 12)
save_figure(figure = plot_scores,             file = "descriptives_scores.svg",             width = 12, height = 7)
save_figure(figure = plot_raters_by_patients, file = "descriptives_raters_by_patients.svg", width = 900 / 96, height = 494 / 96)

set.seed(42)
patient_idx <- sample(data_2_analyze$patient, 3)
tibbie <- data_2_analyze |> filter(patient %in% patient_idx) |> #, item %in% 1:23) |>
  group_by(patient, item, time) |>
  summarise(
    mean_score = mean(score),
    .groups = "keep"
  )# |>
  #mutate(sides = if_else(patient == 1, "bl", "l"))

ggplot(tibbie, aes(x = time, y = mean_score, group = item, color = factor(item), shape = factor(item))) +
  jaspGraphs::scale_JASPcolor_discrete(name = "Item") +
  # scale_color_manual(name = "Item", values = jaspGraphs::JASPcolors(palette = "colorblind", asFunction = TRUE)(23)) +
  scale_shape_manual(name = "Item", values = rep(c(0, 1, 2, 4, 5), length.out = 23)) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(name = "Mean score", breaks = c(1, 5, 10, 15, 18), limits = c(0.8, 18.2), expand = expansion()) +
  xlab("Time") +
  jaspGraphs::geom_rangeframe() +
  facet_grid(cols = vars(patient)) +
  jaspGraphs::themeJaspRaw() +
  theme(
    legend.position = "right",
    strip.text.x = element_blank()
  )
# TODO: in the figure, the y-axis line should not repeat!


table(data_violence$violent_after, data_violence$patient_age_group)
table(data_violence$violent_after, data_violence$diagnosis_group)
table(data_violence$violent_after, data_violence$crime_group)
table(data_violence$violent_after, data_violence$violent_between)
table(data_violence$violent_after, data_violence$violent_before)

vars <- c("patient_age_group", "diagnosis_group", "crime_group", "violent_between", "violent_before", "violent_after")
ddd <- data_violence[
  complete.cases(data_violence[, vars]), vars
]
anyNA(ddd)

table(ddd$violent_after)
ggg <- glm(violent_after ~ 1 + patient_age_group + diagnosis_group + crime_group + violent_between + violent_before, data = ddd, family = binomial())
summary(ggg)


?write.csv(ddd, file = "data/violence_data_test.csv")

trainingFit <- rpart::rpart(
  violent_after ~ 1 + patient_age_group + diagnosis_group + crime_group + violent_between + violent_before, data = ddd,
  method = "class", x = TRUE, y = TRUE
)
plot(trainingFit)

cat(c("\U262E", "\U1F92C"))
ttb <- table(data_violence$violent_after, data_violence$patient_age_group)
rownames(ttb) <- c("\U262E", "\U1F92C")


renv::install("partykit")
renv::install("ggparty")

plotData <- partykit::as.party(trainingFit)
p <- ggparty::ggparty(plotData)
x <- trainingFit
frame <- x$frame
ylevel <- attr(x, "ylevels")
digits <- 3
tfun <- (x$functions)$print
if (!is.null(tfun)) {
  if (is.null(frame$yval2)) {
    yval <- tfun(frame$yval, ylevel, digits, nsmall = 20)
  } else {
    yval <- tfun(frame$yval2, ylevel, digits, nsmall = 20)
  }
} else {
  yval <- format(signif(frame$yval, digits))
}
leafs <- which(x$frame$var == "<leaf>")
labels <- yval[leafs]
# if (purpose == "classification") {
  labels <- strsplit(labels, split = " ")
  labels <- unlist(lapply(labels, `[[`, 1))
# }
nodeNames <- p$data$splitvar
nodeNames[is.na(nodeNames)] <- labels
p$data$info <- paste0(nodeNames, "\nn = ", p$data$nodesize)
# p <-
  p + ggparty::geom_edge() +
  ggparty::geom_edge_label(mapping = ggplot2::aes(label = paste(substr(breaks_label, start = 1, stop = 15))), fill = NA, parse = FALSE,
                           nudge_y = 0.02) +
  ggparty::geom_node_splitvar(mapping = ggplot2::aes(size = max(3, nodesize) / 2, label = info), fill = "white", col = "black") +
  ggparty::geom_node_label(mapping = ggplot2::aes(label = info, size = max(3, nodesize) / 2), ids = "terminal", fill = "white", col = "black") +
  ggplot2::scale_x_continuous(name = NULL, limits = c(min(p$data$x) - abs(0.1 * min(p$data$x)), max(p$data$x) * 1.1)) +
  ggplot2::scale_y_continuous(name = NULL, limits = c(min(p$data$y) - abs(0.1 * min(p$data$y)), max(p$data$y) * 1.1)) +
  jaspGraphs::geom_rangeframe(sides = "") +
  jaspGraphs::themeJaspRaw() +
  ggplot2::theme(
    axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )
p
