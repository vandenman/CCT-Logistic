rm(list = ls())
library(CCTLogistic)
library(dplyr)
library(ggplot2)


all_data <- read_long_data()
data_2_analyze <- all_data |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
  select(-c(patient_age_group, violent_before, violent_between, violent_after,
                 treatement_duration_group, diagnosis_group, crime_group)) |>
  mutate(
    across(c(patient, item, rater, time, rater_group), \(x) normalize_factor(factor(x))),
    # patient     = normalize_factor(patient),
    # item        = normalize_factor(item),
    # rater       = normalize_factor(rater),
    # time        = normalize_factor(time),
    # rater_group = normalize_factor(rater_group),
    score       = score + 1 # transform score from 0 -17 to 1 - 18
  ) |>
  arrange(rater_group, patient, item, rater, time)

data_violence <- all_data |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis_group) & !is.na(crime_group)) |>
  select(c(patient, patient_age_group, violent_before, violent_between, violent_after, treatement_duration_group, diagnosis_group, crime_group)) |>
  mutate(
    patient = normalize_factor(factor(patient)),
    violent_after = as.integer(violent_after) - 1L
  ) |>
  filter(!duplicated(patient))

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


ttt <- table(all_data$item, all_data$rater, all_data$patient)
dim(ttt)
idx <- which(ttt[1, ,] == max(ttt[1, ,]), arr.ind = TRUE, useNames = TRUE)
levels(all_data$rater)[idx[1]]
levels(all_data$patient)[idx[2]]

# 3 observations!!!
all_data |>
  filter(patient == "1640" & rater == "818") |>
  nrow() / 23

table(fff$item, droplevels(fff$rater))
table(fff$item, droplevels(fff$rater), fff$time)

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

colors <- c("0" = "#FFFFFF", setNames(jaspGraphs::JASPcolors(palette = "colorblind", asFunction = TRUE)(2), 1:2))
grey <- "grey95"
plot_raters_by_patients <- ggplot(tib_raters_by_patients, aes(x = rater, y = patient, color = count, fill = factor(count))) +
  geom_raster() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  # jaspGraphs::scale_JASPfill_discrete() +
  # scale_y_continuous(breaks = range(tib_raters_by_patients$patient), labels = range(tib_raters_by_patients$patient), minor_breaks = seq_along(tib_raters_by_patients$patient)-0.5, expand = expansion(), limits = range(tib_raters_by_patients$patient) + c(-0.9, 0.9)) +
  # scale_x_continuous(breaks = range(tib_raters_by_patients$rater),   labels = range(tib_raters_by_patients$rater),   minor_breaks = seq_along(tib_raters_by_patients$rater)  -0.5, expand = expansion(), limits = range(tib_raters_by_patients$rater) + c(-0.9, 0.9)) +
  scale_y_continuous(breaks = range(tib_raters_by_patients$patient), labels = range(tib_raters_by_patients$patient), minor_breaks = seq_along(tib_raters_by_patients$patient)-0.5, expand = expansion(), limits = range(tib_raters_by_patients$patient) + c(-0.5, 0.5)) +
  scale_x_continuous(breaks = range(tib_raters_by_patients$rater),   labels = range(tib_raters_by_patients$rater),   minor_breaks = seq_along(tib_raters_by_patients$rater)  -0.5, expand = expansion(), limits = range(tib_raters_by_patients$rater) + c(-0.5, 0.5)) +
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

tibbie <- data_2_analyze |> filter(patient %in% 1:3, item %in% 1:23) |>
  group_by(patient, item, time) |>
  summarise(
    mean_score = mean(score),
    .groups = "keep"
  ) |>
  mutate(sides = if_else(patient == 1, "bl", "l"))

ggplot(tibbie, aes(x = factor(time), y = mean_score, group = item, color = factor(item), shape = factor(item))) +
  jaspGraphs::scale_JASPcolor_discrete(name = "Item") +
  # scale_color_manual(name = "Item", values = jaspGraphs::JASPcolors(palette = "colorblind", asFunction = TRUE)(23)) +
  scale_shape_manual(name = "Item", values = rep(c(0, 1, 2, 4, 5), length.out = 23)) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(name = "Mean score", breaks = c(1, 5, 10, 15, 18), limits = c(0.8, 18.2), expand = expansion()) +
  jaspGraphs::geom_rangeframe() +
  facet_grid(cols = vars(patient)) +
  jaspGraphs::themeJaspRaw() +
  theme(
    legend.position = "right",
    strip.text.x = element_blank()
  )
# TODO: in the figure, the y-axis line should not repeat!


# plot_scores <-
#   ggplot(data = tib_scores_by_items, aes(x = factor(unique_scores), y = counts, group = factor(time), fill = factor(time))) +
#   geom_col(position = position_dodge()) +
#   jaspGraphs::geom_rangeframe(sides = "bl") +
#   scale_x_discrete() +
#   # scale_y_continuous(name = "Item", breaks = NULL, limits = yLimits) +
#   # scale_y_continuous(name = NULL, breaks = yLimits, limits = yLimits, position = "right",
#   #                    sec.axis = dup_axis(name = "Item", breaks = NULL, labels = NULL)) +
#   scale_y_continuous(name = "Items", breaks = yLimits, limits = yLimits) +
#   jaspGraphs::scale_JASPfill_discrete() +
#   facet_grid(rows = vars(item), switch = "y") +
#   labs(fill = "Time", x = "Score") +
#   jaspGraphs::themeJaspRaw() +
#   theme(
#     panel.spacing = 2*theme_grey()$panel.spacing,
#     strip.text.y.left = element_text(angle = 0),
#     legend.position = "right",
#     # axis.text.y.right = element_blank(),
#     # axis.ticks.y.right = element_blank()
#     axis.text.y.left = element_blank(),
#     axis.ticks.y.left = element_blank()
#   )
#
# table(data_2_analyze$patient)
# table(data_2_analyze$patient, data_2_analyze$time)
