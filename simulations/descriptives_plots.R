rm(list = ls())
library(data.table)
library(ggplot2)
source(file.path("simulations", "utils.R"))

data_long <- read_long_data()

no_patients <- 5
no_items    <- 5

unique_patients_sorted <- names(sort(table(data_long[["patient"]])))

# selected_patients <- rev(unique_patients_sorted)[1:no_patients]
selected_patients <- unique_patients_sorted[1:no_patients]
selected_items    <- unique(data_long[["item"]])[1:no_items]

subset_data_long <- data_long[patient %in% selected_patients & item %in% selected_items, ]

ggplot(data = subset_data_long, aes(x = factor(time), y = score, group = factor(rater), color = factor(rater))) +
  geom_point() +
  geom_line() +
  facet_grid(factor(patient) ~ item) +
  theme_bw() +
  theme(legend.position = "none")

subset_data_long[, lapply(.SD, mean),
                 .SDcols = c("violent_before", "violent_between", "violent_after"),
                 by = c("diagnosis_group", "crime_group")]

unique(data_long$diagnosis_group)
unique(data_long$crime_group)


ggplot(data = subset_data_long, aes(x = factor(time), y = score, group = factor(rater), color = factor(rater))) +
  geom_point() +
  geom_line() +
  facet_grid(factor(patient) ~ item) +
  theme_bw() +
  theme(legend.position = "none")


table(data_long[, c("violent_before", "violent_between", "violent_after")])

tb_ifbe <- table(data_long$item, data_long$score)
tb_ifbe_long <- data.table(
  item  = factor(rownames(tb_ifbe)),
  count = c(tb_ifbe),
  value = as.integer(rep(colnames(tb_ifbe), each = nrow(tb_ifbe)))
)

ggplot(data = tb_ifbe_long, aes(x = value, y = count)) +
  geom_col() +
  facet_wrap(~item)

names(data_long)
ggplot(data = data_long, aes(x = score, y = ..count.., group = rater_group, fill = rater_group)) +
  geom_histogram(position = position_dodge()) +
  scale_x_continuous(breaks = 0:17) +
  facet_wrap(~item, scales = "free_y")

# voor de inrichting, kans op aggresiviteit? delict categorie?
# zijn de ratings van verschillende groepen onderscheidend?

