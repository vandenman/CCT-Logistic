rm(list = ls())
library(CCTLogistic)
library(assertthat)

library(readr)
library(tibble)
library(dplyr)
library(tidyr)

data_wide2 <- read_csv(file.path("data", "IFBE_BASISSET APRIL2019-Bayes_Definitief-2022_edited.csv"))

# Dutch names
old_names <- c(
  "Patientnummer_Random", "Meetmoment", "Leeftijd_M1_4groepen",
  "Behandelduur_M1_4gr", "Geweld_Voor_dicho", "Geweld_Tussen_dicho",
  "Geweld_Na_dicho", "Diagnose_categorie", "Delict_Categorie",
  "PersoneelsNummer_Random", "Functiegroep"
)

# English names
new_names <- c(
  "patient", "time", "age",
  "treatment_duration", "violent_before", "violent_between",
  "violent_after", "diagnosis", "crime",
  "rater", "rater_group"
)

# the actual items
measure_vars <- paste0("IFBE_", 0:22)

# we don't need these
skipped_vars <- c("Aantal_Patienten", "Extra_Kolom")

assert_that(all(colnames(data_wide2) %in% c(old_names, measure_vars, skipped_vars)))

na_if_vec <- function(x, y) {
  idx <- which(x %in% y)
  x[idx] <- NA
  x
}

data_wide22 <- data_wide2 |>
  select(-skipped_vars) |>
  mutate(
    across(old_names,    \(x) na_if_vec(x, c(99, -99))),
    across(measure_vars, \(x) na_if_vec(x, c(99, -99, 0))),
  ) |>
  rename_with(~ new_names, old_names) |>
  mutate(
    across(new_names, factor)
  )

data_long22 <- data_wide22 |>
  pivot_longer(measure_vars, names_to = "item", values_to = "score") |>
  relocate(patient, item, rater, rater_group, time, score) |>
  mutate(
    item = factor(as.numeric(gsub("IFBE_", "", item)) + 1)
  )

data_violence <- data_long22 |>
  filter(
    # !is.na(score) &
    !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)
  ) |>
  select(c(patient, age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  filter(!duplicated(patient))

data_ifte <- data_long22 |>
  filter(
    # !is.na(score) &
    !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)
  ) |>
  select(-c(age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  arrange(rater_group, patient, item, rater, time)

saveRDS(data_wide22,   file.path("data", "data_wide.rds"))
saveRDS(data_long22,   file.path("data", "data_long.rds"))
saveRDS(data_ifte,     file.path("data", "data_ifte.rds"))
saveRDS(data_violence, file.path("data", "data_violence.rds"))
