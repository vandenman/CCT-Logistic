rm(list = ls())
library(CCTLogistic)
library(assertthat)
# library(data.table)

library(readr)
library(tibble)
library(dplyr)
library(tidyr)

# data_wide <- fread(file.path("data", "IFBE_BASISSET_Bayes_Definitief-1.csv"))
#
# old_names <- c(
#   "Patientnummer_Random", "Meetmoment", "Leeftijd_M1_4groepen",
#   "Behandelduur_M1_4gr", "Geweld_Voor_dicho", "Geweld_Tussen_dicho",
#   "Geweld_Na_dicho", "Diagnose_categorie", "Delict_Categorie",
#   "PersoneelsNummer_Random", "Functiegroep"
# )
#
# # the actual items
# measure_vars <- paste0("IFBE_", 0:22)
#
# # recode missing from -99 to NA
# for (col in c(measure_vars, old_names))
#   set(data_wide, which(data_wide[[col]] %in% c(-99, 99)), col, NA)
#
# # we don't need these
# skipped_vars <- c("Aantal_Patienten")
#
# # rename
# new_names <- c(
#   "patient", "time", "patient_age_group",
#   "treatement_duration_group", "violent_before", "violent_between",
#   "violent_after", "diagnosis_group", "crime_group",
#   "rater", "rater_group"
# )
#
# assert_that(all(colnames(data_wide) %in% c(old_names, measure_vars, skipped_vars)))
#
# setnames(data_wide, old_names, new_names)
#
# # transform categorical to factor from character
# for (col in new_names)
#   set(data_wide, j = col, value = as.factor(data_wide[[col]]))
#
# # reshape wide to long
# data_long <- melt(data = data_wide, id.vars = new_names, measure.vars = measure_vars,
#                  variable.name = "item", value.name = "score")
#
# # reorder for convenience
# order <- c("time", "patient", "rater", "item", "score")
# setcolorder(data_long, c(order, setdiff(colnames(data_long), order)))
#
#
# # small summary of observed values
# data_wide[, summary(.SD), .SDcols = measure_vars] # for inspection
#
# # assert the wide and long data yield the same summary statistics
# sum1 <- lapply(data_wide[, measure_vars, with = FALSE], function(x) summary(x))
# sum2 <- tapply(data_long$score, data_long$item, summary, simplify = FALSE)
# assert_that(all.equal(sum1, sum2, check.attributes = FALSE))
#
# saveRDS(data_wide, file.path("data", "data_wide.rds"))
# saveRDS(data_long, file.path("data", "data_long.rds"))

# new approach
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

data_wide22 <- data_wide2 |>
  select(-skipped_vars) |>
  mutate(
    across(c(measure_vars, old_names), \(x) na_if(x,  99)),
    across(c(measure_vars, old_names), \(x) na_if(x, -99)),
  ) |>
  rename_with(~ new_names, old_names) |>
  mutate(
    across(measure_vars, \(x) x + 1), # transform score from 0-17 to 1 - 18
    across(new_names, factor)
  )

data_long22 <- data_wide22 |>
  pivot_longer(measure_vars, names_to = "item", values_to = "score") |>
  relocate(patient, item, rater, rater_group, time, score) |>
  mutate(
    item = factor(as.numeric(gsub("IFBE_", "", item)) + 1)
  )

data_violence <- data_long22 |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(c(patient, age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  filter(!duplicated(patient))

data_ifte <- data_long22 |>
  filter(!is.na(score) & !is.na(violent_before) & !is.na(diagnosis) & !is.na(crime)) |>
  select(-c(age, violent_before, violent_between, violent_after, treatment_duration, diagnosis, crime)) |>
  arrange(rater_group, patient, item, rater, time)

saveRDS(data_wide22,   file.path("data", "data_wide.rds"))
saveRDS(data_long22,   file.path("data", "data_long.rds"))
saveRDS(data_ifte,     file.path("data", "data_ifte.rds"))
saveRDS(data_violence, file.path("data", "data_violence.rds"))
