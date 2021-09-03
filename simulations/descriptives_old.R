rm(list = ls())
library(assertthat)
library(data.table)

dat_wide <- fread(file.path("data", "IFBE_BASISSET_Bayes_Definitief-1.csv"))

old_names <- c(
  "Patientnummer_Random", "Meetmoment", "Leeftijd_M1_4groepen",
  "Behandelduur_M1_4gr", "Geweld_Voor_dicho", "Geweld_Tussen_dicho",
  "Geweld_Na_dicho", "Diagnose_categorie", "Delict_Categorie",
  "PersoneelsNummer_Random", "Functiegroep"
)

# the actual items
measure_vars <- paste0("IFBE_", 0:22)

# recode missing from -99 to NA
for (col in measure_vars)
  set(dat_wide, which(dat_wide[[col]] == -99), col, NA)

# we don't need these
skipped_vars <- c("Aantal_Patienten")

# rename
new_names <- c(
  "patient", "time", "patient_age_group",
  "treatement_duration_group", "violent_before", "violent_between",
  "violent_after", "diagnosis_group", "crime_group",
  "rater", "rater_group"
)

assert_that(all(colnames(dat_wide) %in% c(old_names, measure_vars, skipped_vars)))

setnames(dat_wide, old_names, new_names)

# reshape wide to long
dat_long <- melt(data = dat_wide, id.vars = new_names, measure.vars = measure_vars,
                 variable.name = "item", value.name = "score")

# reorder for convenience
order <- c("time", "patient", "rater", "item", "score")
setcolorder(dat_long, c(order, setdiff(colnames(dat_long), order)))


# small summary of observed values
dat_wide[, summary(.SD), .SDcols = measure_vars] # for inspection

# assert the wide and long data yield the same summary statistics
sum1 <- lapply(dat_wide[, measure_vars, with = FALSE], function(x) summary(x))
sum2 <- tapply(dat_long$score, dat_long$item, summary, simplify = FALSE)
assert_that(all.equal(sum1, sum2, check.attributes = FALSE))

saveRDS(dat_wide, file.path(data, "data_wide.rds"))
saveRDS(dat_long, file.path(data, "data_long.rds"))

# library(readr)

# library(ggplot2)
# library(GGally)
# library(psych)
#
# dd <- read_csv2(file.path("data", "IFBE_BASISSET_Bayes_Definitief-1.csv"))
#
# cnms <- colnames(dd)
# idx0 <- startsWith(colnames(dd), "IFBE")
# idx <- cnms[idx0]
#
# hist(table(dd$PersoneelsNummer_Random))
#
# tb <- table(dd$Patientnummer_Random)
# hist(tb)
# idxMax <- names(tb)[which.max(tb)]
# idxMin <- names(tb)[which.min(tb)]
#
# matplot(dd[dd$Patientnummer_Random == idxMax, idx], type = 'l', bty = 'n', las = 1, xlim = c(1, 25))
# legend("topright", legend = seq_along(idx), lty = 1:5, col = 1:6, ncol = 3)
#
# colnames(dd)
# psych::alpha(dd[, idx], check.keys = TRUE)
#
# apply(dd[, idx], 2L, function(x) {
#   tapply(x, dd[, c("Patientnummer_Random", "PersoneelsNummer_Random")], var)
# })
#
#
# dd$Patientnummer_Random
#
# id.vars <- c("Patientnummer_Random", "Meetmoment", "Leeftijd_M1_4groepen",
#              "Behandelduur_M1_4gr", "Geweld_Voor_dicho", "Geweld_Tussen_dicho",
#              "Geweld_Na_dicho", "Diagnose_categorie", "Delict_Categorie",
#              "PersoneelsNummer_Random", "Functiegroep")
# measure.vars <- colnames(dd)[idx0]
#
# long <- melt(
#   data          = as.data.table(dd),
#   id.vars       = id.vars,
#   measure.vars  = measure.vars,
#   variable.name = "item",
#   value.name    = "score"
# )
# setcolorder(long, c(order, setdiff(colnames(long), order)))
# order <- c("Meetmoment", "Patientnummer_Random", "PersoneelsNummer_Random", "item", "score")
# long <- long[, , with = FALSE]
#
# long[, score := ifelse(score == -99, NA, score)]
# range(long$score, na.rm = TRUE)
#
# old <- colnames(long)
# new <- c(
#   "time", "patient", "rater", "item", "score", "patient_age_group", "treatement_duration_group",
#   "violent_before", "violent_between", "violent_after", "diagnosis_group", "crime_group",
#   "rater_group"
# )
# setnames(long, old, new)
# long
#
# long[item == "IFBE_0", PersoneelsNummer_Random]
# table(sub[Patientnummer_Random == "1640" & item == "IFBE_0" & Meetmoment == "2", PersoneelsNummer_Random])
#
#
# #
# # source("~/ownCloud/PhD/Cultural Consensus Theory/Papers/specialIssue_Batchelder/Cultural-Consensus-Theory-for-the-Evaluation-of-Patients-Mental-Health-Scores-in-Forensic-Psychiatr/R/simulate.R")
# # nr <- 20
# # ni <- 30
# # nc <- 5
# # np <- 50
# # nl <- 3
# # ncatR <- 2
# # ncatP <- 3
# #
# # set.seed(42)
# # obsList <- simulate(nr, ni, nc, np, nl, ncatP = ncatP, ncatR = ncatR)
# # datLongFormat <- reshapeSim(obsList)
#
# head(datLongFormat, 10)
# head(long, 10)
#
# no_patients <- 5
# no_items    <- 5
#
# unique_patients_sorted <- names(sort(table(long$Patientnummer_Random)))
#
# patients <- rev(unique_patients_sorted)[1:no_patients]
# items    <- unique(long$item)[1:no_items]
#
# sub <- long[Patientnummer_Random %in% patients & item %in% items, ]
# ggplot(data = sub, aes(x = factor(Meetmoment), y = score, group = factor(PersoneelsNummer_Random), color = factor(PersoneelsNummer_Random))) +
#   geom_point() +
#   geom_line() +
#   facet_grid(factor(Patientnummer_Random) ~ item) +
#   theme_bw() +
#   theme(legend.position = "none")
#
# table(sub[Patientnummer_Random == "1640" & item == "IFBE_0" & Meetmoment == "2", PersoneelsNummer_Random])
