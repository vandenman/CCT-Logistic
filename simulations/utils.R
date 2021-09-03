read_data <- function(filename)
  readRDS(file.path("data", filename))

read_long_data <- function()
  read_data("data_long.rds")

read_wide_data <- function()
  read_data("data_wide.rds")

read_wider_data <- function()
  read_data("data_wider.rds")

read_wider_imputed_data <- function()
  read_data("data_wider_imputed.rds")
