#' @export
read_data <- function(filename)
  readRDS(file.path("data", filename))

#' @export
read_long_data <- function()
  read_data("data_long.rds")

#' @export
read_wide_data <- function()
  read_data("data_wide.rds")

#' @export
read_wider_data <- function()
  read_data("data_wider.rds")

#' @export
read_wider_imputed_data <- function()
  read_data("data_wider_imputed.rds")

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

save_figure <- function(figure, file, dir = "figures", extension = c("svg", "pdf"), ...) {

  ext <- tools::file_ext(file)
  if (ext == "") {
    extension <- match.arg(extension)
    file <- paste0(file, extension)
  } else {
    extension <- ext
  }
  file <- file.path(dir, file)

  if (extension == "svg") {
    svglite::svglite(file = file, ...)
    plot(figure)
    dev.off()
  } else if (extension == "pdf") {
    grDevices::cairo_pdf(filename = file, ...)
    plot(figure)
    dev.off()
  } else {
    stop("Invalid extension '", extension, "'")
  }
}
