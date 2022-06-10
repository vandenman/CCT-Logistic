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

#' @export
read_violence_data <- function()
  read_data("data_violence.rds")

#' @export
read_ifte_data <- function()
  read_data("data_ifte.rds")


#' @importFrom rlang %||%

#' @export
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
    redraw_figure(figure)
    dev.off()
  } else if (extension == "pdf") {
    grDevices::cairo_pdf(filename = file, ...)
    redraw_figure(figure)
    dev.off()
  } else {
    stop("Invalid extension '", extension, "'")
  }
}

redraw_figure <- function(figure) {
  if (inherits(figure, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(figure)
  } else {
    plot(figure)
  }
}

#' @export
save_figure_obj <- function(figure, file, dir = "figure_r_objs", ...) {

  ext <- tools::file_ext(file)
  if (ext == "") {
    file <- paste0(file, ".rds")
  } else if (ext != "rds") {
    stop("Expected extension .rds")
  }
  file <- file.path(dir, file)

  saveRDS(figure, file)
}

#' @export
load_figure_obj <- function(file, dir = "figure_r_objs", ...) {
  file <- file.path(dir, file)
  readRDS(file)
}

#' @export
normalize_factor <- function(x) {
  if (is.factor(x))
    as.integer(droplevels(x))
  else if (is.integer(x))
    match(x, unique(x))
}
