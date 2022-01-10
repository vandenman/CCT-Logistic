#TODO:

# this could also be
get_prop <- function(x, nm) UseMethod("get_prop", x)
get_prop.default                  <- function(x, nm) x[[nm]]
get_prop.logistic_regression_ltm_data <- function(x, nm) get_prop(x$ltm_data, nm)

# and then get_np and fr
#' @export
get_np <- function(x) get_prop(x, "np")
#' @export
get_ni <- function(x) get_prop(x, "ni")
#' @export
get_nr <- function(x) get_prop(x, "nr")
#' @export
get_nc <- function(x) get_prop(x, "nc")
#' @export
get_no_rater_groups <- function(x) get_prop(x, "no_rater_groups")
#' @export
get_no_time_points <- function(x) get_prop(x, "no_time_points")

#' @export
get_score <- function(x)  UseMethod("get_score", x)
#' @export
get_score.ltm_data <- function(x) x$df$score
#' @export
get_score.logistic_regression_ltm_data <- function(x) get_no_covariates(x$ltm_data)

#' @export
get_no_covariates <- function(x) UseMethod("get_no_covariates", x)
#' @export
get_no_covariates.logistic_regression_ltm_data <- function(x) { ncol(x$log_reg$design_matrix) - get_ni(x) * get_no_time_points(x) }
#' @export
get_no_covariates.logistic_regression_data <- function(x) { ncol(x$log_reg$design_matrix) }

#' #' @export
#' get_np <- function(x) UseMethod("get_np", x)
#' #' @export
#' get_ni <- function(x) UseMethod("get_ni", x)
#' #' @export
#' get_nr <- function(x) UseMethod("get_nr", x)
#' #' @export
#' get_nc <- function(x) UseMethod("get_nc", x)
#' #' @export
#' get_no_rater_groups <- function(x) UseMethod("get_no_rater_groups", x)
#' #' @export
#' get_no_time_points <- function(x) UseMethod("get_no_time_points", x)
#'
#' #' @export
#' get_np.default <- function(x) x$np
#' #' @export
#' get_ni.default <- function(x) x$ni
#' #' @export
#' get_nr.default <- function(x) x$nr
#' #' @export
#' get_nc.default <- function(x) x$nc
#' #' @export
#' get_no_rater_groups <- function(x) x$no_rater_groups
#' #' @export
#' get_no_time_points.default <- function(x) x$no_time_points
#'
#' #' @export
#' get_np.logistic_regression_data <- function(x) get_np(x$ltm_data)
#' #' @export
#' get_ni.logistic_regression_data <- function(x) get_ni(x$ltm_data)
#' #' @export
#' get_nr.logistic_regression_data <- function(x) get_nr(x$ltm_data)
#' #' @export
#' get_nc.logistic_regression_data <- function(x) get_nc(x$ltm_data)
#' #' @export
#' get_no_rater_groups.logistic_regression_data <- function(x) get_no_rater_groups(x$ltm_data)
#' #' @export
#' get_no_time_points.logistic_regression_data <- function(x) get_no_time_points(x$ltm_data)
#'
