assert_counts <- function (np, ni, nr, no_rater_groups) {
  assertthat::assert_that(
    assertthat::is.count(np),
    assertthat::is.count(ni),
    assertthat::is.count(nr),
    assertthat::is.count(no_rater_groups),
    nr > no_rater_groups
  )
}

validate_parameter_dimensions <- function(np, ni, nr, log_a, b, log_E, lt, log_lambda) {
  assertthat::assert_that(
    all(dim(lt)         == c(np, ni)),
    length(log_lambda)  == ni,
    length(log_a)       == nr,
    length(b)           == nr,
    length(log_E)       == nr
  )
}
