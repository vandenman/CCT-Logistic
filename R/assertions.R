assert_counts <- function (np, ni, nr, no_rater_groups) {
  assert_that(
    is.count(np),
    is.count(ni),
    is.count(nr),
    is.count(no_rater_groups),
    nr > no_rater_groups
  )
}

validate_parameter_dimensions <- function(np, ni, nr, log_a, b, log_E, lt, log_lambda) {
  assert_that(
    all(dim(lt)         == c(np, ni)),
    length(log_lambda)  == ni,
    length(log_a)       == nr,
    length(b)           == nr,
    length(log_E)       == nr
  )
}


is.bool <- function(x) {
  # isTRUE(x) || isFALSE(x) simplifies to
  is.logical(x) && length(x) == 1L && !is.na(x)
}
attr(is.bool, "fail") <- function (call, env) {
  sprintf("%s is not TRUE or FALSE (a length one logical vector).", deparse(call$x))
}

is.positive_int <- function(x, larger_than = 0) {
  if (length(x) != 1)
    return(FALSE)
  if (!is.integerish(x))
    return(FALSE)
  x > larger_than && !is.na(x)
}
attr(is.positive_int, "fail") <- function (call, env) {
  sprintf("%s is not a count (a single positive integer) larger than %d", deparse(call$x), call$larger_than)
}

is.integerish <- function (x) {
  # copy of assertthat:::is.integerish
  is.integer(x) || (is.numeric(x) && all(x == trunc(x)) && !is.na(x))
}
