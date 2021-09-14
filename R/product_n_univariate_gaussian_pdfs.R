product_n_univariate_gaussian_pdfs <- function(mus, sds) {

  assertthat::assert_that(length(mus) == length(sds), all(sds > 0))
  precs <- 1 / sds^2
  variance <- 1 / sum(precs)
  mean     <- sum(mus * precs) * variance

  return(c(mean, sqrt(variance)))

}

# TODO: make the stuff below into a unit test!
# foo <- function(x, mus, sds) {
#   assertthat::assert_that(length(mus) == length(sds))
#   res <- numeric(length(x))
#   for (i in seq_along(x)) {
#     # res[i] <- prod(dnorm(x[i], mus, sds))
#     res[i] <- exp(sum(dnorm(x[i], mus, sds, log = TRUE)))
#   }
#   res
# }
#
# ns <- seq(5, 100, 5)
# len_ns <- length(ns)
# nc <- ceiling(sqrt(length(ns)))
# nr <- floor(sqrt(length(ns)))
# xx <- seq(-5, 5, .01)
#
# layout(matrix(c(seq_along(ns), rep(0L, nr * nc - len_ns)), nr, nc, TRUE))
# for (n in ns) {
#   mus <- rnorm(n)
#   sds <- abs(rnorm(n, 4))
#   yy <- foo(xx, mus, sds)
#   int <- integrate(foo, -Inf, Inf, mus = mus, sds = sds, abs.tol = sqrt(.Machine$double.eps), subdivisions = 3e2L)
#   yy <- yy / int$value
#
#   propto <- product_n_univariate_gaussian_pdfs(mus, sds)
#   yrep <- dnorm(xx, propto[1], propto[2])
#
#   testthat::expect_equal(yy, yrep, tolerance = 1e-3)
#
#   plot(xx, yy, type = 'l', main = paste("n =", n))
#   lines(xx, yy, col = 2, lty = 2)
#
# }

