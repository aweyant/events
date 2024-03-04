test_that("Constructed quantile functions match base quantile functions; typical use case", {
  p = 0.5; rate = 1/10
  qexp_wrapped <- function(...) {
    construct_quantile(args = c(alist(p = ),
                                formals(pexp)[-1]),
                       cdf = pexp,
                       cdf_support = c(0, Inf),
                       cdf_search_interval = c(0, 1000),
                       tol = 1e-6,
                       env = parent.frame())
  }
  qexp_constructed <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE) {
    qexp_wrapped()(p = p, rate = rate, lower.tail = lower.tail, log.p = log.p)
  }
  constructed_result <- qexp_constructed(p = p, rate = rate)
  base_r_result <- stats::qexp(p = p, rate = rate)
  expect_equal(constructed_result, base_r_result)
})

test_that("Constructed quantile functions match base quantile functions; checking lower.tail = FALSE", {
  p = 0.75; rate = 1/10
  qexp_wrapped <- function(...) {
    construct_quantile(args = c(alist(p = ),
                                formals(pexp)[-1]),
                       cdf = pexp,
                       cdf_support = c(0, Inf),
                       cdf_search_interval = c(0, 1000),
                       tol = 1e-6,
                       env = parent.frame())
  }
  qexp_constructed <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE) {
    qexp_wrapped()(p = p, rate = rate, lower.tail = lower.tail, log.p = log.p)}
  constructed_result <- qexp_constructed(p = p, rate = rate, lower.tail = FALSE)
  base_r_result <- stats::qexp(p = p, rate = rate, lower.tail = FALSE)
  expect_equal(constructed_result, base_r_result)
})

test_that("Constructed quantile functions match base quantile functions; typical use case, falsely vectorized p", {
  p = c(0,0.2,0.5,0.75,0.95,0.99); rate = 1/10
  qexp_wrapped <- function(...) {
    construct_quantile(args = c(alist(p = ),
                                formals(pexp)[-1]),
                       cdf = pexp,
                       cdf_support = c(0, Inf),
                       cdf_search_interval = c(0, 1000),
                       tol = 1e-6,
                       env = parent.frame())
  }
  qexp_constructed <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE) {
    qexp_wrapped()(p = p, rate = rate, lower.tail = lower.tail, log.p = log.p)}
  qexp_constructed <- Vectorize(qexp_constructed, vectorize.args = c("p", "rate"))
  constructed_result <- qexp_constructed(p = p, rate = rate)
  base_r_result <- stats::qexp(p = p, rate = rate)
  expect_equal(constructed_result, base_r_result)
})
