test_that("pmf and cdf conform", {
  n <- 10
  p <- 0.5
  q <- 0.7
  expect_equal(sum(dhgeom(x = 0:n, prob_p = p, prob_q = q)),
               phgeom(q = n, prob_p = p, prob_q = q))
})

test_that("quantile is truly inverse of cdf", {
  p = c(0.6, 0.75, 0.99)
  prob_p <- 0.42
  prob_q <- 0.7
  # The quantile x is defined to be minimal x s.t. P(X<=x) >= p
  purported_quantiles <- qhgeom(p = p, prob_p = prob_p, prob_q = prob_q)
  inequality_holds <- as.logical(prod(phgeom(q = purported_quantiles, prob_p = prob_p, prob_q = prob_q) >= p))

  minimality_holds <- length(p) == sum(phgeom(q = purported_quantiles - 1, prob_p = prob_p, prob_q = prob_q) < p)

  expect_identical(inequality_holds & minimality_holds,
                   TRUE)
})
