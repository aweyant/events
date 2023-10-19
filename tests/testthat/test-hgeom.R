test_that("pmf and cdf conform", {
  n <- 10
  p <- 0.5
  q <- 0.7
  expect_equal(sum(dhgeom(x = 0:n, prob_p = p, prob_q = q)),
               phgeom(q = n, prob_p = p, prob_q = q))
})
