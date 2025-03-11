test_that("stopping point of a geometric distribution with param p = 0.1 is correctly determined", {
  expect_equal(determine_stopping_point(event_duration_marginal_pmf = dhgeom,
                                        event_duration_marginal_pmf_args = list(prob_q = 0.1,
                                                                                prob_p = 0.1), max_N = 1000), 160)
})

test_that("stopping point of a geometric distribution with param p = 0.1 is correctly determined, even when a user-provided lower.tail = FALSE argument is given", {
  expect_equal(determine_stopping_point(dhgeom,
                                        list(prob_p = 0.1,
                                             prob_q = 0.1),
                                        max_N = 1000), 160)
})

test_that("the search for the stopping point of a geometric distribution with param p = 0.1 throws an error when max_N is set too low", {
  expect_error(determine_stopping_point(dhgeom,
                                        list(prob_p = 0.1, prob_q = 0.1),
                                        max_N = 20))
})


test_that("stopping point is found for an actual estimated set of parameters", {
  dlomax_alpha = 0.41780147499835; dlomax_prob_p = 0.750413360606945; prob_q = 0.611111111111111
  expect_equal(
    determine_stopping_point(
      dhdiscretelomax,
      list(prob_q = prob_q, dlomax_prob_p = dlomax_prob_p, dlomax_alpha = dlomax_alpha),
      tol = 1e-6,
      max_N = 10*(2^9)
    ),
    640
  )
})

