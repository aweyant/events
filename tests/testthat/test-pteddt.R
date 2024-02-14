test_that("pteddt gives expected output given known parameters; iteration 1", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  prob_p = 0.6; prob_q = 0.7
  event_length_arg_name = "n"

  manual_calc_result <- (phgeom(q = 1,
                                prob_p = prob_p,
                                prob_q = prob_q) *
                           psmp_marginal_x(q = q - 1*threshold,
                                           n = 1,
                                           alpha = alpha,
                                           beta = beta)) + (
                                             phgeom(q = 2,
                                                    prob_p = prob_p,
                                                    prob_q = prob_q) *
                                               psmp_marginal_x(q = q - 2*threshold,
                                                               n = 2,
                                                               alpha = alpha,
                                                               beta = beta)
                                             )

  function_calc_result <- pteddt(q = q,
                                 threshold = threshold,
                                 event_duration_marginal_cdf = phgeom,
                                 event_duration_marginal_cdf_args = list(prob_p = prob_p, prob_q = prob_q),
                                 event_length_arg_name = event_length_arg_name,
                                 event_magnitude_conditional_cdf = psmp_marginal_x,
                                 event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                 lower.tail = TRUE,
                                 tol = 10^(-6),
                                 max_N = 100)

  expect_equal(manual_calc_result, function_calc_result)
})

test_that("pteddt gives expected output given known parameters; iteration 2", {
  q = 200; length = 3; threshold = 20
  alpha = 0.02; beta = 0.04
  prob_p = 0.6; prob_q = 0.7
  event_length_arg_name = "n"

  manual_calc_result <- (phgeom(q = 1,
                                prob_p = prob_p,
                                prob_q = prob_q) *
                           psmp_marginal_x(q = q - 1*threshold,
                                           n = 1,
                                           alpha = alpha,
                                           beta = beta)) + (
                                             phgeom(q = 2,
                                                    prob_p = prob_p,
                                                    prob_q = prob_q) *
                                               psmp_marginal_x(q = q - 2*threshold,
                                                               n = 2,
                                                               alpha = alpha,
                                                               beta = beta)) + (
                                                                 phgeom(q = 3,
                                                                        prob_p = prob_p,
                                                                        prob_q = prob_q) *
                                                                   psmp_marginal_x(q = q - 3*threshold,
                                                                                   n = 2,
                                                                                   alpha = alpha,
                                                                                   beta = beta))


  function_calc_result <- pteddt(q = q,
                                 threshold = threshold,
                                 event_duration_marginal_cdf = phgeom,
                                 event_duration_marginal_cdf_args = list(prob_p = prob_p, prob_q = prob_q),
                                 event_length_arg_name = event_length_arg_name,
                                 event_magnitude_conditional_cdf = psmp_marginal_x,
                                 event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                 lower.tail = TRUE,
                                 tol = 10^(-6),
                                 max_N = 100)

  expect_equal(manual_calc_result, function_calc_result)
})

