test_that("pteddt direct evaluation agrees with promise evaluation", {
  q = 100; length = 2; threshold = 20
  prob_p = 0.6; prob_q = 0.7
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"
  # Still need to test

  promise <- pteddt(event_duration_marginal_pmf = dhgeom,
                    event_length_arg_name = event_length_arg_name,
                    event_magnitude_conditional_cdf = psmp_marginal_x)

  promise_eval <- promise(q = q,
                          threshold = threshold,
                          event_duration_marginal_pmf_args = list(prob_p = prob_p,
                                                                  prob_q = prob_q),
                          event_magnitude_conditional_cdf_args = list(alpha = alpha,
                                                                      beta = beta))

  direct_eval <- pteddt(q = q,
                        threshold = threshold,
                        event_duration_marginal_pmf = dhgeom,
                        event_duration_marginal_pmf_args = list(prob_p = prob_p,
                                                                prob_q = prob_q),
                        event_length_arg_name = event_length_arg_name,
                        event_magnitude_conditional_cdf = psmp_marginal_x,
                        event_magnitude_conditional_cdf_args = list(alpha = alpha,
                                                                    beta = beta))
  expect_equal(promise_eval, direct_eval)
})
