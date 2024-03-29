test_that("psdeddt gives expected output given known parameters", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marginal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  function_calc_result <- psdeddt(q = q, length = length, threshold = threshold,
                                  event_magnitude_conditional_cdf = psmp_marginal_x,
                                  event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                  event_length_arg_name = event_length_arg_name)

  expect_equal(function_calc_result,
               manual_calc_result)
})

test_that("psdeddt gives sensible output when values outside support are given; return_NA_outside_support = TRUE", {
  q = 100; length = 2; threshold = 51
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"

  reconstructed_magnitude = q - length*threshold

  # manual_calc_result <-  psmp_marginal_x(q = reconstructed_magnitude,
  #                                        n = length,
  #                                        alpha = alpha, beta = beta)

  function_calc_result <- psdeddt(q = q, length = length, threshold = threshold,
                                  event_magnitude_conditional_cdf = psmp_marginal_x,
                                  event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                  event_length_arg_name = event_length_arg_name,
                                  return_NA_outside_support = TRUE)

  expect_equal(function_calc_result,
               NaN)
})

test_that("psdeddt gives sensible output when values outside support are given; return_NA_outside_support = FALSE", {
  q = 100; length = 2; threshold = 51
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"

  reconstructed_magnitude = q - length*threshold

  # manual_calc_result <-  psmp_marginal_x(q = reconstructed_magnitude,
  #                                        n = length,
  #                                        alpha = alpha, beta = beta)

  function_calc_result <- psdeddt(q = q, length = length, threshold = threshold,
                                  event_magnitude_conditional_cdf = psmp_marginal_x,
                                  event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                  event_length_arg_name = event_length_arg_name,
                                  return_NA_outside_support = FALSE)

  expect_equal(function_calc_result,
               0)
})

test_that("psdeddt fails when the user provides an unrecognized argument for the event length", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "m"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marginal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  expect_error(psdeddt(q = q, length = length, threshold = threshold,
                       event_magnitude_conditional_cdf = psmp_marginal_x,
                       event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                       event_length_arg_name = event_length_arg_name))
})

test_that("psdeddt fails when the user overwrites an existing argument with their argument for the event length", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "q"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marginal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  expect_error(psdeddt(q = q, length = length, threshold = threshold,
                       event_magnitude_conditional_cdf = psmp_marginal_x,
                       event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                       event_length_arg_name = event_length_arg_name))
})

test_that("psdeddt returns an event total CDF which evaluates correctly", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"

  promise <- psdeddt(q = NULL, length = length, threshold = threshold,
                     event_magnitude_conditional_cdf = psmp_marginal_x,
                     event_length_arg_name = event_length_arg_name)

  promise_eval_result <- promise(q = q, length = length, threshold = threshold,
                                 event_magnitude_conditional_cdf_args = list(alpha = alpha,
                                                                             beta = beta))


  direct_eval_result <- psdeddt(q = q, length = length, threshold = threshold,
                                event_magnitude_conditional_cdf = psmp_marginal_x,
                                event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                                event_length_arg_name = event_length_arg_name)

  expect_equal(promise_eval_result, direct_eval_result)
})

# test_that("psdeddt vectorizes over q when event_magnitude_conditional_cdf is vectorized", {
#   q = seq.int(101,110,1); length = c(1,2,2,3,2,1,1,1,1,2); threshold = 20
#   alpha = 0.02; beta = 0.04
#   event_length_arg_name = "n"
#
#   promise <- psdeddt(q = NULL, length = length, threshold = threshold,
#                      event_magnitude_conditional_cdf = psmp_marginal_x,
#                      event_length_arg_name = event_length_arg_name)
#
#   promise_eval_result <- promise(q = q, length = length, threshold = threshold,
#                                  event_magnitude_conditional_cdf_args = list(alpha = alpha,
#                                                                              beta = beta))
#
#
#   direct_eval_result <- psdeddt(q = q, length = length, threshold = threshold,
#                                 event_magnitude_conditional_cdf = psmp_marginal_x,
#                                 event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
#                                 event_length_arg_name = event_length_arg_name)
#
#   expect_equal(promise_eval_result, direct_eval_result)
# })
