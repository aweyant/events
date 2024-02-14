test_that("psdeddt gives expected output given known parameters", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "n"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marignal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  expect_equal(psdeddt(q = q, length = length, threshold = threshold,
                       event_magnitude_conditional_cdf = psmp_marignal_x,
                       event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                       event_length_arg_name = event_length_arg_name),
               manual_calc_result)
})

test_that("psdeddt fails when the user provides an unrecognized argument for the event length", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "m"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marignal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  expect_error(psdeddt(q = q, length = length, threshold = threshold,
                       event_magnitude_conditional_cdf = psmp_marignal_x,
                       event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                       event_length_arg_name = event_length_arg_name))
})

test_that("psdeddt fails when the user overwrites an existing argument with their argument for the event length", {
  q = 100; length = 2; threshold = 20
  alpha = 0.02; beta = 0.04
  event_length_arg_name = "q"

  reconstructed_magnitude = q - length*threshold

  manual_calc_result <-  psmp_marignal_x(q = reconstructed_magnitude,
                                         n = length,
                                         alpha = alpha, beta = beta)

  expect_error(psdeddt(q = q, length = length, threshold = threshold,
                       event_magnitude_conditional_cdf = psmp_marignal_x,
                       event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
                       event_length_arg_name = event_length_arg_name))
})