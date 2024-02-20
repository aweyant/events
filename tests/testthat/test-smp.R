# test_that("psmp_marginal_x() is properly vectorized over q", {
#   q = seq.int(1,200,1)
#   n = 2
#   alpha = 0.004; beta = 2/100
#
#   proper_vectorized_result <- psmp_marginal_x(q = q,
#                                               n = n,
#                                               alpha = alpha,
#                                               beta = beta)
#   homebrew_vectorized_result <- vapply(X = q,
#                                        FUN = function(cur_q) {
#                                          psmp_marginal_x(q = cur_q,
#                                                          n = n,
#                                                          alpha = alpha,
#                                                          beta = beta)
#                                        },
#                                        FUN.VALUE = numeric(1))
#
#   expect_equal(proper_vectorized_result, homebrew_vectorized_result)
# })
