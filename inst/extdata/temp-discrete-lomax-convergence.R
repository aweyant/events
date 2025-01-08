sample <- rdiscretelomax(n = 200, dlomax_alpha = 0.1, dlomax_prob_p = 0.5)
#trace(log_lik_descretelomax_gr, exit = quote(print(c(returnValue()))))
trace(log_lik_discretelomax, exit = quote(print(c(params,returnValue()))))

optim <- optim(par = c(al = 0.05, p = 0.5),
      fn = log_lik_discretelomax,
      #gr = log_lik_descretelomax_gr,
      x = sample,
      lower = c(0, 10e-20),
      upper = c(Inf, 1 - 10e-3),
      method = "L-BFGS-B")

#untrace(log_lik_descretelomax_gr)
untrace(log_lik_discretelomax)

# The convergence issue is due to optim trying to evaluate log_log_discretelomax
#' with p very near 1. This sometimes results in Inf. If I jankily restrict the
#' domain (with the upper argument), then this is not an issue.

sapply(
  X = 1:100,
  FUN = \(i) {
    try(
      expr = unname(janky_mle_discretelomax(rdiscretelomax(n = 200, dlomax_alpha = 0.01, dlomax_prob_p = 0.5))$estimate[1]),
      silent = TRUE
    )
  })

sample_010_05 <- rdiscretelomax(n = 200, dlomax_alpha = 0.1, dlomax_prob_p = 0.5)
sapply(
  X = 1:100,
  FUN = \(i) {
    try(
      expr = unname(janky_mle_discretelomax(rdiscretelomax(n = 200, dlomax_alpha = 0.1, dlomax_prob_p = 0.5))$estimate[1]),
      silent = TRUE
    )
  })

sapply(
  X = 1:100,
  FUN = \(i) {
    try(
      expr = unname(janky_mle_discretelomax(rdiscretelomax(n = 200, dlomax_alpha = 0.2, dlomax_prob_p = 0.5))$estimate[1]),
      silent = TRUE
      )
    })

# Errors are very likely when alpha is nearly 0.

#' Now there is a distinct issue
#' When alpha is small:
#' [82] "Error in optim(par = c(dlomax_alpha = 0.01, dlomax_sigma = 0.5), fn = log_lik_discretelomax,  :  non-finite finite-difference value [1]"

# Log tests ---------------------------------------------------------------
sample <- rdiscretelomax(n = 200, dlomax_alpha = 0.01, dlomax_prob_p = 0.5)
prod(ddiscretelomax(x = sample,
                    dlomax_alpha = 0.01, dlomax_prob_p = 0.5, log = FALSE))

exp(sum(ddiscretelomax(x = sample, dlomax_alpha = 0.01, dlomax_prob_p = 0.5, log = TRUE)))
