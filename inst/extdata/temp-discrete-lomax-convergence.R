optim <- optim(par = c(dlomax_alpha = 0.05, dlomax_sigma = 0.5),
      fn = log_lik_discretelomax,
      gr = log_lik_descretelomax_gr,
      x = rdiscretelomax(n = 200, dlomax_alpha = 0.1, dlomax_prob_p = 0.5),
      lower = c(0, 10e-20),
      upper = c(Inf, 1 - 10e-20),
      method = "L-BFGS-B")$convergence

# Tends to fail when dlomax_alpha is very small
# Probably due to a log of ~0 resulting in a -Inf
# Err

sapply(
  X = 1:100,
  FUN = \(i) {
    try(
      expr = unname(janky_mle_discretelomax(rdiscretelomax(n = 200, dlomax_alpha = 0.01, dlomax_prob_p = 0.5))$estimate[1]),
      silent = TRUE
    )
  })

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
