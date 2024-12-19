#' Direct maximum likelihood estimation of the parameters of a discrete Lomax
#' distribtuion
#'
#' @param v a vector of integers thought to be drawn from a discrete Lomax
#' distribution
#'
#' @return
#' A list with elements:
#' \itemize{
#' \item estimate - a vector of estimates corresponding to each parameter in the
#' distribution
#' }
#' @export
#'
#' @examples
#' sample <- rdiscretelomax(n = 100, dlomax_alpha = 0.01, dlomax_prob_p = 0.5)
#' janky_mle_discretelomax(sample)
janky_mle_discretelomax <- function(v) {
  raw_pars <- unname(optim(par = c(dlomax_alpha = 0.01, dlomax_sigma = 0.5),
                           fn = log_lik_discretelomax,
                           x = v,
                           lower = c(0,0),
                           upper = c(Inf, 1 - 10e-20),
                           method = "L-BFGS-B")$par)

  return(
    list(
      estimate = c(dlomax_alpha = raw_pars[1],
                   dlomax_prob_p = raw_pars[2])
    )
  )
}

log_lik_discretelomax <- function(params, x) {
  al <- params[1]
  p <- params[2]

  -sum(
    log(
      ddiscretelomax(x = x, dlomax_alpha = al, dlomax_prob_p = p)
      )
    )
}
