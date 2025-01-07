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
  optim(par = c(dlomax_alpha = 0.01, dlomax_sigma = 0.5),
        fn = log_lik_discretelomax,
        x = v,
        lower = c(0, 10e-20),
        upper = c(Inf, 1 - 10e-20),
        method = "L-BFGS-B")$trace
  # raw_pars <- unname(optim(par = c(dlomax_alpha = 0.01, dlomax_sigma = 0.5),
  #                          fn = log_lik_discretelomax,
  #                          x = v,
  #                          lower = c(0, 10e-20),
  #                          upper = c(Inf, 1 - 10e-20),
  #                          method = "L-BFGS-B")$par)
  #
  # return(
  #   list(
  #     estimate = c(dlomax_alpha = raw_pars[1],
  #                  dlomax_prob_p = raw_pars[2])
  #   )
  # )
}

log_lik_discretelomax <- function(params, x) {
  al <- params[1]
  p <- params[2]

  -sum(
    log(
      ddiscretelomax(x = x, dlomax_alpha = al, dlomax_prob_p = p)
      )
    )
  #-prod(ddiscretelomax(x = x, dlomax_alpha = al, dlomax_prob_p = p))
}

log_lik_descretelomax_gr <- function(params, x) {
  # params inputed in our prefered parameterization
  alpha_dp_star <- params[1]
  p_dp_star <- params[2]

  # params converted to original parameterization (as in eq 2.1 in Buddana and Kozubowski 2014)
  al = 1/alpha_dp_star; sig = 1/(-alpha_dp_star*(log(1 - p_dp_star)))

  c(
    sum(
      (log((sig)/(x + sig - 1)) / (1 - ((x + sig - 1)/(x + sig))^al)) +
        -(log(sig / (x + sig)) / (((x + sig) / (x + sig - 1))^al - 1))
    ),
    sum(
      (((x + sig - 1)*(x - 1))/(1 - ((x + sig - 1)/(x + sig))^al)) +
        -(((x + sig)*x)/(((x + sig)/(x + sig - 1))^al - 1))
    )
  )
}
