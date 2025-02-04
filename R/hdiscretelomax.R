#' The Hurdle Discrete Pareto Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the 1-inflated discrete lomax distribution with parameters prob_q,
#' dlomax_alpha, and dlomax_prob_p. This is practically a hurdle geometric
#' distribution when dlomax_alpha is near 0. References to "bernoulli trials"
#' and the like here are not formally valid, but are instructive. Deviations
#' from the actual geometric distribution are written in parentheses.
#'
#' @param x,q vector of quantiles representing the number of failures in a
#' sequence of (non)Bernoulli trials before success occurs.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be
#' the number required.
#' @param prob_p (approximate) probability of success in each trial after the
#' first trial, if the first trial was not a success. \eqn{0 <} prob_p \eqn{\le 1}.
#' @param prob_q probability of success in the first trial. \eqn{0 \le} prob_q \eqn{\le 1}
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#'  hdiscretelomax is a dummy function to keep the hurdle geometric documentation
#'  together.
#'
#'  dhdiscretelomax gives the density, phdiscretelomax gives the distribution
#'  function, qhdiscretelomax, gives the quantile function, and rhdiscretelomax
#'  generates random deviates.
#'
#'  Invalid dlomax_prob_p or prob_q will result in return value NaN, with a warning.
#'
#' @rdname hdiscretelomax
#' @export
#' @usage hgeom(x, q, p, n, prob_p, prob_q, log, log.p, lower.tail)
hdiscretelomax <- function(x, q, p, n, prob_q, dlomax_alpha, dlomax_prob_p,
                           log, log.p, lower.tail) {
  list(x, q, p, n, prob_q, dlomax_alpha, dlomax_prob_p, log, log.p, lower.tail)
}

#' @rdname hdiscretelomax
#'
#' @examples
#' # Comparing hurdle discrete Lomax to hurdle discrete geometric mass functions
#' plot(dhdiscretelomax(x = 1:10, prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0.3))
#' points(dhdiscretelomax(x = 1:10, prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0), col = 'red')
#'
#' plot(dhdiscretelomax(x = 1:10, prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0.3)/dhdiscretelomax(x = 1:10, prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0))
#' @export
dhdiscretelomax <- function(x, prob_q, dlomax_prob_p, dlomax_alpha, log = FALSE) {
  prob_out <- numeric(length(x))
  prob_out[which(x == 1)] <- prob_q
  prob_out[which(x > 1)] <- (1 - prob_q) * ddiscretelomax(
    x = x[which(x > 1)] - 1,
    dlomax_alpha = dlomax_alpha, dlomax_prob_p = dlomax_prob_p)

  if(log) {prob_out <- log(prob_out)}
  prob_out
}

#' @rdname hdiscretelomax
#' @export
rhdiscretelomax <- function(n, prob_q, dlomax_prob_p, dlomax_alpha, log = FALSE) {
  n_out <- 1 + stats::rbinom(n = n, size = 1, prob = 1 - prob_q)
  n_gt_1 <- sum(n_out) - n
  n_out[which(n_out > 1)] <- 1 + rdiscretelomax(n = n_gt_1, dlomax_alpha, dlomax_prob_p)
  n_out
}
