#' The Hurdle-geometric Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the 1-inflated geometric distribution with parameters prob_p and prob_q.
#'
#' @param x,q vector of quantiles representing the number of failures in a
#' sequence of Bernoulli trials before success occurs.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be
#' the number required.
#' @param prob_p probability of success in each trial after the first trial, if
#' the first trial was not a success. \eqn{0 <} prob_p \eqn{\le 1}.
#' @param prob_q probability of success in the first trial. \eqn{0 \le} prob_q \eqn{\le 1}
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#'  hgeom is a dummy function to keep the hurdle geometric documentation
#'  together.
#'
#'  dhgeom gives the density, phgeom gives the distribution function, qhgeom
#'  gives the quantile function, and rhgeom generates random deviates.
#'
#'  Invalid prob_p or prob_q will result in return value NaN, with a warning.
#'
#'  rhgeom returns a vector of type integer unless generated values exceed the
#'  maximum representable integer when double values are returned since R
#'  version 4.0.0.
#'
#' @rdname hgeom
#' @export
#' @usage hgeom(x, q, p, n, prob_p, prob_q, log, log.p, lower.tail)
hgeom <- function(x, q, p, n, prob_p, prob_q, log, log.p, lower.tail) {
  list(x, q, p, n, prob_p, prob_q, log, log.p, lower.tail)
  }

#' @rdname hgeom
#' @examples
#' plot(dhgeom(x = 1:25, prob_p = 0.3, prob_q = 0.6))
#' @export
dhgeom <- function(x, prob_p, prob_q, log = FALSE){
  hgeom_check_args(prob_p, prob_q)
  prob_out <- double(length(x))
  prob_out[which(x < 1)] <- 0
  prob_out[which(x == 1)] <- prob_q
  prob_out[which(x > 1)] <- (1 - prob_q)*prob_p*((1-prob_p)^(x[which(x > 1)]-2))
  if(log) {return(log(prob_out))}
  prob_out
}
#' @rdname hgeom
#' @export
phgeom <- function(q, prob_p, prob_q, lower.tail = TRUE, log.p = FALSE){
  #hgeom_check_args(prob_p, prob_q)
  nep <- double(length(q))
  nep[which(q < 1)] <- 0
  nep[which(q == 1)] <- prob_q
  nep[which(q > 1)] <- prob_q + (1 - prob_q) * (1 - ((1 - prob_p)^(q[which(q > 1)] -1)))
  if(lower.tail) {prob_out <- nep}
  else {prob_out <- 1-nep}

  if(log.p) {return(log(prob_out))}
  prob_out
}
#' @rdname hgeom
#' @export
qhgeom <- function(p, prob_p, prob_q, lower.tail = TRUE, log.p = FALSE){
  #hgeom_check_args(prob_p, prob_q)
  x <- integer(length(p))
  if(log.p) {p <- exp(p)}
  if(!lower.tail) {p <- 1-p}

  x[which(p < prob_q)] <- 1
  x[which(p > prob_q)] <- ceiling(1 + log(1 - (p[which(p > prob_q)]-prob_q)/(1-prob_q))/log(1 - prob_p))
  x
}
#' @rdname hgeom
#' @export
rhgeom <- function(n, prob_p, prob_q){
  hgeom_check_args(prob_p, prob_q)
  1 + (stats::rbinom(n = n, size = 1, prob = 1 - prob_q)) * (1 + stats::rgeom(n = n, prob = prob_p))
}

hgeom_check_args <- function(prob_p, prob_q) {
  prob_p_invalid <- FALSE
  if(prob_p <= 0 | prob_p > 1) {
    warning("Prob_p must be greater than 0 and less than or equal to 1. Returning NaN.")
    prob_p_invalid <- TRUE
  }
  prob_q_invalid <- FALSE
  if(prob_q < 0 | prob_q > 1) {
    warning("Prob_q must be greater than or equal to 0 and less than or equal to 1. Returning NaN.")
    prob_q_invalid <- TRUE
  }
  if(prob_p_invalid | prob_q_invalid) {
    call <- rlang::expr(return(NaN))
    rlang::eval_bare(call, env = parent.frame())
  }
  NULL
}
