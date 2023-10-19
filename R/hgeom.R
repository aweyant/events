#' The Hurdle-geometric Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the 1-inflated geometric distribution with parameters prob_p and prob_q.
#'
#' @param x
#' @param prob_p
#' @param prob_q
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dhgeom <- function(x, prob_p, prob_q, log = FALSE){
  prob_out <- double(length(x))
  prob_out[which(x < 1)] <- 0
  prob_out[which(x == 1)] <- prob_q
  prob_out[which(x > 1)] <- (1 - prob_q)*prob_p*((1-prob_p)^(x[which(x > 1)]-2))
  if(log) {return(log(prob_out))}
  prob_out
}
phgeom <- function(q, prob_p, prob_q, lower.tail = TRUE, log.p = FALSE){
  nep <- double(length(q))
  nep[which(q < 1)] <- 0
  nep[which(q == 1)] <- prob_q
  nep[which(q > 1)] <- prob_q + (1 - prob_q) * (1 - ((1 - prob_p)^(q[which(q > 1)] -1)))
  if(lower.tail) {prob_out <- nep}
  else {prob_out <- 1-nep}

  if(log.p) {return(log(prob_out))}
  prob_out
}
qhgeom <- function(p, prob_p, prob_q, lower.tail = TRUE, log.p = FALSE){

}
rhgeom <- function(n, prob_p, prob_q){

}
