#' Direct maximum likelihood estimation of the parameters of a hurdle-geometric
#' distribution
#'
#' This function is labeled "janky" not because it has not passed tests for all
#' realistic use cases so far, but because the author is not yet sure how output
#' of such a function ought to be formatted. This function can be safely used in
#' scripts, but it will either be superseded or "wrapped" later on.
#'
#' @param v a vector of integers thought to be drawn from a hurdle-geometric
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
#' sample_n <- rhgeom(n = 20, prob_p = 0.3, prob_q = 0.7)
#' janky_mle_hgeom(.Last.value)
janky_mle_hgeom <- function(v) {
  return(list(estimate = c(prob_p = janky_mle_hgeom_prob_p(v),
                           prob_q = janky_mle_hgeom_prob_q(v))))
}

janky_mle_hgeom_prob_q <- function(v) {
  n_1 <- sum(v == 1)
  n_2 <- sum(v != 1)
  q <- n_1/(n_1 + n_2)
  return(q)
}

janky_mle_hgeom_prob_p <- function(v) {
  j <- length(v[v != 1])
  p <- j/(sum(v) - length(v))
  return(ifelse(!is.nan(p),p,1))
}
