#' Title
#'
#' @inheritParams smp
#'
#' @name hgsmp
NULL

#' @rdname hgsmp
#' @export
rhgsmp <- function(n, alpha, beta, prob_p, prob_q, rounding = FALSE) {

}

phgsmp <- function(
    x_lower = (1e-10)/beta,
    x_upper = Inf,
    y_lower = (1e-10)/beta,
    y_upper = Inf,
    n_lower = 1,
    n_upper = Inf,
    alpha,
    beta,
    prob_p,
    prob_q,
    x_lims = NULL,
    y_lims = NULL,
    n_lims = NULL) {
}

phgsmpdt <- function(q, alpha, beta, prob_p, prob_q, lower.tail = TRUE) {

}

cdf_hgsmp <- function(x, y, n,
                      alpha, beta, prob_p, prob_q,
                      tol = 10^(-6),
                      max_N = 100) {
}
