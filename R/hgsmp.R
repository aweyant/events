#' Title
#'
#' @inheritParams smp
#' @inheritParams dhgeom
#'
#' @name hgsmp
NULL

#' @rdname hgsmp
#' @export
rhgsmp <- function(n, alpha, beta, prob_p, prob_q, rounding = FALSE) {
  rsmp(n = n,
       m = rhgeom(n = n, prob_p = prob_p, prob_q = prob_q),
       alpha = alpha,
       beta = beta,
       rounding = rounding)
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

#' @rdname hgsmp
#' @export
cdf_hgsmp <- function(x, y, n,
                      alpha, beta, prob_p, prob_q) {
  cdf_hgsmp_wrapped()(duration = n,
                      x = x, y = y, n = NULL,
                      alpha = alpha, beta = beta,
                      prob_p = prob_p, prob_q = prob_q,
                      log = FALSE)
}

cdf_hgsmp_wrapped <- function(...) {
  construct_cdf_ted(args = c(alist(duration = ),
                             formals(cdf_smp),
                             formals(dhgeom)[-1]),
                    event_bivariate_conditional_cdf = cdf_smp,
                    event_duration_marginal_pmf = dhgeom,
                    conditional_cdf_argument_transformations = list("n = k"))
}
