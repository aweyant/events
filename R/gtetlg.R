#' Title
#'
#' @inheritParams bgge
#'
#' @name gtetlg
NULL

#' @rdname gtetlg
#' @export
rgtetlg <- function(n, beta, prob_p, prob_q, rounding = FALSE) {
  rbgge(n = n,
        m = rhgeom(n = n, prob_p = prob_p, prob_q = prob_q),
        beta = beta,
        rounding = rounding)
}

pgtetlg <- function(
    x_lower = (1e-10)/beta,
    x_upper = Inf,
    y_lower = (1e-10)/beta,
    y_upper = Inf,
    n_lower = 1,
    n_upper = Inf,
    beta,
    prob_p,
    prob_q,
    x_lims = NULL,
    y_lims = NULL,
    n_lims = NULL) {
}

pgtetlgdt <- function(q, beta, prob_p, prob_q, lower.tail = TRUE) {

}

#
# cdf_gtetlg <- function(x, y, n, prob_p, prob_q, beta) {
#   # substitute(do.call(what = cdf_gtetlg_wrapped,
#   #                    args = list(x = x, y = y, n = n,
#   #                                event_duration_marginal_pmf_args = list(prob_p = prob_p, prob_q = prob_q),
#   #                                event_bivariate_conditional_cdf_args = list(beta = beta))))
#   cdf_gtetlg_wrapped(x = x, y = y, n = n,
#                      event_duration_marginal_pmf_args = list(prob_p = prob_p, prob_q = prob_q),
#                      event_bivariate_conditional_cdf_args = list(beta = beta))
# }

# cdf_gtetlg_wrapped <- do.call(what = cdf_ted,
#                               args = list(event_duration_marginal_pmf = dhgeom,
#                                           event_bivariate_conditional_cdf = cdf_bgge))


# dhgeom_clone <- as.function(alist(x = ., prob_p = ., prob_q = .,
#                                   dhgeom(x,prob_p,prob_q)))


# cdf_gtetlg_wrapped <- cdf_ted(event_duration_marginal_pmf = dhgeom,
#                               event_bivariate_conditional_cdf = cdf_bgge)

# cdf_gtetlg_wrapped <- methods(substitute(do.call(what = cdf_ted,
#                                                  args = list(event_duration_marginal_pmf = dhgeom,
#                                                              event_bivariate_conditional_cdf = cdf_bgge))))

