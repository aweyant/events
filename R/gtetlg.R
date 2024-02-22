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

cdf_gtetlg <- function(x, y, n,
                       beta, prob_p, prob_q,
                       tol = 10^(-6),
                       max_N = 100) {
  cdf_ted(x = x, y = y, n = n,
          event_bivariate_conditional_cdf = cdf_bgge)
}
