#' **H**urdle **G**eometric **S**ums & **M**axima of **P**areto II Vectors
#'
#' Probability and random generation functions for the trivariate HGSMP distribution
#'
#' @inheritParams smp
#' @inheritParams dhgeom
#' @param n_lower,n_upper each are individual numbers; limits of summation of
#' for the duration of the a trivariate distribution, where n_lower is the
#' lowest duration included in the summation and n_upper is the largest included
#' in the summation; for example, n_lower = 3, n_upper = 5 sums over {3,4,5},
#' whereas n_lower = 3, n_upper = 3 sums over {3}.
#' @inheritParams create_events_from_vector
#' @param log logical; if TRUE, probabilities p are given as log(p).
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

#' @rdname hgsmp
#' @export
phgsmp <- function(
    x_lower = 0,
    x_upper = Inf,
    y_lower = 0,
    y_upper = Inf,
    n_lower = 1,
    n_upper = Inf,
    alpha,
    beta,
    prob_p,
    prob_q) {
  phgsmp_wrapped()(x_lower = x_lower,
                   x_upper = x_upper,
                   y_lower = y_lower,
                   y_upper = y_upper,
                   n_lower = n_lower,
                   n_upper = n_upper,
                   alpha = alpha,
                   beta = beta,
                   prob_p = prob_p,
                   prob_q = prob_q)
}

phgsmp_wrapped <- function(...) {
  construct_pted(args = c(alist(x_lower = ,
                                x_upper = ,
                                y_lower = ,
                                y_upper = ,
                                n_lower = ,
                                n_upper = ,
                                tol = ,
                                max_N = ),
                          formals(cdf_smp),
                          formals(dhgeom)[-1]),
                 pbed = psmp,
                 event_duration_marginal_pmf = dhgeom,
                 pbed_argument_transformations = list("n = k"),
                 tol = 10^(-6),
                 max_N = 160)
}

#' @rdname hgsmp
#' @export
phgsmpdt <- function(q, threshold, alpha, beta, prob_p, prob_q, lower.tail = TRUE, log = FALSE) {
  phgsmpdt_wrapped()(q = q,
                     threshold = threshold,
                     alpha = alpha, beta = beta,
                     prob_p = prob_p, prob_q = prob_q, lower.tail = lower.tail, log = log)
}

phgsmpdt_wrapped <- function(...) {
  construct_pteddt(args = c(alist(threshold = ),
                            formals(psmp_marginal_x),
                            formals(dhgeom)[-1]),
                   event_magnitude_conditional_cdf = psmp_marginal_x,
                   event_duration_marginal_pmf = dhgeom,
                   event_magnitude_conditional_cdf_argument_transformations = list("n = k"))
}

#' @rdname hgsmp
#' @export
qhgsmpdt <- function(p, threshold,
                     alpha, beta,
                     prob_p, prob_q,
                     lower.tail = TRUE,
                     log = FALSE) {
  qhgsmpdt_wrapped()(p = p,
                     threshold = threshold,
                     alpha = alpha,
                     beta = beta,
                     prob_p = prob_p,
                     prob_q = prob_q,
                     lower.tail = lower.tail,
                     log = log)}

qhgsmpdt_wrapped <- function(...) {
  construct_quantile(args = c(alist(p = ),
                              formals(phgsmpdt)[-1]),
                     cdf = phgsmpdt,
                     cdf_support = c(0, Inf),
                     cdf_search_interval = c(0, 10000),
                     tol = 1e-6,
                     env = parent.frame())
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
