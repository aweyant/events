#' **H**urdle **D**iscrete **P**areto **S**ums & **M**axima of **P**areto II Vectors
#'
#' Probability and random generation functions for the trivariate HDPSMP distribution
#'
#' @inheritParams smp
#' @inheritParams dhdiscretelomax
#' @param n_lower,n_upper each are individual numbers; limits of summation of
#' for the duration of the a trivariate distribution, where n_lower is the
#' lowest duration included in the summation and n_upper is the largest included
#' in the summation; for example, n_lower = 3, n_upper = 5 sums over {3,4,5},
#' whereas n_lower = 3, n_upper = 3 sums over {3}.
#' @inheritParams create_events_from_vector
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @name hdpsmp
NULL

#' @rdname hdpsmp
#' @export
rhdpsmp <- function(n, alpha, beta, prob_q, dlomax_prob_p, dlomax_alpha, rounding = FALSE) {
  rsmp(n = n,
       m = rhdiscretelomax(n = n, prob_q = prob_q,
                            dlomax_prob_p = dlomax_prob_p,
                            dlomax_alpha = dlomax_alpha),
       alpha = alpha, beta = beta,
       rounding = rounding)
}

#' @rdname hdpsmp
#' @export
phdpsmp <- function(
    x_lower = 0, x_upper = Inf,
    y_lower = 0, y_upper = Inf,
    n_lower = 1, n_upper = Inf,
    alpha, beta,
    prob_q, dlomax_prob_p, dlomax_alpha) {
  phdpsmp_wrapped()(x_lower = x_lower,
                   x_upper = x_upper,
                   y_lower = y_lower,
                   y_upper = y_upper,
                   n_lower = n_lower,
                   n_upper = n_upper,
                   alpha = alpha,
                   beta = beta,
                   prob_q = prob_q,
                   dlomax_prob_p = dlomax_prob_p,
                   dlomax_alpha = dlomax_alpha)
}

phdpsmp_wrapped <- function(...) {
  construct_pted(args = c(alist(x_lower = ,
                                x_upper = ,
                                y_lower = ,
                                y_upper = ,
                                n_lower = ,
                                n_upper = ,
                                tol = ,
                                max_N = ),
                          formals(cdf_smp),
                          formals(dhdiscretelomax)[-1]),
                 pbed = psmp,
                 event_duration_marginal_pmf = dhdiscretelomax,
                 pbed_argument_transformations = list("n = k"),
                 tol = 10^(-10),
                 max_N = 10*(2^8))
}

#' @rdname hdpsmp
#' @examples
#' 1/(1-phdpsmpdt(q = 12, threshold = 1, alpha = 0.2, beta = 1/0.4,
#' prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0.3))
#'
#' 1/(1-phdpsmpdt(q = 12, threshold = 1, alpha = 0.2, beta = 1/0.4,
#' prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0))
#' @export
phdpsmpdt <- function(q, threshold, alpha, beta, prob_q,
                      dlomax_prob_p, dlomax_alpha,
                      lower.tail = TRUE, log = FALSE) {
  phdpsmpdt_wrapped()(q = q,
                     threshold = threshold,
                     alpha = alpha, beta = beta,
                     prob_q = prob_q,
                     dlomax_prob_p = dlomax_prob_p,
                     dlomax_alpha = dlomax_alpha,
                     lower.tail = lower.tail, log = log)
}

phdpsmpdt_wrapped <- function(...) {
  construct_pteddt(args = c(alist(threshold = ),
                            formals(psmp_marginal_x),
                            formals(dhdiscretelomax)[-1]),
                   event_magnitude_conditional_cdf = psmp_marginal_x,
                   event_duration_marginal_pmf = dhdiscretelomax,
                   event_magnitude_conditional_cdf_argument_transformations = list("n = k"))
}

#' @rdname hdpsmp
#' @examples
#' mean = 2; p_0 = 1/mean; alpha_0 = 0
#' alpha_1 = 0.1; p_1 = 0.538632
#' alpha_2 = 0.3; p_2 = 0.6348
#' qhdpsmpdt(p = 0.999,threshold = 1, alpha = 0.2, beta = 1,
#' prob_q = 0.7, dlomax_prob_p = p_0, dlomax_alpha = alpha_0)
#' qhdpsmpdt(p = 0.999,threshold = 1, alpha = 0.2, beta = 1,
#' prob_q = 0.7, dlomax_prob_p = p_1, dlomax_alpha = alpha_1)
#' qhdpsmpdt(p = 0.999,threshold = 1, alpha = 0.2, beta = 1,
#' prob_q = 0.7, dlomax_prob_p = p_2, dlomax_alpha = alpha_2)
#' @export
qhdpsmpdt <- function(p, threshold,
                     alpha, beta,
                     prob_q,
                     dlomax_prob_p, dlomax_alpha,
                     lower.tail = TRUE,
                     log = FALSE) {
  qhdpsmpdt_wrapped()(p = p,
                     threshold = threshold,
                     alpha = alpha,
                     beta = beta,
                     prob_q = prob_q,
                     dlomax_prob_p = dlomax_prob_p,
                     dlomax_alpha = dlomax_alpha,
                     lower.tail = lower.tail,
                     log = log)}

qhdpsmpdt_wrapped <- function(...) {
  construct_quantile(args = c(alist(p = ),
                              formals(phdpsmpdt)[-1]),
                     cdf = phdpsmpdt,
                     cdf_support = c(0, Inf),
                     cdf_search_interval = c(0, 10000),
                     tol = 1e-10,
                     env = parent.frame())
}

#' @rdname hdpsmp
#' @example
#' cdf_hdpsmp(x = 20, y = 5, n = 8, alpha = 0.2, beta = 1/0.4, prob_q = 0.7, dlomax_prob_p = 0.6, dlomax_alpha = 0.3)
#' @export
cdf_hdpsmp <- function(x, y, n,
                       alpha, beta,
                       prob_q, dlomax_prob_p, dlomax_alpha) {
  cdf_hdpsmp_wrapped()(duration = n,
                      x = x, y = y, n = NULL,
                      alpha = alpha,
                      beta = beta,
                      prob_q = prob_q,
                      dlomax_prob_p = dlomax_prob_p,
                      dlomax_alpha = dlomax_alpha,
                      log = FALSE)
}

cdf_hdpsmp_wrapped <- function(...) {
  construct_cdf_ted(args = c(alist(duration = ),
                             formals(cdf_smp),
                             formals(dhdiscretelomax)[-1]),
                    event_bivariate_conditional_cdf = cdf_smp,
                    event_duration_marginal_pmf = dhdiscretelomax,
                    conditional_cdf_argument_transformations = list("n = k"))
}

#' #' @rdname hdpsmp
#' #' @export
#' qhdpsmpdt <- function(p, threshold,
#'                      alpha, beta,
#'                      prob_q, dlomax_prob_p, dlomax_alpha,
#'                      lower.tail = TRUE,
#'                      log = FALSE) {
#'   qhdpsmpdt_wrapped()(p = p,
#'                      threshold = threshold,
#'                      alpha = alpha,
#'                      beta = beta,
#'                      prob_q = prob_q,
#'                      dlomax_prob_p = dlomax_prob_p,
#'                      dlomax_alpha = dlomax_alpha,
#'                      lower.tail = lower.tail,
#'                      log = log)}
#'
#' qhspsmpdt_wrapped <- function(...) {
#'   construct_quantile(args = c(alist(p = ),
#'                               formals(phdpsmpdt)[-1]),
#'                      cdf = phgsmpdt,
#'                      cdf_support = c(0, Inf),
#'                      cdf_search_interval = c(0, 10000),
#'                      tol = 1e-8,
#'                      env = parent.frame())
#' }
