#' Sums & Maxima of Multivariate Pareto II Vectors
#'
#' @param n number of observations.
#' @param alpha,beta tail/shape parameter and scale parameters of a Multivariate
#' Pareto II distribution (see details for the parameterization used).
#' @param m the number of elements in the vectors. This is expected to either be
#' a single number (all generated vectors have the same number of elements) or a
#' vector of length n (in which case the number of elements in each vector is
#' specified)
#' @param rounding logical; if TRUE, generated Pareto II values are rounded to
#' the nearest whole number. This feature is present in case this is important
#' for a comparison of actual observations to some theoretical output.
#'
#' @return
#' a data.frame() with three columns:
#' \itemize{
#' \item n the number of elements in the Pareto II vector
#' \item x the sum of the elements in the Pareto II vector
#' \item y the maximal value in the Pareto II vector
#' }
#' @export
#'
#' @examples
#' rsmp(n = 3, m = c(3,4,2), alpha = 0.03, beta = 1/100)
#' \dontrun{
#'   n        x        y
#' 1 3 241.6733 157.2649
#' 2 4 357.7778 187.8859
#' 3 2 219.1081 140.8779
#' }
rsmp <- function(n, alpha, beta, m, rounding = FALSE) {
  event_vectors <- rmvlomax(n = n,
                            m = m,
                            alpha = alpha,
                            beta = beta,
                            rounding = rounding)
  data.frame(n = vapply(X = event_vectors, FUN = length, FUN.VALUE = 1L),
             x = vapply(X = event_vectors, FUN = sum, FUN.VALUE = 1.05),
             y = vapply(X = event_vectors, FUN = max, FUN.VALUE = 1.05))
}
