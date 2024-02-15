#' **B**ivariate Distribution with **G**amma and **G**eneralized
#' **E**xponential Marginals
#'
#' rbgge generates sums and maxima of IID exponential vectors of specified
#' lengths
#'
#' pbgge evaulates the marginal distribution of the sum of the exponential
#' components, called "magnitudes"; This marginal distribution is merely an
#' Erlang (Gamma) distribution with shape n and rate beta
#'
#' @inheritParams smp
#' @param beta numeric > 0; the rate parameter of the exponential distribution
#' which describes individual exceedances. It is the reciprocal of the mean.
#' @return
#' rbgge generates random values of the distribution and returns a them
#' data.frame() with one obervation per row and three columns:
#' \itemize{
#' \item n the number of elements in the IID Expoential vector
#' \item x the sum of the elements in the vector
#' \item y the maximal value in the vector
#' }
#'
#' @name bgge
NULL

#' @rdname bgge
#' @examples
#' rbgge(n = 3, m = c(3,4,2), beta = 1/100)
#' @export
rbgge <- function(n, m, beta, rounding = FALSE){
  event_vectors <- rmvlomax(n = n,
                            m = m,
                            alpha = 0,
                            beta = beta,
                            rounding = rounding)
  data.frame(n = vapply(X = event_vectors, FUN = length, FUN.VALUE = 1L),
             x = vapply(X = event_vectors, FUN = sum, FUN.VALUE = 1.05),
             y = vapply(X = event_vectors, FUN = max, FUN.VALUE = 1.05))
}

#' @rdname bgge
#' @export
pbgge_marginal_x <- function(q, n, beta, lower.tail = TRUE) {
  stats::pgamma(q = q, shape = n, rate = beta, lower.tail = lower.tail)
}
