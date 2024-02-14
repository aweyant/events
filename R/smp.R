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


psmp <- function(x_lower = 0,
                 x_upper = Inf,
                 y_lower = 0,
                 y_upper = Inf,
                 n = 1,
                 alpha,
                 beta,
                 x_lims = NULL,
                 y_lims = NULL) {

  if(!is.null(x_lims)) {
    x_lower = x_lims[1]
    x_upper = x_lims[2]
  }
  if(!is.null(y_lims)) {
    y_lower = y_lims[1]
    y_upper = y_lims[2]
  }
  if(n == 1) {
    return(cdf_smp(x = x_upper, y = x_upper, n = n, alpha = alpha, beta = beta) +
             - cdf_smp(x = x_lower, y = x_lower, n = n, alpha = alpha, beta = beta))
  }
  cdf_smp(x = x_upper, y = y_upper, n = n, alpha = alpha, beta = beta) +
    -cdf_smp(x = x_lower, y = y_upper, n = n, alpha = alpha, beta = beta) +
    #-cdf_smp(x = x_lower, y = min(c(x_lower,y_upper)), n = n, alpha = alpha, beta = beta) +
    -cdf_smp(x = x_upper, y = y_lower, n = n, alpha = alpha, beta = beta) +
    #cdf_smp(x = x_lower, y = min(c(x_lower,y_lower)), n = n, alpha = alpha, beta = beta)
    cdf_smp(x = x_lower, y = y_lower, n = n, alpha = alpha, beta = beta)
}

#' Title
#'
#' @param q
#' @param n
#' @param alpha
#' @param beta
#' @param lower.tail
#'
#' @return
#' @export
#'
#' @examples
psmp_marginal_x <- function(q, n, alpha, beta, lower.tail = TRUE) {
  if(lower.tail) {
    return(psmp(x_lower = 0,
                x_upper = q,
                y_lower = 0,
                y_upper = Inf,
                n = n,
                alpha,
                beta,
                x_lims = NULL,
                y_lims = NULL))
  }
  else{
    return(psmp(x_lower = q,
                x_upper = Inf,
                y_lower = 0,
                y_upper = Inf,
                n = n,
                alpha,
                beta,
                x_lims = NULL,
                y_lims = NULL))
  }
}

#' Title
#'
#' @param q
#' @param n
#' @param alpha
#' @param beta
#' @param lower.tail
#'
#' @return
#' @export
#'
#' @examples
psmp_marginal_y <- function(q, n, alpha, beta, lower.tail = TRUE) {
  if(lower.tail) {
    return(psmp(x_lower = 0,
                x_upper = Inf,
                y_lower = 0,
                y_upper = q,
                n = n,
                alpha,
                beta,
                x_lims = NULL,
                y_lims = NULL))
  }
  else{
    return(psmp(x_lower = 0,
                x_upper = Inf,
                y_lower = q,
                y_upper = Inf,
                n = n,
                alpha,
                beta,
                x_lims = NULL,
                y_lims = NULL))
  }
}

cdf_smp <- function(x, y, n, alpha, beta) {
  if(is.infinite(x) & x >= 0 &
     is.infinite(y) & y >= 0) {
    return(1)
  }
  if(n == 1) {
    if(x <= 0) {return(0)}
    if(is.infinite(x)) {return(1)}
    return(1 - (1 + alpha*beta*x)^(-1/alpha))
  }
  if(x < 0 | y < 0) {
    # WRITTEN AND TESTED
    # print("case 1")
    return(0)
  }
  else if(0 <= x & x <= y) {
    # WRITTEN, UNTESTED
    # print("case 2")
    return(1 - sum(sapply(X = 0:(n-1),
                          FUN = function(i) {
                            exp( i * log(alpha*beta*x) - lfactorial(i) +
                                   lgamma(i + 1/alpha) - lgamma(1/alpha) +
                                   (-1/alpha - i) * log(1 + alpha*beta*x))
                          })))
  }
  else if(0 <= n*y & n*y <= x) {
    # WRITTEN, UNTESTED
    # print("case 4")
    return(sum(sapply(X = (0:n),
                      FUN = function(i) {
                        ((-1)^i) * exp(lchoose(n,i) + (-1/alpha)*log(1 + i*alpha*beta*y))
                      })))
  }
  else{
    # WRITTEN, UNTESTED
    # print("case 3")
    k <- numeric()
    for(i in seq.int(from = 1, to = n-1, by = 1)) {
      if(0 <= x/(i+1) & x/(i+1) <= y & y <= x/i) {
        k <- i
      }
    }
    if(is.infinite(x)) {
      return(sum(sapply(X = 0:n,
                        FUN = function(i){
                          ((-1)^i) * exp(lchoose(n,i) + (-1/alpha)*log(1 + i*alpha*beta*y))
                        })))
    }
    return((fun_b(x,y,n,alpha,beta,k) +
              fun_c(x,y,n,alpha,beta,k)) -
             fun_d(x,y,n,alpha,beta,k))
  }
}


# Internal Functions ------------------------------------------------------
fun_b <- function(x,y,n,alpha,beta,k) {
  # INSPECTED, SEEMS GOOD
  factorial(n) * sum(sapply(X = 1:(k-1),
                            FUN = function(s) {
                              (((-1)^(s+1))/(factorial(s)*factorial(n-s))) *
                                ((1-s/k)^(n-1)-(1 + alpha*beta*s*y)^(-1/alpha))
                            })) +
    factorial(n) * sum(sapply(X = 1:(k-1),
                              FUN = function(s){
                                (((-1)^(s+1))/(factorial(s)*factorial(n-s))) *
                                  sum(sapply(X = 0:(n-2),
                                             FUN = function(m){
                                               (((alpha*beta*y*k)^m)/factorial(m))*
                                                 (gamma(m + 1/alpha)/gamma(1/alpha)) *
                                                 (((1 - s/k)^m) - (1 - s/k)^(n-1)) *
                                                 ((1 + alpha * beta * k * y)^(-1/alpha - m))
                                             }))
                              }))
}

fun_c <- function(x,y,n,alpha,beta,k) {
  if(k == 1) {
    return((1 - sum(sapply(X = (0:(n-1)),
                           FUN = function(i){
                             (((alpha*beta*x)^i)/factorial(i))*
                               (gamma(i + 1/alpha)/gamma(1/alpha)) *
                               ((1 + alpha * beta * x)^(-1/alpha - i))
                           }))) *
             (factorial(n) * sum(sapply(X = (1:(n-1)),
                                        FUN = function(s){
                                          ((-1)^(s+1))*((1-s/n)^(n-1))/(factorial(s)*factorial(n-s))
                                        })))
    )
  }
  (1 - sum(sapply(X = (0:(n-1)),
                  FUN = function(i){
                    (((alpha*beta*x)^i)/factorial(i))*
                      (gamma(i + 1/alpha)/gamma(1/alpha)) *
                      ((1 + alpha * beta * x)^(-1/alpha - i))
                  }))) *
    (factorial(n) * sum(sapply(X = (1:(n-1)),
                               FUN = function(s){
                                 ((-1)^(s+1))*((1-s/n)^(n-1))/(factorial(s)*factorial(n-s))
                               })) -
       factorial(n) * sum(sapply(X = (1:(k-1)),
                                 FUN = function(s){
                                   ((-1)^(s+1))*((1-s/k)^(n-1))/(factorial(s)*factorial(n-s))
                                 })))
}

fun_d <- function(x,y,n,alpha,beta,k) {
  # INSPECTED, LOOKS GOOD
  factorial(n) * sum(sapply(X = (1:k),
                            FUN = function(s){
                              (((-1)^(s+1))/(factorial(s)*factorial(n-s))) *
                                sum(sapply(X = (0:(n-1)),
                                           FUN = function(m){
                                             (((alpha*beta)^m)/factorial(m))*
                                               (gamma(1/alpha + m)/gamma(1/alpha)) *
                                               (((y^m)*
                                                   ((k-s)^m - (((1-s/k)^(n-1))*(k^m)))) *
                                                  ((1 + alpha*beta*k*y)^(-1/alpha - m)) +
                                                  ((x^m)*((1-s/k)^(n-1)) - (x-s*y)^m) *
                                                  ((1 + alpha*beta*x)^(-1/alpha - m)))
                                           }))
                            }))
}
