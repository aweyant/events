#' **B**ivariate Distribution with **G**amma and **G**eneralized
#' **E**xponential Marginals (BGGE)
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

#' @rdname bgge
#' @export
pbgge <- function(x_lower = 0,
                  x_upper = Inf,
                  y_lower = 0,
                  y_upper = Inf,
                  n = 1,
                  beta) {
  psmp(x_lower = x_lower,
       x_upper = x_upper,
       y_lower = y_lower,
       y_upper = y_upper,
       n = n,
       alpha = 0,
       beta = beta)
}

#' @rdname bgge
#' @export
cdf_bgge <- function(x, y, n, beta) {
  if(is.infinite(x) & x >= 0 &
     is.infinite(y) & y >= 0) {
    return(1)
  }
  if(n == 1) {
    x <- y
    if(x <= 0) {return(0)}
    if(is.infinite(x)) {return(1)}
    return(stats::pexp(q = x, rate = beta))
  }
  if(x <= 0 | y <= 0) {
    # WRITTEN AND TESTED
    # print("case 1")
    return(0)
  }
  else if(0 <= x & x <= y) {
    # WRITTEN, UNTESTED
    # print("case 2")
    index <- seq.int(0,n-1,1)
    return(1 - sum(exp(index * log(beta*x) - beta*x - lfactorial(index))))
  }
  else if(0 <= n*y & n*y <= x) {
    # WRITTEN, UNTESTED
    # print("case 4")
    return((1-exp(-beta*y))^n)
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
    # print(k)
    return((fun_A(x,y,n,beta,k) +
              fun_E(x,y,n,beta,k)) +
             -fun_D(x,y,n,beta,k))
  }
}


# Helper functions --------------------------------------------------------
fun_A <- function(x,y,n,beta,k) {
  if(k == 1) {return(0)}
  index_s <- seq.int(1,k-1,1)
  index_m <- seq.int(0,n-2,1)

  factorial(n) * sum( (((-1)^(index_s+1))/(factorial(index_s)*factorial(n-index_s))) * ((1-index_s/k)^(n-1) - exp(-index_s*beta*y)) ) +
    factorial(n) * exp(-k*beta*y) * sum(vapply(X = index_s,
                                               FUN = function(s) {
                                                 (((-1)^(s + 1))/(factorial(s) * factorial(n - s))) *
                                                   sum(
                                                     (((y*beta*k)^index_m)/factorial(index_m)) * ( ((1-s/k)^index_m) - ((1-s/k)^(n-1)) )
                                                     )
                                               },
                                               FUN.VALUE = numeric(1)))
}
fun_E <- function(x,y,n,beta,k) {
  index_m = seq.int(0,n-1,1)
  index_s1 = seq.int(1,n-1,1)
  if(k == 1) {
    return((1 - exp(-beta*x)*sum(
      ((beta*x)^index_m)/factorial(index_m))) *
        (factorial(n) * sum((((-1)^(index_s1 + 1)) * (1 - index_s1/n)^(n-1))/(factorial(index_s1)*factorial(n-index_s1))) +
           -0
         )
      )
  }
  index_s2 = seq.int(1,k-1,1)

  (1 - exp(-beta*x)*sum(
    ((beta*x)^index_m)/factorial(index_m))) *
    (factorial(n) * sum((((-1)^(index_s1 + 1)) * (1 - index_s1/n)^(n-1))/(factorial(index_s1)*factorial(n-index_s1))) +
       -factorial(n) * sum((((-1)^(index_s2+1)) * ((1 - index_s2/k)^(n-1)))/(factorial(index_s2)*factorial(n-index_s2)))
     )
}
fun_D <- function(x,y,n,beta,k) {
  index_s <- seq.int(1,k,1)
  index_m <- seq.int(0,n-1,1)

  factorial(n) * sum(vapply(X = index_s,
                            FUN = function(s) {
                              (((-1)^(s+1))/(factorial(s)*factorial(n-s))) * sum(
                                ((beta^index_m)/factorial(index_m)) * (
                                  (((y^index_m)*exp(-beta*k*y))*(((k-s)^index_m) - (((1-s/k)^(n-1))*(k^index_m)))) +
                                    ( exp(-beta*x) * ((((x^index_m) * ((1-s/k)^(n-1)))) - (x-s*y)^index_m) )
                                  )
                              )
                            },
                            FUN.VALUE = numeric(1)))
}

