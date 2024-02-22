#' Multivariate Pareto II Vectors
#'
#' @inheritParams rsmp
#'
#' @return
#' A list of random Pareto II vectors
#' @export
#'
#' @examples
#' rmvlomax(n = 3, m = c(3,4,2), alpha = 0.03, beta = 1/100)
#' \dontrun{
#' [[1]]
#' [1] 170.56225 120.49477  14.40258
#'
#' [[2]]
#' [1] 192.361464 274.643827  64.666721   9.729169
#'
#' [[3]]
#' [1] 108.0837 255.7757
#' }
rmvlomax <- function(n,
                     m, # length of each vector
                     alpha,
                     beta,
                     rounding = FALSE) {
  if(length(m) == 0) {
    warning("No vector lengths provided - assuming that each observations is
            meant to be scalar.")
    m <- 1
  }
  else if(length(m) == 1) {
    message("A single length m was provided - assuming each observation is meant
    to be this length.")
  }
  else if(n != length(m) & length(m) != 1) {
    warning("Warning: Requested sample size is not equal to the length of the
            vector which specifies the length of each observation.")
    n <- length(m)
  }

  lapply(X = seq_len(n),
         FUN = function(i) {
           stats::rexp(n = m[min(length(m),i)], rate = beta) %>%
             {if(alpha > 0) ./stats::rgamma(n = 1, shape = 1/alpha, scale = alpha) else{.}} %>%
             {if(rounding){round(.)} else {.}}
         })
}
