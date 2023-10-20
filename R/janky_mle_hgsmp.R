#' Maximum Likelihood Estimator Calculation for the Hurdle-Geometric Sums &
#' Maxima of Pareto Vectors distribution
#'
#' @param magnitudes sums of Pareto II vectors
#' @param durations number of elements in Pareto II vectors
#'
#' @return
#' A list with elements:
#' \itemize{
#' \item estimate raw ML estimates of the four HGSMP parameters}
#' @export
#'
#' @examples
#' random_events_df <- rsmp(n = 100, m = c(rep(1,70), rep(2,22), rep(3,7), 4),
#' alpha = 0.03, beta = 1/100)
#' janky_mle_hgsmp(random_events_df$x, random_events_df$n)
janky_mle_hgsmp <- function(magnitudes, durations) {
  # ESTIMATE p & q
  p = janky_mle_hgeom_prob_p(durations)
  q = janky_mle_hgeom_prob_q(durations)

  # check to see if GTETLG is immediately the best choice
  if(fun_Sm(magnitudes,durations) < 0){
    # return(list(p = p, q = q, beta = sum(durations)/sum(magnitudes), alpha = 0))
    return(list(estimate = c(prob_p = p, prob_q = q,
                             alpha = 0, beta = sum(durations)/sum(magnitudes))))
  }

  # numerically estimate general hgsmp params
  else {
    alpha = alpha_est(magnitudes, durations)
    # print(alpha)
    beta = beta_est(magnitudes, durations, alpha = alpha)
    # return(list(p = p, q = q, beta = beta, alpha = alpha))
    return(list(estimate = c(prob_p = p, prob_q = q,
                             alpha = alpha,beta = beta)))
  }
}

#' @describeIn janky_mle_hgsmp Helper function to estimate alpha
#' @param num_search logical; whether or not to initiate a numerical search for
#' alpha. This is mostly for internal use.
#' @export
alpha_est <- function(magnitudes, durations, num_search = TRUE) {
  if(!num_search) {
    return(0)
  }
  else{
    exp_sample <- FALSE
    tryCatch({upper = stats::uniroot(f = fun_v_minus_v_at_0,
                              lower = 0.01,
                              upper = 1,
                              extendInt = "yes",
                              magnitudes = magnitudes, durations = durations)$root * 1.1},
             error = function(e){exp_sample <<- TRUE},
             warning = function(w){exp_sample <<- TRUE})
    # print(exp_sample)
    if(exp_sample || upper < 0) {
      return(0)
    }

    # print(upper)
    stats::optimize(f = fun_v,
             lower = 0,
             upper = upper,
             maximum = TRUE,
             magnitudes = magnitudes, durations = durations)$maximum
  }
}

#' @describeIn janky_mle_hgsmp Helper function to estimate beta
#'
#' @param alpha optimal alpha used to calculate beta. This is only used by beta_est()
#' because beta's optimal value depends on the optimal value of for alpha. If you
#' estimate both alpha and beta, there is no reason to repeat the first step.
#' @export
beta_est <- function(magnitudes, durations,
                     alpha = NULL) {
  if(is.null(alpha)) {
    alpha = alpha_est(magnitudes, durations, num_search = TRUE)
  }
  if(alpha == 0) {
    return(sum(durations)/sum(magnitudes))
  }
  else {
    exp_sample <- FALSE
    tryCatch({beta_hat = stats::uniroot(f = fun_h_minus_1,
                                 lower = .Machine$double.eps^0.25,
                                 upper = 2 * sum(durations)/sum(magnitudes),
                                 extendInt = "upX",
                                 alpha = alpha,
                                 magnitudes = magnitudes,
                                 durations = durations)$root},
             error = function(e){exp_sample <<- TRUE},
             warning = function(w){exp_sample <<- TRUE})
    # print(exp_sample)
    if(exp_sample || beta_hat < 0) {
      return(sum(durations)/sum(magnitudes))
    }
    return(beta_hat)
  }
}

fun_Sm <- function(magnitudes, durations) {
  ((sum(magnitudes^2) / ((sum(magnitudes))^2)) +
     (sum(durations^2)/((sum(durations))^2)) -
     1/sum(durations) -
     2 * sum(magnitudes*durations)/(sum(durations)*sum(magnitudes)))
}

fun_v_minus_v_at_0 <- function(alpha, magnitudes, durations) {
  if(alpha == 0) {
    0
  }
  else{
    beta_argmax = stats::uniroot(f = fun_h_minus_1,
                          lower = .Machine$double.eps^0.25, upper = 1000,
                          extendInt = "upX",
                          alpha = alpha, magnitudes = magnitudes,
                          durations = durations)$root
    fun_g(alpha, beta_argmax, magnitudes, durations) + mean(durations)*(log(mean(magnitudes)/mean(durations)) + 1)

  }
}

fun_h_minus_1 <- function(beta, alpha, magnitudes, durations) {
  (1/(sum(durations))) * sum((1 + alpha*durations)*(beta*magnitudes)/(1 + alpha*beta*magnitudes)) - 1
}

fun_g <- function(alpha, beta, magnitudes, durations) {
  mean(sapply(X = durations,
              FUN = function(N) {
                #print(sum(log(1 + alpha * (0:(N-1)))))
                sum(log(1 + alpha * (0:(N-1))))
              })) - mean(
                sapply(X = seq_along(durations),
                       FUN = function(i){
                         #print((durations[i] + 1/alpha) * log(1 + alpha * beta * magnitudes[i]))
                         (durations[i] + 1/alpha) * log(1 + alpha * beta * magnitudes[i])
                       })
              ) + mean(durations)*log(beta)
}

fun_v <- function(alpha, magnitudes, durations) {
  if(alpha <= .Machine$double.eps^0.25) {
    -mean(durations) * (log(mean(magnitudes)/mean(durations)) + 1)
  }
  else {
    # beta_argmax = (newtonRaphson(fun = fun_h_minus_1,
    #                              x0 = sum(durations)/sum(magnitudes), # first guess
    #                              dfun = fun_dh_dbeta,
    #                              alpha = alpha, magnitudes = magnitudes, durations = durations))$root
    #print(beta_argmax)
    tryCatch({beta_argmax = stats::uniroot(f = fun_h_minus_1,
                                    lower = .Machine$double.eps^0.25,
                                    upper = 1,
                                    extendInt = "upX",
                                    alpha = alpha, magnitudes = magnitudes, durations = durations)$root},
             error = function(e) {
               beta_argmax = mean(durations)/mean(magnitudes)
             },
             warning = function(w) {
               beta_argmax = mean(durations)/mean(magnitudes)
             })

    fun_g(alpha, beta_argmax, magnitudes, durations)
  }
}
