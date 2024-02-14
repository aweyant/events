psdeddt <- function(q, # prob set duration event distribution derived total
                    length = 1,
                    threshold,
                    event_magnitude_conditional_cdf,
                    event_magnitude_conditional_cdf_args,
                    event_length_arg_name = "n",
                    lower.tail = TRUE,
                    return_NA_outside_support = FALSE) {
  # length is passed to the proper arg for event_magnitude_conditional_cdf, as
  #' its name cannot be anticipated
  if(!is.null(eval(str2expression(paste0("event_magnitude_conditional_cdf_args$", event_length_arg_name))))) {
    message("Your provided argment name for length conflicts with another
    argument of the event magnitude marginal CDF. The argument for the numerical
    argument for event length should be provided directly to psdeddt, as should
    its name as it is given to the event magnitude marginal CDF/. If this error
    is not due to a mistake on your part, you might need to write a wrapper
    function for your magitude marginal CDF.")
    stop()
  }
  eval(str2expression(paste0("event_magnitude_conditional_cdf_args$", event_length_arg_name,"<-",length)))

  # make sure that given lower.tail arguments are non-contradictory
  if(!is.null(event_magnitude_conditional_cdf_args$lower.tail)){
    message(paste0("lower.tail argument is provided for ",
                   as.character(substitute(test_fun)),
                   ". This argument will be ignored.\n If you want to use lower.tail, pass it to psdeddt() directly"))
    event_magnitude_conditional_cdf_args <- event_magnitude_conditional_cdf_args[event_magnitude_conditional_cdf_args != "lower.tail"]
  }

  outside_support <- FALSE
  # Given an actual event total, duration, and threshold, calculate
  #' corresponding event magnitude
  event_derived_magnitude <- q - (length * threshold)

  if(event_derived_magnitude < 0) { # check if event total is possible
    outside_support <- TRUE
  }

  nep <- do.call(what = event_magnitude_conditional_cdf,
                 args = c(q = event_derived_magnitude,
                          lower.tail = TRUE,
                          event_magnitude_conditional_cdf_args))
  if(outside_support) {
    if(return_NA_outside_support) {return(NaN)}
    else {return(0)}
  }
  else if(lower.tail) {
    return(nep)
  }
  else{
    return(1 - nep)
  }
}

pteddt <- function(q,
                   threshold,
                   event_duration_marginal_cdf,
                   event_duration_marginal_cdf_args,
                   event_length_arg_name = "n",
                   event_magnitude_conditional_cdf,
                   event_magnitude_conditional_cdf_args,
                   lower.tail = TRUE,
                   tol = 10^(-6),
                   max_N = 100) {
  N <- determine_stopping_point(event_duration_marginal_cdf,
                                event_duration_marginal_cdf_args,
                                tol = tol,
                                max_N = max_N)
  nep <- sum(vapply(X = seq(1,N,1),
                    FUN = function(n) {
                      event_duration_marginal_cdf_args$q <- n
                      psdeddt(q = q,
                              length = n,
                              threshold = threshold,
                              event_magnitude_conditional_cdf = event_magnitude_conditional_cdf,
                              event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
                              event_length_arg_name = event_length_arg_name,
                              lower.tail = TRUE,
                              return_NA_outside_support = FALSE) *
                        do.call(what = event_duration_marginal_cdf,
                                args = event_duration_marginal_cdf_args)
                    },
                    FUN.VALUE = numeric(1)))
  if(lower.tail) {return(nep)}
  else {return(1-nep)}
}

#' Find a quantile of a discrete distribution to achieve a low exceedance probability
#'
#' This is a utility function used when calculating marginal distributions of
#' magnitudes and maxima of trivariate events. These marginal distributions are
#' generally calculated as sums of CDFs conditioned on duration multiplied by
#' probabilities of durations. Distributions of durations of all events we
#' anticipate seeing will be decreasing after a certain duration. This function
#' finds the duration we need to calculate the marginal probabilities to an
#' arbitrary precision (the "tol").
#'
#' @param event_duration_marginal_cdf a function; the marginal CDF for event durations
#' @param event_duration_marginal_cdf_args a list; arguments for the marginal duration CDF
#' @param tol a very small number; minimal exceedance probability desired
#' @param max_N a number; highest N to search
#'
#' @return
#' an integer
#' @export
#'
#' @examples
#' determine_stopping_point(pgeom, list(prob = 0.1), max_N = 200)
determine_stopping_point <- function(event_duration_marginal_cdf,
                                     event_duration_marginal_cdf_args,
                                     tol = 10^(-6),
                                     max_N = 160) {
  # make sure that given lower.tail arguments are non-contradictory
  if(!is.null(event_duration_marginal_cdf_args$lower.tail)){
    message(paste0("lower.tail argument is provided for ",
                   as.character(substitute(test_fun)),
                   ". This argument will be ignored.\n If you want to use lower.tail, pass it to psdeddt() directly"))
    event_duration_marginal_cdf_args <- event_duration_marginal_cdf_args[event_duration_marginal_cdf_args != "lower.tail"]
  }
  max_N_exceeded <- FALSE
  N <- 10
  event_duration_marginal_cdf_args$q <- N
  event_duration_marginal_cdf_args$lower.tail <- TRUE
  while(1 - do.call(what = event_duration_marginal_cdf, args = event_duration_marginal_cdf_args) > tol) {
    if(N > max_N) {
      max_N_exceeded <- TRUE
      break
    }
    N <- 2*N
    event_duration_marginal_cdf_args$q <- N
  }
  if(max_N_exceeded) {
    message("max_N exceeded. Check your event duration marginal CDF or set max_N to a higher value.")
    stop()
  }
  else{
    return(N)
  }
}
