#' **P**robability of **T**rivariate **E**vent **D**istribution **D**erived
#' **T**otal (P TED DT) and **P**robability of **S**et **D**uration **E**vent
#' **D**istribution **D**erived **T**otal
#'
#' This evaluates the CDF of the total accumulation of an amount over an event.
#'
#' pteddt() by way of example:
#' a 3 day precipitation event might have a magnitude of 430m. If
#' the threshold for defining the event was 75mm/day, then the actual total
#' precipitation over this event was 655mm \[430 + 3*75\].
#'
#' Note that durations and daily accumulations are, in general, both random
#' variables, so if we are interested in probabilties (return intervals) of
#' uninterrupted precipitation accumulations, we must construct an event total
#' CDF from (at least) two components: the (discrete) PMF describing event
#' durations and the CDF for magnitudes, conditioned upon a duration. The Law of
#' Total Probability lets us create a CDF for the event total, which is quite
#' meaningful if we are studying precipitation.
#'
#' psdeddt() evaluates the CDF for event totals of a given duration (hence, "set
#' duration").
#'
#' @inheritParams determine_stopping_point
#' @param q a number; the event total, including the below-threshold portion
#' @param threshold a number; the minimum threshold used to define events
#' @param event_length_arg_name a string; the name of the argument used for the
#' duration of an event in the CDF of the magnitude conditioned duration. For
#' example, this would be "n" for psmp_marginal_x()
#' @param event_magnitude_conditional_cdf a function; the CDF of event magnitudes
#' conditioned on durations, e.g. psmp_marginal_x()
#' @param event_magnitude_conditional_cdf_args a list; other arguments, such as
#' parameters of the distribution
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @param length integer; the duration of the event
#' @param return_NA_outside_support logical; should an unsupported event cause
#' a return of NA or 0? Decide here.
#'
#' @return a number in the interval \[0,1\], the probability of not exceeding or
#' exceeding q
#'
#' @name ted
NULL

#' @rdname ted
#' @examples
#' q = 100; length = 2; threshold = 20
#' alpha = 0.02; beta = 0.04
#' prob_p = 0.6; prob_q = 0.7
#' event_length_arg_name = "n"
#'
#' pteddt(q = q,
#' threshold = threshold,
#' event_duration_marginal_pmf = dhgeom,
#' event_duration_marginal_pmf_args = list(prob_p = prob_p, prob_q = prob_q),
#' event_length_arg_name = event_length_arg_name,
#' event_magnitude_conditional_cdf = psmp_marginal_x,
#' event_magnitude_conditional_cdf_args = list(alpha = alpha, beta = beta),
#' lower.tail = TRUE,
#' tol = 10^(-6),
#' max_N = 100)
#' @export
pteddt <- function(q,
                   threshold,
                   event_duration_marginal_pmf,
                   event_duration_marginal_pmf_args,
                   event_length_arg_name = "n",
                   event_magnitude_conditional_cdf,
                   event_magnitude_conditional_cdf_args,
                   lower.tail = TRUE,
                   tol = 10^(-6),
                   max_N = 100) {
  threshold <- unname(threshold)
  N <- determine_stopping_point(event_duration_marginal_pmf,
                                event_duration_marginal_pmf_args,
                                tol = tol,
                                max_N = max_N)
  nep <- sum(vapply(X = seq(1,N,1),
                    FUN = function(n) {
                      event_duration_marginal_pmf_args$x <- n
                      psdeddt(q = q,
                              length = n,
                              threshold = threshold,
                              event_magnitude_conditional_cdf = event_magnitude_conditional_cdf,
                              event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
                              event_length_arg_name = event_length_arg_name,
                              lower.tail = TRUE,
                              return_NA_outside_support = FALSE) *
                        do.call(what = event_duration_marginal_pmf,
                                args = event_duration_marginal_pmf_args)
                    },
                    FUN.VALUE = numeric(1)))
  if(lower.tail) {return(nep)}
  else {return(1-nep)}
}

#' @rdname ted
#' @examples
#' rsmp(n = 3, m = c(3,4,2), alpha = 0.03, beta = 1/100)
#' \dontrun{
#'   n        x        y
#' 1 3 241.6733 157.2649
#' 2 4 357.7778 187.8859
#' 3 2 219.1081 140.8779
#' }
#' @export
psdeddt <- function(q,
                    length = 1,
                    threshold,
                    event_magnitude_conditional_cdf,
                    event_magnitude_conditional_cdf_args,
                    event_length_arg_name = "n",
                    lower.tail = TRUE,
                    return_NA_outside_support = FALSE) {
  threshold <- unname(threshold)
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
  # Given an actual event total, duration, and threshold, calculate corresponding event magnitude
  event_derived_magnitude <- q - (length * threshold)

  if(event_derived_magnitude <= 0) { # check if event total is possible
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
#' @param event_duration_marginal_pmf a function; the marginal CDF for event durations
#' @param event_duration_marginal_pmf_args a list; arguments for the marginal duration CDF
#' @param tol a very small number; minimal exceedance probability desired
#' @param max_N a number; highest N to search
#'
#' @return
#' an integer
#'
#' @examples
#' determine_stopping_point(pgeom, list(prob = 0.1), max_N = 200)
#' @export
determine_stopping_point <- function(event_duration_marginal_pmf,
                                     event_duration_marginal_pmf_args,
                                     tol = 10^(-6),
                                     max_N = 160) {
  # make sure that given lower.tail arguments are non-contradictory
  # if(!is.null(event_duration_marginal_pmf_args$lower.tail)){
  #   message(paste0("lower.tail argument is provided for ",
  #                  as.character(substitute(test_fun)),
  #                  ". This argument will be ignored.\n If you want to use lower.tail, pass it to psdeddt() directly"))
  #   event_duration_marginal_pmf_args <- event_duration_marginal_pmf_args[event_duration_marginal_pmf_args != "lower.tail"]
  # }
  max_N_exceeded <- FALSE
  N <- 10

  event_duration_marginal_pmf_args <- within(event_duration_marginal_pmf_args, rm("x"))
  event_duration_marginal_pmf_args <- within(event_duration_marginal_pmf_args, rm("log"))

  while(1 - sum(
    eval(
      parse(
        text = paste0("vapply(X = seq(1,N,1),FUN = event_duration_marginal_pmf,FUN.VALUE = numeric(1),",paste(names(event_duration_marginal_pmf_args), "=", unlist(event_duration_marginal_pmf_args), collapse = ","), ")")
        )
      )
    ) > tol) {
    if(N > max_N) {
      max_N_exceeded <- TRUE
      break
    }
    N <- 2*N
  }
  if(max_N_exceeded) {
    message("max_N exceeded. Check your event duration marginal PMF or set max_N to a higher value.")
    stop()
  }
  else{
    return(N)
  }
}
