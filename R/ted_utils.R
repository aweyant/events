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
#' @param event_magnitude_conditional_cdf a function; the CDF of event magnitudes
#' conditioned on durations, e.g. psmp_marginal_x()
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
pteddt <- function(q = NULL,
                   threshold,
                   event_duration_marginal_pmf,
                   event_duration_marginal_pmf_args = NULL,
                   event_length_arg_name = "n",
                   event_magnitude_conditional_cdf,
                   event_magnitude_conditional_cdf_args,
                   lower.tail = TRUE,
                   tol = 1e-6,
                   max_N = 10*(2^9)) {
  evaluate_cdf <- FALSE
  if(!is.null(q)) {evaluate_cdf <- TRUE}

  if(!is.null(event_duration_marginal_pmf_args) && !is.null(q) && !is.null(event_magnitude_conditional_cdf)) {
    max_N <- determine_stopping_point(event_duration_marginal_pmf,
                                      event_duration_marginal_pmf_args,
                                      tol = tol,
                                      max_N = max_N)
  }

  derived_total_cdf <- function(q,
                                threshold,
                                event_duration_marginal_pmf_args,
                                event_magnitude_conditional_cdf_args,
                                lower.tail = TRUE) {
    threshold <- unname(threshold)

    set_duration_cdfs_l <- lapply(X = seq(1,max_N,1),
                                  FUN = function(n) {
                                    event_duration_marginal_pmf_args$x <- n
                                    psdeddt(q = q,
                                            length = n,
                                            threshold = threshold,
                                            event_magnitude_conditional_cdf = event_magnitude_conditional_cdf,
                                            event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
                                            event_length_arg_name = event_length_arg_name,
                                            lower.tail = lower.tail,
                                            return_NA_outside_support = FALSE)
                                  })

    nep <- sum(vapply(X = seq(1,max_N,1),
                      FUN = function(n) {
                        event_duration_marginal_pmf_args$x <- n
                        set_duration_cdfs_l[[n]] * do.call(what = event_duration_marginal_pmf,
                                                           args = event_duration_marginal_pmf_args)
                      },
                      FUN.VALUE = numeric(1)))

    # nep <- sum(vapply(X = seq(1,N,1),
    #                   FUN = function(n) {
    #                     event_duration_marginal_pmf_args$x <- n
    #                     psdeddt(q = q,
    #                             length = n,
    #                             threshold = threshold,
    #                             event_magnitude_conditional_cdf = event_magnitude_conditional_cdf,
    #                             event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
    #                             event_length_arg_name = event_length_arg_name,
    #                             lower.tail = TRUE,
    #                             return_NA_outside_support = FALSE) *
    #                       do.call(what = event_duration_marginal_pmf,
    #                               args = event_duration_marginal_pmf_args)
    #                   },
    #                   FUN.VALUE = numeric(1)))
    if(lower.tail) {return(nep)}
    else {return(1-nep)}
  }
  if(evaluate_cdf) {
    return(derived_total_cdf(q = q,
                             threshold = threshold,
                             event_duration_marginal_pmf_args = event_duration_marginal_pmf_args,
                             event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
                             lower.tail = lower.tail))
  }
  return(derived_total_cdf)
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
psdeddt <- function(q = NULL,
                    length = 1,
                    threshold,
                    event_magnitude_conditional_cdf,
                    event_magnitude_conditional_cdf_args,
                    event_length_arg_name = "n",
                    lower.tail = TRUE,
                    return_NA_outside_support = FALSE) {
  evaluate_cdf <- FALSE
  if(!is.null(q)) {evaluate_cdf <- TRUE}

  # Given an actual event total, duration, and threshold, calculate corresponding event magnitude
  derived_total_cdf <- function(q = NULL,
                                length = 1,
                                threshold,
                                event_magnitude_conditional_cdf_args,
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

    event_derived_magnitude <- q - (length * threshold)
    outside_support <- FALSE
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
  if(evaluate_cdf) {
    return(derived_total_cdf(q = q,
                             length = length,
                             threshold = threshold,
                             event_magnitude_conditional_cdf_args = event_magnitude_conditional_cdf_args,
                             lower.tail = lower.tail,
                             return_NA_outside_support = return_NA_outside_support))
  }
  else {
    return(derived_total_cdf)
  }
}


cdf_ted <- function(x = NULL, y = NULL, n = NULL,
                    event_duration_marginal_pmf,
                    event_duration_marginal_pmf_args = NULL,
                    # event_length_arg_name = "n",
                    event_bivariate_conditional_cdf,
                    event_bivariate_conditional_cdf_args = NULL) {
  evaluate_cdf <- FALSE
  if(!is.null(x) && !is.null(y) && !is.null(n)) {evaluate_cdf <- TRUE}

  trivariate_cdf <- function(x = NULL, y = NULL, n = NULL,
                             event_duration_marginal_pmf_args,
                             event_bivariate_conditional_cdf_args) {
    event_bivariate_conditional_cdf_args$x <- x
    event_bivariate_conditional_cdf_args$y <- y
    nep <- sum(vapply(X = seq(1, n, 1),
                      FUN = function(k) {
                        event_bivariate_conditional_cdf_args$n <- k
                        event_duration_marginal_pmf_args$x <- k

                        do.call(what = event_bivariate_conditional_cdf,
                                args = event_bivariate_conditional_cdf_args) *
                          do.call(what = event_duration_marginal_pmf,
                                  args = event_duration_marginal_pmf_args)
                      },
                      FUN.VALUE = numeric(1)))
    return(nep)
  }

  if(evaluate_cdf) {
    return(trivariate_cdf(x = x, y = y, n = n,
                          event_duration_marginal_pmf_args = event_duration_marginal_pmf_args,
                          event_bivariate_conditional_cdf_args = event_bivariate_conditional_cdf_args))
  }
  return(trivariate_cdf)
}

#' @rdname ted
construct_cdf_ted <- function(args,
                              event_bivariate_conditional_cdf,
                              event_duration_marginal_pmf,
                              conditional_cdf_argument_transformations,
                              env = parent.frame()) {
  body = substitute({
    arguments = c(as.list(environment()))
    #print(arguments)
    #print(names(arguments) %in% names(formals(duration_pmf)))
    pmf_arguments = arguments[names(arguments) %in% names(formals(event_duration_marginal_pmf))]
    #print(pmf_arguments)
    cdf_arguments = arguments[names(arguments) %in% names(formals(event_bivariate_conditional_cdf))]
    #print("...")
    #print(names(arguments) %in% names(formals(conditional_cdf)))
    #print(cdf_arguments)

    if(duration == 0) {
      return(0)
    }

    sum(vapply(X = seq.int(1,duration,1),
               FUN = function(k) {

                 temp_str = paste(conditional_cdf_argument_transformations, collapse = ";")
                 transformed_cdf_arguments = within(data = cdf_arguments,
                                                    expr = {
                                                      eval(parse(text = temp_str))
                                                    })

                 pmf_arguments$x <- k

                 #print(pmf_arguments)
                 #print(".....")
                 #print(transformed_cdf_arguments)

                 do.call(what = event_duration_marginal_pmf,
                         args = pmf_arguments) *
                   do.call(what = event_bivariate_conditional_cdf,
                           transformed_cdf_arguments)
               },
               FUN.VALUE = numeric(1)))
  })
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
}

#' @rdname ted
construct_pted <- function(args,
                           pbed,
                           event_duration_marginal_pmf,
                           pbed_argument_transformations,
                           env = parent.frame(),
                           tol = 1e-6,
                           max_N = 10*(2^9)) {
  body = substitute({
    arguments = c(as.list(environment()))
    #print(arguments)
    #print(names(arguments) %in% names(formals(duration_pmf)))
    pmf_arguments = arguments[names(arguments) %in% names(formals(event_duration_marginal_pmf))]
    #print(pmf_arguments)
    pbed_arguments = arguments[names(arguments) %in% names(formals(pbed))]
    #print("...")
    #print(names(arguments) %in% names(formals(conditional_cdf)))
    #print(pbed_arguments)
    if(n_lower < 1){
      n_lower <- 1
    }
    if(is.infinite(n_upper)) {
      n_upper <- determine_stopping_point(event_duration_marginal_pmf = event_duration_marginal_pmf,
                                          event_duration_marginal_pmf_args = pmf_arguments,
                                          tol = tol,
                                          max_N = max_N)
    }

    sum(vapply(X = seq.int(n_lower, n_upper, 1),
               FUN = function(k) {

                 temp_str = paste(pbed_argument_transformations, collapse = ";")
                 transformed_pbed_arguments = within(data = pbed_arguments,
                                                    expr = {
                                                      eval(parse(text = temp_str))
                                                    })

                 pmf_arguments$x <- k

                 do.call(what = event_duration_marginal_pmf,
                         args = pmf_arguments) *
                   do.call(what = pbed,
                           transformed_pbed_arguments)
               },
               FUN.VALUE = numeric(1)))
  })
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
}

construct_pteddt <- function(args, # q, threshold, lower.tail, log.p
                             event_magnitude_conditional_cdf,
                             event_duration_marginal_pmf,
                             event_magnitude_conditional_cdf_argument_transformations,
                             env = parent.frame(),
                             tol = 1e-6,
                             max_N = 10*(2^9)) {
  body = substitute({
    arguments = c(as.list(environment()))
    #print(arguments)
    #print(names(arguments) %in% names(formals(duration_pmf)))
    pmf_arguments = arguments[names(arguments) %in% names(formals(event_duration_marginal_pmf))]
    #print(pmf_arguments)
    event_magnitude_conditional_cdf_arguments = arguments[names(arguments) %in% names(formals(event_magnitude_conditional_cdf))]
    #print("...")
    #print(names(arguments) %in% names(formals(conditional_cdf)))
    #print(event_magnitude_conditional_cdf_arguments)

    N <- determine_stopping_point(event_duration_marginal_pmf = event_duration_marginal_pmf,
                                  event_duration_marginal_pmf_args = pmf_arguments,
                                  tol = tol,
                                  max_N = max_N)

    nep <- sum(vapply(X = seq.int(1, N, 1),
                      FUN = function(k) {

                        temp_str = paste(event_magnitude_conditional_cdf_argument_transformations, collapse = ";")
                        transformed_event_magnitude_conditional_cdf_arguments = within(data = event_magnitude_conditional_cdf_arguments,
                                                                                       expr = {
                                                                                         eval(parse(text = temp_str))
                                                                                         q = q - (k*threshold)
                                                                                         lower.tail = TRUE
                                                                                       })
                        pmf_arguments$log <- FALSE
                        pmf_arguments$x <- k
                        if(transformed_event_magnitude_conditional_cdf_arguments$q <= 0) {return(0)}
                        else{
                          do.call(what = event_duration_marginal_pmf,
                                  args = pmf_arguments) *
                            do.call(what = event_magnitude_conditional_cdf,
                                    transformed_event_magnitude_conditional_cdf_arguments)
                        }
                      },
                      FUN.VALUE = numeric(1)))
    return_prob = nep
    if(!lower.tail) {return_prob <- 1 - return_prob}
    if(log) {return_prob <- log(return_prob)}
    return(return_prob)
  })
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
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
#' determine_stopping_point(dhgeom, list(prob_q = 0.5, prob_p = 0.5), max_N = 200)
#' @export
determine_stopping_point <- function(event_duration_marginal_pmf,
                                     event_duration_marginal_pmf_args,
                                     tol = 1e-6,
                                     max_N = 10*(2^9)) {
  # make sure that given lower.tail arguments are non-contradictory
  # if(!is.null(event_duration_marginal_pmf_args$lower.tail)){
  #   message(paste0("lower.tail argument is provided for ",
  #                  as.character(substitute(test_fun)),
  #                  ". This argument will be ignored.\n If you want to use lower.tail, pass it to psdeddt() directly"))
  #   event_duration_marginal_pmf_args <- event_duration_marginal_pmf_args[event_duration_marginal_pmf_args != "lower.tail"]
  # }
  max_N_exceeded <- FALSE
  N <- 10

  suppressWarnings({
    event_duration_marginal_pmf_args <- within(event_duration_marginal_pmf_args, rm("x"))
    event_duration_marginal_pmf_args <- within(event_duration_marginal_pmf_args, x <- 1:N)
    event_duration_marginal_pmf_args <- within(event_duration_marginal_pmf_args, rm("log"))
  })

  while(1 - sum(
    do.call(what = event_duration_marginal_pmf,
            args = event_duration_marginal_pmf_args)
    # eval(
    #   parse(
    #     text = paste0("vapply(X = seq(1,N,1),FUN = event_duration_marginal_pmf,FUN.VALUE = numeric(1),",paste(names(event_duration_marginal_pmf_args), "=", unlist(event_duration_marginal_pmf_args), collapse = ","), ")")
    #     )
    #   )
    ) > tol) {
    if(N > max_N) {
      max_N_exceeded <- TRUE
      break
    }
    N <- 2*N
    #print(paste0("Here at N = ", N))
    event_duration_marginal_pmf_args$x <- 1:N
    #print(paste0("pmf_args are now ", event_duration_marginal_pmf_args))
  }
  if(max_N_exceeded) {
    message("max_N exceeded. Check your event duration marginal PMF or set max_N to a higher value.")
    stop()
  }
  else{
    return(N)
  }
}
