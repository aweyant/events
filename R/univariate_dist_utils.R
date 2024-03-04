construct_quantile <- function(args,
                               cdf,
                               cdf_support,
                               cdf_search_interval,
                               tol = 1e-6,
                               env = parent.frame()) {
  body = substitute({
    arguments = c(as.list(environment()))
    #params_and_other_cdf_args = arguments[which(names(arguments) != "p")]
    #print(params_and_other_cdf_args)

    # check if probabiltiy given is valid
    if(p < 0) {
      warning("Nans produced")
      return(NaN)
    }
    if(p > 1) {
      warning("Nans produced")
      return(NaN)
    }
    # set up search interval extension argument
    if(min(cdf_support) < min(cdf_search_interval) && max(cdf_support) > max(cdf_search_interval)) {
      extendInt = "yes"
    }
    else if(min(cdf_support) < min(cdf_search_interval)) {
       extendInt = "downX"
    }
    else if(max(cdf_support) > max(cdf_search_interval)) {
      extendInt = "upX"
    }
    else {
      extendInt = "no"
    }


    # should probably also check if the tolerance is suitable


    do.call(what = uniroot,
            args = list(f = function(q, ...) {
              all_args = list(...)[[1]]
              #print(all_args)
              other_cdf_args = all_args[which(names(all_args) != "p")]
              #print(other_cdf_args)
              arg_p = all_args[which(names(all_args) == "p")]
              #print("...")
              #print(arg_p)
              #print(c(q = q, other_cdf_args))
              do.call(what = cdf, args = c(q = q, other_cdf_args)) - arg_p$p
              },
              interval = cdf_search_interval,
              tol = tol,
              extendInt = "yes",
              arguments))$root

    # uniroot(f = function(p, ...) {
    #   cdf(...) - p},
    #   interval = cdf_support,
    #   tol = tol,
    #   extendInt = "yes",
    #   ... = ...)$root

    # uniroot(f = phgsmpdt_offset,
    #         interval = c(0,100),
    #         tol = 0.001,
    #         extendInt = "up",
    #         p = p,
    #         threshold = threshold,
    #         prob_p = prob_p, prob_q = prob_q,
    #         alpha = alpha, beta = beta)$root
  })
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
}
