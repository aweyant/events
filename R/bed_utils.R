construct_pbed <- function(args, cdf_bed, env = parent.frame()) {
  body = substitute({
    arguments <- c(as.list(environment()))
    #print(arguments)
    cdf_param_args <- arguments[!(names(arguments) %in% c("x_lower", "y_lower",
                                                          "x_upper", "y_upper" , "n"))]
    if(n <= 0) {
      return(0)
    }
    if(x_upper <= y_lower) {return(0)}
    if(n == 1) {
      if(y_upper <= x_lower) { #1 *
        # args1 <- c(cdf_param_args,
        #            x = arguments$x_upper,
        #            y = arguments$x_upper,
        #            n = n)
        # args2 <- c(cdf_param_args,
        #            x = arguments$x_lower,
        #            y = arguments$x_lower,
        #            n = n)
        #print("case 1")
        return(0)
      }
      if(x_lower < y_upper && y_upper <= x_upper) { #2 *
        args1 <- c(cdf_param_args,
                   x = arguments$y_upper,
                   y = arguments$y_upper,
                   n = n)
        args2 <- c(cdf_param_args,
                   x = arguments$x_lower,
                   y = arguments$x_lower,
                   n = n)
      }
      if(y_lower <= x_lower && x_upper <= y_upper) { #3 *
        args1 <- c(cdf_param_args,
                   x = arguments$x_upper,
                   y = arguments$x_upper,
                   n = n)
        args2 <- c(cdf_param_args,
                   x = arguments$x_lower,
                   y = arguments$x_lower,
                   n = n)
      }
      if(x_lower <= y_lower && y_upper <= x_upper) { #4 *
        args1 <- c(cdf_param_args,
                   x = arguments$y_upper,
                   y = arguments$y_upper,
                   n = n)
        args2 <- c(cdf_param_args,
                   x = arguments$y_lower,
                   y = arguments$y_lower,
                   n = n)
      }
      if(x_lower <= y_lower && x_upper <= y_upper) { #5 *
        args1 <- c(cdf_param_args,
                   x = arguments$x_upper,
                   y = arguments$x_upper,
                   n = n)
        args2 <- c(cdf_param_args,
                   x = arguments$y_lower,
                   y = arguments$y_lower,
                   n = n)
      }
      # args1 <- c(cdf_param_args,
      #            x = arguments$x_upper,
      #            y = arguments$y_upper,
      #            n = n)
      # args2 <- c(cdf_param_args,
      #            x = arguments$x_lower,
      #            y = arguments$y_lower,
      #            n = n)
      #print(args1)
      return(
        do.call(what = cdf_bed, args = args1) +
          -do.call(what = cdf_bed, args = args2)
      )
    }
    else {
      args1 <- c(cdf_param_args,
                 x = arguments$x_upper,
                 y = arguments$y_upper,
                 n = n)
      args2 <- c(cdf_param_args,
                 x = arguments$x_lower,
                 y = arguments$y_upper,
                 n = n)
      args3 <- c(cdf_param_args,
                 x = arguments$x_upper,
                 y = arguments$y_lower,
                 n = n)
      args4 <- c(cdf_param_args,
                 x = arguments$x_lower,
                 y = arguments$y_lower,
                 n = n)
      return(
        do.call(cdf_bed, args1) +
          -do.call(cdf_bed, args2) +
          -do.call(cdf_bed, args3) +
          do.call(cdf_bed, args4)
      )
    }
  })
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
}
