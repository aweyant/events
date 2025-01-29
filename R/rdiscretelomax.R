#' Discrete Lomax (Pareto II) Distribution
#'
#' @name discretelomax
#'
#' @param x
#' @param dlomax_alpha
#' @param dlomax_prob_p
#'
#' @details
#' In the output vector, it is the prop_p which vary more quickly. If you
#' provide 2 dlomax_alphas and 3 dlomax_prob_ps, the first 3 numbers in the output are
#' associated with the first dlomax_alpha.
#'
#'
#' @examples
#' \dontrun{
#' ddiscretelomax(x = 1:10, dlomax_alpha = c(0,0.2,0.4),
#' prob = c(0.6), return_type = "data.frame") %>%
#' ggplot2::ggplot(
#' ggplot2::aes(x = x, y = prob_mass, group = dlomax_alpha, fill = dlomax_alpha)) +
#' ggplot2::geom_col(position = "dodge") +
#' ggplot2::scale_fill_binned() +
#' ggplot2::scale_y_continuous(trans = "log")
#'
#' pdiscretelomax(q = seq(1:15),
#' dlomax_alpha = c(0,0.05,0.1,0.2),
#' dlomax_prob_p = c(0.5),
#' return_type = "data.frame", lower.tail = F)  %>%
#' ggplot2::ggplot(ggplot2::aes(x = q,
#' y = non_exceedance_prob, group = dlomax_alpha, color = dlomax_alpha)) +
#' ggplot2::geom_path() +
#' ggplot2::scale_x_continuous(breaks = 1:15) +
#' ggplot2::scale_y_continuous(trans = "log10")
#'
#' qdiscretelomax(p = seq(0.54,0.99, by = 0.05),
#' dlomax_alpha = c(0,0.1,0.3),
#' dlomax_prob_p = 0.5,
#' return_type = "data.frame") %>%
#' ggplot2::ggplot(ggplot2::aes(x = lte_prob,
#' y= q, color = dlomax_alpha)) +
#' ggplot2::geom_point() +
#' ggplot2::lims(y = c(0,NA))
#' }
NULL

#' @rdname discretelomax
#' @export
ddiscretelomax <- function(x, dlomax_alpha, dlomax_prob_p, return_type = "vector", log = FALSE) {
  # PRESERVE ORIGINAL Xs IN CASE BAD ARGUMENTS ARE GIVEN BUT WE STILL WANT TO
  # RETURN A DATAFRAME WITH DETAILS
  if(return_type == "data.frame") {
    x_grid_df <- rep(x, times = length(dlomax_prob_p) * length(dlomax_alpha))
    dlomax_alpha_grid_df <- rep(dlomax_alpha, each = length(x) * length(dlomax_prob_p))
    dlomax_prob_p_grid_df <- rep(dlomax_prob_p, times = length(x) * length(dlomax_alpha))
  }

  # CHECK X INPUT
  discretelomax_check_x()
  # CHECK PARAMETERS
  discretelomax_check_param_args()

  # Expand combinations of arguments into a grid
  x_grid <- rep(x, times = length(dlomax_prob_p) * length(dlomax_alpha))
  dlomax_alpha_grid <- rep(dlomax_alpha, each = length(x) * length(dlomax_prob_p))
  dlomax_prob_p_grid <- rep(dlomax_prob_p, times = length(x) * length(dlomax_alpha))

  # Create vector of probabilities to return
  prob_out <- numeric(length(dlomax_alpha_grid))

  # Calculate probs for the cases alpha = 0 and alpha > 0. In practice, I need
  # to fudge it a bit by writing the condition as alpha < small number or
  # alpha >= small number.
  prob_out[which(dlomax_alpha_grid <= 10e-5)] <- dgeom(
    x = x - 1,
    prob = dlomax_prob_p_grid[which(dlomax_alpha_grid <= 10e-5)])
  prob_out[which(dlomax_alpha_grid > 10e-5)] <- {
    (1 - dlomax_alpha_grid[which(dlomax_alpha_grid > 10e-5)]*(x-1)*log(1-dlomax_prob_p_grid[which(dlomax_alpha_grid > 10e-5)])
    )^(-1/dlomax_alpha_grid[which(dlomax_alpha_grid > 10e-5)]) - (
      1 - dlomax_alpha_grid[which(dlomax_alpha_grid > 10e-5)]*x*log(1-dlomax_prob_p_grid[which(dlomax_alpha_grid > 10e-5)])
    )^(-1/dlomax_alpha_grid[which(dlomax_alpha_grid > 10e-5)])
  }
  if(log) {prob_out <- log(prob_out)}
  if(return_type == "data.frame") {
    return(data.frame(
      x = x_grid_df,
      dlomax_prob_p = dlomax_prob_p_grid_df,
      dlomax_alpha = dlomax_alpha_grid_df,
      prob_mass = prob_out
    ))
  }
  else {
    return(prob_out)
  }
}

discretelomax_check_x <- function() {
  # Check x
  call_x_gt_0 <- rlang::expr(x[which(x <= 0)] <- NaN)
  rlang::eval_bare(call_x_gt_0, env = parent.frame())

  call_x_warning <- rlang::expr({
    if(any(is.nan(x))) {
      warning("x is supported only on positive integers. NaN(s) returned.")
    }
  })

  rlang::eval_bare(call_x_warning, env = parent.frame())
  NULL
}

discretelomax_check_q <- function() {
  # Check x
  call_q_gt_0 <- rlang::expr(
    q[FALSE#which(q <= 0)
      ] <- NaN)
  rlang::eval_bare(call_q_gt_0, env = parent.frame())

  call_q_warning <- rlang::expr({
    if(any(is.nan(q))) {
      warning("q is supported only on positive integers. NaN(s) returned.")
    }
  })

  rlang::eval_bare(call_q_warning, env = parent.frame())
  NULL
}

discretelomax_check_p <- function() {
  # Check p
  call_p_gte_0 <- rlang::expr(p[which(p < 0)] <- NaN)
  rlang::eval_bare(call_p_gte_0, env = parent.frame())

  call_p_lt_1 <- rlang::expr(p[which(p >= 1)] <- NaN)
  rlang::eval_bare(call_p_lt_1, env = parent.frame())

  call_p_warning <- rlang::expr({
    if(any(is.nan(p))) {
      warning("p is supported on [0,1). NaN(s) returned.")
    }
  })

  rlang::eval_bare(call_p_warning, env = parent.frame())
  NULL
}

discretelomax_check_param_args <- function() {
  # Check dlomax_alpha
  call_dlomax_alpha_gt_0 <- rlang::expr(dlomax_alpha[which(dlomax_alpha < 0)] <- NaN)
  rlang::eval_bare(call_dlomax_alpha_gt_0, env = parent.frame())

  # Check prob p
  call_dlomax_prob_p_gt_0 <- rlang::expr(dlomax_prob_p[which(dlomax_prob_p <= 0)] <- NaN)
  call_dlomax_prob_p_lt_1 <- rlang::expr(dlomax_prob_p[which(dlomax_prob_p > 1)] <- NaN)

  rlang::eval_bare(call_dlomax_prob_p_gt_0, env = parent.frame())
  rlang::eval_bare(call_dlomax_prob_p_lt_1, env = parent.frame())

  # Throw warnings, if necessary
  call_dlomax_alpha_warning <- rlang::expr({
    if(any(is.nan(dlomax_alpha))) {
      warning("dlomax_alpha must be >= 0. NaN(s) returned.")
    }
  })
  call_dlomax_prob_p_warning <- rlang::expr({
    if(any(is.nan(dlomax_prob_p))) {
      warning("dlomax_prob_p must be >0 and <1. NaN(s) returned.")
    }
  })

  rlang::eval_bare(call_dlomax_alpha_warning, env = parent.frame())
  rlang::eval_bare(call_dlomax_prob_p_warning, env = parent.frame())
  NULL
}

#' @rdname discretelomax
#' @export
rdiscretelomax <- function(n, dlomax_alpha, dlomax_prob_p) {

  discretelomax_check_param_args()

  dlomax_alpha_grid = rep_len(dlomax_alpha, length.out = n)
  dlomax_prob_p_grid = rep(dlomax_prob_p, each = length(dlomax_alpha), length.out = n)

  count_out <- numeric(n)

  count_out[dlomax_alpha_grid == 0] <- 1 + rgeom(
    n = sum(dlomax_alpha_grid == 0),
    prob = dlomax_prob_p_grid[which(dlomax_alpha_grid == 0)]
    )

  count_out[dlomax_alpha_grid > 0] <- ceiling(
    stats::rexp(n = sum(dlomax_alpha_grid > 0)) / rgamma(
      n = sum(dlomax_alpha_grid > 0),
      shape = 1/dlomax_alpha_grid[which(dlomax_alpha_grid > 0)],
      rate = 1 / (-dlomax_alpha_grid[which(dlomax_alpha_grid > 0)] * log(1 - dlomax_prob_p_grid[which(dlomax_alpha_grid > 0)])) )
  )

  count_out
}

#' @rdname discretelomax
#' @export
pdiscretelomax <- function(q, dlomax_alpha, dlomax_prob_p,
                           lower.tail = TRUE, log.p = FALSE, return_type = "vector") {
  # PRESERVE ORIGINAL Xs IN CASE BAD ARGUMENTS ARE GIVEN BUT WE STILL WANT TO
  # RETURN A DATAFRAME WITH DETAILS
  if(return_type == "data.frame") {
    q_grid_df <- rep(q, times = length(dlomax_prob_p) * length(dlomax_alpha))
    dlomax_alpha_grid_df <- rep(dlomax_alpha, each = length(q) * length(dlomax_prob_p))
    dlomax_prob_p_grid_df <- rep(dlomax_prob_p, times = length(q) * length(dlomax_alpha))
  }

  # CHECK X INPUT
  discretelomax_check_q()
  # CHECK PARAMETERS
  discretelomax_check_param_args()

  # Expand combinations of arguments into a grid
  q_grid <- rep(q, times = length(dlomax_prob_p) * length(dlomax_alpha))
  dlomax_alpha_grid <- rep(dlomax_alpha, each = length(q) * length(dlomax_prob_p))
  dlomax_prob_p_grid <- rep(dlomax_prob_p, times = length(q) * length(dlomax_alpha))

  # Create vector of probabilities to return
  prob_out <- numeric(length(dlomax_alpha_grid))

  prob_out[which(dlomax_alpha_grid == 0)] <- stats::pgeom(
    q = floor(q_grid - 1),
    prob = dlomax_prob_p_grid[which(dlomax_alpha_grid == 0)])
  prob_out[which(dlomax_alpha_grid > 0)] <- {
    1 - (1 - dlomax_alpha_grid[which(dlomax_alpha_grid > 0)]*floor(q)*log(1-dlomax_prob_p_grid[which(dlomax_alpha_grid > 0)])
    )^(-1/dlomax_alpha_grid[which(dlomax_alpha_grid > 0)])
  }

  prob_out[which(q_grid < 1)] <- 0

  if(log.p) {
    prob_out = log(prob_out)
  }
  if(!lower.tail) {
    prob_out <- 1 - prob_out
  }

  if(return_type == "data.frame") {
    return(data.frame(
      q = q_grid_df,
      dlomax_prob_p = dlomax_prob_p_grid_df,
      dlomax_alpha = dlomax_alpha_grid_df,
      non_exceedance_prob = prob_out
    ))
  }
  else {
    return(prob_out)
  }
}
#' @rdname discretelomax
#' @export
qdiscretelomax <- function(p, dlomax_alpha, dlomax_prob_p, lower.tail = TRUE, return_type = "vector") {
  if(return_type == "data.frame") {
    p_grid_df <- rep(p, times = length(dlomax_prob_p) * length(dlomax_alpha))
    dlomax_alpha_grid_df <- rep(dlomax_alpha, each = length(p) * length(dlomax_prob_p))
    dlomax_prob_p_grid_df <- rep(dlomax_prob_p, times = length(p) * length(dlomax_alpha))
  }

  # CHECK p INPUT
  if(!lower.tail) {p <- 1 - p}
  discretelomax_check_p()
  # CHECK PARAMETERS
  discretelomax_check_param_args()

  # Expand combinations of arguments into a grid
  p_grid <- rep(p, times = length(dlomax_prob_p) * length(dlomax_alpha))
  dlomax_alpha_grid <- rep(dlomax_alpha, each = length(p) * length(dlomax_prob_p))
  dlomax_prob_p_grid <- rep(dlomax_prob_p, times = length(p) * length(dlomax_alpha))

  # Create vector of probabilities to return
  quant_out <- numeric(length(dlomax_alpha_grid))

  # quant_out[which(dlomax_alpha_grid == 0)] <- 1 + qgeom(
  #   p = p_grid[which(dlomax_alpha_grid == 0)],
  #   prob = dlomax_prob_p_grid[which(dlomax_alpha_grid == 0)])
  quant_out[which(dlomax_alpha_grid == 0)] <- ceiling(
    log(1-p_grid[which(dlomax_alpha_grid == 0)])/log(1 - dlomax_prob_p_grid[which(dlomax_alpha_grid == 0)]))
  quant_out[which(dlomax_alpha_grid > 0)] <- {
    ceiling(
      -(1/dlomax_alpha_grid[which(dlomax_alpha_grid > 0)]) *
        (1/log(1-dlomax_prob_p_grid[which(dlomax_alpha_grid > 0)])) *
        (1 - (1 - p_grid[which(dlomax_alpha_grid > 0)])^dlomax_alpha_grid[which(dlomax_alpha_grid > 0)]) / ((1 - p_grid[which(dlomax_alpha_grid > 0)])^(dlomax_alpha_grid[which(dlomax_alpha_grid > 0)]))
    )
  }
  quant_out[which(p_grid == 0)] <- 1

  if(return_type == "data.frame") {
    return(data.frame(
      q = quant_out,
      dlomax_prob_p = dlomax_prob_p_grid_df,
      dlomax_alpha = dlomax_alpha_grid_df,
      lte_prob = p_grid_df
    ))
  }
  else {
    return(quant_out)
  }
}

