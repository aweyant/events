#' Title
#'
#' @param x
#' @param alpha
#' @param prob_p
#'
#' @return
#' @export
#'
#' @details
#' In the output vector, it is the prop_p which vary more quickly. If you
#' provide 2 alphas and 3 prob_ps, the first 3 numbers in the output are
#' associated with the first alpha.
#'
#'
#' @examples
#' \dontrun{
#' ddiscretelomax(x = 1:10, alpha = c(0,0.2,0.4),
#' prob = c(0.6), return_type = "data.frame") %>%
#' ggplot2::ggplot(
#' ggplot2::aes(x = x, y = prob_mass, group = alpha, fill = alpha)) +
#' ggplot2::geom_col(position = "dodge") +
#' ggplot2::scale_fill_binned() +
#' ggplot2::scale_y_continuous(trans = "log")
#' }
ddiscretelomax <- function(x, alpha, prob_p, return_type = "vector") {
  # PRESERVE ORIGINAL Xs IN CASE BAD ARGUMENTS ARE GIVEN BUT WE STILL WANT TO
  # RETURN A DATAFRAME WITH DETAILS
  if(return_type == "data.frame") {
    x_grid_df <- rep(x, times = length(prob_p) * length(alpha))
    alpha_grid_df <- rep(alpha, each = length(x) * length(prob_p))
    prob_p_grid_df <- rep(prob_p, times = length(x) * length(alpha))
  }

  # CHECK PARAMETERS
  discretelomax_check_args()

  # Expand combinations of arguments into a grid
  x_grid <- rep(x, times = length(prob_p) * length(alpha))
  alpha_grid <- rep(alpha, each = length(x) * length(prob_p))
  prob_p_grid <- rep(prob_p, times = length(x) * length(alpha))

  # Create vector of probabilities to return
  prob_out <- numeric(length(alpha_grid))

  prob_out[which(alpha_grid == 0)] <- dgeom(
    x = x - 1,
    prob = prob_p_grid[which(alpha_grid == 0)])
  prob_out[which(alpha_grid > 0)] <- {
    (1 - alpha_grid[which(alpha_grid > 0)]*(x-1)*log(1-prob_p_grid[which(alpha_grid > 0)])
    )^(-1/alpha_grid[which(alpha_grid > 0)]) - (
      1 - alpha_grid[which(alpha_grid > 0)]*x*log(1-prob_p_grid[which(alpha_grid > 0)])
    )^(-1/alpha_grid[which(alpha_grid > 0)])
  }
  if(return_type == "data.frame") {
    return(data.frame(
      x = x_grid_df,
      prob_p = prob_p_grid_df,
      alpha = alpha_grid_df,
      prob_mass = prob_out
    ))
  }
  else {
    return(prob_out)
  }
}


discretelomax_check_args <- function() {
  # Check x
  call_x_gt_0 <- rlang::expr(x[which(x <= 0)] <- NaN)
  rlang::eval_bare(call_x_gt_0, env = parent.frame())

  # Check alpha
  call_alpha_gt_0 <- rlang::expr(alpha[which(alpha < 0)] <- NaN)
  rlang::eval_bare(call_alpha_gt_0, env = parent.frame())

  # Check prob p
  call_prob_p_gt_0 <- rlang::expr(prob_p[which(prob_p <= 0)] <- NaN)
  call_prob_p_lt_1 <- rlang::expr(prob_p[which(prob_p > 1)] <- NaN)

  rlang::eval_bare(call_prob_p_gt_0, env = parent.frame())
  rlang::eval_bare(call_prob_p_lt_1, env = parent.frame())

  # Throw warnings, if necessary
  call_alpha_warning <- rlang::expr({
    if(any(is.nan(alpha))) {
      warning("Alpha must be >= 0. NaN(s) returned.")
    }
  })
  call_prob_p_warning <- rlang::expr({
    if(any(is.nan(prob_p))) {
      warning("prob_p must be >0 and <1. NaN(s) returned.")
    }
  })
  call_x_warning <- rlang::expr({
    if(any(is.nan(x))) {
      warning("x is supported only on positive integers. NaN(s) returned.")
    }
  })
  rlang::eval_bare(call_alpha_warning, env = parent.frame())
  rlang::eval_bare(call_prob_p_warning, env = parent.frame())
  rlang::eval_bare(call_x_warning, env = parent.frame())
  NULL
}
