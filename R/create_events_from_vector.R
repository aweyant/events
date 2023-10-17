#' Create events from a regularly-sampled time series which is represented as a
#' vector
#'
#' @param v A vector (typically numeric) representing some variable of interest
#' which has been measured at regular time intervals.
#' @param threshold A numeric threshold which defines the minimum value for
#' which an observation would be considered part of an event
#'
#' @return unknown for now...
#' @export
#'
#' @examples
#' v <- c(333,587,0,201,82,225) # perhaps precipitation measured in 10ths of mm
#' create_events_from_vector(v, 100)
#' # ... some output ...
create_events_from_vector <- function(v, threshold) {
  # . . .
  return(v)
}
