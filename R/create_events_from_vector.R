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
  temp_rle <- rle(ifelse(v > threshold,1,0))
  end_indices <- cumsum(temp_rle$lengths)[which(temp_rle$values != 0)] # end
  #' indices of runs
  lengths <- temp_rle$lengths[which(temp_rle$values == 1)] # durations
  start_indices <- (end_indices - lengths) + 1 # start indices of runs

  return(cbind(data.frame(start_index = start_indices,
                          length = lengths),
               do.call(rbind,
                       lapply(X = seq_along(start_indices),
                              FUN = function(i){
                                vals = (v-threshold)[start_indices[i]:end_indices[i]]
                                data.frame(magnitude = sum(vals),
                                           maximum = max(vals))
                              }))
  ))
}

# scratch
v <- c(333,587,0,201,82,225,70,62,47,101,187)
threshold <- 100
create_events_from_vector(v, threshold)



