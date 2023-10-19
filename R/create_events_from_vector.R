#' Create events from a regularly-sampled time series which is represented as a
#' vector
#'
#' @param v A vector (typically numeric) representing some variable of interest
#' which has been measured at regular time intervals.
#' @param threshold A numeric threshold which defines the minimum value for
#' which an observation would be considered part of an event
#'
#' @return A list containing:
#' \itemize{
#' \item event_summaries_df, a data.frame() -
#' Each row provides a basic summary of events, including their starting
#' indices, durations, magnitudes, maxima, and whether or not the duration,
#' magnitude, and maximum are to be taken at face value}
#' @export
#'
#' @examples
#' v <- c(333,587,0,201,82,225,70,62,47,101,187,NA,NA,126,222,86,NA,125)
#' # ^^ perhaps precipitation measured in 10ths of mm, with some missing values
#' #' in the mix
#' create_events_from_vector(v, 100)
#'
#' \dontrun{
#' $event_summaries_df
#'   start_index length magnitude maximum length_uncertain
#' 1           1      2       487     487             TRUE
#' 2           4      1       101     101            FALSE
#' 3           6      1       125     125            FALSE
#' 4          10      2        88      87             TRUE
#' 5          14      2       148     122             TRUE
#' 6          18      1        25      25             TRUE
#' }
create_events_from_vector <- function(v, threshold) {
  temp_rle <- rle(ifelse(v > threshold,1,0))
  end_indices <- cumsum(temp_rle$lengths)[which(temp_rle$values != 0)] # end
  # indices of runs
  lengths <- temp_rle$lengths[which(temp_rle$values == 1)] # durations
  start_indices <- (end_indices - lengths) + 1 # start indices of runs
  event_summaries_df <- cbind(data.frame(start_index = start_indices,
                                         length = lengths),
                              do.call(rbind,
                                      lapply(X = seq_along(start_indices),
                                             FUN = function(i){
                                               # vals = (v-threshold)[start_indices[i]:end_indices[i]]
                                               # data.frame(magnitude = sum(vals),
                                               #            maximum = max(vals))

                                               vals = (v-threshold)[(start_indices[i]-1):(end_indices[i]+1)]
                                               #length_uncertain = ifelse(length(vals[1]) == 0 | is.na(vals[1]) | is.na(vals[length(vals)]), TRUE, FALSE)
                                               data.frame(magnitude = sum(vals[c(-1,-length(vals))]),
                                                          maximum = max(vals[c(-1,-length(vals))]),
                                                          length_uncertain = ifelse(start_indices[i] == 1 | is.na(vals[1]) | is.na(vals[length(vals)]), TRUE, FALSE))
                                             }))
  )
  return(list(event_summaries_df = event_summaries_df))
}





