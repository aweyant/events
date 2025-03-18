#' Create events from observed time series and a threshold
#'
#' @name create_events
NULL

#' @rdname create_events
#' @param v A vector (assumed to be numeric) representing some variable of
#' interest which has been measured at regular time intervals; Missing values
#' or "gaps in observed series" must be explicitly filled with NAs.
#' @param threshold A numeric threshold which defines the minimum value for
#' which an observation would be considered part of an event
#'
#' @return A list containing:
#' \itemize{
#' \item threshold, the vector above which the series was truncated in order to
#' define events
#' \item event_summaries_df, a data.frame() -
#' Each row provides a basic summary of events, including their starting
#' indices, durations, magnitudes, maxima, and whether or not the duration,
#' magnitude, and maximum are to be taken at face value}
#' @examples
#' v <- c(333,587,0,201,82,225,70,62,47,101,187,NA,NA,126,222,86,NA,125)
#' # ^^ perhaps precipitation measured in 10ths of mm, with some missing values
#' #' in the mix
#' create_events_from_vector(v, 100)
#'
#' df <- data.frame(
#' precip = v,
#' tmax = rnorm(n = 18, mean = 10, sd = 2),
#' date = seq.Date(from = as.Date("2025-01-01"), by = "day", length.out = 19)[-10])
#'
#' create_events_from_table(df = df, var = "precip", aux_vars = "tmax")
#'
#'
#' \dontrun{
#' $threshold
#' [1] 100
#'
#' $event_summaries_df
#'   start_index length magnitude maximum length_uncertain
#' 1           1      2       487     487             TRUE
#' 2           4      1       101     101            FALSE
#' 3           6      1       125     125            FALSE
#' 4          10      2        88      87             TRUE
#' 5          14      2       148     122             TRUE
#' 6          18      1        25      25             TRUE
#' }
#' @export
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
  return(list(threshold = threshold,
              event_summaries_df = event_summaries_df))
}

#' @rdname create_events
#' @param df A data.frame of observations recorded at regular intervals (
#' restricted, for now, to days). A column called 'date' must exist as well as
#' var, specified by the user-provided argument var.
#' @param var (character) the variable within df for which events are defined,
#' e.g. "precip"
#' @param time_colname (character) NOT YET IMPLEMENTED the name of the column
#' with time information; "date" by default
#' @param threshold_colname (character) NOT YET IMPLEMENTED the name of the
#' column containing a pre-calcuated threshold; This option is provided for
#' users who wish to experiment with thresholds beyond what threshold_desc
#' allows for.
#' @param threshold_desc (character) NOT YET IMPLEMENTED a string of a
#' particular format (see )
#' @param aux_vars (character vector) NOT YET IMPLEMENTED other variables within
#' df which you want to summarize alongside events. For example, "tmin" might be
#' of interest during precipitation events
#' @param aux_funs (list of functions) NOT YET IMPLEMENTED summary functions to
#' to apply to aux_vars
#' @export
create_events_from_table <- function(df, var,
                                     time_colname = "date",
                                     threshold_colname = NULL,
                                     threshold_desc = "75p_above_0.3",
                                     aux_vars = NULL,
                                     aux_funs = list(min,max,mean)) {
  # ADD EXPLICIT NAs WHENEVER OBS ARE MISSING
  df <- df_complete_time_column(df)

  threshold = threshold_given_desc(v = df[[var]], threshold_desc = threshold_desc)

  out_list <- create_events_from_vector(v = df[[var]], threshold = threshold)

  out_list$event_summaries_df$start_time <- df[[time_colname]][out_list$event_summaries_df$start_index]

  out_list
}





