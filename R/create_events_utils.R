#' Data prep utils used by create_events\* functions
#' @name create_events_utils
NULL

#' @name create_events_utils
#'
#' @param df data.frame() with a column called time_colname with Date or POSIXt
#' entries
#' @param time_colname (character) the name of the date or time column
#' @param by (character) the increment of the sequence, e.g. "day" or "hour";
#' See base::seq.Date or base::seq.POSIXt for details
#'
#' @examples
#' df <- data.frame(date = as.Date(c("2025-01-01", "2025-01-03")), precip = c(0,1.2))
#' df_complete_time_column(df)
#'
#' @export
df_complete_time_column <- function(df, time_colname = "date", by = "day") {

  time_df <- data.frame(time_col = seq.Date(
    from = min(df[[time_colname]]),
    to = max(df[[time_colname]]),
    by = by))

  names(time_df)[1] <- time_colname

  merge(x = time_df, y = df, by = time_colname, all.x = TRUE)
}

#' @rdname create_events_utils
#' @export
threshold_given_desc <- function(v, threshold_desc) {
  # Calculates a threshold given a vector and a description of it
  threshold_fun <- parse_threshold_desc(threshold_desc)
  threshold_fun(v)
}

#' @rdname create_events_utils
#' @export
parse_threshold_desc <- function(threshold_desc) {
  # given a description such as "75p_above_1", this returns the function used
  # to calculate the threshold
  #
  # For now, a default placeholder
  function(x) {
    x <- stats::na.omit(x)
    stats::quantile(x[which(x > 1)], probs = 0.75)
    }
}
