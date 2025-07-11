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
#' # in the mix
#' create_events_from_vector(v, 100)
#'
#' df <- data.frame(
#' precip = v,
#' tmax = rnorm(n = 18, mean = 10, sd = 2),
#' threshold = 100,
#' date = seq.Date(from = as.Date("2025-01-01"), by = "day", length.out = 19)[-10])
#'
#' event_var = "precip"
#' thresh_col = "threshold"
#' time_col = "date"
#'
#' aux_vars = c("tmax", "tmin")
#' aux_funs = list(max = max, min = min, sum = sum, mean = mean)
#'
#' create_events_from_df(df = df, event_var = "precip",
#' aux_vars = aux_vars, aux_funs = aux_funs,
#' thresh_col = thresh_col, time_col = time_col)
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
#' restricted, for now, to days).
#' @param event_var (character) the variable within df for which events are
#' defined, e.g. "precip"
#' @param time_col (character) the name of the column with time information;
#' "date" by default, but it could also be something like "date_time_pacific"
#' if finer-res data are to be considered
#' @param thresh_col (character)the name of the column containing a
#' pre-calcuated threshold
#' @param aux_vars (character vector) other variables within df which you want
#' to summarize alongside events. For example, "tmin" or some measure of
#' stability might be of interest during precipitation events
#' @param aux_funs (list of functions) summary functions to to apply to,
#' aux_vars such as min or max; Internally, lists of functions are applied to
#' multiple aux_vars by means of the dplyr package's mutate/summarize-across
#' construct. See documentation for dplyr::across for details.
#' @export
create_events_from_df <- function(
    df, event_var = "precip", thresh_col = "threshold", time_col = "date",
    aux_vars = c("tmax", "tmin"),
    aux_funs = list(max = max, min = min, sum = sum, mean = mean)) {
  # Event processing --------------------------------------------------------
  unagg_events_df <- df %>%
    dplyr::mutate(run_id = data.table::rleid(df[,event_var] > df[,thresh_col])) %>%
    dplyr::group_by(`run_id`) %>%
    dplyr::mutate(active_event =
                    dplyr::if_else(
                      condition = max(.data[[event_var]] - .data[[thresh_col]]) > 0,
                      true = TRUE, false = FALSE, missing = FALSE)) %>%
    dplyr::ungroup()

  agg_events_df <- unagg_events_df %>%
    dplyr::group_by(`run_id`) %>%
    dplyr::mutate(
      start_time = min(.data[[time_col]]),
      active_event = dplyr::first(.data$active_event),
      duration = dplyr::n(),
      magnitude = sum(.data[[event_var]] - .data[[thresh_col]]),
      maximum = max(.data[[event_var]] - .data[[thresh_col]]),
      total = sum(.data[[event_var]]),
      abs_maximum = max(.data[[event_var]]),
      dplyr::across(.cols = dplyr::any_of(aux_vars), .fns = aux_funs),
      .before = 1
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      duration_uncertain = dplyr::case_when(
        dplyr::row_number() == 1 | dplyr::row_number() == dplyr::n()
        ~ TRUE,
        .data$active_event &
          .data$run_id != dplyr::lag(.data$run_id) &
          is.na(dplyr::lag(.data[[event_var]]))
        ~ TRUE,
        .data$active_event &
          .data$run_id != dplyr::lead(.data$run_id) &
          is.na(dplyr::lead(.data[[event_var]]))
        ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::mutate(duration_between_events = dplyr::lag(.data$duration), .before = 1) %>%
    dplyr::filter(.data$active_event) %>%
    dplyr::group_by(.data$run_id) %>%
    dplyr::mutate(duration_uncertain = any(.data$duration_uncertain), .keep = "all") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(event_number = dplyr::row_number(), .before = 1) %>%
    dplyr::select(-any_of(c("run_id", "active_event", {{event_var}}, {{aux_vars}}, {{time_col}})))

  return(list(unagg_events_df = unagg_events_df, agg_events_df = agg_events_df))
}

# create_events_from_table <- function(df, var,
#                                      time_colname = "date",
#                                      threshold_colname = NULL,
#                                      threshold_desc = "75p_above_0.3",
#                                      aux_vars = NULL,
#                                      aux_funs = list(min,max,mean)) {
#   # ADD EXPLICIT NAs WHENEVER OBS ARE MISSING
#   df <- df_complete_time_column(df)
#
#   threshold = threshold_given_desc(v = df[[var]], threshold_desc = threshold_desc)
#
#   out_list <- create_events_from_vector(v = df[[var]], threshold = threshold)
#
#   out_list$event_summaries_df$start_time <- df[[time_colname]][out_list$event_summaries_df$start_index]
#
#   out_list
# }
