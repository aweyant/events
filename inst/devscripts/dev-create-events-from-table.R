# Inputs ------------------------------------------------------------------
v <- c(333,587,0,201,82,225,70,62,47,101,187,NA,NA,126,222,86,NA,125)
# ^^ perhaps precipitation measured in 10ths of mm, with some missing values
#' in the mix
create_events_from_vector(v, 100)

df <- data.frame(
precip = v,
threshold = 100,
tmax = rnorm(n = 18, mean = 10, sd = 2),
date = seq.Date(from = as.Date("2025-01-01"), by = "day", length.out = 19)[-10])


# Other function arguments ------------------------------------------------
event_var = "precip"
thresh_col = "threshold"
time_col = "date"

aux_vars = c("tmax", "tmin")
aux_funs = list(max = max, min = min, sum = sum, mean = mean)

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
      dplyr::across(.cols = any_of(aux_vars), .fns = aux_funs),
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
}

#dplyr::select(start_time, time_since_last_event, precip, active_event, run_id, duration_uncertain, any_of(event_var))

temp_rle <- rle(ifelse(df[,event_var] > df[,thresh_col],1,0))

end_indices <- cumsum(temp_rle$lengths)[which(temp_rle$values != 0)] # end
# indices of runs
lengths <- temp_rle$lengths[which(temp_rle$values == 1)] # durations
start_indices <- (end_indices - lengths) + 1 # start indices of runs
event_summaries_df <- cbind(
  data.frame(start_index = start_indices,
             length = lengths),
  do.call(rbind,
          lapply(X = seq_along(start_indices),
                 FUN = function(i){
                   # vals = (v-threshold)[start_indices[i]:end_indices[i]]
                   # data.frame(magnitude = sum(vals),
                   #            maximum = max(vals))
                   vals = (df[,event_var] - df[,thresh_col])[(start_indices[i]-1):(end_indices[i]+1)]
                   #length_uncertain = ifelse(length(vals[1]) == 0 | is.na(vals[1]) | is.na(vals[length(vals)]), TRUE, FALSE)
                   data.frame(magnitude = sum(vals[c(-1,-length(vals))]),
                              maximum = max(vals[c(-1,-length(vals))]),
                              length_uncertain = ifelse(
                                start_indices[i] == 1 | is.na(vals[1]) |is.na(vals[length(vals)]), TRUE, FALSE))
                   }))
)
