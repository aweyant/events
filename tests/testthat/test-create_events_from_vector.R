test_that("create_events_from_vector() outputs expected event summary data.frame() given an imperfect numeric vector input", {
  v <- c(333,587,0,201,82,225,70,62,47,101,187,NA,NA,126,222,86,NA,125)
  threshold <- 100
    expect_equivalent(create_events_from_vector(v, threshold)$event_summaries_df,
               data.frame(start_index = c(1,4,6,10,14,18),
                          length = c(2,1,1,2,2,1),
                          magnitude = c(487,101,125,88,148,25),
                          maximum = c(487,101,125,87,122,25),
                          length_uncertain = c(TRUE,FALSE,FALSE,TRUE,TRUE,TRUE)))
})
