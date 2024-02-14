test_that("stopping point of a geometric distribution with param p = 0.1 is correctly determined", {
  expect_equal(determine_stopping_point(pgeom, list(prob = 0.1), max_N = 200), 160)
})

test_that("stopping point of a geometric distribution with param p = 0.1 is correctly determined, even when a user-provided lower.tail = FALSE argument is given", {
  expect_equal(determine_stopping_point(pgeom, list(prob = 0.1, lower.tail = FALSE), max_N = 200), 160)
})

test_that("the search for the stopping point of a geometric distribution with param p = 0.1 throws an error when max_N is set too low", {
  expect_error(determine_stopping_point(pgeom, list(prob = 0.1), max_N = 20))
})
