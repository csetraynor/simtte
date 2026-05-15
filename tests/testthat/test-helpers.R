test_that("helper .get_tte returns correct index", {
  # Survival probabilities decreasing
  p <- c(0.9, 0.7, 0.5, 0.3, 0.1)
  # U = 0.6 means first time p < U is at index 3 (p=0.5)

  idx <- simtte:::.get_tte(0.6, p)
  expect_true(idx > 0)

  # U very small, event at first step

  idx2 <- simtte:::.get_tte(0.95, p)
  expect_equal(idx2, 1L)
})

test_that("helper .get_tte returns -99 when no event", {
  p <- c(0.99, 0.98, 0.97, 0.96)
  idx <- simtte:::.get_tte(0.01, p)
  expect_equal(idx, -99L)
})

test_that(".get_time and .get_max_time work correctly", {
  df <- data.frame(ID = 1, time = c(1, 2, 3, 4, 5))
  expect_equal(simtte:::.get_time(df, 3), 3)
  expect_equal(simtte:::.get_max_time(df), 5)
})
