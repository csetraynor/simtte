test_that("sim_tte validates inputs", {
  expect_error(sim_tte(pi = "abc"), "'pi' must be numeric")
  expect_error(sim_tte(pi = 1, mu = c(1, 2)), "'mu' must be a numeric scalar")
  expect_error(sim_tte(pi = 1, coefs = "a"), "'coefs' must be numeric")
  expect_error(sim_tte(pi = 1, type = "badtype"), "'arg' should be one of")
  expect_error(sim_tte(pi = 1, coefs = c(1, 2), type = "weibull"),
    "Only 1 coefficient")
  expect_error(sim_tte(pi = 1, coefs = -1, type = "weibull"),
    "shape parameter must be positive")
  expect_error(sim_tte(pi = 1, type = "ms", basis = NULL),
    "'basis' must be provided")
})

test_that("sim_tte runs Weibull simulation", {
  skip_on_cran()
  skip_if_not_installed("mrgsolve")
  set.seed(42)
  lp <- matrix(runif(10), nrow = 10)
  result <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
    time = seq(0.1, 50, by = 0.1), type = "weibull", end_time = 50)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("sim_time", "sim_status", "ID", "lp") %in% names(result)))
  expect_equal(nrow(result), 10)
  expect_true(all(result$sim_time > 0))
  expect_true(all(result$sim_status %in% c(0, 1)))
})

test_that("sim_tte is reproducible with set.seed", {

  skip_on_cran()
  skip_if_not_installed("mrgsolve")
  lp <- matrix(c(0.1, 0.2, 0.3), nrow = 3)
  set.seed(123)
  r1 <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
    time = seq(0.1, 20, by = 0.1), type = "weibull", end_time = 20)
  set.seed(123)
  r2 <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
    time = seq(0.1, 20, by = 0.1), type = "weibull", end_time = 20)
  expect_identical(r1, r2)
})

test_that("sim_tte M-splines simulation runs", {
  skip_on_cran()
  skip_if_not_installed("mrgsolve")
  data("ms_data", package = "simtte")
  set.seed(42)
  lp <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))
  result <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
    coefs = ms_data$coefs, time = ms_data$time, type = "ms")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("sim_time", "sim_status", "ID") %in% names(result)))
})
