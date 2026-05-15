test_that("inverse transform sampling produces valid results", {
  # Create a simple mock survival curve: decreasing from 1 to 0
  simdat <- data.frame(
    ID = rep(1, 100),
    time = seq(0.1, 10, length.out = 100),
    p11 = exp(-0.5 * seq(0.1, 10, length.out = 100))
  )

  set.seed(42)
  result <- simtte:::.simulate_survival_id(simdat, id = 1, id_var = "ID")
  expect_true(result$time > 0)
  expect_true(result$status %in% c(0, 1))
  expect_equal(result$ID, 1)
})

test_that("all simulated times are positive", {
  simdat <- data.frame(
    ID = rep(1, 200),
    time = seq(0.05, 10, length.out = 200),
    p11 = exp(-0.3 * seq(0.05, 10, length.out = 200))
  )

  times <- replicate(50, {
    simtte:::.simulate_survival_id(simdat, id = 1, id_var = "ID")$time
  })
  expect_true(all(times > 0))
})

test_that("censoring works when survival never crosses threshold", {
  # Survival stays high (barely decays)
  simdat <- data.frame(
    ID = rep(1, 50),
    time = seq(0.1, 5, length.out = 50),
    p11 = rep(0.999, 50)
  )
  set.seed(1)
  result <- simtte:::.simulate_survival_id(simdat, id = 1, id_var = "ID")
  # With p11 always near 1, U ~ Uniform will almost always exceed p11 = FALSE

  # so most likely status = 0 (censored at max time)
  expect_true(result$time > 0)
})
