# Tests for sim_tte_df bugfix

test_that("default behavior works", {
  set.seed(1)
  dat <- data.frame(
    ID = rep(1:3, each = 20),
    time = rep(seq(0.5, 10, by = 0.5), 3),
    p11 = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 3)
  )
  result <- sim_tte_df(dat)
  expect_s3_class(result, "data.frame")

  expect_true(all(c("sim_time", "sim_status", "ID") %in% names(result)))
  expect_equal(nrow(result), 3)
})

test_that("reordered columns do not change results", {
  dat <- data.frame(
    ID = rep(1:3, each = 20),
    time = rep(seq(0.5, 10, by = 0.5), 3),
    p11 = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 3)
  )
  dat_reordered <- dat[, c("time", "p11", "ID")]

  set.seed(42)
  r1 <- sim_tte_df(dat)
  set.seed(42)

  r2 <- sim_tte_df(dat_reordered)
  expect_identical(r1, r2)
})

test_that("custom variable names work", {
  set.seed(1)
  dat <- data.frame(
    subj = rep(c(101, 205), each = 20),
    tgrid = rep(seq(0.5, 10, by = 0.5), 2),
    surv = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 2)
  )
  result <- sim_tte_df(dat, id_var = "subj", time_var = "tgrid",
    surv_var = "surv")
  expect_true(all(result$ID %in% c(101, 205)))
})

test_that("original IDs are preserved", {
  set.seed(1)
  dat <- data.frame(
    ID = rep(c(101, 205), each = 20),
    time = rep(seq(0.5, 10, by = 0.5), 2),
    p11 = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 2)
  )
  result <- sim_tte_df(dat)
  expect_equal(sort(result$ID), c(101, 205))
})

test_that("missing required columns fail clearly", {
  dat <- data.frame(ID = 1, time = 1, p11 = 0.5)
  expect_error(sim_tte_df(dat, time_var = "t"), "Column 't' not found")
  expect_error(sim_tte_df(dat, id_var = "subj"), "Column 'subj' not found")
  expect_error(sim_tte_df(dat, surv_var = "S"), "Column 'S' not found")
})

test_that("bad input types fail clearly", {
  dat <- data.frame(ID = 1:2, time = c("a", "b"), p11 = c(0.9, 0.8),
    stringsAsFactors = FALSE)
  expect_error(sim_tte_df(dat), "must be numeric")

  dat2 <- data.frame(ID = 1:2, time = c(1, 2), p11 = c("a", "b"),
    stringsAsFactors = FALSE)
  expect_error(sim_tte_df(dat2), "must be numeric")
})

test_that("duplicated column names fail", {
  dat <- data.frame(ID = 1, time = 1, p11 = 0.5, p11 = 0.5,
    check.names = FALSE)
  expect_error(sim_tte_df(dat), "duplicated column names")
})

test_that("xdata join works by ID", {
  set.seed(1)
  dat <- data.frame(
    ID = rep(c(10, 20), each = 20),
    time = rep(seq(0.5, 10, by = 0.5), 2),
    p11 = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 2)
  )
  xdata <- data.frame(ID = c(10, 20), trt = c("A", "B"))
  result <- sim_tte_df(dat, xdata = xdata)
  expect_true("trt" %in% names(result))
  expect_equal(nrow(result), 2)

  # xdata missing id_var
  xdata_bad <- data.frame(subj = 1, trt = "A")
  expect_error(sim_tte_df(dat, xdata = xdata_bad), "not found in 'xdata'")

  # xdata duplicated IDs
  xdata_dup <- data.frame(ID = c(10, 10), trt = c("A", "B"))
  expect_error(sim_tte_df(dat, xdata = xdata_dup), "duplicated IDs")
})

test_that(".get_tte uses correct boundary rule", {
  # S(t) = 0.9, 0.7, 0.4; U = 0.7

  # First index where S <= U is index 2
  pcurr <- c(0.9, 0.7, 0.4)
  expect_equal(simtte:::.get_tte(0.7, pcurr), 2L)

  # No event case
  pcurr2 <- c(0.9, 0.8, 0.7)
  expect_equal(simtte:::.get_tte(0.5, pcurr2), -99L)
})
