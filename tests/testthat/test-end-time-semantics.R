# Pre-Phase-B audit, Phase 6: end_time boundary semantics.
#
# Weibull: end_time >= 0 (closed-form hazard defined for all t >= 0);
#   end_time = 0 is a valid degenerate zero-duration follow-up.
# M-spline: min(time) <= end_time <= max(time) (the baseline hazard is
#   only defined on the supplied grid; no extrapolation in either
#   direction).

test_that("M-spline: end_time < min(time) is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(2, 3), basehaz = matrix(c(1, 2),
                ncol = 1), end_time = 1),
        "earlier than min\\(time\\)"
    )
})

test_that("M-spline: end_time == min(time) is accepted (one-point grid)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(2, 3), basehaz = matrix(c(1, 2), ncol = 1),
        end_time = 2)
    expect_equal(out$time, 2)
    expect_equal(out$p11, 1) # no interval has elapsed yet
})

test_that("M-spline: end_time == max(time) is accepted (default case)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 2)
    expect_equal(max(out$time), 2)
})

test_that("M-spline: end_time > max(time) is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
                ncol = 1), end_time = 5),
        "exceeds max\\(time\\)"
    )
})

test_that("Weibull: end_time = 0 is accepted and censors every subject at 0", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(0, 0.5, -0.5), nrow = 3)
    result <- sim_tte(pi = lp, mu = -1, coefs = 1.2,
        time = seq(0, 5, by = 0.5), end_time = 0)
    expect_true(all(result$sim_time == 0))
    expect_true(all(result$sim_status == 0))
})

test_that("Weibull: end_time = 0 works even when 'time' does not include 0", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    result <- sim_tte(pi = 0, mu = -1, coefs = 1.2,
        time = seq(0.1, 5, by = 0.1), end_time = 0)
    expect_equal(result$sim_time, 0)
    expect_equal(result$sim_status, 0)
})

test_that("negative end_time is still rejected for both model types", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte(pi = 0, mu = -1, coefs = 1, time = c(0, 1),
        end_time = -1), "non-negative")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(0, 1), basehaz = matrix(c(1, 2),
                ncol = 1), end_time = -1),
        "non-negative"
    )
})

test_that("M-spline: end_time = 0 is accepted when 0 is within [min(time), max(time)]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 0)
    expect_equal(out$time, 0)
    expect_equal(out$p11, 1)
})

test_that("M-spline: end_time = 0 is rejected when the grid does not start at 0", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(0.2, 1, 2), basehaz = matrix(c(1, 2, 4),
                ncol = 1), end_time = 0),
        "earlier than min\\(time\\)"
    )
})
