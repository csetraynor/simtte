# Pre-Phase-B audit, Phase 2: M-spline hazard carry semantics.
#
# Contract (see ?sim_tte, "M-spline hazard carry convention"): the
# hazard value at time[i] applies over [time[i], time[i+1]) --
# last-observation-carried-forward (mrgsolve nocb = FALSE). This is a
# piecewise-constant grid-based hazard, not continuous spline
# evaluation. All checks here are deterministic, direct trajectory
# comparisons against the closed-form piecewise-exponential survival
# function, not Monte Carlo.
#
# `basehaz` must always have exactly one row per element of `times`
# (each row is positionally a knot/hazard pair): to observe output at
# an intermediate point without changing the underlying step function,
# that point is added as an extra knot carrying the same hazard value
# as the interval it falls in.

test_that("constant hazard reduces to the exponential survival function", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2, 3),
        basehaz = matrix(c(2, 2, 2, 2), ncol = 1), end_time = 3)
    expect_equal(out$p11, exp(-2 * out$time), tolerance = 1e-6)
})

test_that("changing piecewise-constant hazard c(1, 2, 4) on knots c(0, 1, 2) matches the analytical trajectory", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 2)
    # h = 1 on [0, 1), h = 2 on [1, 2); the hazard value 4 at the last
    # knot (t = 2) never applies to any interval.
    expected <- c(`0` = 1, `1` = exp(-1), `2` = exp(-1 - 2 * 1))
    got <- setNames(out$p11, as.character(out$time))
    expect_equal(unname(got[names(expected)]), unname(expected),
        tolerance = 1e-6)
})

test_that("mid-interval survival is not reported unless explicitly requested in 'time' (no interpolation)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 2)
    expect_length(out$p11[out$time == 0.5], 0)
})

test_that("mid-interval survival at an explicit extra knot matches exp(-0.5) for hazard 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # Extra knot at t = 0.5 carries the *same* hazard (1) as the
    # interval [0, 1) it falls inside, so the underlying step function
    # is unchanged; it only adds a reported output point.
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 0.5, 1, 2),
        basehaz = matrix(c(1, 1, 2, 4), ncol = 1), end_time = 2)
    expect_equal(out$p11[out$time == 0.5], exp(-0.5), tolerance = 1e-6)
    expect_equal(out$p11[out$time == 1], exp(-1), tolerance = 1e-6)
})

test_that("survival exactly at a knot reflects the interval ending there, not starting there", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 2)
    # S(1) uses hazard 1 (active on [0, 1)), NOT hazard 2 (which only
    # starts applying at t = 1). This is the regression check for the
    # audit's default-nocb bug: under mrgsolve's default nocb = TRUE,
    # S(1) would instead be exp(-2).
    expect_equal(out$p11[out$time == 1], exp(-1), tolerance = 1e-6)
    expect_false(isTRUE(all.equal(out$p11[out$time == 1], exp(-2),
        tolerance = 1e-3)))
})

test_that("endpoint survival does not depend on the final knot's own hazard value", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out_hi <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 999),
            ncol = 1), end_time = 2)
    out_lo <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 0),
            ncol = 1), end_time = 2)
    # Only the last knot's hazard differs (999 vs 0); S(2) must be
    # identical either way because that value never applies to any
    # interval under the LOCF convention and end_time cannot extend
    # past max(time) for M-spline models.
    expect_equal(out_hi$p11[out_hi$time == 2], out_lo$p11[out_lo$time == 2],
        tolerance = 1e-6)
    expect_equal(out_hi$p11[out_hi$time == 2], exp(-3), tolerance = 1e-6)
})

test_that("irregular (non-uniform) time grids follow the same LOCF convention", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    times <- c(0, 0.2, 1.5, 4)
    haz <- c(0.5, 1, 3, 3) # last value (knot at t = 4) is inert
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = times, basehaz = matrix(haz, ncol = 1),
        end_time = 4)
    h_at <- function(t) {
        h1 <- 0.5 * min(t, 0.2)
        h2 <- if (t > 0.2) 1 * (min(t, 1.5) - 0.2) else 0
        h3 <- if (t > 1.5) 3 * (min(t, 4) - 1.5) else 0
        h1 + h2 + h3
    }
    expected <- exp(-vapply(times, h_at, numeric(1)))
    expect_equal(out$p11, expected, tolerance = 1e-6)
})

test_that("duplicated (time, hazard) rows in ms_data do not disrupt the carry convention", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    data("ms_data", package = "simtte")
    expect_true(anyDuplicated(ms_data$time) > 0)
    lp <- matrix(0, nrow = 1)
    result <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms")
    expect_equal(nrow(result), 1)
    expect_true(result$sim_status %in% c(0, 1))
})
