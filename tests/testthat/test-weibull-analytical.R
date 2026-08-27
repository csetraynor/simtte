# Phase A1: Weibull correctness and output-grid semantics.
#
# These tests validate the survival trajectory produced by the built-in
# Weibull mrgsolve model (inst/models/weibull.cpp) directly against the
# closed-form analytical survival function
#
#   S(t) = exp[-exp(mu + lp) * t^shape]
#
# using .sim_surv_df(), which is the internal function that drives both
# sim_tte() and explore_pi_tq_surv(). A tolerance of 1e-8 is used
# throughout: the corrected model computes p11 as a closed-form
# expression (see inst/models/weibull.cpp), so the only expected
# discrepancy from the analytical formula is double-precision rounding,
# not solver/discretization error. This tolerance is not tied to a
# particular machine or mrgsolve version.

analytical_survival <- function(time, mu, lp, shape) {
    eta <- exp(mu + lp)
    exp(-eta * time^shape)
}

test_that("Weibull p11 matches the analytical survival function (shape < 1)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    grid <- c(seq(0.001, 0.099, by = 0.007), seq(0.1, 5, by = 0.3))
    out <- simtte:::.sim_surv_df(log_hr = c(0.2, -0.4), mu = -1,
        shape = 0.5, type = "weibull", times = grid)
    expect_equal(out$p11, analytical_survival(out$time, -1, out$lp, 0.5),
        tolerance = 1e-8)
})

test_that("Weibull p11 matches the analytical survival function (shape = 1)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    grid <- c(seq(0.001, 0.099, by = 0.007), seq(0.1, 5, by = 0.3))
    out <- simtte:::.sim_surv_df(log_hr = c(0.2, -0.4), mu = -1,
        shape = 1, type = "weibull", times = grid)
    expect_equal(out$p11, analytical_survival(out$time, -1, out$lp, 1),
        tolerance = 1e-8)
    # shape = 1 is the exponential model: hazard is constant.
    expect_equal(out$p11, exp(-exp(-1 + out$lp) * out$time), tolerance = 1e-8)
})

test_that("Weibull p11 matches the analytical survival function (shape > 1)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    grid <- c(seq(0.001, 0.099, by = 0.007), seq(0.1, 5, by = 0.3))
    out <- simtte:::.sim_surv_df(log_hr = c(0.2, -0.4), mu = -1,
        shape = 2.5, type = "weibull", times = grid)
    expect_equal(out$p11, analytical_survival(out$time, -1, out$lp, 2.5),
        tolerance = 1e-8)
})

test_that("Weibull p11 is exact near zero, including t < 0.1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # This is precisely the region the old near-zero workaround got
    # wrong (SOLVERTIME < 0.1 was approximated by a shape-independent
    # exponential model). shape = 0.3 has a divergent instantaneous
    # hazard as t -> 0, the most stringent case.
    near_zero <- c(0.001, 0.01, 0.02, 0.05, 0.08, 0.099)
    for (shape in c(0.3, 0.7, 1.3, 3)) {
        out <- simtte:::.sim_surv_df(log_hr = 0.1, mu = -0.5,
            shape = shape, type = "weibull", times = near_zero)
        expect_equal(out$p11,
            analytical_survival(out$time, -0.5, 0.1, shape),
            tolerance = 1e-8,
            info = paste("shape =", shape))
    }
})

test_that("regression: old near-zero exponential approximation would fail this", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # Under the old implementation, for SOLVERTIME < 0.1 the ODE used
    # dxdt_p11 = -p11 * eta regardless of shape, i.e. treated the
    # hazard as if shape == 1. For shape far from 1, this produces a
    # survival value very different from the true Weibull formula at
    # t = 0.05. The corrected model must match the true formula to
    # near machine precision; the old (shape = 1) approximation would
    # miss by a wide, easily-detectable margin.
    shape <- 3
    mu <- -1
    lp <- 0.2
    t <- 0.05
    out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
        type = "weibull", times = t)
    correct <- analytical_survival(t, mu, lp, shape)
    old_wrong <- exp(-exp(mu + lp) * t) # shape treated as 1

    expect_equal(out$p11, correct, tolerance = 1e-8)
    # The old (incorrect) value differs from the correct one by more
    # than can be explained by numerical noise.
    expect_gt(abs(old_wrong - correct), 1e-3)
    expect_gt(abs(out$p11 - old_wrong), 1e-3)
})

test_that("Weibull handles multiple linear predictors simultaneously", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- c(-1, -0.5, 0, 0.5, 1, 2)
    grid <- seq(0.1, 5, by = 0.5)
    out <- simtte:::.sim_surv_df(log_hr = lp, mu = -0.7, shape = 1.4,
        type = "weibull", times = grid)
    expect_equal(sort(unique(out$lp)), sort(lp))
    expect_equal(out$p11, analytical_survival(out$time, -0.7, out$lp, 1.4),
        tolerance = 1e-8)
})

test_that("the public 'time' argument controls the mrgsolve output grid", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    grid <- seq(0.1, 1, by = 0.1)
    out <- simtte:::.sim_surv_df(log_hr = 0.1, mu = -1, shape = 1.1,
        type = "weibull", times = grid)
    # This is the regression check for the original bug report: with
    # the old implementation, mrgsolve fell back to its own default
    # output schedule (start = 0, end = 24, delta = 1) and the
    # requested 0.1-spaced grid was silently ignored.
    expect_equal(sort(unique(out$time)), grid)
    expect_equal(nrow(out), length(grid))
})

test_that("different follow-up grids give distinct, deterministic output grids", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    coarse <- seq(1, 10, by = 1)
    fine <- seq(1, 10, by = 0.5)
    out_coarse <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1.2,
        type = "weibull", times = coarse)
    out_fine <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1.2,
        type = "weibull", times = fine)
    expect_equal(sort(unique(out_coarse$time)), coarse)
    expect_equal(sort(unique(out_fine$time)), fine)
    expect_true(nrow(out_fine) > nrow(out_coarse))
})

test_that("sim_tte() administratively censors exactly at end_time", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # An extremely small hazard (very negative mu) makes the survival
    # probability stay overwhelmingly close to 1 across the whole
    # grid, so censoring at end_time is the only realistic outcome;
    # this is checked structurally (rather than relying on it) via the
    # resolved output grid.
    lp <- matrix(rep(0, 5), nrow = 5)
    result <- sim_tte(pi = lp, mu = -50, coefs = 1.2,
        time = seq(0.1, 5, by = 0.1), type = "weibull", end_time = 5)
    expect_true(all(result$sim_status == 0))
    expect_true(all(result$sim_time == 5))
})

test_that("end_time truncates/extends the grid and controls censoring (regression)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # Previously, sim_tte()'s `end_time` argument was accepted but
    # never forwarded to the internal simulation call, so it had no
    # effect whenever it differed from max(time). This test fails
    # under that old behaviour: with the same `time` grid but two
    # different `end_time` values (below and above max(time)), the
    # resolved output grid -- and hence the deterministic censoring
    # time -- must differ.
    lp <- matrix(0, nrow = 3)
    grid <- seq(0.1, 5, by = 0.1)

    short <- sim_tte(pi = lp, mu = -50, coefs = 1.2, time = grid,
        type = "weibull", end_time = 2)
    long <- sim_tte(pi = lp, mu = -50, coefs = 1.2, time = grid,
        type = "weibull", end_time = 8)

    expect_true(all(short$sim_time == 2))
    expect_true(all(long$sim_time == 8))
})

test_that("end_time may extend the Weibull grid beyond max(time)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = seq(0.1, 1, by = 0.1), end_time = 5)
    expect_equal(max(out$time), 5)
    expect_true(5 %in% out$time)
})
