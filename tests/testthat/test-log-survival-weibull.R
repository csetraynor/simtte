# Phase B: "log_survival" against the Weibull closed form.
#
# The analytical Weibull event time T = (-log(U) / eta)^(1/shape),
# eta = exp(mu + lp), is used ONLY as a validation reference here. The
# package implementation (.interpolate_log_survival()) never evaluates
# this formula and remains fully model-agnostic; these tests exercise it
# against Weibull-generated trajectories purely to characterize expected
# exactness/approximation behavior (see PHASE_B_DESIGN_AUDIT.md §7-8).

weibull_T <- function(U, mu, lp, shape) {
    eta <- exp(mu + lp)
    (-log(U) / eta)^(1 / shape)
}

# Given a fixed U and a Weibull trajectory (times, p11) from
# .sim_surv_df(), locate the bracketing pair and interpolate.
interpolate_from_trajectory <- function(times, p11, U) {
    i <- max(which(p11 > U))
    ip1 <- i + 1L
    stopifnot(ip1 <= length(times), p11[ip1] <= U)
    simtte:::.interpolate_log_survival(t_i = times[i], t_ip1 = times[ip1],
        s_i = p11[i], s_ip1 = p11[ip1], u = U)
}

test_that("shape = 1: log_survival interpolation exactly matches the analytical Weibull time", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- -1; lp <- 0.2; shape <- 1
    grid <- seq(0, 5, by = 0.5)
    out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
        type = "weibull", times = grid)

    for (U in c(0.8, 0.5, 0.2, 0.05)) {
        T_true <- weibull_T(U, mu, lp, shape)
        if (T_true <= max(grid)) {
            got <- interpolate_from_trajectory(out$time, out$p11, U)
            expect_equal(got, T_true, tolerance = 1e-6,
                info = paste("U =", U))
        }
    }
})

test_that("shape != 1: log_survival interpolation approximates the analytical Weibull time", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    for (shape in c(0.6, 1.8)) {
        mu <- -1; lp <- 0.2
        grid <- seq(0, 5, by = 0.5)
        out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
            type = "weibull", times = grid)
        U <- 0.5
        T_true <- weibull_T(U, mu, lp, shape)
        if (T_true <= max(grid)) {
            got <- interpolate_from_trajectory(out$time, out$p11, U)
            # Approximation, not exact: allow a grid-spacing-scale
            # tolerance, not floating-point tolerance.
            expect_equal(got, T_true, tolerance = 0.5,
                info = paste("shape =", shape))
            # ...but it should not be *coincidentally* exact either --
            # confirm this really is an approximation for a coarse grid.
            expect_gt(abs(got - T_true), 1e-6)
        }
    }
})

test_that("refining the grid reduces the log_survival interpolation error (shape != 1)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- -1; lp <- 0.2; shape <- 2; U <- 0.5
    T_true <- weibull_T(U, mu, lp, shape)

    error_for_grid <- function(grid) {
        out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
            type = "weibull", times = grid)
        got <- interpolate_from_trajectory(out$time, out$p11, U)
        abs(got - T_true)
    }

    coarse <- error_for_grid(seq(0, 5, by = 1))
    medium <- error_for_grid(seq(0, 5, by = 0.25))
    fine <- error_for_grid(seq(0, 5, by = 0.02))

    # Robust monotone-decreasing check (not a strict convergence-rate
    # assertion, which would be fragile): each refinement must not make
    # things worse, and the finest grid must be substantially better
    # than the coarsest.
    expect_lte(medium, coarse)
    expect_lte(fine, medium)
    expect_lt(fine, coarse / 10)
})

test_that("shape = 1 interpolation error does not shrink further with refinement (already exact)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- -1; lp <- 0.2; shape <- 1; U <- 0.5
    T_true <- weibull_T(U, mu, lp, shape)

    error_for_grid <- function(grid) {
        out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
            type = "weibull", times = grid)
        got <- interpolate_from_trajectory(out$time, out$p11, U)
        abs(got - T_true)
    }

    coarse <- error_for_grid(seq(0, 5, by = 1))
    fine <- error_for_grid(seq(0, 5, by = 0.02))
    expect_lt(coarse, 1e-6)
    expect_lt(fine, 1e-6)
})

test_that("public sim_tte() pipeline runs with log_survival for shape < 1, = 1, > 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(-0.3, 0, 0.3), nrow = 3)
    for (shape in c(0.6, 1, 1.8)) {
        result <- sim_tte(pi = lp, mu = -1, coefs = shape,
            time = seq(0.1, 10, by = 0.1), type = "weibull", end_time = 10,
            event_time_method = "log_survival")
        expect_true(all(is.finite(result$sim_time)))
        expect_equal(nrow(result), 3)
    }
})

test_that("coarse vs fine grid through the public pipeline both run cleanly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(0, nrow = 1)
    coarse <- sim_tte(pi = lp, mu = -1, coefs = 1.5, time = seq(0, 10, by = 2),
        type = "weibull", end_time = 10, event_time_method = "log_survival")
    fine <- sim_tte(pi = lp, mu = -1, coefs = 1.5, time = seq(0, 10, by = 0.05),
        type = "weibull", end_time = 10, event_time_method = "log_survival")
    expect_true(is.finite(coarse$sim_time))
    expect_true(is.finite(fine$sim_time))
})
