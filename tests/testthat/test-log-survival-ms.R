# Phase B: "log_survival" against the M-spline piecewise-constant-hazard
# representation established by the nocb = FALSE hazard-carry convention
# (see ?sim_tte "M-spline hazard carry convention" and
# test-ms-hazard-carry.R). Because the hazard genuinely is
# piecewise-constant on the reported grid, linear-in-H interpolation is
# expected to recover the exact event time implied by that discretized
# hazard (not just an approximation), verified here against
# hand-derived formulas.

test_that("constant hazard: interpolated event time matches the exponential formula exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    haz <- 2
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2, 3),
        basehaz = matrix(rep(haz, 4), ncol = 1), end_time = 3)

    for (U in c(0.9, 0.5, 0.1)) {
        T_true <- -log(U) / haz
        i <- max(which(out$p11 > U))
        got <- simtte:::.interpolate_log_survival(t_i = out$time[i],
            t_ip1 = out$time[i + 1L], s_i = out$p11[i],
            s_ip1 = out$p11[i + 1L], u = U)
        expect_equal(got, T_true, tolerance = 1e-6, info = paste("U =", U))
    }
})

test_that("changing piecewise-constant hazard c(1, 2, 4) on c(0, 1, 2): exact agreement with hand-derived event times", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = matrix(c(1, 2, 4),
            ncol = 1), end_time = 2)
    # h = 1 on [0, 1), h = 2 on [1, 2); S(0) = 1, S(1) = exp(-1),
    # S(2) = exp(-3) (established in test-ms-hazard-carry.R).

    # Crossing in [0, 1): H(t) = 1 * t, so T = -log(U).
    U1 <- (1 + exp(-1)) / 2 # strictly between S(1) and S(0)
    T1_true <- -log(U1)
    got1 <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = 1, s_ip1 = exp(-1), u = U1)
    expect_equal(got1, T1_true, tolerance = 1e-9)
    expect_true(got1 > 0 && got1 < 1)

    # Crossing in [1, 2): H(t) = 1 + 2 * (t - 1), so
    # T = 1 + (-log(U) - 1) / 2.
    U2 <- (exp(-1) + exp(-3)) / 2 # strictly between S(2) and S(1)
    T2_true <- 1 + (-log(U2) - 1) / 2
    got2 <- simtte:::.interpolate_log_survival(t_i = 1, t_ip1 = 2,
        s_i = exp(-1), s_ip1 = exp(-3), u = U2)
    expect_equal(got2, T2_true, tolerance = 1e-9)
    expect_true(got2 > 1 && got2 < 2)
})

test_that("irregular grid: exact agreement with hand-derived piecewise-constant-hazard event time", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    times <- c(0, 0.2, 1.5, 4)
    haz <- c(0.5, 1, 3, 3) # last value inert, per the carry convention
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = times, basehaz = matrix(haz, ncol = 1),
        end_time = 4)

    h_at <- function(t) {
        h1 <- 0.5 * min(t, 0.2)
        h2 <- if (t > 0.2) 1 * (min(t, 1.5) - 0.2) else 0
        h3 <- if (t > 1.5) 3 * (min(t, 4) - 1.5) else 0
        h1 + h2 + h3
    }
    S_at <- function(t) exp(-h_at(t))

    # Crossing strictly inside [0.2, 1.5): H(t) = H(0.2) + 1*(t-0.2).
    U <- (S_at(0.2) + S_at(1.5)) / 2
    T_true <- 0.2 + ((-log(U)) - h_at(0.2)) / 1
    got <- simtte:::.interpolate_log_survival(t_i = 0.2, t_ip1 = 1.5,
        s_i = S_at(0.2), s_ip1 = S_at(1.5), u = U)
    expect_equal(got, T_true, tolerance = 1e-9)
    expect_true(got > 0.2 && got < 1.5)
})

test_that("full sim_tte() pipeline with the ms model and log_survival runs cleanly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    data("ms_data", package = "simtte")
    lp <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))
    result <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms",
        event_time_method = "log_survival")
    expect_true(all(is.finite(result$sim_time)))
    expect_equal(nrow(result), nrow(ms_data$basis))
})

test_that("ms grid + log_survival classification matches grid method for the same seed", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    data("ms_data", package = "simtte")
    lp <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))

    set.seed(555)
    r_grid <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms",
        event_time_method = "grid")
    set.seed(555)
    r_log <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms",
        event_time_method = "log_survival")

    expect_identical(r_grid$sim_status, r_log$sim_status)
})
