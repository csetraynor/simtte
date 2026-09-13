# Phase G: time-varying lp(t) for the built-in Weibull and M-spline
# models (sim_tte()'s `lp_data` argument). See the "Time-varying lp(t)"
# section of ?sim_tte and PHASE_G_REPORT.md for the full contract.
#
# Categories covered in this file:
#   1. lp_data = NULL reproduces pre-Phase-G behavior exactly
#   2. population-level lp_data (no ID column)
#   3. subject-specific lp_data (ID column)
#   4. Weibull piecewise-constant analytical agreement
#   5. M-spline + time-varying-lp analytical agreement
#   6. validation of malformed lp_data
#   7. both event_time_method values with lp_data
#   8. lp_data cannot leak through '...'

# ---- 1. lp_data = NULL: byte-for-byte unchanged behavior ----

test_that("omitting lp_data and passing lp_data = NULL give identical .sim_surv_df() output", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    grid <- seq(0.1, 2, by = 0.1)
    out_default <- simtte:::.sim_surv_df(log_hr = c(-0.3, 0.4), mu = -1,
        shape = 1.3, type = "weibull", times = grid, end_time = 2)
    out_explicit <- simtte:::.sim_surv_df(log_hr = c(-0.3, 0.4), mu = -1,
        shape = 1.3, type = "weibull", times = grid, end_time = 2,
        lp_data = NULL)
    expect_identical(out_default, out_explicit)
})

test_that("sim_tte() with lp_data = NULL reproduces the pre-Phase-G Weibull trajectory exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(-0.2, 0, 0.2), nrow = 3)
    set.seed(42)
    r1 <- sim_tte(pi = lp, mu = -1, coefs = 1.2, time = seq(0.1, 5, by = 0.1),
        type = "weibull", end_time = 5)
    set.seed(42)
    r2 <- sim_tte(pi = lp, mu = -1, coefs = 1.2, time = seq(0.1, 5, by = 0.1),
        type = "weibull", end_time = 5, lp_data = NULL)
    expect_identical(r1, r2)
})

test_that("sim_tte() with lp_data = NULL reproduces the pre-Phase-G M-spline trajectory exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    data("ms_data", package = "simtte")
    lp <- matrix(0, nrow = 2)
    set.seed(7)
    r1 <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms")
    set.seed(7)
    r2 <- sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
        coefs = ms_data$coefs, time = ms_data$time, type = "ms",
        lp_data = NULL)
    expect_identical(r1, r2)
})

# ---- 2. Population-level lp_data (no ID column) ----

test_that(".canonicalize_lp_data() recycles a population-level trajectory to every subject", {
    lp_data <- data.frame(time = c(0, 1, 2), lp = c(0, 0.5, -0.5))
    out <- simtte:::.canonicalize_lp_data(lp_data, n_subjects = 3)
    expect_equal(sort(unique(out$ID)), 1:3)
    for (id in 1:3) {
        sub <- out[out$ID == id, ]
        expect_equal(sub$time, lp_data$time)
        expect_equal(sub$lp, lp_data$lp)
    }
})

test_that("population-level lp_data gives identical p11 trajectories across subjects with equal pi", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(time = c(0, 1, 2), lp = c(0, 1, -1))
    lp_canonical <- simtte:::.canonicalize_lp_data(lp_data, n_subjects = 2)
    # log_hr already folded into lp, per sim_tte()'s own convention; use
    # equal baseline offsets so both subjects share one trajectory.
    lp_canonical$lp <- lp_canonical$lp + 0
    out <- simtte:::.sim_surv_df(log_hr = c(0, 0), mu = -1, shape = 1.2,
        type = "weibull", times = c(0.5, 1, 1.5, 2), end_time = 2,
        lp_data = lp_canonical)
    p1 <- out$p11[out$ID == 1]
    p2 <- out$p11[out$ID == 2]
    expect_equal(p1, p2, tolerance = 1e-10)
})

test_that("sim_tte() accepts population-level lp_data end to end for both models", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(time = c(0, 1, 2, 3), lp = c(0, 0.5, 1, -0.3))
    lp <- matrix(c(-0.1, 0, 0.1, 0.2), nrow = 4)
    r_wb <- sim_tte(pi = lp, mu = -1, coefs = 1.1, time = seq(0.1, 3, by = 0.1),
        type = "weibull", end_time = 3, lp_data = lp_data)
    expect_equal(nrow(r_wb), 4)
    expect_true(all(r_wb$sim_status %in% c(0, 1)))

    basis <- matrix(c(1, 1, 1, 1), ncol = 1)
    r_ms <- sim_tte(pi = lp, mu = 0, basis = basis, coefs = 1,
        time = c(0, 1, 2, 3), type = "ms", end_time = 3,
        lp_data = lp_data)
    expect_equal(nrow(r_ms), 4)
    expect_true(all(r_ms$sim_status %in% c(0, 1)))
})

# ---- 3. Subject-specific lp_data (ID column) ----

test_that(".canonicalize_lp_data() preserves subject-specific trajectories", {
    lp_data <- data.frame(ID = c(1, 1, 2, 2), time = c(0, 1, 0, 1),
        lp = c(0, 1, 0, -1))
    out <- simtte:::.canonicalize_lp_data(lp_data, n_subjects = 2)
    expect_equal(out$ID, lp_data$ID)
    expect_equal(out$time, lp_data$time)
    expect_equal(out$lp, lp_data$lp)
})

test_that("subject-specific lp_data produces distinct p11 trajectories per subject", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(ID = c(1, 1, 2, 2), time = c(0, 1, 0, 1),
        lp = c(0, 2, 0, -2))
    lp_canonical <- simtte:::.canonicalize_lp_data(lp_data, n_subjects = 2)
    out <- simtte:::.sim_surv_df(log_hr = c(0, 0), mu = -1, shape = 1,
        type = "weibull", times = c(0.5, 1, 2), end_time = 2,
        lp_data = lp_canonical)
    p1 <- out$p11[out$ID == 1 & out$time == 2]
    p2 <- out$p11[out$ID == 2 & out$time == 2]
    expect_false(isTRUE(all.equal(p1, p2)))
    # Subject 1's lp jumps up (higher hazard) -> lower survival than
    # subject 2's, whose lp jumps down.
    expect_true(p1 < p2)
})

test_that("sim_tte() combines subject-specific lp_data with each subject's own pi offset", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(ID = c(1, 2), time = c(0, 0), lp = c(0, 0))
    lp <- matrix(c(-5, 5), nrow = 2)
    # With mu very negative and lp_data constant at 0 for both subjects,
    # the effective linear predictor is dominated by pi: subject 1
    # (pi = -5) should essentially never have an event; subject 2
    # (pi = 5) should essentially always have one, well before end_time.
    r <- sim_tte(pi = lp, mu = -3, coefs = 1, time = seq(0, 10, by = 0.1),
        type = "weibull", end_time = 10, lp_data = lp_data)
    expect_equal(r$sim_status[r$ID == 1], 0)
    expect_equal(r$sim_status[r$ID == 2], 1)
})

test_that("mismatched subject-specific IDs raise an informative error", {
    lp_data <- data.frame(ID = c(1, 3), time = c(0, 0), lp = c(0, 0))
    expect_error(
        simtte:::.canonicalize_lp_data(lp_data, n_subjects = 2),
        "one trajectory per"
    )
})

# ---- 4. Weibull piecewise-constant analytical agreement ----

test_that("weibull_tv matches the hand-derived piecewise-constant-lp cumulative hazard", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- -1
    shape <- 1.5
    # lp(t): 0 on [0, 1), 1 on [1, 2), -0.5 on [2, 3).
    lp_data <- data.frame(ID = 1, time = c(0, 1, 2, 3),
        lp = c(0, 1, -0.5, -0.5))
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = mu, shape = shape,
        type = "weibull", times = c(1, 2, 3), end_time = 3,
        lp_data = lp_data)

    H1 <- exp(mu + 0) * (1^shape - 0^shape)
    H2 <- H1 + exp(mu + 1) * (2^shape - 1^shape)
    H3 <- H2 + exp(mu - 0.5) * (3^shape - 2^shape)
    ref <- setNames(exp(-c(H1, H2, H3)), c("1", "2", "3"))

    got <- setNames(out$p11, as.character(out$time))
    expect_equal(unname(got[names(ref)]), unname(ref), tolerance = 1e-8)
})

test_that("weibull_tv reduces to the closed-form constant-lp baseline when lp(t) is constant", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- -0.7
    shape <- 0.9
    lp_const <- 0.3
    grid <- seq(0.1, 4, by = 0.3)
    lp_data <- data.frame(ID = 1, time = c(0, grid), lp = lp_const)

    out_tv <- simtte:::.sim_surv_df(log_hr = 0, mu = mu, shape = shape,
        type = "weibull", times = grid, end_time = max(grid),
        lp_data = lp_data)
    out_baseline <- simtte:::.sim_surv_df(log_hr = lp_const, mu = mu,
        shape = shape, type = "weibull", times = grid)

    got <- setNames(out_tv$p11, as.character(out_tv$time))
    ref <- setNames(out_baseline$p11, as.character(out_baseline$time))
    expect_equal(unname(got[as.character(grid)]),
        unname(ref[as.character(grid)]), tolerance = 1e-8)
})

test_that("weibull_tv handles t = 0 explicitly with p11 = 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(ID = 1, time = c(0, 1), lp = c(0, 1))
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = c(0, 1), end_time = 1,
        lp_data = lp_data)
    expect_equal(out$p11[out$time == 0], 1)
})

test_that("weibull_tv's internal grid need not match the requested output grid", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # lp_data knots (0, 0.3, 1.7) are decoupled from the requested output
    # grid (0.5, 1, 1.5, 2): the returned rows must still be exactly the
    # requested grid, and the segment ending at each still reflects
    # every knot that falls before it.
    lp_data <- data.frame(ID = 1, time = c(0, 0.3, 1.7), lp = c(0, 1, -1))
    grid <- c(0.5, 1, 1.5, 2)
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = grid, end_time = 2, lp_data = lp_data)
    expect_equal(sort(unique(out$time)), grid)

    mu <- -1
    H <- function(t) {
        segs <- c(0, 0.3, 1.7, t)
        lps <- c(0, 1, -1)
        total <- 0
        for (j in seq_along(lps)) {
            lo <- segs[j]; hi <- min(segs[j + 1], t)
            if (hi > lo) total <- total + exp(mu + lps[j]) * (hi - lo)
        }
        total
    }
    ref <- exp(-vapply(grid, H, numeric(1)))
    got <- setNames(out$p11, as.character(out$time))
    expect_equal(unname(got[as.character(grid)]), ref, tolerance = 1e-8)
})

# ---- 5. M-spline + time-varying-lp analytical agreement ----

test_that("ms + lp_data matches h(t) = basehaz(t) * exp(mu + lp(t))", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- 0
    basehaz <- matrix(c(1, 2, 4), ncol = 1)
    lp_data <- data.frame(ID = 1, time = c(0, 1, 2), lp = c(0, 1, 1))
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = mu, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = basehaz, end_time = 2,
        lp_data = lp_data)

    H1 <- 1 * exp(mu + 0) * 1
    H2 <- H1 + 2 * exp(mu + 1) * 1
    ref <- setNames(c(1, exp(-H1), exp(-H2)), c("0", "1", "2"))
    got <- setNames(out$p11, as.character(out$time))
    expect_equal(unname(got[names(ref)]), unname(ref), tolerance = 1e-6)
})

test_that("ms + lp_data reduces to the constant-lp baseline when lp(t) is constant", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    basehaz <- matrix(c(1, 2, 4), ncol = 1)
    lp_const <- 0.4
    lp_data <- data.frame(ID = 1, time = c(0, 1, 2), lp = lp_const)
    out_tv <- simtte:::.sim_surv_df(log_hr = 0, mu = -0.2, shape = NULL,
        type = "ms", times = c(0, 1, 2), basehaz = basehaz, end_time = 2,
        lp_data = lp_data)
    out_baseline <- simtte:::.sim_surv_df(log_hr = lp_const, mu = -0.2,
        shape = NULL, type = "ms", times = c(0, 1, 2), basehaz = basehaz,
        end_time = 2)
    expect_equal(out_tv$p11, out_baseline$p11, tolerance = 1e-6)
})

test_that("ms + lp_data on an independently timed lp grid is LOCF-merged correctly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mu <- 0
    basehaz <- matrix(c(1, 1), ncol = 1) # constant hazard = 1 on [0, 2)
    times <- c(0, 2)
    # lp(t) changes at t = 0.5 and t = 1.5, knots not shared with `times`.
    lp_data <- data.frame(ID = 1, time = c(0, 0.5, 1.5), lp = c(0, 1, -1))
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = mu, shape = NULL,
        type = "ms", times = times, basehaz = basehaz, end_time = 2,
        lp_data = lp_data)

    H2 <- 1 * exp(0) * 0.5 + 1 * exp(1) * 1 + 1 * exp(-1) * 0.5
    expect_equal(out$p11[out$time == 0], 1, tolerance = 1e-6)
    expect_equal(out$p11[out$time == 2], exp(-H2), tolerance = 1e-6)
})

# ---- 6. Validation of malformed lp_data ----

test_that("duplicate lp_data times with conflicting lp values are rejected", {
    lp_data <- data.frame(ID = 1, time = c(0, 1, 1), lp = c(0, 1, 2))
    expect_error(
        simtte:::.validate_lp_data_trajectories(lp_data),
        "conflicting lp values"
    )
})

test_that("duplicate lp_data times with identical lp values are accepted", {
    lp_data <- data.frame(ID = 1, time = c(0, 1, 1), lp = c(0, 1, 1))
    expect_true(simtte:::.validate_lp_data_trajectories(lp_data))
})

test_that("unsorted lp_data times are rejected, not silently sorted", {
    lp_data <- data.frame(ID = 1, time = c(0, 2, 1), lp = c(0, 1, 2))
    expect_error(
        simtte:::.validate_lp_data_trajectories(lp_data),
        "sorted in ascending order"
    )
})

test_that("lp_data missing an observation at time = 0 is rejected", {
    lp_data <- data.frame(ID = 1, time = c(0.5, 1), lp = c(0, 1))
    expect_error(
        simtte:::.check_lp_data_coverage(lp_data, end_time = 1,
            strict_coverage = FALSE),
        "time = 0"
    )
})

test_that("non-finite lp_data values are rejected", {
    expect_error(
        simtte:::.canonicalize_lp_data(
            data.frame(time = c(0, 1), lp = c(0, NA_real_)), n_subjects = 1),
        "NA, NaN, Inf, or -Inf"
    )
    expect_error(
        simtte:::.canonicalize_lp_data(
            data.frame(time = c(0, Inf), lp = c(0, 1)), n_subjects = 1),
        "NA, NaN, Inf, or -Inf"
    )
})

test_that("negative lp values are allowed (unlike basehaz)", {
    out <- simtte:::.canonicalize_lp_data(
        data.frame(time = c(0, 1), lp = c(-5, -0.001)), n_subjects = 1)
    expect_equal(out$lp, c(-5, -0.001))
})

test_that("insufficient M-spline coverage (lp_data does not reach end_time) is rejected", {
    lp_data <- data.frame(ID = 1, time = c(0, 1), lp = c(0, 1))
    expect_error(
        simtte:::.check_lp_data_coverage(lp_data, end_time = 2,
            strict_coverage = TRUE),
        "follow-up horizon"
    )
    # The same trajectory is fine for Weibull, which carries forward.
    expect_true(
        simtte:::.check_lp_data_coverage(lp_data, end_time = 2,
            strict_coverage = FALSE)
    )
})

# ---------------------------------------------------------------------
# Phase 2.5 (reports/04_author_decisions.md "After the test runbook /
# Phase 2.5", decision 3): `.check_lp_data_coverage()`'s `type` string
# was replaced by an explicit `strict_coverage` boolean naming the
# actual behaviour. Direct unit test of both values, independent of
# which model happens to map to which value at the sim_tte()/
# sim_tte_ode() call sites (those call sites are exercised above and in
# test-sim-tte-ode-covariates.R).
# ---------------------------------------------------------------------
test_that(".check_lp_data_coverage()'s strict_coverage argument controls only the end_time-reach requirement", {
    lp_data <- data.frame(ID = 1, time = c(0, 1), lp = c(0, 1))
    # strict_coverage = FALSE: short trajectories are fine (LOCF).
    expect_true(simtte:::.check_lp_data_coverage(lp_data, end_time = 100,
        strict_coverage = FALSE))
    # strict_coverage = TRUE: same trajectory, now rejected.
    expect_error(simtte:::.check_lp_data_coverage(lp_data, end_time = 100,
        strict_coverage = TRUE), "follow-up horizon")
    # Both values still enforce the time = 0 requirement identically.
    no_zero <- data.frame(ID = 1, time = c(0.5, 1), lp = c(0, 1))
    expect_error(simtte:::.check_lp_data_coverage(no_zero, end_time = 1,
        strict_coverage = TRUE), "time = 0")
    expect_error(simtte:::.check_lp_data_coverage(no_zero, end_time = 1,
        strict_coverage = FALSE), "time = 0")
})

test_that("sim_tte() propagates lp_data validation errors for both model types", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(0, nrow = 2)
    bad <- data.frame(time = c(1, 0), lp = c(0, 1)) # unsorted
    expect_error(
        sim_tte(pi = lp, mu = -1, coefs = 1, time = seq(0, 2, by = 0.5),
            type = "weibull", end_time = 2, lp_data = bad),
        "sorted in ascending order"
    )
    basis <- matrix(c(1, 1, 1), ncol = 1)
    short <- data.frame(time = c(0, 1), lp = c(0, 1)) # doesn't reach end_time
    expect_error(
        sim_tte(pi = lp, mu = 0, basis = basis, coefs = 1,
            time = c(0, 1, 2), type = "ms", end_time = 2, lp_data = short),
        "follow-up horizon"
    )
})

# ---- 7. Both event_time_method values with lp_data ----

test_that("event_time_method = 'grid' and 'log_survival' agree on event/censoring status with lp_data (Weibull)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(time = c(0, 1, 2), lp = c(0, 1, -0.5))
    lp <- matrix(rep(0, 5), nrow = 5)
    set.seed(1)
    r_grid <- sim_tte(pi = lp, mu = -1, coefs = 1.5,
        time = seq(0.1, 5, by = 0.1), type = "weibull", end_time = 5,
        lp_data = lp_data, event_time_method = "grid")
    set.seed(1)
    r_log <- sim_tte(pi = lp, mu = -1, coefs = 1.5,
        time = seq(0.1, 5, by = 0.1), type = "weibull", end_time = 5,
        lp_data = lp_data, event_time_method = "log_survival")
    expect_identical(r_grid$sim_status, r_log$sim_status)
    # log_survival must only ever refine within-grid crossings, never
    # move a status-defining boundary.
    events <- r_grid$sim_status == 1
    expect_true(all(r_log$sim_time[events] <= r_grid$sim_time[events]))
    expect_true(all(r_log$sim_time[events] >
        r_grid$sim_time[events] - diff(range(seq(0.1, 5, by = 0.1)[1:2]))))
})

test_that("event_time_method = 'grid' and 'log_survival' agree on event/censoring status with lp_data (M-spline)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(time = c(0, 1, 2), lp = c(0, 1, 1))
    basis <- matrix(c(1, 1, 1), ncol = 1)
    lp <- matrix(rep(0, 4), nrow = 4)
    set.seed(2)
    r_grid <- sim_tte(pi = lp, mu = 0, basis = basis, coefs = 1,
        time = c(0, 1, 2), type = "ms", end_time = 2, lp_data = lp_data,
        event_time_method = "grid")
    set.seed(2)
    r_log <- sim_tte(pi = lp, mu = 0, basis = basis, coefs = 1,
        time = c(0, 1, 2), type = "ms", end_time = 2, lp_data = lp_data,
        event_time_method = "log_survival")
    expect_identical(r_grid$sim_status, r_log$sim_status)
})

# ---- 8. lp_data cannot leak through '...' ----

test_that("lp_data is a formal argument and cannot be supplied a second time via '...'", {
    lp <- matrix(0, nrow = 1)
    expect_error(
        do.call(sim_tte, list(pi = lp, mu = -1, coefs = 1, time = c(0, 1),
            end_time = 1, lp_data = data.frame(time = 0, lp = 0),
            lp_data = data.frame(time = 0, lp = 1))),
        "formal argument"
    )
})

test_that("mrgsim()-reserved arguments are still rejected when lp_data is also supplied", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(time = 0, lp = 0)
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = c(0, 1), end_time = 1,
            lp_data = data.frame(ID = 1, time = 0, lp = 0), tgrid = c(0, 1)),
        "'tgrid'"
    )
})

test_that("legitimate '...' arguments still reach mrgsim() when lp_data is supplied", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp_data <- data.frame(ID = 1, time = 0, lp = 0)
    result <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = c(0, 1), end_time = 1,
        lp_data = lp_data, recsort = 1)
    expect_s3_class(result, "data.frame")
    expect_equal(sort(unique(result$time)), c(0, 1))
})
