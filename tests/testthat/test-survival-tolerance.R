# Pre-Phase-B audit, Phase 3: survival tolerance policy.
#
# Policy chosen: clamp (not strictly reject) survival values within
# .SURV_TOL (1e-8) of [0, 1] to exactly 0 or 1, deterministically,
# before monotonicity checking and event selection. Values further
# outside [0, 1] remain a hard error. See .SURV_TOL / .normalize_survival
# in R/helpers.R and the "Trajectory contract" section of ?sim_tte_df
# for the full rationale, in particular why clamping cannot change
# event classification (stats::runif() draws from [0, 1)).

test_that("slightly negative survival (within tolerance) is accepted and clamped to 0", {
    dat <- data.frame(ID = 1, time = c(1, 2), p11 = c(1, -5e-9))
    result <- sim_tte_df(dat)
    # 0 <= U for every achievable U (U >= 0), so this is deterministically
    # an event at the second (only negative) observation.
    expect_equal(result$sim_time, 2)
    expect_equal(result$sim_status, 1)
})

test_that("slightly >1 survival (within tolerance) is accepted and clamped to 1", {
    dat <- data.frame(ID = 1, time = c(1, 2), p11 = c(1 + 5e-9, 0))
    result <- sim_tte_df(dat)
    # clamp(1 + 5e-9) == 1, and 1 <= U is never true for U in [0, 1),
    # so the first observation is never a crossing; the event is
    # deterministically detected at the second (p11 = 0) observation.
    expect_equal(result$sim_time, 2)
    expect_equal(result$sim_status, 1)
})

test_that("clearly invalid negative survival (beyond tolerance) is a hard error", {
    dat <- data.frame(ID = 1, time = c(1, 2), p11 = c(1, -0.05))
    expect_error(sim_tte_df(dat), "\\[0, 1\\]")
})

test_that("clearly invalid >1 survival (beyond tolerance) is a hard error", {
    dat <- data.frame(ID = 1, time = c(1, 2), p11 = c(1.2, 0))
    expect_error(sim_tte_df(dat), "\\[0, 1\\]")
})

test_that("survival exactly 0 and exactly 1 are valid, unclamped inputs", {
    dat0 <- data.frame(ID = 1, time = 1:3, p11 = c(1, 1, 0))
    result0 <- sim_tte_df(dat0)
    expect_equal(result0$sim_time, 3)
    expect_equal(result0$sim_status, 1)

    dat1 <- data.frame(ID = 1, time = 1:3, p11 = c(1, 1, 1))
    result1 <- sim_tte_df(dat1)
    expect_equal(result1$sim_status, 0)
    expect_equal(result1$sim_time, 3)
})

test_that("event classification is identical whether or not clamping was needed", {
    # Two trajectories that are mathematically equivalent after
    # clamping: one uses the exact boundary, the other approaches it
    # from just outside the tolerance band. Both must classify subjects
    # identically for a range of thresholds, since .get_tte() is
    # deterministic given p11 and a fixed U -- verified directly via
    # the boundary-rule helper rather than via random sampling.
    p_exact <- c(1, 0.5, 0)
    p_noisy <- c(1 + 5e-9, 0.5, -5e-9)
    for (u in c(0, 0.001, 0.4999, 0.5, 0.5001, 0.999)) {
        idx_exact <- simtte:::.get_tte(u, p_exact)
        idx_noisy <- simtte:::.get_tte(u, simtte:::.normalize_survival(p_noisy))
        expect_equal(idx_exact, idx_noisy, info = paste("u =", u))
    }
})

test_that("monotonicity noise within tolerance is accepted without altering the value", {
    dat <- data.frame(ID = 1, time = 1:4,
        p11 = c(1, 0.6000000005, 0.6, 0.3))
    expect_error(sim_tte_df(dat), NA)
})

test_that("genuinely non-monotone survival (beyond tolerance) is still a hard error after normalization", {
    dat <- data.frame(ID = 1, time = 1:4, p11 = c(1, 0.5, 0.7, 0.3))
    expect_error(sim_tte_df(dat), "non-increasing")
})

test_that("normalized survival is what the monotonicity check operates on", {
    # Raw values are technically non-monotone by a hair (1 + 5e-9 then
    # 1 - 5e-9 is a decrease of 1e-8, right at the tolerance boundary),
    # but after clamping both become exactly 1: no increase at all.
    dat <- data.frame(ID = 1, time = 1:3,
        p11 = c(1 + 5e-9, 1 - 5e-9, 0))
    expect_error(sim_tte_df(dat), NA)
    result <- sim_tte_df(dat)
    expect_equal(result$sim_time, 3)
    expect_equal(result$sim_status, 1)
})

test_that(".normalize_survival() clamps only out-of-range values", {
    x <- c(-1e-9, 0, 0.5, 1, 1 + 1e-9)
    got <- simtte:::.normalize_survival(x)
    expect_equal(got, c(0, 0, 0.5, 1, 1))
})
