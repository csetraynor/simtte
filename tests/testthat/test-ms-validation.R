# Pre-Phase-B audit, Phase 4: M-spline hazard input validation.
#
# Validated before mrgsolve is invoked: `basis` numeric matrix, finite;
# nrow(basis) == length(time); coefs numeric, finite; resulting
# basehaz = basis %*% coefs finite and non-negative; duplicate times
# accepted only when the resulting hazard values agree exactly.

test_that("nrow(basis) != length(time) is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, 1),
            basis = matrix(1, nrow = 2, ncol = 2), time = c(0, 1, 2),
            type = "ms"),
        "one row per element"
    )
})

test_that("ncol(basis) != length(coefs) is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, 1, 1),
            basis = matrix(1, nrow = 2, ncol = 2), time = c(0, 1),
            type = "ms"),
        "same length"
    )
})

test_that("NA in basis is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(1, NA, 0.5, 0.5), nrow = 2)
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, 1), basis = b,
            time = c(0, 1), type = "ms"),
        "NA, NaN, Inf"
    )
})

test_that("Inf in basis is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(1, Inf, 0.5, 0.5), nrow = 2)
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, 1), basis = b,
            time = c(0, 1), type = "ms"),
        "NA, NaN, Inf"
    )
})

test_that("NA coefficient is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(1, 1, 0.5, 0.5), nrow = 2)
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, NA), basis = b,
            time = c(0, 1), type = "ms"),
        "NA, NaN, Inf"
    )
})

test_that("Inf coefficient is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(1, 1, 0.5, 0.5), nrow = 2)
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = c(1, Inf), basis = b,
            time = c(0, 1), type = "ms"),
        "NA, NaN, Inf"
    )
})

test_that("negative resulting hazard is rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(1, -1), ncol = 1)
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = 1, basis = b, time = c(0, 1),
            type = "ms"),
        "non-negative"
    )
})

test_that("zero hazard is a valid input (not rejected)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    b <- matrix(c(0, 0), ncol = 1)
    result <- sim_tte(pi = 0, mu = -1, coefs = 1, basis = b,
        time = c(0, 1), type = "ms")
    # zero hazard throughout: certain censoring, S(t) == 1 everywhere.
    expect_equal(result$sim_status, 0)
    expect_equal(result$sim_time, 1)
})

test_that("duplicate times with identical hazard values are accepted", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
        type = "ms", times = c(0, 0, 1, 2),
        basehaz = matrix(c(1, 1, 2, 4), ncol = 1), end_time = 2)
    expect_s3_class(out, "data.frame")
    expect_equal(sort(unique(out$time)), c(0, 1, 2))
})

test_that("duplicate times with conflicting hazard values are rejected", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(0, 0, 1),
            basehaz = matrix(c(1, 2, 3), ncol = 1), end_time = 1),
        "conflicting"
    )
})

test_that("basis must be a numeric matrix", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = 1, basis = c(1, 2),
            time = c(0, 1), type = "ms"),
        "numeric matrix"
    )
})

test_that("basehaz must be a numeric matrix (explore_pi_tq_surv path)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        explore_pi_tq_surv(pi = 0, mu = -1, type = "ms",
            times = c(0, 1), basehaz = c(1, 2), end_time = 1),
        "numeric matrix"
    )
})

test_that("bundled ms_data passes all validation as-is", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    data("ms_data", package = "simtte")
    expect_true(is.matrix(ms_data$basis))
    expect_true(is.numeric(ms_data$basis))
    expect_equal(nrow(ms_data$basis), length(ms_data$time))
    expect_true(all(is.finite(ms_data$basis)))
    expect_true(all(is.finite(ms_data$coefs)))
    basehaz <- as.numeric(ms_data$basis %*% ms_data$coefs)
    expect_true(all(basehaz >= 0))
    lp <- matrix(0, nrow = 1)
    expect_error(
        sim_tte(pi = lp, mu = ms_data$mu, basis = ms_data$basis,
            coefs = ms_data$coefs, time = ms_data$time, type = "ms"),
        NA
    )
})
