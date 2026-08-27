# Pre-Phase-B audit, Phase 5: Weibull numerical robustness and
# parameter validation. The closed-form calculation itself is unchanged
# (still exact); this hardens it against overflow for extreme but
# finite inputs, and adds upfront validation so invalid input fails
# with a clear message rather than reaching mrgsolve.

analytical_survival <- function(time, mu, lp, shape) {
    eta <- exp(mu + lp)
    exp(-eta * time^shape)
}

test_that("time = 0 always gives S(0) = 1, for any shape", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    for (shape in c(0.2, 1, 3)) {
        out <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = shape,
            type = "weibull", times = 0)
        expect_equal(out$p11, 1, tolerance = 1e-12)
    }
})

test_that("time > 0 matches the analytical formula for shape < 1, = 1, > 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    for (shape in c(0.5, 1, 2)) {
        out <- simtte:::.sim_surv_df(log_hr = 0.3, mu = -1, shape = shape,
            type = "weibull", times = c(0.1, 1, 5))
        expect_equal(out$p11, analytical_survival(out$time, -1, 0.3, shape),
            tolerance = 1e-8)
    }
})

test_that("very small positive shape does not crash or produce NaN", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1e-4,
        type = "weibull", times = c(0, 0.001, 1, 100))
    expect_true(all(is.finite(out$p11)))
    expect_equal(out$p11, analytical_survival(out$time, -1, 0, 1e-4),
        tolerance = 1e-8)
})

test_that("large but finite mu no longer produces NaN or Inf (regression)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # mu = 710 alone overflows exp(mu + lp) to Inf in the naive
    # implementation, producing Inf * 0 = NaN at t = 0.
    result <- sim_tte(pi = 0, mu = 710, coefs = 1, time = c(0, 1))
    expect_true(all(is.finite(result$sim_time)))
    expect_equal(result$sim_status, 1)
    expect_equal(result$sim_time, 1) # certain, immediate event by t = 1

    out <- simtte:::.sim_surv_df(log_hr = 0, mu = 710, shape = 1,
        type = "weibull", times = c(0, 1))
    expect_equal(out$p11[out$time == 0], 1)
    expect_equal(out$p11[out$time == 1], 0)
})

test_that("large positive linear predictor behaves like large mu (regression)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = 710, mu = 0, shape = 1,
        type = "weibull", times = c(0, 1))
    expect_true(all(is.finite(out$p11)))
    expect_equal(out$p11[out$time == 0], 1)
    expect_equal(out$p11[out$time == 1], 0)
})

test_that("very negative linear predictor gives survival close to 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    out <- simtte:::.sim_surv_df(log_hr = -700, mu = 0, shape = 1,
        type = "weibull", times = c(0, 1, 1000))
    expect_true(all(is.finite(out$p11)))
    expect_true(all(out$p11 >= 1 - 1e-6))
})

test_that("NA/Inf/-Inf 'mu' is rejected with a clear message", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte(pi = 0, mu = NA_real_, coefs = 1, time = c(0, 1)),
        "finite")
    expect_error(sim_tte(pi = 0, mu = Inf, coefs = 1, time = c(0, 1)),
        "finite")
    expect_error(sim_tte(pi = 0, mu = -Inf, coefs = 1, time = c(0, 1)),
        "finite")
})

test_that("NA/Inf/-Inf 'pi' is rejected with a clear message", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte(pi = NA_real_, mu = -1, coefs = 1, time = c(0, 1)),
        "NA, NaN, Inf")
    expect_error(sim_tte(pi = Inf, mu = -1, coefs = 1, time = c(0, 1)),
        "NA, NaN, Inf")
})

test_that("empty cohort ('pi' of length 0) is rejected with a clear message", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte(pi = numeric(0), mu = -1, coefs = 1,
        time = c(0, 1)), "at least one individual")
})

test_that("invalid shape (NA, Inf, zero, negative) is rejected with a clear message", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte(pi = 0, mu = -1, coefs = NA_real_,
        time = c(0, 1)), "NA, NaN, Inf")
    expect_error(sim_tte(pi = 0, mu = -1, coefs = Inf,
        time = c(0, 1)), "NA, NaN, Inf")
    expect_error(sim_tte(pi = 0, mu = -1, coefs = 0, time = c(0, 1)),
        "positive")
    expect_error(sim_tte(pi = 0, mu = -1, coefs = -1, time = c(0, 1)),
        "positive")
})

test_that("shape < 1, = 1, > 1 all run end to end without NaN/Inf", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(-0.5, 0, 0.5), nrow = 3)
    for (shape in c(0.6, 1, 1.8)) {
        result <- sim_tte(pi = lp, mu = -1, coefs = shape,
            time = seq(0.1, 5, by = 0.1), type = "weibull", end_time = 5)
        expect_true(all(is.finite(result$sim_time)))
        expect_equal(nrow(result), 3)
    }
})
