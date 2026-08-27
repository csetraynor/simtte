# Pre-Phase-B audit, Phase 1: package-controlled mrgsim() arguments must
# not be overridable via `...`, through either the public sim_tte() /
# explore_pi_tq_surv() entry points or the shared internal
# .sim_surv_df() funnel.

test_that(".sim_surv_df() rejects a conflicting 'tgrid' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = 1, basehaz = NULL, end_time = 1,
            tgrid = c(0, 0.5, 1)),
        "'tgrid'"
    )
})

test_that(".sim_surv_df() rejects a conflicting 'obsonly' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = 1, basehaz = NULL, end_time = 1,
            obsonly = FALSE),
        "'obsonly'"
    )
})

test_that(".sim_surv_df() rejects a conflicting 'nocb' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = 0, shape = NULL,
            type = "ms", times = c(0, 1, 2),
            basehaz = matrix(c(1, 2, 4), ncol = 1), end_time = 2,
            nocb = TRUE),
        "'nocb'"
    )
})

test_that(".sim_surv_df() rejects a conflicting 'carry_out' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = 1, basehaz = NULL, end_time = 1,
            carry_out = "foo"),
        "'carry_out'"
    )
})

test_that(".sim_surv_df() rejects a conflicting 'carry.out' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = 1, basehaz = NULL, end_time = 1,
            carry.out = "foo"),
        "'carry.out'"
    )
})

test_that(".sim_surv_df() rejects a conflicting 'data' passed via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
            type = "weibull", times = 1, basehaz = NULL, end_time = 1,
            data = data.frame(ID = 99)),
        "'data'"
    )
})

test_that("sim_tte() propagates the reserved-argument error from '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = 1, time = c(0, 1), tgrid = c(0, 1)),
        "'tgrid'"
    )
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = 1, time = c(0, 1),
            carry_out = "foo"),
        "'carry_out'"
    )
})

test_that("explore_pi_tq_surv() propagates the reserved-argument error from '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        explore_pi_tq_surv(pi = seq(-1, 1, by = 1), mu = -1, shape = 1.1,
            end_time = 5, type = "weibull", obsonly = FALSE),
        "'obsonly'"
    )
})

test_that("legitimate '...' arguments are still forwarded to mrgsim()", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # recsort is a genuine mrgsim() record-sorting option that does not
    # conflict with the package's grid/trajectory contract.
    result <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = c(0, 1), basehaz = NULL, end_time = 1,
        recsort = 1)
    expect_s3_class(result, "data.frame")
    expect_equal(sort(unique(result$time)), c(0, 1))

    # etasrc= (mrgsolve ETA-source option) is also unreserved.
    result2 <- simtte:::.sim_surv_df(log_hr = 0, mu = -1, shape = 1,
        type = "weibull", times = c(0, 1), basehaz = NULL, end_time = 1,
        etasrc = "omega")
    expect_s3_class(result2, "data.frame")
})

test_that("sim_tte() still works end to end with no extra '...' arguments", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(-0.2, 0, 0.2), nrow = 3)
    result <- sim_tte(pi = lp, mu = -1, coefs = 1.2,
        time = seq(0.1, 5, by = 0.1), type = "weibull", end_time = 5)
    expect_equal(nrow(result), 3)
})
