# Phase 6: right censoring (reports/16_censoring_design.md, option A).
# add_censoring()/censoring_rate_for() work on any events-shaped data
# frame (ID/sim_time/sim_status), so most tests here build a small
# fabricated events data frame directly -- no mrgsolve compile needed,
# keeps this file fast. The sim_tte_ode(censoring = ...) integration
# tests do need mrgsolve and are guarded with skip_if_not_installed(),
# same convention as every other ODE test file.

# .fake_events() is in helper-events.R (shared with
# test-interval-censoring.R and test-visit-process.R).

# ---------------------------------------------------------------------
# 1. Malformed spec errors [fast].
# ---------------------------------------------------------------------
test_that("add_censoring() rejects a censoring spec with no 'dist'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(rate = 0.1), end = 10),
        "'dist' element")
    expect_error(add_censoring(ev, censoring = "exponential", end = 10),
        "'dist' element")
})

test_that("add_censoring() rejects an unknown dist", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "cauchy"),
        end = 10), "one of")
})

test_that("add_censoring() validates exponential's 'rate'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "exponential"),
        end = 10), "rate")
    expect_error(add_censoring(ev, censoring = list(dist = "exponential",
        rate = -1), end = 10), "rate")
    expect_error(add_censoring(ev, censoring = list(dist = "exponential",
        rate = c(1, 2)), end = 10), "rate")
})

test_that("add_censoring() validates weibull's 'shape'/'scale'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "weibull",
        scale = 5), end = 10), "shape.*scale")
    expect_error(add_censoring(ev, censoring = list(dist = "weibull",
        shape = 1.5, scale = -5), end = 10), "shape.*scale")
})

test_that("add_censoring() validates uniform's 'min'/'max'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "uniform",
        min = 5, max = 5), end = 10), "min.*max")
    expect_error(add_censoring(ev, censoring = list(dist = "uniform",
        min = -1, max = 5), end = 10), "min.*max")
})

test_that("add_censoring() validates lognormal's 'meanlog'/'sdlog'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "lognormal",
        sdlog = 0.5), end = 10), "meanlog.*sdlog")
    expect_error(add_censoring(ev, censoring = list(dist = "lognormal",
        meanlog = 2, sdlog = -1), end = 10), "meanlog.*sdlog")
})

test_that("add_censoring() validates gamma's 'shape'/'rate', and rejects 'scale'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "gamma",
        rate = 0.1), end = 10), "shape.*rate")
    expect_error(add_censoring(ev, censoring = list(dist = "gamma",
        shape = 2, rate = -0.1), end = 10), "shape.*rate")
    expect_error(add_censoring(ev, censoring = list(dist = "gamma",
        shape = 2, scale = 10), end = 10), "scale.*not accepted")
})

test_that("add_censoring() validates a user dist's 'fun'", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    expect_error(add_censoring(ev, censoring = list(dist = "user"),
        end = 10), "function\\(n\\)")
    expect_error(add_censoring(ev, censoring = list(dist = "user",
        fun = function(n) 1:(n - 1)), end = 10), "length")
    expect_error(add_censoring(ev, censoring = list(dist = "user",
        fun = function(n) rep(NA_real_, n)), end = 10),
        "non-finite or negative")
    expect_error(add_censoring(ev, censoring = list(dist = "user",
        fun = function(n) rep(-1, n)), end = 10),
        "non-finite or negative")
})

test_that("add_censoring() validates 'end' and required columns", {
    ev <- .fake_events(c(1, 2), c(1, 1))
    spec <- list(dist = "exponential", rate = 0.1)
    expect_error(add_censoring(ev, censoring = spec, end = -1), "end")
    expect_error(add_censoring(ev, censoring = spec, end = c(1, 2)), "end")
    expect_error(add_censoring(ev[, "ID", drop = FALSE], censoring = spec,
        end = 10), "sim_time")
})

# ---------------------------------------------------------------------
# 2. Each distribution runs and produces a well-formed result [fast].
# ---------------------------------------------------------------------
test_that("each dist option draws successfully and returns valid output", {
    ev <- .fake_events(sim_time = seq(1, 20, length.out = 20),
        sim_status = rep(1L, 20))
    specs <- list(
        exponential = list(dist = "exponential", rate = 0.05),
        weibull = list(dist = "weibull", shape = 1.5, scale = 15),
        uniform = list(dist = "uniform", min = 2, max = 25),
        lognormal = list(dist = "lognormal", meanlog = 2, sdlog = 0.5),
        gamma = list(dist = "gamma", shape = 2, rate = 0.1),
        user = list(dist = "user", fun = function(n) stats::rgamma(n, 2, 0.2))
    )
    for (nm in names(specs)) {
        out <- add_censoring(ev, censoring = specs[[nm]], end = 20, seed = 1)
        expect_true(all(out$sim_status %in% c(0L, 1L)))
        expect_true(all(out$sim_time >= 0 & out$sim_time <= 20))
        expect_true(all(out$sim_reason %in%
            c("event", "censored", "administrative")))
        expect_equal(nrow(out), 20)
    }
})

# ---------------------------------------------------------------------
# 3. 'end' always wins [fast].
# ---------------------------------------------------------------------
test_that("'end' caps every observed time regardless of the censoring draw", {
    ev <- .fake_events(sim_time = rep(5, 10), sim_status = rep(1L, 10))
    out <- add_censoring(ev, censoring = list(dist = "user",
        fun = function(n) rep(1e6, n)), end = 8, seed = 1)
    expect_true(all(out$sim_time <= 8))
    # Every subject's real event time (5) beat both the huge draw and
    # end (8), so every subject should still be an event.
    expect_true(all(out$sim_status == 1L))
    expect_true(all(out$sim_reason == "event"))
})

test_that("an administratively-censored subject can still be pulled earlier by C", {
    # sim_status = 0 at sim_time = end (the sim_tte_ode() convention):
    # true T_i is unknown but >= end, so a C_i < end always determines
    # the new observed time and status stays 0.
    ev <- .fake_events(sim_time = c(20, 20), sim_status = c(0L, 0L))
    out <- add_censoring(ev, censoring = list(dist = "user",
        fun = function(n) c(5, 25)), end = 20, seed = 1)
    expect_equal(out$sim_time, c(5, 20))
    expect_equal(out$sim_status, c(0L, 0L))
    expect_equal(out$sim_reason, c("censored", "administrative"))
})

# ---------------------------------------------------------------------
# 4. status/reason consistency [fast].
# ---------------------------------------------------------------------
test_that("sim_status == 1 iff sim_reason == 'event', across all dists", {
    set.seed(42)
    n <- 200
    ev <- .fake_events(sim_time = stats::runif(n, 0, 30),
        sim_status = rbinom(n, 1, 0.7))
    for (spec in list(list(dist = "exponential", rate = 0.03),
        list(dist = "weibull", shape = 0.8, scale = 20),
        list(dist = "uniform", min = 0, max = 30),
        list(dist = "lognormal", meanlog = 2.5, sdlog = 0.8),
        list(dist = "gamma", shape = 2, rate = 0.08))) {
        out <- add_censoring(ev, censoring = spec, end = 30, seed = 1)
        expect_equal(out$sim_status == 1L, out$sim_reason == "event")
        non_event <- out$sim_reason != "event"
        expect_equal(out$sim_reason[non_event] == "administrative",
            out$sim_time[non_event] == 30)
    }
})

# ---------------------------------------------------------------------
# 5. Reproducibility [fast].
# ---------------------------------------------------------------------
test_that("add_censoring() with the same seed gives identical results", {
    ev <- .fake_events(sim_time = seq(1, 30, length.out = 30),
        sim_status = rep(1L, 30))
    spec <- list(dist = "exponential", rate = 0.04)
    out1 <- add_censoring(ev, censoring = spec, end = 25, seed = 7)
    out2 <- add_censoring(ev, censoring = spec, end = 25, seed = 7)
    expect_identical(out1, out2)
})

test_that("sim_tte_ode(censoring = ...) is reproducible via its own seed", {
    skip_if_not_installed("mrgsolve")
    spec <- list(dist = "exponential", rate = 0.05)
    s1 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 20, end = 30, delta = 2, censoring = spec, seed = 3)
    s2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 20, end = 30, delta = 2, censoring = spec, seed = 3)
    expect_identical(s1$events, s2$events)
    # And the censoring draw genuinely does something relative to the
    # uncensored call with the same seed (not a no-op).
    s0 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 20, end = 30, delta = 2, seed = 3)
    expect_false(identical(s0$events, s1$events))
})

# ---------------------------------------------------------------------
# 6. sim_tte_ode() integration: $events always has sim_reason; censoring
#    updates it; other columns untouched [fast].
# ---------------------------------------------------------------------
test_that("sim_tte_ode() $events always carries sim_reason, even with no censoring spec", {
    skip_if_not_installed("mrgsolve")
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.2),
        n = 30, end = 10, delta = 1, seed = 1)
    expect_true("sim_reason" %in% names(sim$events))
    expect_equal(sim$events$sim_status == 1L, sim$events$sim_reason == "event")
    expect_true(all(sim$events$sim_reason[sim$events$sim_status == 0L] ==
        "administrative"))
})

test_that("sim_tte_ode(censoring = ...) can only pull events earlier, never remove censoring", {
    skip_if_not_installed("mrgsolve")
    spec <- list(dist = "exponential", rate = 0.1)
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.02),
        n = 50, end = 30, delta = 2, censoring = spec, seed = 5)
    expect_true(any(sim$events$sim_reason == "censored"))
    expect_true(all(sim$events$sim_time <= 30))
})

# ---------------------------------------------------------------------
# 7. Legacy sim_tte()/sim_tte_df() output goes through add_censoring()
#    with other columns unchanged [fast].
# ---------------------------------------------------------------------
test_that("add_censoring() on sim_tte() output leaves other columns untouched", {
    set.seed(1)
    lp <- matrix(rnorm(20, 0, 0.5), nrow = 20)
    legacy <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
        time = seq(0, 20, by = 1), type = "weibull", end_time = 20)
    out <- add_censoring(legacy, censoring = list(dist = "exponential",
        rate = 0.05), end = 20, seed = 2)
    expect_identical(out$lp, legacy$lp)
    expect_identical(out$ID, legacy$ID)
    expect_true("sim_reason" %in% names(out))
    expect_false("sim_reason" %in% names(legacy))
})

test_that("no-spec sim_tte_ode() output is identical apart from the new sim_reason column", {
    skip_if_not_installed("mrgsolve")
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 15, end = 10, delta = 1, seed = 1)
    without_reason <- sim$events
    without_reason$sim_reason <- NULL
    expect_identical(without_reason,
        sim$events[, c("ID", "sim_time", "sim_status")])
})

# ---------------------------------------------------------------------
# 8. Slow: censoring_rate_for() hits the target fraction; distributional
#    check of censored times.
# ---------------------------------------------------------------------
test_that("censoring_rate_for() solves an exponential rate hitting the target fraction [slow]", {
    skip_if_not_slow()
    set.seed(1)
    times <- stats::rexp(3000, rate = 0.03) # stand-in for an uncensored run's event times
    target <- 0.25
    rate <- censoring_rate_for(times, target = target, dist = "exponential")
    ev <- .fake_events(sim_time = times, sim_status = rep(1L, length(times)))
    out <- add_censoring(ev, censoring = list(dist = "exponential",
        rate = rate), end = max(times) + 1, seed = 2)
    observed <- mean(out$sim_status == 0L)
    expect_equal(observed, target, tolerance = 0.06)
})

test_that("censoring_rate_for() solves a weibull scale hitting the target fraction [slow]", {
    skip_if_not_slow()
    set.seed(1)
    times <- stats::rweibull(3000, shape = 1.3, scale = 25)
    target <- 0.15
    scale <- censoring_rate_for(times, target = target, dist = "weibull",
        shape = 2)
    ev <- .fake_events(sim_time = times, sim_status = rep(1L, length(times)))
    out <- add_censoring(ev, censoring = list(dist = "weibull", shape = 2,
        scale = scale), end = max(times) + 1, seed = 2)
    observed <- mean(out$sim_status == 0L)
    expect_equal(observed, target, tolerance = 0.06)
})

test_that("censoring_rate_for() solves a lognormal meanlog hitting the target fraction [slow]", {
    skip_if_not_slow()
    set.seed(1)
    times <- stats::rlnorm(3000, meanlog = 2.5, sdlog = 0.8)
    target <- 0.2
    meanlog <- censoring_rate_for(times, target = target, dist = "lognormal",
        sdlog = 0.6)
    ev <- .fake_events(sim_time = times, sim_status = rep(1L, length(times)))
    out <- add_censoring(ev, censoring = list(dist = "lognormal",
        meanlog = meanlog, sdlog = 0.6), end = max(times) + 1, seed = 2)
    observed <- mean(out$sim_status == 0L)
    expect_equal(observed, target, tolerance = 0.06)
})

test_that("censoring_rate_for() solves a gamma rate hitting the target fraction [slow]", {
    skip_if_not_slow()
    set.seed(1)
    times <- stats::rgamma(3000, shape = 2, rate = 0.05)
    target <- 0.3
    rate <- censoring_rate_for(times, target = target, dist = "gamma",
        shape = 2)
    ev <- .fake_events(sim_time = times, sim_status = rep(1L, length(times)))
    out <- add_censoring(ev, censoring = list(dist = "gamma", shape = 2,
        rate = rate), end = max(times) + 1, seed = 2)
    observed <- mean(out$sim_status == 0L)
    expect_equal(observed, target, tolerance = 0.06)
})

test_that("censored times are distributionally consistent with the requested exponential dist [slow]", {
    skip_if_not_slow()
    set.seed(3)
    n <- 4000
    rate <- 0.04
    # Event times far beyond any realistic censoring horizon so almost
    # every subject's observed time IS the censoring draw -- isolates
    # the censoring-time distribution itself for a clean KS check.
    ev <- .fake_events(sim_time = rep(1e6, n), sim_status = rep(1L, n))
    out <- add_censoring(ev, censoring = list(dist = "exponential",
        rate = rate), end = 1e7, seed = 4)
    expect_true(all(out$sim_status == 0L))
    ks <- suppressWarnings(stats::ks.test(out$sim_time, "pexp", rate = rate))
    expect_gt(ks$p.value, 0.01)
})
