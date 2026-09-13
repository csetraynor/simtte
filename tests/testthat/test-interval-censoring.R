# Phase 6b: interval censoring (reports/18_interval_censoring_design.md,
# option A). add_interval_censoring()/visit_schedule() work on any
# events-shaped data frame -- most tests build a small fabricated
# events data frame directly, no mrgsolve compile needed, keeping this
# file fast. The sim_tte_ode(visits = ...) integration tests need
# mrgsolve and are guarded with skip_if_not_installed(), same
# convention as every other ODE test file.

# .fake_events() is in helper-events.R (shared with test-censoring.R
# and test-visit-process.R).

# ---------------------------------------------------------------------
# 1. Interval conventions, one row each [fast].
# ---------------------------------------------------------------------
test_that("event mid-schedule gives (L, R]", {
    ev <- .fake_events(9, 1L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12, 16))
    expect_equal(out$sim_time_left, 8)
    expect_equal(out$sim_time_right, 12)
})

test_that("event before the first post-baseline visit gives (0, R_1]", {
    ev <- .fake_events(2, 1L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    expect_equal(out$sim_time_left, 0)
    expect_equal(out$sim_time_right, 4)
})

test_that("event exactly at a visit time gives (L_prev, R]", {
    ev <- .fake_events(8, 1L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    expect_equal(out$sim_time_left, 4)
    expect_equal(out$sim_time_right, 8)
})

test_that("event-free at the last visit (administratively censored) gives (L_last, Inf)", {
    ev <- .fake_events(24, 0L) # sim_time == end, the sim_tte_ode() convention
    out <- add_interval_censoring(ev, visits = c(0, 8, 16, 24))
    expect_equal(out$sim_time_left, 24)
    expect_true(is.infinite(out$sim_time_right))
})

test_that("right-censored between visits gives (L, Inf), L = last visit before censoring", {
    ev <- .fake_events(10, 0L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12, 16))
    expect_equal(out$sim_time_left, 8)
    expect_true(is.infinite(out$sim_time_right))
})

test_that("censored exactly at a visit time gives (that visit, Inf)", {
    ev <- .fake_events(8, 0L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    expect_equal(out$sim_time_left, 8)
    expect_true(is.infinite(out$sim_time_right))
})

test_that("an event that escapes the visit schedule (T > max(visits)) gives (max(visits), Inf)", {
    ev <- .fake_events(15, 1L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    expect_equal(out$sim_time_left, 12)
    expect_true(is.infinite(out$sim_time_right))
    # sim_status/sim_time are the exact ground truth, untouched.
    expect_equal(out$sim_status, 1L)
    expect_equal(out$sim_time, 15)
})

test_that("a t = 0 event gives the degenerate interval (0, 0]", {
    ev <- .fake_events(0, 1L)
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    expect_equal(out$sim_time_left, 0)
    expect_equal(out$sim_time_right, 0)
})

test_that("sim_time/sim_status/sim_reason are never modified by interval censoring", {
    ev <- .fake_events(c(3, 24), c(1L, 0L))
    ev$sim_reason <- c("event", "administrative")
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12, 16, 20, 24))
    expect_identical(out$sim_time, ev$sim_time)
    expect_identical(out$sim_status, ev$sim_status)
    expect_identical(out$sim_reason, ev$sim_reason)
})

# ---------------------------------------------------------------------
# 2. Common vs. per-subject schedules [fast].
# ---------------------------------------------------------------------
test_that("a common vector schedule applies identically to every subject", {
    ev <- .fake_events(c(2, 9, 15), c(1L, 1L, 1L))
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12, 16))
    expect_equal(out$sim_time_left, c(0, 8, 12))
    expect_equal(out$sim_time_right, c(4, 12, 16))
})

test_that("a per-subject data frame schedule gives each subject its own bounds", {
    ev <- .fake_events(c(2, 9), c(1L, 1L))
    visits <- data.frame(ID = c(1, 1, 1, 2, 2, 2),
        time = c(0, 4, 8, 0, 10, 20))
    out <- add_interval_censoring(ev, visits = visits)
    expect_equal(out$sim_time_left, c(0, 0))
    expect_equal(out$sim_time_right, c(4, 10))
})

# ---------------------------------------------------------------------
# 3. Schedule validation errors [fast].
# ---------------------------------------------------------------------
test_that("a vector schedule not starting at 0 is an error", {
    ev <- .fake_events(2, 1L)
    expect_error(add_interval_censoring(ev, visits = c(4, 8, 12)),
        "start at 0")
})

test_that("an unsorted or duplicated vector schedule is an error", {
    ev <- .fake_events(2, 1L)
    expect_error(add_interval_censoring(ev, visits = c(0, 8, 4, 12)),
        "strictly increasing")
    expect_error(add_interval_censoring(ev, visits = c(0, 4, 4, 12)),
        "strictly increasing")
})

test_that("a negative visit time is an error", {
    ev <- .fake_events(2, 1L)
    expect_error(add_interval_censoring(ev, visits = c(-1, 0, 4)),
        "non-negative")
})

test_that("a per-subject schedule missing an events ID, or with an unknown ID, errors", {
    ev <- .fake_events(c(2, 9), c(1L, 1L))
    expect_error(add_interval_censoring(ev,
        visits = data.frame(ID = 1, time = c(0, 4))),
        "missing a schedule.*2")
    expect_error(add_interval_censoring(ev,
        visits = data.frame(ID = c(1, 2, 3), time = c(0, 0, 0))),
        "not present in 'events'.*3")
})

test_that("a malformed 'visits' argument (not numeric, not a data frame) errors", {
    ev <- .fake_events(2, 1L)
    expect_error(add_interval_censoring(ev, visits = "not a schedule"),
        "numeric vector.*data frame")
    expect_error(add_interval_censoring(ev, visits = list(every = 4)),
        "numeric vector.*data frame")
})

# ---------------------------------------------------------------------
# 4. visit_schedule(): seeded, bounded jitter, never crosses order
#    [fast]. jitter/every bound: uniform needs jitter < every / 2;
#    normal (jitter_trunc default 2) needs 2 * jitter < every / 2.
#    every = 4 below -> bound is 2 for uniform, 1 for normal at the
#    default jitter_trunc; 0.9 clears both.
# ---------------------------------------------------------------------
test_that("visit_schedule() is reproducible with the same seed, for each jitter_dist", {
    for (dist in c("uniform", "normal")) {
        a <- visit_schedule(n = 10, every = 4, end = 20, jitter = 0.9,
            jitter_dist = dist, seed = 1)
        b <- visit_schedule(n = 10, every = 4, end = 20, jitter = 0.9,
            jitter_dist = dist, seed = 1)
        expect_identical(a, b)
    }
})

test_that("jitter_dist = 'uniform' is visit_schedule()'s unchanged default", {
    a <- visit_schedule(n = 10, every = 4, end = 20, jitter = 0.9, seed = 1)
    b <- visit_schedule(n = 10, every = 4, end = 20, jitter = 0.9,
        jitter_dist = "uniform", seed = 1)
    expect_identical(a, b)
})

test_that("'uniform' jitter output is byte-identical to the pre-truncation-bound implementation", {
    # Captured from the implementation immediately before this session's
    # jitter/every bound was added (reports/23_phase6c_report.md): since
    # jitter = 1 < every / 2 = 2 here, the old pmax(0, .)/sort(unique(.))
    # calls were already no-ops, so this session's simplification must
    # reproduce the exact same draws for the same seed.
    out <- visit_schedule(n = 5, every = 4, end = 20, jitter = 1,
        jitter_dist = "uniform", seed = 1)
    ref <- c(0, 3.5310173262842, 7.74424779927358, 12.1457067267038,
        16.8164155799896, 19.4033638620749, 0, 4.79677936993539,
        8.8893505372107, 12.3215955849737, 16.2582280877978,
        19.1235725409351, 0, 3.4119491497986, 7.35311350505799,
        12.3740456933156, 15.7682074364275, 20.5396828399971, 0,
        3.99539848417044, 8.43523701652884, 12.983812189661,
        15.7600703588687, 20.5548904426396, 0, 4.86941046221182,
        7.42428504256532, 12.3033475321718, 15.2511101919226,
        19.5344413374551)
    expect_equal(out$time, ref)
})

test_that("visit_schedule() never produces a negative or unsorted time, for each jitter_dist, at the largest allowed jitter", {
    for (dist in c("uniform", "normal")) {
        max_jitter <- if (dist == "uniform") 1.999 else 0.999
        sched <- visit_schedule(n = 50, every = 4, end = 20,
            jitter = max_jitter, jitter_dist = dist, seed = 2)
        expect_true(all(sched$time >= 0))
        for (id in unique(sched$ID)) {
            expect_false(is.unsorted(sched$time[sched$ID == id],
                strictly = TRUE))
        }
    }
})

test_that("a jitter/jitter_trunc combination too large for 'every' errors before any draw, for each jitter_dist", {
    expect_error(
        visit_schedule(n = 5, every = 4, end = 20, jitter = 2,
            jitter_dist = "uniform"),
        "too large relative to 'every'")
    expect_error(
        visit_schedule(n = 5, every = 4, end = 20, jitter = 1,
            jitter_dist = "normal"), # jitter_trunc = 2 default -> max_dev = 2
        "too large relative to 'every'")
    expect_error(
        visit_schedule(n = 5, every = 4, end = 20, jitter = 0.7,
            jitter_dist = "normal", jitter_trunc = 3), # max_dev = 2.1 >= every/2 = 2
        "too large relative to 'every'")
    # A jitter comfortably inside the bound never errors, for either dist.
    expect_no_error(visit_schedule(n = 20, every = 4, end = 20,
        jitter = 0.1, seed = 3))
    expect_no_error(visit_schedule(n = 20, every = 4, end = 20,
        jitter = 0.1, jitter_dist = "normal", seed = 3))
})

test_that("the reordering invariant is unreachable at the maximum allowed jitter, for each jitter_dist, at scale", {
    for (dist in c("uniform", "normal")) {
        max_jitter <- if (dist == "uniform") 1.999 else 0.999
        expect_no_error(expect_no_warning(
            visit_schedule(n = 5000, every = 4, end = 20,
                jitter = max_jitter, jitter_dist = dist, seed = 7)))
    }
})

test_that("visit_schedule() keeps the baseline visit at exactly 0 and stays ordered per subject", {
    sched <- visit_schedule(n = 20, every = 3, end = 15, jitter = 1.4, seed = 3)
    for (id in unique(sched$ID)) {
        v <- sched$time[sched$ID == id]
        expect_equal(v[1], 0)
        expect_false(is.unsorted(v, strictly = TRUE))
    }
})

test_that("visit_schedule() output is directly usable as add_interval_censoring()'s visits", {
    sched <- visit_schedule(n = 3, every = 4, end = 12, jitter = 0, seed = 1)
    ev <- .fake_events(c(1, 5, 9), c(1L, 1L, 1L))
    out <- add_interval_censoring(ev, visits = sched)
    expect_true(all(is.finite(out$sim_time_left)))
})

test_that("visit_schedule() validates its own arguments", {
    expect_error(visit_schedule(n = 0, every = 4, end = 20), "n")
    expect_error(visit_schedule(n = 5, every = -1, end = 20), "every")
    expect_error(visit_schedule(n = 5, every = 4, end = -1), "end")
    expect_error(visit_schedule(n = 5, every = 4, end = 20, jitter = -1),
        "jitter")
    expect_error(visit_schedule(n = 5, every = 4, end = 20,
        jitter_dist = "gaussian"), "should be one of")
    expect_error(visit_schedule(n = 5, every = 4, end = 20,
        jitter_trunc = -1), "jitter_trunc")
})

# ---------------------------------------------------------------------
# 5. survival::Surv(type = "interval2") round-trip [fast].
# ---------------------------------------------------------------------
test_that("sim_time_right = Inf converts to NA and round-trips through Surv(interval2)", {
    skip_if_not_installed("survival")
    ev <- .fake_events(c(9, 10, 0), c(1L, 0L, 1L))
    out <- add_interval_censoring(ev, visits = c(0, 4, 8, 12))
    right_for_surv <- ifelse(is.infinite(out$sim_time_right), NA,
        out$sim_time_right)
    s <- survival::Surv(time = out$sim_time_left, time2 = right_for_surv,
        type = "interval2")
    printed <- trimws(format(s))
    expect_equal(printed[1], "[8, 12]") # genuine interval
    expect_equal(printed[2], "8+")      # open upper bound -> right-censored notation
    expect_equal(printed[3], "0")       # degenerate t = 0 event -> exact
})

# ---------------------------------------------------------------------
# 6. Ordering with add_censoring(): right order vs. wrong order [fast].
# ---------------------------------------------------------------------
test_that("interval mapping before right censoring gives bounds inconsistent with the final outcome", {
    ev <- .fake_events(25, 1L) # a real event, well before end = 30
    spec <- list(dist = "user", fun = function(n) 10) # deterministic pull to 10

    # Correct order: censor first, map second.
    right_order <- add_censoring(ev, censoring = spec, end = 30)
    right_order <- add_interval_censoring(right_order, visits = c(0, 10, 20, 30))
    expect_true(right_order$sim_time_left <= right_order$sim_time)
    expect_true(is.infinite(right_order$sim_time_right))

    # Wrong order: map first (against the pre-censoring truth), censor
    # second -- the interval columns are never touched by add_censoring(),
    # so they go stale relative to the now-different final sim_time.
    wrong_order <- add_interval_censoring(ev, visits = c(0, 10, 20, 30))
    wrong_order <- add_censoring(wrong_order, censoring = spec, end = 30)
    expect_false(wrong_order$sim_time_left <= wrong_order$sim_time &&
        wrong_order$sim_time <= wrong_order$sim_time_right)
})

# ---------------------------------------------------------------------
# 7. Legacy sim_tte() output goes through unchanged in other columns
#    [fast].
# ---------------------------------------------------------------------
test_that("add_interval_censoring() on sim_tte() output leaves other columns untouched", {
    set.seed(1)
    lp <- matrix(rnorm(20, 0, 0.5), nrow = 20)
    legacy <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
        time = seq(0, 20, by = 1), type = "weibull", end_time = 20)
    out <- add_interval_censoring(legacy, visits = c(0, 5, 10, 15, 20))
    expect_identical(out$lp, legacy$lp)
    expect_identical(out$ID, legacy$ID)
    expect_true(all(c("sim_time_left", "sim_time_right") %in% names(out)))
    expect_false("sim_time_left" %in% names(legacy))
})

# ---------------------------------------------------------------------
# 8. sim_tte_ode() convenience argument [fast].
# ---------------------------------------------------------------------
test_that("sim_tte_ode(visits = <schedule>) equals the standalone call", {
    skip_if_not_installed("mrgsolve")
    visits <- c(0, 4, 8, 12, 16, 20, 24)
    via_arg <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 20, end = 24, delta = 2, visits = visits, seed = 1)
    manual <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 20, end = 24, delta = 2, seed = 1)
    manual_events <- add_interval_censoring(manual$events, visits = visits)
    expect_identical(via_arg$events, manual_events)
})

test_that("sim_tte_ode(visits = <schedule>) applies after censoring = , not before", {
    skip_if_not_installed("mrgsolve")
    spec <- list(dist = "exponential", rate = 0.1)
    visits <- c(0, 4, 8, 12, 16, 20, 24)
    via_arg <- sim_tte_ode(model = "exponential", param = list(H0 = 0.02),
        n = 30, end = 24, delta = 2, censoring = spec, visits = visits,
        seed = 5)
    manual <- sim_tte_ode(model = "exponential", param = list(H0 = 0.02),
        n = 30, end = 24, delta = 2, censoring = spec, seed = 5)
    manual_events <- add_interval_censoring(manual$events, visits = visits)
    expect_identical(via_arg$events, manual_events)
})

test_that("sim_tte_ode(visits = list(every=, jitter=)) builds a schedule inside the seeded call", {
    skip_if_not_installed("mrgsolve")
    s1 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2, visits = list(every = 4, jitter = 1),
        seed = 9)
    s2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2, visits = list(every = 4, jitter = 1),
        seed = 9)
    expect_identical(s1$events, s2$events)
    expect_true(all(c("sim_time_left", "sim_time_right") %in% names(s1$events)))
})

test_that("sim_tte_ode(visits = list(..., jitter_dist = 'normal')) is passed through", {
    skip_if_not_installed("mrgsolve")
    # every = 4 -> bound is every / 2 = 2; jitter = 0.4 clears both the
    # uniform bound (0.4 < 2) and the default jitter_trunc = 2 normal
    # bound (2 * 0.4 = 0.8 < 2).
    s1 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2,
        visits = list(every = 4, jitter = 0.4, jitter_dist = "normal"),
        seed = 9)
    s2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2,
        visits = list(every = 4, jitter = 0.4, jitter_dist = "normal"),
        seed = 9)
    expect_identical(s1$events, s2$events)
    # A different jitter_dist gives a different schedule (not a no-op).
    s_uniform <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2,
        visits = list(every = 4, jitter = 0.4, jitter_dist = "uniform"),
        seed = 9)
    expect_false(identical(s1$events, s_uniform$events))
})

test_that("sim_tte_ode(visits = list(..., jitter_trunc = )) is passed through", {
    skip_if_not_installed("mrgsolve")
    # jitter_trunc = 4 widens the normal truncation bound (max_dev = 4 *
    # 0.4 = 1.6 < every / 2 = 2), still valid, and changes the draw.
    s_default <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2,
        visits = list(every = 4, jitter = 0.4, jitter_dist = "normal"),
        seed = 9)
    s_wide <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 15, end = 20, delta = 2,
        visits = list(every = 4, jitter = 0.4, jitter_dist = "normal",
            jitter_trunc = 4),
        seed = 9)
    expect_false(identical(s_default$events, s_wide$events))
})

# ---------------------------------------------------------------------
# 9. Slow: bounds always bracket the truth at scale; a nonparametric
#    interval-censored fit recovers the known survival curve.
# ---------------------------------------------------------------------
test_that("sim_time always lies in (sim_time_left, sim_time_right] at scale [slow]", {
    skip_if_not_slow()
    skip_if_not_installed("mrgsolve")
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.08),
        n = 3000, end = 30, delta = 2, seed = 1)
    out <- add_interval_censoring(sim$events, visits = seq(0, 30, by = 3))
    expect_true(all(out$sim_time <= out$sim_time_right))
    expect_true(all(out$sim_time >= out$sim_time_left))
})

test_that("a Turnbull NPMLE fit on interval-censored data recovers the known exponential survival curve [slow]", {
    skip_if_not_slow()
    skip_if_not_installed("mrgsolve")
    skip_if_not_installed("survival")
    rate <- 0.08
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = rate),
        n = 3000, end = 30, delta = 2, seed = 1)
    visit_times <- seq(0, 30, by = 3)
    out <- add_interval_censoring(sim$events, visits = visit_times)
    right_for_surv <- ifelse(is.infinite(out$sim_time_right), NA,
        out$sim_time_right)
    fit <- survival::survfit(survival::Surv(time = sim_time_left,
        time2 = right_for_surv, type = "interval2") ~ 1,
        data = data.frame(sim_time_left = out$sim_time_left,
            right_for_surv = right_for_surv))
    fitted_surv <- summary(fit, times = visit_times, extend = TRUE)$surv
    analytical_surv <- exp(-rate * visit_times)
    expect_equal(fitted_surv, analytical_surv, tolerance = 0.05)
})
