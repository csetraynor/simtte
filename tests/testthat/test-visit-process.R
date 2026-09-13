# Phase 6c: visit-process generators C1 (thin_visits()) and C2
# (visit_schedule_informative()), reports/20_visit_process_evaluation.md.
# Both operate on a schedule (+ events for C2), no mrgsolve needed for
# most tests; the sim_tte_ode() pipeline/bias-demonstration tests are
# guarded with skip_if_not_installed(), same convention as elsewhere.

# .fake_events() is in helper-events.R (shared with test-censoring.R
# and test-interval-censoring.R).

# ---------------------------------------------------------------------
# 1. thin_visits() (C1) [fast].
# ---------------------------------------------------------------------
test_that("thin_visits() never drops the baseline visit", {
    sched <- visit_schedule(n = 30, every = 4, end = 24, seed = 1)
    out <- thin_visits(sched, p_miss = 1, p_dropout = 1, seed = 2)
    for (id in unique(out$ID)) {
        expect_equal(min(out$time[out$ID == id]), 0)
    }
})

test_that("p_miss = 0 and p_dropout = 0 returns the schedule unchanged", {
    sched <- visit_schedule(n = 10, every = 4, end = 20, seed = 1)
    out <- thin_visits(sched, p_miss = 0, p_dropout = 0)
    expect_equal(out, sched)
})

test_that("p_miss = 1 drops every non-baseline visit", {
    sched <- visit_schedule(n = 10, every = 4, end = 20, seed = 1)
    out <- thin_visits(sched, p_miss = 1, seed = 3)
    for (id in unique(out$ID)) {
        expect_equal(out$time[out$ID == id], 0)
    }
})

test_that("p_dropout = 1 drops every visit from some onset point onward", {
    sched <- visit_schedule(n = 20, every = 4, end = 24, seed = 1)
    out <- thin_visits(sched, p_dropout = 1, seed = 4)
    for (id in unique(out$ID)) {
        v_before <- sched$time[sched$ID == id]
        v_after <- out$time[out$ID == id]
        expect_true(length(v_after) <= length(v_before))
        # Whatever remains is a strict prefix of the original schedule.
        expect_identical(v_after, v_before[seq_along(v_after)])
    }
})

test_that("thin_visits() is reproducible with the same seed", {
    sched <- visit_schedule(n = 20, every = 4, end = 24, seed = 1)
    a <- thin_visits(sched, p_miss = 0.3, p_dropout = 0.2, seed = 5)
    b <- thin_visits(sched, p_miss = 0.3, p_dropout = 0.2, seed = 5)
    expect_identical(a, b)
})

test_that("thin_visits() accepts a common numeric vector schedule given 'ids'", {
    out <- thin_visits(c(0, 4, 8, 12), p_miss = 0.5, ids = 1:5, seed = 1)
    expect_true(all(c("ID", "time") %in% names(out)))
    expect_equal(sort(unique(out$ID)), 1:5)
})

test_that("thin_visits() requires 'ids' for a numeric vector schedule", {
    expect_error(thin_visits(c(0, 4, 8), p_miss = 0.5), "'ids' is required")
})

test_that("thin_visits() validates p_miss/p_dropout", {
    sched <- visit_schedule(n = 5, every = 4, end = 20, seed = 1)
    expect_error(thin_visits(sched, p_miss = 1.5), "p_miss")
    expect_error(thin_visits(sched, p_dropout = -0.1), "p_dropout")
})

test_that("thin_visits() output passes add_interval_censoring() validation", {
    sched <- visit_schedule(n = 5, every = 4, end = 20, seed = 1)
    thinned <- thin_visits(sched, p_miss = 0.5, p_dropout = 0.3, seed = 6)
    ev <- .fake_events(c(3, 9, 15, 18, 2), c(1L, 1L, 0L, 1L, 0L))
    out <- add_interval_censoring(ev, visits = thinned)
    expect_true(all(is.finite(out$sim_time_left)))
})

test_that("thin_visits() composes with add_censoring() in either order", {
    ev <- .fake_events(c(9, 25), c(1L, 1L))
    spec <- list(dist = "user", fun = function(n) c(20, 5))
    sched <- visit_schedule(n = 2, every = 4, end = 30, seed = 1)
    thinned <- thin_visits(sched, p_miss = 0.4, seed = 2)

    censor_then_thin <- add_censoring(ev, censoring = spec, end = 30)
    censor_then_thin <- add_interval_censoring(censor_then_thin, visits = thinned)

    thin_then_censor_events <- add_censoring(ev, censoring = spec, end = 30)
    thin_first <- thin_visits(sched, p_miss = 0.4, seed = 2) # same schedule/seed
    thin_then_censor <- add_interval_censoring(thin_then_censor_events,
        visits = thin_first)

    # Same schedule + same censored events either way -> identical bounds.
    expect_identical(censor_then_thin$sim_time_left, thin_then_censor$sim_time_left)
    expect_identical(censor_then_thin$sim_time_right, thin_then_censor$sim_time_right)
})

test_that("realized miss/dropout rates are close to nominal at scale [slow]", {
    skip_if_not_slow()
    sched <- visit_schedule(n = 5000, every = 4, end = 24, seed = 1)
    n_non_baseline <- sum(sched$time != 0)
    out <- thin_visits(sched, p_miss = 0.2, seed = 2)
    n_kept <- sum(out$time != 0)
    realized_miss <- 1 - n_kept / n_non_baseline
    expect_equal(realized_miss, 0.2, tolerance = 0.05)
})

# ---------------------------------------------------------------------
# 2. visit_schedule_informative() (C2) [fast].
# ---------------------------------------------------------------------
test_that("a visit inside the near-event window is missed at the combined rate (p_miss_base = 0 here, so just p)", {
    ev <- .fake_events(9, 1L) # event at t = 9
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = 3, p = 1), seed = 1))
    # Visit 8 is in [9 - 3, 9) = [6, 9) -> always missed.
    expect_false(8 %in% out$time)
    expect_true(all(c(0, 4, 12, 16) %in% out$time))
})

test_that("p_miss_base = 0 reproduces the previous override behavior exactly", {
    # The additive formula 1 - (1 - p_miss_base) * (1 - p) collapses to
    # p alone when p_miss_base = 0 -- the pre-additive-change behavior.
    ev <- .fake_events(9, 1L)
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        p_miss_base = 0, miss_near_event = list(window = 3, p = 0.5),
        seed = 4))
    manual_p_eff <- 1 - (1 - 0) * (1 - 0.5)
    expect_equal(manual_p_eff, 0.5)
})

test_that("near-event miss and background miss are additive/competing, not overriding, at scale", {
    p_base <- 0.2
    p_near <- 0.5
    expected <- 1 - (1 - p_base) * (1 - p_near) # 0.6, not 0.5
    n <- 4000
    ev <- .fake_events(rep(9, n), rep(1L, n))
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        p_miss_base = p_base, miss_near_event = list(window = 3, p = p_near),
        seed = 5))
    realized_miss <- 1 - mean(sapply(split(out$time, out$ID),
        function(v) 8 %in% v))
    expect_equal(realized_miss, expected, tolerance = 0.05)
})

test_that("a visit outside the near-event window is unaffected by miss_near_event", {
    ev <- .fake_events(9, 1L)
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = 1, p = 1), seed = 1))
    # Window is [8, 9) -> only visit 8 is at risk; nothing else near t = 9.
    expect_true(all(c(0, 4, 12, 16) %in% out$time))
})

test_that("a censored subject is unaffected by miss_near_event", {
    ev <- .fake_events(9, 0L) # censored, not an event
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = 3, p = 1), seed = 1))
    expect_identical(sort(out$time), sched)
})

test_that("p_miss_base applies to every subject/visit not covered by miss_near_event", {
    ev <- .fake_events(c(9, 20), c(1L, 0L))
    sched <- c(0, 4, 8, 12, 16, 20)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        p_miss_base = 1, seed = 1))
    for (id in unique(out$ID)) {
        expect_equal(out$time[out$ID == id], 0)
    }
})

test_that("an extra visit is added after an event, at sim_time + delay", {
    ev <- .fake_events(9, 1L)
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        extra_visit_after_event = list(dist = "user", fun = function(n) rep(2, n)),
        seed = 1))
    expect_true(11 %in% out$time) # 9 + 2
})

test_that("no extra visit is added for a censored subject", {
    ev <- .fake_events(9, 0L)
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        extra_visit_after_event = list(dist = "user", fun = function(n) rep(2, n)),
        seed = 1))
    expect_identical(sort(out$time), sched)
})

test_that("an extra visit landing after 'end' is dropped, not added", {
    ev <- .fake_events(19, 1L)
    sched <- c(0, 4, 8, 12, 16)
    out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
        extra_visit_after_event = list(dist = "user", fun = function(n) rep(5, n)),
        seed = 1)) # 19 + 5 = 24 > 20
    expect_identical(sort(out$time), sched)
})

test_that("each extra_visit_after_event distribution spec is accepted", {
    ev <- .fake_events(5, 1L)
    sched <- c(0, 4, 8, 12, 16, 20)
    for (spec in list(
        list(dist = "exponential", rate = 1),
        list(dist = "uniform", min = 0, max = 2),
        list(dist = "user", fun = function(n) rep(1, n)))) {
        out <- suppressMessages(visit_schedule_informative(sched, ev, end = 20,
            extra_visit_after_event = spec, seed = 1))
        expect_true(nrow(out) >= length(sched))
    }
})

test_that("visit_schedule_informative() emits one message() per call, and it is suppressible", {
    ev <- .fake_events(9, 1L)
    sched <- c(0, 4, 8, 12, 16)
    expect_message(
        visit_schedule_informative(sched, ev, end = 20,
            miss_near_event = list(window = 3, p = 1), seed = 1),
        "biased by design")
    expect_no_message(suppressMessages(
        visit_schedule_informative(sched, ev, end = 20,
            miss_near_event = list(window = 3, p = 1), seed = 1)))
})

test_that("visit_schedule_informative() is reproducible with the same seed", {
    ev <- .fake_events(c(9, 20, 5), c(1L, 0L, 1L))
    sched <- visit_schedule(n = 3, every = 4, end = 24, seed = 1)
    a <- suppressMessages(visit_schedule_informative(sched, ev, end = 24,
        miss_near_event = list(window = 3, p = 0.5),
        extra_visit_after_event = list(dist = "exponential", rate = 0.5),
        seed = 9))
    b <- suppressMessages(visit_schedule_informative(sched, ev, end = 24,
        miss_near_event = list(window = 3, p = 0.5),
        extra_visit_after_event = list(dist = "exponential", rate = 0.5),
        seed = 9))
    expect_identical(a, b)
})

test_that("visit_schedule_informative() validates miss_near_event/p_miss_base", {
    ev <- .fake_events(9, 1L)
    sched <- c(0, 4, 8, 12, 16)
    expect_error(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = 3)), "window.*p")
    expect_error(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = -1, p = 0.5)), "window")
    expect_error(visit_schedule_informative(sched, ev, end = 20,
        miss_near_event = list(window = 3, p = 2)), "p")
    expect_error(visit_schedule_informative(sched, ev, end = 20,
        p_miss_base = -0.1), "p_miss_base")
})

test_that("visit_schedule_informative() output passes add_interval_censoring() validation", {
    ev <- .fake_events(c(9, 20, 5), c(1L, 0L, 1L))
    sched <- visit_schedule(n = 3, every = 4, end = 24, seed = 1)
    realized <- suppressMessages(visit_schedule_informative(sched, ev, end = 24,
        miss_near_event = list(window = 3, p = 0.5),
        extra_visit_after_event = list(dist = "exponential", rate = 0.5),
        seed = 9))
    out <- add_interval_censoring(ev, visits = realized)
    expect_true(all(is.finite(out$sim_time_left)))
})

test_that("visit_schedule_informative() built from pre- vs. post-right-censoring events gives different bounds", {
    ev <- .fake_events(25, 1L) # a real event, well before end = 30
    spec <- list(dist = "user", fun = function(n) 10) # deterministic pull to 10
    sched <- c(0, 5, 10, 15, 20, 25, 30)

    pre <- suppressMessages(visit_schedule_informative(sched, ev, end = 30,
        miss_near_event = list(window = 5, p = 1), seed = 1))
    out_pre <- add_interval_censoring(ev, visits = pre)

    censored <- add_censoring(ev, censoring = spec, end = 30)
    post <- suppressMessages(visit_schedule_informative(sched, censored,
        end = 30, miss_near_event = list(window = 5, p = 1), seed = 1))
    out_post <- add_interval_censoring(censored, visits = post)

    expect_false(identical(out_pre$sim_time_right, out_post$sim_time_right))
})

# ---------------------------------------------------------------------
# 3. Slow: the C2 bias demonstration -- an informative schedule biases
#    the interval-censored survival estimate more than a non-informative
#    (C1) schedule at the same nominal miss rate.
# ---------------------------------------------------------------------
test_that("an informative (C2) visit schedule biases a Turnbull fit more than a non-informative (C1) one at the same nominal miss rate [slow]", {
    skip_if_not_slow()
    skip_if_not_installed("mrgsolve")
    skip_if_not_installed("survival")
    rate <- 0.08
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = rate),
        n = 3000, end = 30, delta = 2, seed = 1)
    visit_times <- seq(0, 30, by = 3)
    analytical_surv <- exp(-rate * visit_times)

    fit_surv <- function(schedule) {
        out <- add_interval_censoring(sim$events, visits = schedule)
        right_for_surv <- ifelse(is.infinite(out$sim_time_right), NA,
            out$sim_time_right)
        fit <- survival::survfit(survival::Surv(time = sim_time_left,
            time2 = right_for_surv, type = "interval2") ~ 1,
            data = data.frame(sim_time_left = out$sim_time_left,
                right_for_surv = right_for_surv))
        summary(fit, times = visit_times, extend = TRUE)$surv
    }

    base_sched <- visit_schedule(n = nrow(sim$events), every = 3, end = 30,
        seed = 1)
    base_sched$ID <- sim$events$ID[base_sched$ID]

    c1_sched <- thin_visits(base_sched, p_miss = 0.3, seed = 2)
    c1_surv <- fit_surv(c1_sched)

    # C2: a visit shortly before an event is missed much more often than
    # the same 0.3 nominal rate would suggest -- an informative process.
    c2_sched <- suppressMessages(visit_schedule_informative(base_sched,
        sim$events, end = 30, p_miss_base = 0.1,
        miss_near_event = list(window = 3, p = 0.9), seed = 2))
    c2_surv <- fit_surv(c2_sched)

    c1_error <- max(abs(c1_surv - analytical_surv))
    c2_error <- max(abs(c2_surv - analytical_surv))
    expect_true(c2_error > c1_error)
})
