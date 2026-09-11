# Phase 2: sim_tte_ode() Weibull library model (weibull_ode.cpp). See
# reports/03_implementation_plan.md Phase 2, reports/02_technical_design.md
# section 2, and reports/06_phase2_report.md (risk R1: Weibull-ODE
# shape-range stability) for the design this implements.
#
# Fast/slow split (interlude session, reports/07_test_runbook.md): the
# heaviest tests in this file (analytical/cross-method agreement, the
# full R1 shape sweep, the 2000-subject boundary-guard check) are gated
# behind skip_if_not_slow() and paired with a smaller-n/fewer-shapes
# counterpart that always runs -- same assertion shape, same `delta`
# (grid fineness, not just `n`, matters for accuracy here; see the note
# in section 1), automatically widened tolerance via binom_tol()'s own
# `n` dependence.

# ---------------------------------------------------------------------
# 1. Analytical agreement (4-SE binomial tolerance), several shapes
#    including one < 1. A reasonably fine `delta` is used deliberately,
#    at both scales: the *refinement* step (event-time-refinement,
#    ?sim_tte_ode) assumes a constant hazard within each *reported*
#    interval, which is only exact at shape = 1 -- for shape != 1 the
#    approximation error shrinks as the reported grid is refined
#    (confirmed directly in this session, reports/06_phase2_report.md
#    section "R1"): a coarse `delta` here would conflate that
#    already-documented, expected refinement approximation with a
#    genuine model defect. Only `n` (not `delta`) differs between the
#    fast and slow variants, for exactly this reason.
# ---------------------------------------------------------------------
check_weibull_distribution <- function(shape, mu = -1, n = 3000, end = 20,
    delta = 0.25, seed = 20260910) {
    sim <- sim_tte_ode(model = "weibull", param = list(mu = mu, shape = shape),
        n = n, end = end, delta = delta, seed = seed)
    eta <- exp(mu)
    analytic_S <- function(t) exp(-eta * t^shape)

    p_cens_analytic <- analytic_S(end)
    p_cens_empirical <- mean(sim$events$sim_status == 0)
    expect_lt(abs(p_cens_empirical - p_cens_analytic),
        binom_tol(p_cens_analytic, n), label = paste0("shape=", shape))

    # Checkpoints chosen as analytical quantiles (event probability
    # exactly 0.2/0.4/0.6/0.8 by construction), not fixed absolute
    # times: a fixed-time checkpoint can saturate to (empirical and
    # analytic both) exactly 0 or 1 for a steep hazard (e.g. shape = 2
    # at t = 15), giving a zero-width binomial tolerance band that any
    # floating-point-scale noise trivially "fails" -- not a real
    # disagreement. Quantile-based checkpoints keep every comparison at
    # an informative, non-degenerate probability regardless of shape.
    for (p_target in c(0.2, 0.4, 0.6, 0.8)) {
        t_j <- (-log(1 - p_target) / eta)^(1 / shape)
        p_event_empirical <- mean(sim$events$sim_time <= t_j &
            sim$events$sim_status == 1)
        expect_lt(abs(p_event_empirical - p_target),
            binom_tol(p_target, n),
            label = paste0("shape=", shape, " t_j=", round(t_j, 3)))
    }
}

test_that("sim_tte_ode() Weibull matches the closed-form S(t) [fast]: shape = 1 (exponential-equivalent)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_weibull_distribution(shape = 1, n = 500)
})
test_that("sim_tte_ode() Weibull matches the closed-form S(t) [slow]: shape = 1 (exponential-equivalent)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_weibull_distribution(shape = 1, n = 3000)
})

test_that("sim_tte_ode() Weibull matches the closed-form S(t) [fast]: shape > 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_weibull_distribution(shape = 2, n = 500)
})
test_that("sim_tte_ode() Weibull matches the closed-form S(t) [slow]: shape > 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_weibull_distribution(shape = 2, n = 3000)
})

test_that("sim_tte_ode() Weibull matches the closed-form S(t) [fast]: shape < 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_weibull_distribution(shape = 0.5, n = 500)
})
test_that("sim_tte_ode() Weibull matches the closed-form S(t) [slow]: shape < 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_weibull_distribution(shape = 0.5, n = 3000)
})

# ---------------------------------------------------------------------
# 2. Risk R1: shape-range pass/fail (no crash, no NaN, no non-monotone
#    p11). The full 10-shape sweep is a *permanent* check (slow); the
#    fast counterpart covers the two extremes and shape = 1 -- enough to
#    catch a regression that reintroduces the pre-fix crash, not a
#    substitute for the full sweep. See reports/06_phase2_report.md "R1"
#    for the experiment that determined this range and the T_FLOOR fix.
# ---------------------------------------------------------------------
check_weibull_shape_stability <- function(shapes) {
    for (shape in shapes) {
        # suppressMessages(): a steep shape (e.g. 10) at this coarse
        # delta can trigger the Phase 2.5 refinement-fallback message
        # (see section 5 below) -- expected and irrelevant to what this
        # loop checks (p11 finiteness/monotonicity at reported rows).
        sim <- suppressMessages(sim_tte_ode(model = "weibull",
            param = list(mu = -1, shape = shape), n = 200, end = 20,
            delta = 2, keep_trajectory = TRUE, seed = 1))
        traj <- sim$trajectory
        expect_true(all(is.finite(traj$p11)), label = paste0("shape=", shape))
        by_id <- split(traj$p11, traj$ID)
        nonmono <- vapply(by_id, function(p) any(diff(p) > 1e-8), logical(1))
        expect_true(all(!nonmono), label = paste0("shape=", shape))
    }
}

test_that("weibull_ode.cpp runs without error/NaN/non-monotone p11 [fast: shape 0.05, 1, 10]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_weibull_shape_stability(c(0.05, 1, 10))
})

test_that("weibull_ode.cpp runs without error/NaN/non-monotone p11 [slow: full shape 0.05-10 sweep]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_weibull_shape_stability(c(0.05, 0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 5, 10))
})

test_that("weibull_ode.cpp documented shape range has no hard lower/upper bound (no error/warning for any shape > 0)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # Per ?sim_tte_ode's "Weibull shape support" section: this model
    # always runs (no error, no warning) for any shape > 0; there is no
    # "unsupported range" that refuses or warns -- accuracy degrades
    # gracefully instead (documented, not enforced at runtime). This
    # test asserts exactly that documented policy, including at a shape
    # well below the evidenced-accurate range (0.01), so a future change
    # that silently turns this into an error/warning is caught. Fast
    # already (n = 20); not gated.
    # (mrgsolve's own routine "Building .../Loading model from cache."
    # message is expected on every call and is not itself evidence of a
    # problem -- suppressed here so it isn't conflated with a genuine
    # warning, which is what this test actually checks for.)
    warned <- FALSE
    withCallingHandlers(
        suppressMessages(
            sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 0.01),
                n = 20, end = 10, delta = 2, seed = 1)),
        warning = function(w) {
            warned <<- TRUE
            invokeRestart("muffleWarning")
        })
    expect_false(warned)
})

# ---------------------------------------------------------------------
# 3. Cross-method agreement vs sim_tte(type = "weibull") on identical
#    parameters/seed: quantile comparison (median, IQR, 90th), not just
#    a single-number tolerance, since the two methods' error concentrates
#    differently (design report section 6 item 2).
# ---------------------------------------------------------------------
check_weibull_cross_method <- function(n, quantile_tol) {
    mu <- -1
    shape <- 1.8
    end <- 15
    seed <- 5050

    set.seed(seed)
    lp <- matrix(rep(0, n), nrow = n)
    ref_grid <- sim_tte(pi = lp, mu = mu, coefs = shape,
        time = seq(0.05, end, by = 0.05), type = "weibull", end_time = end,
        event_time_method = "grid")
    set.seed(seed)
    ref_log <- sim_tte(pi = lp, mu = mu, coefs = shape,
        time = seq(0.05, end, by = 0.05), type = "weibull", end_time = end,
        event_time_method = "log_survival")

    sim <- sim_tte_ode(model = "weibull", param = list(mu = mu, shape = shape),
        n = n, end = end, delta = 0.25, seed = seed)

    ode_events <- sim$events$sim_time[sim$events$sim_status == 1]
    grid_events <- ref_grid$sim_time[ref_grid$sim_status == 1]
    log_events <- ref_log$sim_time[ref_log$sim_status == 1]

    probs <- c(0.25, 0.5, 0.75, 0.9)
    q_ode <- quantile(ode_events, probs)
    q_grid <- quantile(grid_events, probs)
    q_log <- quantile(log_events, probs)

    eta <- exp(mu)
    q_analytic <- (-log(1 - probs) / eta)^(1 / shape)

    expect_equal(unname(q_ode), q_analytic, tolerance = quantile_tol)
    expect_equal(unname(q_grid), q_analytic, tolerance = quantile_tol)
    expect_equal(unname(q_log), q_analytic, tolerance = quantile_tol * 0.8)
}

test_that("sim_tte_ode() Weibull quantiles agree with sim_tte(type = 'weibull') [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_weibull_cross_method(n = 500, quantile_tol = 0.3)
})

test_that("sim_tte_ode() Weibull quantiles agree with sim_tte(type = 'weibull') [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_weibull_cross_method(n = 3000, quantile_tol = 0.15)
})

# ---------------------------------------------------------------------
# 4. Boundary-guard regression (shared helper, helper-ode-models.R).
#    40 subjects is already fast; 2000 is gated.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() Weibull never reports sim_time > end [fast: 40 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    events <- .run_ode_boundary_guard_check("weibull",
        list(mu = -1, shape = 1.5), n = 40, seed = 2041)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})

test_that("sim_tte_ode() Weibull never reports sim_time > end [slow: 2000 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    events <- .run_ode_boundary_guard_check("weibull",
        list(mu = -1, shape = 1.5), n = 2000, seed = 4000)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})

# ---------------------------------------------------------------------
# 5. Phase 2.5 (reports/04_author_decisions.md "After the test runbook /
#    Phase 2.5"): grid-free in-solver refinement and the Weibull
#    shape < 0.05 guardrail message. check_ode_bracket_containment()/
#    check_ode_delta_independence() live in helper-ode-models.R,
#    parameterized over model/param like the other shared checks. All
#    cheap (n <= 300); not gated.
# ---------------------------------------------------------------------
test_that("refined sim_time stays within its subject's [T_PRE, TEVT] in-solver bracket", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_bracket_containment("weibull", list(mu = -1, shape = 1.8))
    check_ode_bracket_containment("weibull", list(mu = -1, shape = 0.3))
})

test_that("refined sim_time is (near-)independent of delta: shape > 1 (tol from inst/validation/10_ode_grid_free_refinement.R)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_delta_independence("weibull", list(mu = -1, shape = 2),
        seed = 2026, tol = 0.02)
})

test_that("refined sim_time is (near-)independent of delta: shape < 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_delta_independence("weibull", list(mu = -1, shape = 0.5),
        seed = 2026, tol = 0.01)
})

test_that("sim_tte_ode() messages once when Weibull shape < 0.05, and not at shape = 0.05", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # Below threshold: exactly one message, and the returned events are
    # identical to the same call with the message suppressed (the
    # guardrail is informational only, never changes what runs).
    expect_message(
        below <- sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 0.03),
            n = 20, end = 10, delta = 2, seed = 1),
        "shape' = 0.03 is below 0.05")
    below_quiet <- suppressMessages(sim_tte_ode(model = "weibull",
        param = list(mu = -1, shape = 0.03), n = 20, end = 10, delta = 2,
        seed = 1))
    expect_identical(below$events, below_quiet$events)

    # At and above threshold: no such message. All messages (including
    # mrgsolve's own routine "Building.../Loading model from cache.")
    # are collected and muffled here, then filtered by pattern --
    # suppressMessages() would swallow the one message this test needs
    # to inspect, so it is not used (unlike the "no warning" test above,
    # which only needs to detect *warnings*, unaffected by
    # suppressMessages()).
    expect_no_shape_message <- function(shape) {
        msgs <- character(0)
        withCallingHandlers(
            sim_tte_ode(model = "weibull", param = list(mu = -1, shape = shape),
                n = 20, end = 10, delta = 2, seed = 1),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            })
        expect_false(any(grepl("below 0.05", msgs)), label = paste0("shape=", shape))
    }
    expect_no_shape_message(0.05)
    expect_no_shape_message(1.5)
})
