# Phase 3: sim_tte_ode() M-spline library model (mspline_ode_k{3,5,7}.cpp).
# See reports/02_technical_design.md section 4 and reports/09_phase3_report.md
# for the design/convention this implements. Fast/slow split (as in the
# other sim-tte-ode-*.R files): only `n` differs unless noted.
#
# Shared default fixture: .MSPLINE_TEST_DEFAULT_ARGS (helper-ode-models.R),
# the K = 3 variant, knots = c(5, 10, 15), coefs = rep(1, 6), mu = -1,
# boundary_knots defaulting to c(0, end).

.mspline_analytic_S <- function(t, mu, lp, knots, boundary_knots, coefs) {
    eta <- exp(mu + lp)
    I <- splines2::iSpline(t, knots = knots, Boundary.knots = boundary_knots,
        degree = simtte:::.MSPLINE_DEGREE, intercept = TRUE)
    exp(-eta * as.numeric(I %*% coefs))
}

# ---------------------------------------------------------------------
# 1. Convention equivalence: the in-model C++ basis (captured via HAZ)
#    must match splines2::mSpline(..., degree = 2, intercept = TRUE)
#    to ~1e-8, for every shipped knot-count variant -- the specific
#    claim reports/09_phase3_report.md section 1 makes. Deterministic
#    (no seed dependence); cheap; not gated.
# ---------------------------------------------------------------------
check_mspline_convention <- function(K, knots, boundary_knots = c(0, 20)) {
    skip_if_not_installed("splines2")
    skip_if_not_installed("mrgsolve")
    M <- K + simtte:::.MSPLINE_DEGREE + 1L
    coefs <- seq(0.2, 1.5, length.out = M)   # varied, not all-1s, so a
                                              # basis-vs-coefficient mixup
                                              # would not go unnoticed
    mu <- -0.5
    # HAZ is captured every reported row regardless of event status, so
    # no idata/U override is needed to keep the trajectory dense.
    sim <- sim_tte_ode(model = "mspline", knots = knots,
        boundary_knots = boundary_knots, coefs = coefs, param = list(mu = mu),
        n = 1, end = boundary_knots[2], delta = 1.3, keep_trajectory = TRUE)
    traj <- sim$trajectory
    eta <- exp(mu)
    ref_basis <- splines2::mSpline(traj$time, knots = knots,
        Boundary.knots = boundary_knots, degree = simtte:::.MSPLINE_DEGREE,
        intercept = TRUE)
    ref_haz <- eta * as.numeric(ref_basis %*% coefs)
    expect_equal(traj$HAZ, ref_haz, tolerance = 1e-8, label = paste0("K=", K))
}

test_that("mspline_ode_k3.cpp basis matches splines2::mSpline() (K = 3)", {
    skip_on_cran()
    check_mspline_convention(3, c(5, 10, 15))
})
test_that("mspline_ode_k5.cpp basis matches splines2::mSpline() (K = 5)", {
    skip_on_cran()
    check_mspline_convention(5, c(3.33, 6.67, 10, 13.33, 16.67))
})
test_that("mspline_ode_k7.cpp basis matches splines2::mSpline() (K = 7)", {
    skip_on_cran()
    check_mspline_convention(7, c(2.5, 5, 7.5, 10, 12.5, 15, 17.5))
})

# ---------------------------------------------------------------------
# 2. Analytical agreement (4-SE binomial tolerance), quantile checkpoints
#    against the closed-form I-spline cumulative hazard (H(t) = eta *
#    sum(c * I(t)), S(t) = exp(-H(t))).
# ---------------------------------------------------------------------
check_mspline_distribution <- function(n) {
    skip_if_not_installed("splines2")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    end <- 20
    sim <- sim_tte_ode(model = "mspline", knots = args$knots,
        coefs = args$coefs, param = args$param, n = n, end = end,
        delta = 1, seed = 20260910)

    analytic_S <- function(t) .mspline_analytic_S(t, mu = args$param$mu,
        lp = 0, knots = args$knots, boundary_knots = c(0, end),
        coefs = args$coefs)

    p_cens_analytic <- analytic_S(end)
    p_cens_empirical <- mean(sim$events$sim_status == 0)
    expect_lt(abs(p_cens_empirical - p_cens_analytic), binom_tol(p_cens_analytic, n))

    for (t_j in c(2, 6, 10, 14, 18)) {
        p_event_analytic <- 1 - analytic_S(t_j)
        p_event_empirical <- mean(sim$events$sim_time <= t_j &
            sim$events$sim_status == 1)
        expect_lt(abs(p_event_empirical - p_event_analytic),
            binom_tol(p_event_analytic, n), label = paste0("t_j=", t_j))
    }
}

test_that("sim_tte_ode() M-spline matches the closed-form I-spline S(t) [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_mspline_distribution(n = 500)
})
test_that("sim_tte_ode() M-spline matches the closed-form I-spline S(t) [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_mspline_distribution(n = 3000)
})

# ---------------------------------------------------------------------
# 3. Bracket containment and delta-independence (shared helpers,
#    helper-ode-models.R, generalized with `...` this session to forward
#    knots/coefs/boundary_knots).
# ---------------------------------------------------------------------
test_that("refined sim_time stays within its subject's [T_PRE, TEVT] in-solver bracket", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    check_ode_bracket_containment("mspline", args$param,
        knots = args$knots, coefs = args$coefs)
})

test_that("refined sim_time is (near-)independent of delta (tol from inst/validation/11_ode_mspline_validation.R)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    check_ode_delta_independence("mspline", args$param, seed = 2026,
        tol = 0.01, knots = args$knots, coefs = args$coefs)
})

# ---------------------------------------------------------------------
# 4. Risk-R1-style stability sweep: a range of coefficient magnitudes,
#    including a sharply peaked one (reports/09_phase3_report.md's
#    steepness stress test) -- reported p11 stays finite/monotone
#    regardless (the dense-output interpolant is well-behaved even when
#    an internal evaluation is momentarily pathological, exactly as
#    found for Weibull shape >= 5 in Phase 2.5); the refinement fallback
#    message may fire and is suppressed here, since this test is only
#    about p11, not about which refinement method produced sim_time.
# ---------------------------------------------------------------------
test_that("mspline_ode_k7.cpp runs without error/NaN/non-monotone p11 across a range of coefficient magnitudes, including a sharply peaked one", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    knots <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5)
    for (peak in c(1, 5, 50, 200)) {
        coefs <- rep(0.01, 10)
        coefs[6] <- peak
        sim <- suppressMessages(sim_tte_ode(model = "mspline", knots = knots,
            coefs = coefs, param = list(mu = -1), n = 100, end = 20,
            delta = 1, keep_trajectory = TRUE, seed = 1))
        traj <- sim$trajectory
        expect_true(all(is.finite(traj$p11)), label = paste0("peak=", peak))
        by_id <- split(traj$p11, traj$ID)
        nonmono <- vapply(by_id, function(p) any(diff(p) > 1e-8), logical(1))
        expect_true(all(!nonmono), label = paste0("peak=", peak))
    }
})

# ---------------------------------------------------------------------
# 5. Boundary-guard regression (shared helper, helper-ode-models.R,
#    generalized with `...` this session).
# ---------------------------------------------------------------------
test_that("sim_tte_ode() M-spline never reports sim_time > end [fast: 40 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    events <- .run_ode_boundary_guard_check("mspline", args$param, n = 40,
        seed = 7041, end = 18, knots = args$knots, coefs = args$coefs,
        boundary_knots = c(0, 20))
    expect_true(all(events$sim_time <= 18 + 1e-9))
})
test_that("sim_tte_ode() M-spline never reports sim_time > end [slow: 2000 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    events <- .run_ode_boundary_guard_check("mspline", args$param, n = 2000,
        seed = 7042, end = 18, knots = args$knots, coefs = args$coefs,
        boundary_knots = c(0, 20))
    expect_true(all(events$sim_time <= 18 + 1e-9))
})

# ---------------------------------------------------------------------
# 6. Input validation: unsupported knot count, negative coefficients,
#    wrong coefficient length, end exceeding boundary_knots[2],
#    knots/coefs/boundary_knots supplied for a non-mspline model,
#    knots not strictly inside the boundary interval.
# ---------------------------------------------------------------------
test_that("model = 'mspline' requires both 'knots' and 'coefs'", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "mspline", n = 5, end = 10),
        "requires both 'knots' and 'coefs'")
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10, 15),
        n = 5, end = 10), "requires both 'knots' and 'coefs'")
})

test_that("an unsupported interior-knot count is rejected, naming the supported counts", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10),
        coefs = rep(1, 5), n = 5, end = 10),
        "3, 5, 7")
})

test_that("negative coefficients are rejected", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10, 15),
        coefs = c(-1, 1, 1, 1, 1, 1), param = list(mu = -1), n = 5,
        end = 20), "non-negative")
})

test_that("a wrong-length coefficient vector is rejected", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10, 15),
        coefs = rep(1, 5), param = list(mu = -1), n = 5, end = 20),
        "length\\(knots\\)")
})

test_that("'end' exceeding boundary_knots[2] is rejected, not silently extrapolated", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10, 15),
        coefs = rep(1, 6), boundary_knots = c(0, 18), param = list(mu = -1),
        n = 5, end = 20), "exceeds boundary_knots")
})

test_that("knots not lying strictly inside boundary_knots are rejected", {
    skip_on_cran()
    # Last knot equal to bk_hi (default boundary_knots = c(0, end)) --
    # must be strictly less, not merely <=.
    expect_error(sim_tte_ode(model = "mspline", knots = c(5, 10, 20),
        coefs = rep(1, 6), param = list(mu = -1), n = 5, end = 20),
        "strictly inside")
})

test_that("'knots'/'boundary_knots'/'coefs' are rejected for non-mspline models", {
    skip_on_cran()
    expect_error(sim_tte_ode(model = "weibull", knots = c(5, 10, 15),
        coefs = rep(1, 6), param = list(mu = -1, shape = 1), n = 5,
        end = 10), "only used when model")
})

# ---------------------------------------------------------------------
# 7. Time-varying covariates: a constant (single time = 0 row)
#    population-level covariate must be exactly equivalent to folding
#    the same value into a fixed 'lp' offset -- the cleanest possible
#    direct check that Phase 2's covariates/beta mechanism, unmodified,
#    works correctly for this new model (reports/09_phase3_report.md
#    section "Time-varying covariates").
# ---------------------------------------------------------------------
test_that("time-varying covariates combine correctly for the M-spline model", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    cov <- data.frame(time = 0, x = 1)
    s1 <- sim_tte_ode(model = "mspline", knots = args$knots,
        coefs = args$coefs, param = args$param, covariates = cov,
        beta = c(x = 0.7), n = 150, end = 20, delta = 1, seed = 55)
    s2 <- sim_tte_ode(model = "mspline", knots = args$knots,
        coefs = args$coefs, param = c(args$param, list(lp = 0.7)),
        n = 150, end = 20, delta = 1, seed = 55)
    expect_identical(s1$events, s2$events)
})

# ---------------------------------------------------------------------
# 8. Reproducibility.
# ---------------------------------------------------------------------
test_that("same seed gives identical() $events for the M-spline model", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    args <- .MSPLINE_TEST_DEFAULT_ARGS
    s1 <- sim_tte_ode(model = "mspline", knots = args$knots,
        coefs = args$coefs, param = args$param, n = 100, end = 20,
        delta = 1, seed = 4242)
    s2 <- sim_tte_ode(model = "mspline", knots = args$knots,
        coefs = args$coefs, param = args$param, n = 100, end = 20,
        delta = 1, seed = 4242)
    expect_identical(s1$events, s2$events)
})
