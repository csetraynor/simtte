# Phase 2: sim_tte_ode() Gompertz library model (gompertz_ode.cpp). See
# reports/03_implementation_plan.md Phase 2 and
# reports/02_technical_design.md section 2 for the design this
# implements. Unlike Weibull, the Gompertz hazard h(t) = eta*exp(gamma*t)
# has no t -> 0+ singularity for any gamma (including gamma < 0, a
# decreasing hazard) -- there is no risk-R1-style shape-range concern
# here, confirmed directly below.
#
# Fast/slow split (interlude session, reports/07_test_runbook.md): the
# heaviest tests here (analytical agreement, the gamma = 0 vs.
# exponential comparison, the 10-checkpoint empirical-survival check,
# the 2000-subject boundary guard) are gated behind skip_if_not_slow()
# with a smaller-n counterpart that always runs. The gamma sweep (unlike
# Weibull's R1 shape sweep) is already fast (~0.4s) and is not one of
# this session's three named slow categories, so it is left ungated.

analytic_gompertz_S <- function(t, mu, gamma) {
    eta <- exp(mu)
    exp(-(eta / gamma) * (exp(gamma * t) - 1))
}

# ---------------------------------------------------------------------
# 1. Analytical agreement (4-SE binomial tolerance): an increasing
#    hazard (gamma > 0) and a decreasing hazard (gamma < 0).
# ---------------------------------------------------------------------
check_gompertz_distribution <- function(gamma, mu = -2, n = 3000, end = 20,
    delta = 2, seed = 20260910) {
    sim <- sim_tte_ode(model = "gompertz", param = list(mu = mu, gamma = gamma),
        n = n, end = end, delta = delta, seed = seed)
    analytic_S <- function(t) analytic_gompertz_S(t, mu, gamma)

    p_cens_analytic <- analytic_S(end)
    p_cens_empirical <- mean(sim$events$sim_status == 0)
    expect_lt(abs(p_cens_empirical - p_cens_analytic),
        binom_tol(p_cens_analytic, n), label = paste0("gamma=", gamma))

    for (t_j in c(3, 8, 13, 18)) {
        p_event_analytic <- 1 - analytic_S(t_j)
        p_event_empirical <- mean(sim$events$sim_time <= t_j &
            sim$events$sim_status == 1)
        expect_lt(abs(p_event_empirical - p_event_analytic),
            binom_tol(p_event_analytic, n),
            label = paste0("gamma=", gamma, " t_j=", t_j))
    }
}

test_that("sim_tte_ode() Gompertz matches the closed-form S(t) [fast]: increasing hazard (gamma > 0)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_gompertz_distribution(gamma = 0.15, n = 500)
})
test_that("sim_tte_ode() Gompertz matches the closed-form S(t) [slow]: increasing hazard (gamma > 0)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_gompertz_distribution(gamma = 0.15, n = 3000)
})

test_that("sim_tte_ode() Gompertz matches the closed-form S(t) [fast]: decreasing hazard (gamma < 0)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_gompertz_distribution(gamma = -0.1, mu = -1, n = 500)
})
test_that("sim_tte_ode() Gompertz matches the closed-form S(t) [slow]: decreasing hazard (gamma < 0)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_gompertz_distribution(gamma = -0.1, mu = -1, n = 3000)
})

# ---------------------------------------------------------------------
# 2. No t -> 0+ singularity: unlike Weibull, no shape-range sweep is
#    needed, but confirm directly (no error/NaN/non-monotone p11) across
#    a comparably wide gamma range, including gamma = 0 (reduces to
#    exponential) as an edge case. Already fast; not gated.
# ---------------------------------------------------------------------
test_that("gompertz_ode.cpp runs without error/NaN/non-monotone p11 across a wide gamma range, including gamma = 0", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    gammas <- c(-0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, 1)
    for (gamma in gammas) {
        sim <- sim_tte_ode(model = "gompertz",
            param = list(mu = -1, gamma = gamma), n = 200, end = 20,
            delta = 2, keep_trajectory = TRUE, seed = 1)
        traj <- sim$trajectory
        expect_true(all(is.finite(traj$p11)), label = paste0("gamma=", gamma))
        by_id <- split(traj$p11, traj$ID)
        nonmono <- vapply(by_id, function(p) any(diff(p) > 1e-8), logical(1))
        expect_true(all(!nonmono), label = paste0("gamma=", gamma))
    }
})

check_gompertz_gamma0_vs_exponential <- function(n) {
    mu <- -1
    end <- 20
    sim_g <- sim_tte_ode(model = "gompertz", param = list(mu = mu, gamma = 0),
        n = n, end = end, delta = 2, seed = 77)
    sim_e <- sim_tte_ode(model = "exponential", param = list(H0 = exp(mu)),
        n = n, end = end, delta = 2, seed = 77)
    expect_identical(sim_g$events$sim_status, sim_e$events$sim_status)
    expect_equal(sim_g$events$sim_time, sim_e$events$sim_time,
        tolerance = 1e-6)
}

test_that("gompertz_ode.cpp at gamma = 0 matches the exponential model [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_gompertz_gamma0_vs_exponential(n = 400)
})
test_that("gompertz_ode.cpp at gamma = 0 matches the exponential model [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_gompertz_gamma0_vs_exponential(n = 2000)
})

# ---------------------------------------------------------------------
# 3. Cross-method agreement: no built-in sim_tte(type = "gompertz")
#    exists, so this checks internal consistency instead -- the
#    empirical Kaplan-Meier-style event probability from sim_tte_ode()
#    vs. the closed-form Gompertz S(t), at more checkpoints than the
#    basic analytical-agreement test above, as the closest available
#    analogue to the Weibull cross-method check.
# ---------------------------------------------------------------------
check_gompertz_checkpoints <- function(n) {
    mu <- -2
    gamma <- 0.12
    end <- 20
    sim <- sim_tte_ode(model = "gompertz", param = list(mu = mu, gamma = gamma),
        n = n, end = end, delta = 1, seed = 909)
    checkpoints <- seq(2, 18, by = 2)
    for (t_j in checkpoints) {
        p_a <- 1 - analytic_gompertz_S(t_j, mu, gamma)
        p_e <- mean(sim$events$sim_time <= t_j & sim$events$sim_status == 1)
        expect_lt(abs(p_a - p_e), binom_tol(p_a, n), label = paste0("t=", t_j))
    }
}

test_that("sim_tte_ode() Gompertz empirical survival matches closed-form S(t) at 10 checkpoints [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_gompertz_checkpoints(n = 600)
})
test_that("sim_tte_ode() Gompertz empirical survival matches closed-form S(t) at 10 checkpoints [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_gompertz_checkpoints(n = 4000)
})

# ---------------------------------------------------------------------
# 4. Boundary-guard regression (shared helper, helper-ode-models.R).
#    40 subjects is already fast; 2000 is gated.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() Gompertz never reports sim_time > end [fast: 40 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    events <- .run_ode_boundary_guard_check("gompertz",
        list(mu = -1, gamma = 0.1), n = 40, seed = 3041)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})
test_that("sim_tte_ode() Gompertz never reports sim_time > end [slow: 2000 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    events <- .run_ode_boundary_guard_check("gompertz",
        list(mu = -1, gamma = 0.1), n = 2000, seed = 5000)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})

# ---------------------------------------------------------------------
# 5. Phase 2.5 (reports/04_author_decisions.md "After the test runbook /
#    Phase 2.5"): grid-free in-solver refinement. See
#    test-sim-tte-ode-weibull.R section 5 for the shared helpers
#    (helper-ode-models.R) and rationale; not repeated here. Cheap
#    (n <= 300); not gated.
# ---------------------------------------------------------------------
test_that("refined sim_time stays within its subject's [T_PRE, TEVT] in-solver bracket", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_bracket_containment("gompertz", list(mu = -1, gamma = 0.1))
    check_ode_bracket_containment("gompertz", list(mu = -1, gamma = -0.05))
})

test_that("refined sim_time is (near-)independent of delta: increasing hazard (tol from inst/validation/10_ode_grid_free_refinement.R)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_delta_independence("gompertz", list(mu = -1, gamma = 0.1),
        seed = 2026, tol = 0.001)
})

test_that("refined sim_time is (near-)independent of delta: decreasing hazard", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_ode_delta_independence("gompertz", list(mu = -1, gamma = -0.05),
        seed = 2026, tol = 0.001)
})
