# Shared fixtures for the sim_tte_ode() test files (Phase 1/2). testthat
# sources every tests/testthat/helper-*.R file before any test-*.R file,
# in every run mode (devtools::test(), R CMD check, testthat::test_dir()),
# so this is the correct place for cross-file constants -- not a
# top-level binding in one test file relied on by another.

# Default $PARAM overrides sufficient to run each library model, used by
# every test parameterized "for every library model" across
# test-sim-tte-ode-{exponential,weibull,gompertz}.R, rather than
# hardcoding per-model literals at each call site (Phase 2,
# reports/03_implementation_plan.md: "parameterise the existing tests
# over the model names rather than copying them").
.PHASE2_TEST_DEFAULT_PARAM <- list(
    exponential = list(H0 = 0.1),
    weibull = list(mu = -1, shape = 1.5),
    gompertz = list(mu = -1, gamma = 0.1)
)

# Boundary-guard regression scenario, generalized (Phase 2) over
# `model`/`param` rather than copied per model: a covariate-update row
# 0.05 time units before `end`, updating `lp` for every subject --
# forces the solver to restart its step-size search close to the
# administrative horizon (design report section 2.7), the mechanism
# PHASE_E_THRESHOLD_TRACKING_REPORT.md's own M-spline-structure test
# used to force violations. Originally introduced (Phase 1) hardcoded
# to the exponential model only.
.run_ode_boundary_guard_check <- function(model, param, n, seed, end = 10) {
    set.seed(seed)
    data <- data.frame(ID = seq_len(n), time = end - 0.05, lp = 0.5,
        evid = 1, amt = 0, cmt = 1)
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 1, data = data, seed = seed)
    sim$events
}

# 4-SE binomial tolerance for empirical-vs-analytical event-probability
# comparisons; same methodology as test-weibull-distribution.R /
# inst/validation/01_weibull_validation.R, shared across every
# analytical-agreement test in the sim_tte_ode() test files. Widens
# automatically for a smaller `n` (see skip_if_not_slow()'s fast
# counterparts below), so a fast test needs no separate tolerance
# constant -- only a smaller `n` passed through consistently to both the
# simulation and this function.
binom_tol <- function(p, n, z = 4) {
    z * sqrt(p * (1 - p) / n)
}

# ---------------------------------------------------------------------
# Fast/slow test split (interlude session, reports/07_test_runbook.md).
#
# The 2000-subject boundary-guard runs, the 3000-subject analytical/
# cross-method checks, and the R1 Weibull shape sweep are what make a
# targeted `sim-tte-ode` test run slow (measured: the four heaviest
# Weibull analytical/cross-method tests alone total ~20s of a ~37s
# sim-tte-ode-only run). Each such test is gated behind
# skip_if_not_slow() and paired with a smaller-n/coarser-grid
# counterpart (weaker tolerance via the same binom_tol() formula, same
# assertion shape) that always runs, so `devtools::test()`'s default
# invocation still exercises every code path -- only the *scale* differs.
#
# `R CMD check`/CI runs the slow set (SIMTTE_SLOW_TESTS=true is set by
# `dev/run-tests.R --check` and the Makefile's `check` target), so
# nothing is lost at check time; only interactive/targeted runs default
# to fast.
skip_if_not_slow <- function() {
    if (!identical(Sys.getenv("SIMTTE_SLOW_TESTS"), "true")) {
        testthat::skip("slow test (set SIMTTE_SLOW_TESTS=true to run; see reports/07_test_runbook.md)")
    }
}

# ---------------------------------------------------------------------
# Phase 2.5 (reports/04_author_decisions.md "After the test runbook /
# Phase 2.5"): grid-free in-solver refinement. Shared, parameterized
# over `model`/`param` like the fixtures above, rather than copied into
# each of test-sim-tte-ode-{weibull,gompertz}.R.
# ---------------------------------------------------------------------

# Every refined event time must lie within its own subject's
# [T_PRE, TEVT] bracket (the interval .refine_ode_event_time_insolver()
# interpolates within) -- a direct structural check, independent of
# analytical accuracy. Cheap (small n); not one of the three gated
# categories (boundary guard / analytical-cross-method / R1 sweep), so
# always runs.
check_ode_bracket_containment <- function(model, param, n = 300,
    end = 20, seed = 1) {
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 2, keep_trajectory = TRUE, seed = seed)
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(sim$events$ID, last$ID), ]
    is_event <- sim$events$sim_status == 1L
    testthat::expect_true(all(sim$events$sim_time[is_event] >=
        last$T_PRE[is_event] - 1e-8))
    testthat::expect_true(all(sim$events$sim_time[is_event] <=
        last$TEVT[is_event] + 1e-8))
}

# The grid-free bracket does not depend on the reported output grid, so
# a fixed-seed run's refined event times should barely move between a
# coarse and a fine `delta` -- unlike the pre-Phase-2.5 reported-grid
# method (reports/06_phase2_report.md section 4). `tol` is set per
# caller from the corresponding accuracy scale measured in
# inst/validation/10_ode_grid_free_refinement.R, not an arbitrary
# constant.
#
# Event/censoring status itself is compared too, but only for a small
# allowed mismatch rate, not required to match exactly: the raw TEVT
# (which decides event-vs-censored) is itself a solver-quantized time
# and, for a subject whose true event time sits within one internal
# solver step of `end`, can fall on either side of the `end` boundary
# depending on `delta` (mrgsolve/lsoda does not take an internal step
# past the next requested output time) -- a real, pre-existing property
# of the boundary rule, not a refinement defect, and not what this
# check is targeting.
check_ode_delta_independence <- function(model, param, n = 300, end = 20,
    seed = 1, delta_coarse = 4, delta_fine = 0.25, tol,
    max_mismatch_frac = 0.05) {
    sim_coarse <- sim_tte_ode(model = model, param = param, n = n,
        end = end, delta = delta_coarse, seed = seed)
    sim_fine <- sim_tte_ode(model = model, param = param, n = n,
        end = end, delta = delta_fine, seed = seed)
    ev <- sim_coarse$events
    ev_fine <- sim_fine$events[match(ev$ID, sim_fine$events$ID), ]
    n_events <- sum(ev$sim_status == 1L)
    mismatch <- sum(ev$sim_status != ev_fine$sim_status)
    testthat::expect_lt(mismatch, max(1, max_mismatch_frac * n_events))
    both_event <- ev$sim_status == 1L & ev_fine$sim_status == 1L
    diff <- abs(ev$sim_time[both_event] - ev_fine$sim_time[both_event])
    testthat::expect_lt(max(diff), tol)
}
