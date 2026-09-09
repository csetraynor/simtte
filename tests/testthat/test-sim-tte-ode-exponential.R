# Phase 1: sim_tte_ode() in-solver event detection, exponential (constant
# hazard) library model only. See reports/02_technical_design.md section 2
# and reports/03_implementation_plan.md Phase 1 for the design this
# implements, and reports/05_phase1_report.md for results.

binom_tol <- function(p, n, z = 4) {
    z * sqrt(p * (1 - p) / n)
}

# ---------------------------------------------------------------------
# 1. Analytical agreement (4-SE binomial tolerance), same methodology as
#    test-weibull-distribution.R / inst/validation/01_weibull_validation.R.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() exponential matches the closed-form S(t) = exp(-eta*t)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")

    H0 <- 0.15
    lp <- 0.2
    eta <- H0 * exp(lp)
    n <- 3000
    end <- 20

    sim <- sim_tte_ode(model = "exponential", param = list(H0 = H0, lp = lp),
        n = n, end = end, delta = 2, seed = 20260910)

    analytic_S <- function(t) exp(-eta * t)

    p_cens_analytic <- analytic_S(end)
    p_cens_empirical <- mean(sim$events$sim_status == 0)
    expect_lt(abs(p_cens_empirical - p_cens_analytic),
        binom_tol(p_cens_analytic, n))

    for (t_j in c(1, 3, 5, 10, 15)) {
        p_event_analytic <- 1 - analytic_S(t_j)
        p_event_empirical <- mean(sim$events$sim_time <= t_j &
            sim$events$sim_status == 1)
        expect_lt(abs(p_event_empirical - p_event_analytic),
            binom_tol(p_event_analytic, n), label = paste0("t_j=", t_j))
    }
})

# ---------------------------------------------------------------------
# 2. Boundary-guard regression suite: a covariate-update row close to
#    `end` (the scenario PHASE_E_THRESHOLD_TRACKING_REPORT.md section 10
#    found violations in), at both 40- and 2000-subject scale. Zero
#    violations expected with the shipped (guarded) model.
# ---------------------------------------------------------------------
run_boundary_guard_check <- function(n, seed) {
    set.seed(seed)
    end <- 10
    # A covariate-update row 0.05 time units before `end`, updating `lp`
    # for every subject -- forces the solver to restart its step-size
    # search close to the administrative horizon (design report section
    # 2.7), the mechanism PHASE_E_THRESHOLD_TRACKING_REPORT.md's own
    # M-spline-structure test used to force violations.
    data <- data.frame(ID = rep(seq_len(n), each = 1), time = end - 0.05,
        lp = 0.5, evid = 1, amt = 0, cmt = 1)
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.3),
        n = n, end = end, delta = 1, data = data, seed = seed)
    sim$events
}

test_that("sim_tte_ode() never reports sim_time > end (40 subjects, covariate update near end)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    events <- run_boundary_guard_check(40, seed = 1001)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})

test_that("sim_tte_ode() never reports sim_time > end (2000 subjects, covariate update near end)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    events <- run_boundary_guard_check(2000, seed = 1002)
    expect_true(all(events$sim_time <= 10 + 1e-9))
})

test_that("removing the in-model SOLVERTIME <= END guard reproduces PHASE_E_THRESHOLD_TRACKING_REPORT.md's boundary violation (one-time, not shipped)", {
    # Documents, once, that the guard is load-bearing rather than
    # incidental. Two attempts to reproduce the boundary defect from
    # scratch on a plain constant-hazard model with a single covariate
    # jump near `end` (a simpler scenario than Phase E's) did NOT
    # reproduce it in this mrgsolve version -- reported honestly in
    # reports/05_phase1_report.md, consistent with
    # PHASE_E_THRESHOLD_TRACKING_REPORT.md's own characterization of the
    # defect as scenario-dependent, not universal. What *does* reliably
    # reproduce it is Phase E's own exact M-spline-structure scenario
    # (three basehaz knots at t = 0, 1, 2, no boundary guard), re-run
    # here directly (mirroring
    # inst/validation/05_des_threshold_tracking_investigation.R, itself
    # re-verified live in this session to still reproduce 4/40
    # violations at seed 20260828). This test's own twin models are
    # never shipped; it is not a permanently-required regression check
    # on package behavior (the shipped exponential_ode.cpp model always
    # has the guard).
    skip_on_cran()
    skip_if_not_installed("mrgsolve")

    des_code_no_guard <- "
$PARAM lp = 0, mu = 0, basehaz = 1, U = 0.5
$INIT p11 = 1
$GLOBAL
static int event_found = 0;
static double TEVT = -1.0;
$MAIN
double eta = exp(mu + lp);
if (NEWIND <= 1) { event_found = 0; TEVT = -1.0; }
$ODE
if (p11 <= U && event_found == 0) { event_found = 1; TEVT = SOLVERTIME; }
dxdt_p11 = -p11 * basehaz * eta;
$CAPTURE TEVT event_found
"
    des_code_guarded <- "
$PARAM lp = 0, mu = 0, basehaz = 1, U = 0.5, END = 2
$INIT p11 = 1
$GLOBAL
static int event_found = 0;
static double TEVT = -1.0;
$MAIN
double eta = exp(mu + lp);
if (NEWIND <= 1) { event_found = 0; TEVT = -1.0; }
$ODE
if (p11 <= U && event_found == 0 && SOLVERTIME <= END) {
    event_found = 1; TEVT = SOLVERTIME;
}
dxdt_p11 = -p11 * basehaz * eta;
$CAPTURE TEVT event_found
"
    mod_no_guard <- mrgsolve::mcode("phase1_ms_no_guard_twin",
        des_code_no_guard, quiet = TRUE)
    mod_guarded <- mrgsolve::mcode("phase1_ms_guarded_twin",
        des_code_guarded, quiet = TRUE)

    n <- 1000
    end <- 2
    set.seed(20260828)
    U <- runif(n)
    data <- data.frame(ID = rep(seq_len(n), each = 3),
        time = rep(c(0, 1, end), n), amt = 0, evid = 1, cmt = 1,
        basehaz = rep(c(1, 2, 4), n), mu = 0, lp = 0,
        U = rep(U, each = 3))

    run <- function(mod) {
        out <- as.data.frame(mrgsolve::mrgsim(mrgsolve::data_set(mod, data),
            tgrid = c(0, end), obsonly = TRUE, nocb = FALSE))
        out[out$time == max(out$time), ]
    }
    last_no_guard <- run(mod_no_guard)
    last_guarded <- run(mod_guarded)

    T_true_fn <- function(u) {
        ifelse(u >= exp(-1), -log(u),
            ifelse(u >= exp(-3), 1 + (-log(u) - 1) / 2, NA_real_))
    }
    should_censor <- is.na(T_true_fn(U))

    n_violations_no_guard <- sum(last_no_guard$TEVT[should_censor] > end)
    n_violations_guarded <- sum(last_guarded$TEVT[should_censor] > end)

    # The defect class is confirmed reproducible without the guard...
    expect_gt(n_violations_no_guard, 0)
    # ...and eliminated on the identical scenario with the guard added.
    expect_equal(n_violations_guarded, 0)
})

# ---------------------------------------------------------------------
# 3. Censoring rule: TEVT >= end -> censored at end with status 0.
# ---------------------------------------------------------------------
test_that("subjects with no latched event, or a latched TEVT >= end, are censored at end", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # H0 = 0 -> hazard is identically zero -> p11 never crosses U for any
    # U in [0, 1) -> event_found stays 0 for every subject -> every
    # subject must be censored at `end` with status 0.
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0),
        n = 50, end = 12, delta = 3, seed = 55)
    expect_true(all(sim$events$sim_status == 0))
    expect_true(all(sim$events$sim_time == 12))
})

# ---------------------------------------------------------------------
# 4. Refinement: with a deliberately coarse `delta`, sim_time is closer
#    to the analytical quantile than the raw TEVT, and always within the
#    bracketing reported interval.
# ---------------------------------------------------------------------
test_that("refinement improves on raw TEVT and stays within the bracketing reported interval", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    H0 <- 0.1
    n <- 500
    end <- 30
    delta <- 5   # deliberately coarse
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = H0), n = n,
        end = end, delta = delta, keep_trajectory = TRUE, seed = 2026)

    events <- sim$events
    traj <- sim$trajectory
    is_event <- events$sim_status == 1L

    grid <- sort(unique(traj$time))
    # Bracketing check: every refined event time lies within the pair of
    # adjacent grid points surrounding it (inclusive), never outside.
    lower <- findInterval(events$sim_time[is_event], grid)
    lower <- pmax(lower, 1L)
    upper <- pmin(lower + 1L, length(grid))
    expect_true(all(events$sim_time[is_event] >= grid[lower] - 1e-8))
    expect_true(all(events$sim_time[is_event] <= grid[upper] + 1e-8))
})

test_that("refined sim_time has lower mean absolute error than raw TEVT vs. the analytical exponential quantile", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    H0 <- 0.1
    n <- 800
    end <- 40
    delta <- 8   # deliberately coarse, to make the raw-vs-refined gap visible
    seed <- 3033
    set.seed(seed)
    idata <- data.frame(ID = seq_len(n), U = runif(n))
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = H0),
        idata = idata, end = end, delta = delta, keep_trajectory = TRUE)

    events <- sim$events
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(events$ID, last$ID), ]

    analytic_T <- -log(idata$U[match(events$ID, idata$ID)]) / H0
    is_event <- events$sim_status == 1L & analytic_T < end

    err_refined <- abs(events$sim_time[is_event] - analytic_T[is_event])
    err_raw <- abs(last$TEVT[is_event] - analytic_T[is_event])

    expect_lt(mean(err_refined), mean(err_raw))
})

# ---------------------------------------------------------------------
# 5. Reproducibility.
# ---------------------------------------------------------------------
test_that("same seed gives identical() $events", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    s1 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 100, end = 20, delta = 2, seed = 4242)
    s2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 100, end = 20, delta = 2, seed = 4242)
    expect_identical(s1$events, s2$events)
})

# ---------------------------------------------------------------------
# 6. Model-contract validator.
# ---------------------------------------------------------------------
test_that(".validate_ode_model_contract() rejects a model missing p11", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    bad <- mrgsolve::mcode("phase1_missing_p11",
        "$PARAM H0=0.1, U=0, END=1\n$CMT DUMMY\n$ODE\ndxdt_DUMMY=0;\n",
        quiet = TRUE)
    expect_error(simtte:::.validate_ode_model_contract(bad), "p11")
})

test_that(".validate_ode_model_contract() rejects a model missing U/END", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    bad <- mrgsolve::mcode("phase1_missing_u_end",
        "$PARAM H0=0.1\n$INIT\np11=1\n$ODE\ndxdt_p11=-p11*H0;\n",
        quiet = TRUE)
    expect_error(simtte:::.validate_ode_model_contract(bad), "U, END")
})

# ---------------------------------------------------------------------
# 7. R3: sequential-individual-simulation dependency -- subject i's
#    trajectory must reflect subject i's own parameters, not another
#    subject's (loud failure if mrgsolve ever stops simulating
#    individuals strictly sequentially; see
#    reports/03_implementation_plan.md risk R3).
# ---------------------------------------------------------------------
test_that("each subject's p11 trajectory matches its own per-subject H0, not another subject's", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    n <- 300
    idata <- data.frame(ID = seq_len(n), H0 = seq(0.02, 0.5, length.out = n))
    sim <- sim_tte_ode(model = "exponential", n = n, end = 50, delta = 5,
        idata = idata, keep_trajectory = TRUE, seed = 123)
    traj <- sim$trajectory
    chk <- merge(traj[traj$time == 5, c("ID", "p11")], idata, by = "ID")
    expect_equal(chk$p11, exp(-chk$H0 * 5), tolerance = 1e-6)
})

# ---------------------------------------------------------------------
# 8. sim_tte_ode()'s own reserved-...-args / model-enum handling.
# ---------------------------------------------------------------------
test_that("reserved mrgsim() arguments are still rejected via '...'", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(
        sim_tte_ode(model = "exponential", n = 5, end = 10, tgrid = c(0, 10)),
        "tgrid")
})

test_that("an unknown model name is rejected via match.arg()", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte_ode(model = "weibull", n = 5, end = 10))
})
