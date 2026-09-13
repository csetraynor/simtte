# Phase 4: sim_tte_ode() PK/PD-linked hazard library (pk_hazard,
# irm1_hazard-irm4_hazard, tmdd_hazard). See reports/10_phase4_report.md
# for what was built and reports/04_author_decisions.md "After Phase 3"
# decision 2 for the CRAN-run-time policy this file follows: unlike
# every other sim_tte_ode() test file (Phase 1-3, all skip_on_cran()'d
# throughout), the compile-and-smoke test per model below -- plus the
# beta=0 equivalence check, the reserved-args check, and the (cheap,
# pure-R, no compilation) differentiated-fallback-message tests -- run
# on CRAN; everything else here is skip_on_cran() + skip_if_not_slow()
# gated, per decision 2's "Slow set" list.

# ---------------------------------------------------------------------
# 1. Compile-and-smoke test per model [CRAN].
# ---------------------------------------------------------------------
for (.m in .PKPD_MODELS) {
    local({
        model <- .m
        test_that(paste0("sim_tte_ode() model = '", model,
            "' compiles and a 20-subject run returns well-formed $events [CRAN]"), {
            skip_if_not_installed("mrgsolve")
            sim <- sim_tte_ode(model = model,
                param = .PKPD_TEST_DEFAULT_PARAM[[model]], n = 20, end = 20,
                delta = 2, data = .pkpd_dose_data(20), seed = 1)
            expect_s3_class(sim, "simtte_ode_sim")
            expect_named(sim$events,
                c("ID", "sim_time", "sim_status", "sim_reason"))
            expect_equal(nrow(sim$events), 20)
            expect_true(all(sim$events$sim_status %in% c(0L, 1L)))
            expect_true(all(sim$events$sim_time >= 0 &
                sim$events$sim_time <= 20 + 1e-9))
        })
    })
}
rm(.m)

# ---------------------------------------------------------------------
# 2. beta = 0 exponential-equivalence, one model [CRAN].
# ---------------------------------------------------------------------
test_that("sim_tte_ode() pk_hazard with beta_cp = 0 matches the exponential model with hazard H0 [CRAN]", {
    skip_if_not_installed("mrgsolve")
    sim_pk <- sim_tte_ode(model = "pk_hazard",
        param = list(H0 = 0.02, beta_cp = 0), n = 40, end = 20, delta = 1,
        data = .pkpd_dose_data(40), seed = 42)
    sim_exp <- sim_tte_ode(model = "exponential", param = list(H0 = 0.02),
        n = 40, end = 20, delta = 1, seed = 42)
    expect_identical(sim_pk$events$sim_status, sim_exp$events$sim_status)
    expect_equal(sim_pk$events$sim_time, sim_exp$events$sim_time,
        tolerance = 1e-6)
})

# ---------------------------------------------------------------------
# 3. Reserved-args behaviour with a dosing data frame [CRAN].
# ---------------------------------------------------------------------
test_that("sim_tte_ode() accepts a dosing 'data' frame for pk_hazard while reserved dots still error [CRAN]", {
    skip_if_not_installed("mrgsolve")
    dose <- .pkpd_dose_data(10)
    expect_no_error(sim_tte_ode(model = "pk_hazard",
        param = list(H0 = 0.02, beta_cp = 0.05), n = 10, end = 20,
        delta = 2, data = dose, seed = 1))
    expect_error(sim_tte_ode(model = "pk_hazard",
        param = list(H0 = 0.02, beta_cp = 0.05), n = 10, end = 20,
        delta = 2, data = dose, seed = 1, tgrid = seq(0, 20, 5)),
        "controlled internally")
})

# ---------------------------------------------------------------------
# 4. Differentiated fallback message (decision 1): .classify_ode_fallback()
#    directly, plus .resolve_ode_events() end-to-end with synthetic
#    bracket columns for all three signatures at once. Pure R, no model
#    compilation -- cheap; runs on CRAN like section 1-3 above.
# ---------------------------------------------------------------------
test_that(".classify_ode_fallback() identifies the weibull_type signature (P_POST outside [0, 1])", {
    expect_identical(simtte:::.classify_ode_fallback(t_pre = 4, p_pre = 0.9,
        tevt = 5, p_post = -0.05), "weibull_type")
    expect_identical(simtte:::.classify_ode_fallback(t_pre = 4, p_pre = 0.9,
        tevt = 5, p_post = 1.2), "weibull_type")
})

test_that(".classify_ode_fallback() identifies the mspline_type signature (T_PRE == TEVT, legitimate probabilities)", {
    expect_identical(simtte:::.classify_ode_fallback(t_pre = 5, p_pre = 0.6,
        tevt = 5, p_post = 0.55), "mspline_type")
})

test_that(".classify_ode_fallback() falls back to unclassified otherwise", {
    expect_identical(simtte:::.classify_ode_fallback(t_pre = NaN,
        p_pre = 0.6, tevt = 6, p_post = 0.5), "unclassified")
    # A degenerate bracket (t_pre < tevt would not be degenerate; here
    # t_pre > tevt, which .refine_ode_event_time_insolver() also treats
    # as degenerate) with legitimate probabilities but no T_PRE == TEVT:
    expect_identical(simtte:::.classify_ode_fallback(t_pre = 6, p_pre = 0.5,
        tevt = 5, p_post = 0.55), "unclassified")
})

test_that(".resolve_ode_events() emits one differentiated message naming all three fallback signatures", {
    traj <- rbind(
        # Subject 1 (weibull_type): P_POST outside [0, 1].
        data.frame(ID = 1, time = c(0, 10), p11 = c(1, 0.3), U = 0.5,
            END = 10, event_found = 1, TEVT = 5, T_PRE = 4, P_PRE = 0.9,
            P_POST = -0.05),
        # Subject 2 (mspline_type): T_PRE == TEVT, legitimate probabilities.
        data.frame(ID = 2, time = c(0, 10), p11 = c(1, 0.4), U = 0.5,
            END = 10, event_found = 1, TEVT = 5, T_PRE = 5, P_PRE = 0.6,
            P_POST = 0.55),
        # Subject 3 (unclassified): non-finite T_PRE.
        data.frame(ID = 3, time = c(0, 10), p11 = c(1, 0.4), U = 0.5,
            END = 10, event_found = 1, TEVT = 6, T_PRE = NaN, P_PRE = 0.6,
            P_POST = 0.5),
        # Subject 4 (control): a clean, non-degenerate bracket -- must
        # NOT be counted in the fallback message at all.
        data.frame(ID = 4, time = c(0, 10), p11 = c(1, 0.4), U = 0.5,
            END = 10, event_found = 1, TEVT = 4, T_PRE = 3, P_PRE = 0.8,
            P_POST = 0.5)
    )
    # log(negative P_POST) inside .interpolate_log_survival() produces
    # an expected NaN-with-warning for the synthetic weibull_type row
    # (the same mechanism a real unphysical solver overshoot triggers);
    # not the thing under test here.
    msgs <- suppressWarnings(capture_messages(
        out <- simtte:::.resolve_ode_events(traj)))
    expect_length(msgs, 1L)
    expect_match(msgs, "3 subject\\(s\\) used the reported-grid")
    expect_match(msgs, "1 subject\\(s\\) had an unphysical P_POST")
    expect_match(msgs, "1 subject\\(s\\) had .*T_PRE == TEVT")
    expect_match(msgs, "1 subject\\(s\\) fell back for an unclassified reason")
    expect_match(msgs, "finer 'delta'")
    expect_match(msgs, "tighter 'rtol'/'atol' resolves this, not 'delta'")
    # Subject 4 used the in-solver bracket directly, not the fallback.
    expect_false(out$sim_time[out$ID == 4] %in%
        c(out$sim_time[out$ID == 1], out$sim_time[out$ID == 2]))
    expect_identical(names(out),
        c("ID", "sim_time", "sim_status", "sim_reason"))
})

test_that(".resolve_ode_events() emits no message when no subject falls back", {
    traj <- data.frame(ID = 1, time = c(0, 10), p11 = c(1, 0.4), U = 0.5,
        END = 10, event_found = 1, TEVT = 4, T_PRE = 3, P_PRE = 0.8,
        P_POST = 0.5)
    expect_length(capture_messages(simtte:::.resolve_ode_events(traj)), 0L)
})

# ---------------------------------------------------------------------
# 5. Quantile agreement vs. a fine-grid sim_tte_df() reference [slow].
#    No closed form exists for a PK/PD-linked hazard in general
#    (reports/10_phase4_report.md section on validation), so this
#    mirrors inst/validation/12_ode_pkpd_validation.R's own reference
#    method on a smaller scale, checking internal consistency between
#    sim_tte_ode()'s in-solver mechanism and the package's independent
#    grid-based mechanism rather than a closed form.
# ---------------------------------------------------------------------
check_pkpd_vs_sim_tte_df <- function(model, param, n, seed = 909, end = 30,
    dose_amt = 100) {
    dose <- .pkpd_dose_data(n, amt = dose_amt)
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 1, data = dose, seed = seed)

    mod <- simtte:::.read_ode_library_model(model)
    mod <- mrgsolve::param(mod, param)
    idata <- data.frame(ID = seq_len(n))
    fine_grid <- seq(0, end, by = 0.001)
    traj <- as.data.frame(mrgsolve::mrgsim(mod, idata = idata, data = dose,
        tgrid = fine_grid, obsonly = TRUE))
    # sim_tte_df() treats each subject's own last reported time as their
    # censoring horizon (?sim_tte_df); every subject shares the same
    # fine_grid here, so that horizon is 'end' for all of them, the same
    # administrative censoring sim_tte_ode() applies.
    set.seed(seed)
    ref <- sim_tte_df(traj, event_time_method = "log_survival")

    ev <- sim$events[order(sim$events$ID), ]
    rf <- ref[order(ref$ID), ]
    n_events <- sum(ev$sim_status == 1L)
    both_event <- ev$sim_status == 1L & rf$sim_status == 1L
    expect_gt(sum(both_event), 0.5 * n_events)
    diff <- abs(ev$sim_time[both_event] - rf$sim_time[both_event])
    # Loose tolerance: sim_tte_df() draws its own independent U (a
    # different RNG stream than sim_tte_ode()'s idata$U), so this is a
    # distributional/quantile check, not a per-subject one.
    expect_lt(stats::median(diff), 1)
}

for (.m in .PKPD_MODELS) {
    local({
        model <- .m
        test_that(paste0("sim_tte_ode() '", model,
            "' broadly agrees with an independent fine-grid reference [slow]"), {
            skip_on_cran()
            skip_if_not_installed("mrgsolve")
            skip_if_not_slow()
            check_pkpd_vs_sim_tte_df(model, .PKPD_TEST_DEFAULT_PARAM[[model]],
                n = 300)
        })
    })
}
rm(.m)

# ---------------------------------------------------------------------
# 6. Boundary guard [slow].
# ---------------------------------------------------------------------
for (.m in .PKPD_MODELS) {
    local({
        model <- .m
        test_that(paste0("sim_tte_ode() '", model,
            "' never reports sim_time > end [slow: 2000 subjects]"), {
            skip_on_cran()
            skip_if_not_installed("mrgsolve")
            skip_if_not_slow()
            events <- .run_ode_boundary_guard_check(model,
                .PKPD_TEST_DEFAULT_PARAM[[model]], n = 2000, seed = 5001,
                end = 10)
            expect_true(all(events$sim_time <= 10 + 1e-9))
        })
    })
}
rm(.m)
test_that("sim_tte_ode() PK/PD models never report sim_time > end [fast: 40 subjects]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    for (model in .PKPD_MODELS) {
        events <- .run_ode_boundary_guard_check(model,
            .PKPD_TEST_DEFAULT_PARAM[[model]], n = 40, seed = 3042, end = 10)
        expect_true(all(events$sim_time <= 10 + 1e-9),
            label = paste0("model=", model))
    }
})

# ---------------------------------------------------------------------
# 7. Bracket containment [slow].
# ---------------------------------------------------------------------
test_that("refined sim_time stays within its subject's [T_PRE, TEVT] in-solver bracket, every PK/PD model [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    for (model in .PKPD_MODELS) {
        check_ode_bracket_containment(model, .PKPD_TEST_DEFAULT_PARAM[[model]],
            n = 300, end = 20, data = .pkpd_dose_data(300))
    }
})

# ---------------------------------------------------------------------
# 8. Mechanistic direction checks [slow]: increasing dose/beta shifts
#    event times the way the link says it should.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() pk_hazard: a higher dose (beta_cp > 0) shortens event times [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    n <- 800
    sim_lo <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.05), n = n, end = 30, delta = 1,
        data = .pkpd_dose_data(n, amt = 20), seed = 1)
    sim_hi <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.05), n = n, end = 30, delta = 1,
        data = .pkpd_dose_data(n, amt = 200), seed = 1)
    p_event_lo <- mean(sim_lo$events$sim_status == 1)
    p_event_hi <- mean(sim_hi$events$sim_status == 1)
    expect_gt(p_event_hi, p_event_lo)
})

test_that("sim_tte_ode() irm3_hazard: sign of beta_r sets the hazard direction [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    # irm3 is a stimulation-of-production model: dosing robustly pushes
    # RESP *above* RESP0 (unlike irm1/irm2's inhibition models, where
    # dosing pushes RESP below baseline -- see irm1_hazard.cpp's own
    # doc: "positive makes a response *above* baseline harmful", which
    # for an inhibition-type backbone means beta_r > 0 is the
    # *protective* direction instead, since RESP falls below baseline
    # there; irm3 is the unambiguous case for this direction check).
    n <- 800
    dose <- .pkpd_dose_data(n, amt = 100)
    sim_harm <- sim_tte_ode(model = "irm3_hazard",
        param = list(H0 = 0.01, beta_r = 3), n = n, end = 30, delta = 1,
        data = dose, seed = 1)
    sim_prot <- sim_tte_ode(model = "irm3_hazard",
        param = list(H0 = 0.01, beta_r = -3), n = n, end = 30, delta = 1,
        data = dose, seed = 1)
    p_event_harm <- mean(sim_harm$events$sim_status == 1)
    p_event_prot <- mean(sim_prot$events$sim_status == 1)
    expect_gt(p_event_harm, p_event_prot)
})

# ---------------------------------------------------------------------
# 9. BSV propagation [CRAN-tier, cheap: no model compilation past what
#    section 1 already exercises]. UPDATED after the BSV review/
#    implementation session (reports/11_bsv_review.md,
#    reports/12_bsv_implementation_report.md): Phase 4's finding that
#    'omega=' could not add BSV to any of the six models (the original
#    version of this test asserted it errors even on pk_hazard) turned
#    out to be fixable, cheaply, via per-subject idata columns rather
#    than a declared $OMEGA block -- see test-sim-tte-ode-bsv.R for the
#    full BSV test suite. 'sigma=' was reviewed too and correctly found
#    to need no such fix (no model's hazard depends on a residual-
#    error-perturbed quantity) -- that half of this test is unchanged.
#    What's still tested here: a model with neither a declared block
#    nor a BSV registry entry (any of the four pre-PK/PD models) still
#    gets the informative error .apply_ode_matlist() added in Phase 4.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() 'omega'/'sigma' raise an informative error on a model with neither a declared block nor a BSV registry entry", {
    skip_if_not_installed("mrgsolve")
    om <- matrix(c(0.3, 0, 0, 0.3), 2, 2)
    expect_error(sim_tte_ode(model = "exponential",
        param = list(H0 = 0.1), omega = om[1, 1, drop = FALSE], n = 5,
        end = 10, delta = 2, seed = 1),
        "Could not apply 'omega'.*does not declare a matching")
    expect_error(sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        sigma = matrix(0.1, 1, 1), n = 5, end = 10, delta = 2, seed = 1),
        "Could not apply 'sigma'.*does not declare a matching")
})

# ---------------------------------------------------------------------
# 10. Fallback classification, real-world integration [slow]: a
#     genuinely steep run (tmdd_hazard, TMDD's own stiff binding
#     kinetics, decision-1's "steepness stress test" case) must not
#     error, and .resolve_ode_events()'s differentiated message (if any
#     subject falls back) must be one of the three known signatures.
# ---------------------------------------------------------------------
test_that("sim_tte_ode() tmdd_hazard under a coarse delta does not error, and any fallback is classified [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    n <- 300
    dose <- .pkpd_dose_data(n, amt = 500)
    msgs <- capture_messages(sim <- sim_tte_ode(model = "tmdd_hazard",
        param = list(H0 = 0.02, beta_rc = 0.3), n = n, end = 30, delta = 4,
        data = dose, seed = 1))
    expect_equal(nrow(sim$events), n)
    fallback_msgs <- grep("used the reported-grid", msgs, value = TRUE)
    if (length(fallback_msgs)) {
        expect_match(fallback_msgs[1], "unphysical P_POST|T_PRE == TEVT|unclassified")
    }
})
