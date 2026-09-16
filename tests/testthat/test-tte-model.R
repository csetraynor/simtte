# Phase 5a: tte_model(), the PK/PD-model converter, and sim_tte_ode()'s
# acceptance of its output. See reports/13_converter_design.md for the
# design and reports/14_phase5a_report.md for what was built. Every
# mrgsolve mechanic this file exercises (autodec/HAZ, $DES as an $ODE
# synonym, [BLOCK]/$BLOCK, singleton GLOBAL/MAIN, compile = FALSE
# introspection) was verified directly under reports/experiments/
# before being relied on here -- see the design report's citations.

# ---------------------------------------------------------------------
# 1. Hard precondition: $ODE/$DES required [fast, no compile].
# ---------------------------------------------------------------------
test_that("tte_model() rejects a $PKMODEL-only model, naming an ODE-form equivalent", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20",
        "$PKMODEL ncmt = 1")
    expect_error(tte_model(code, hazard = "H0", params = list(H0 = 0.1)),
        "PKMODEL.*ODE-form equivalent|ODE-form equivalent")
})

test_that("tte_model() rejects a model with no $ODE/$DES/$PKMODEL/$PRED block at all", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT")
    expect_error(tte_model(code, hazard = "H0", params = list(H0 = 0.1)),
        "no \\$ODE/\\$DES block")
})

test_that("tte_model() accepts a $DES block exactly like $ODE", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$DES", "dxdt_CENT = -(CL/V)*CENT;")
    tm <- tte_model(code, hazard = "0.05", params = list())
    expect_s3_class(tm, "simtte_model")
    expect_true("p11" %in% names(mrgsolve::init(tm$mod)))
})

test_that("tte_model() accepts a model mixing $BLOCK and [BLOCK] header syntax", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "[CMT] CENT",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    tm <- tte_model(code, hazard = "0.05", params = list())
    expect_s3_class(tm, "simtte_model")
})

# ---------------------------------------------------------------------
# 2. GLOBAL/MAIN: created when absent, edited in place when present
#    [fast].
# ---------------------------------------------------------------------
test_that("tte_model() creates new $GLOBAL/$MAIN blocks when the model has neither", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    tm <- tte_model(code, hazard = "0.05", params = list())
    expect_true(any(grepl("^\\$GLOBAL", tm$code)))
    expect_true(any(grepl("^\\$MAIN", tm$code)))
    expect_equal(sum(grepl("^\\$GLOBAL", tm$code)), 1L)
    expect_equal(sum(grepl("^\\$MAIN", tm$code)), 1L)
})

test_that("tte_model() edits an existing $GLOBAL/$MAIN in place, not a second block", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$GLOBAL", "#define TWO_CL (2*CL)",
        "$MAIN", "double dummy_local = CL;",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    tm <- tte_model(code, hazard = "0.05", params = list())
    expect_equal(sum(grepl("^\\$GLOBAL", tm$code)), 1L)
    expect_equal(sum(grepl("^\\$MAIN", tm$code)), 1L)
    expect_true(any(grepl("TWO_CL", tm$code)))
    expect_true(any(grepl("dummy_local", tm$code)))
    expect_true(any(grepl("event_found", tm$code)))
})

# ---------------------------------------------------------------------
# 3. Name collisions [fast, no compile: caught by the compile = FALSE
#    introspection before any C++ build is attempted].
# ---------------------------------------------------------------------
test_that("tte_model() rejects a model that already declares a reserved scaffold name", {
    skip_if_not_installed("mrgsolve")
    for (nm in c("p11", "U", "END", "HAZ", "TEVT", "event_found",
        "T_PRE", "P_PRE", "P_POST", "lp")) {
        code <- c(paste0("$PARAM ", nm, " = 1, V = 20"), "$CMT CENT",
            "$ODE", "dxdt_CENT = -(1/V)*CENT;")
        expect_error(tte_model(code, hazard = "0.05", params = list()),
            "already declares", info = nm)
    }
})

test_that("tte_model() rejects a 'params' name that clashes with an existing parameter", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    expect_error(tte_model(code, hazard = "H0", params = list(CL = 1)),
        "'params' names CL")
})

test_that("tte_model() rejects 'params' naming lp/U/END (always added automatically)", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    for (nm in c("lp", "U", "END")) {
        params <- setNames(list(1), nm)
        expect_error(tte_model(code, hazard = "H0", params = params),
            "always added automatically", info = nm)
    }
})

# ---------------------------------------------------------------------
# 4. Round-trip [fast].
# ---------------------------------------------------------------------
test_that("tte_model() refuses to convert an already-converted simtte_model", {
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    expect_error(tte_model(tm, hazard = "0.05", params = list()),
        "already a simtte_model")
})

# ---------------------------------------------------------------------
# 5. Compile failure surfaces the converted code [fast: a tiny model].
# ---------------------------------------------------------------------
test_that("tte_model() attaches the converted code to a hazard-expression compile failure", {
    skip_if_not_installed("mrgsolve")
    code <- c("$PARAM CL = 1, V = 20", "$CMT CENT",
        "$ODE", "dxdt_CENT = -(CL/V)*CENT;")
    expect_error(tte_model(code, hazard = "NOT_A_REAL_IDENTIFIER",
        params = list()), "converted model code")
})

# ---------------------------------------------------------------------
# 6. print()/code accessor/contract validation on a real converted
#    model [fast: one compile, cached and reused by later tests too].
# ---------------------------------------------------------------------
test_that("a converted model prints its hazard/parameters/BSV route", {
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    out <- capture.output(print(tm))
    expect_true(any(grepl("HAZ = H0", out)))
    expect_true(any(grepl("H0 = ", out)))
    expect_true(any(grepl("between-subject variability", out)))
})

test_that("a converted model's $code is directly printable and contract-valid", {
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    expect_type(tm$code, "character")
    expect_true(any(grepl("^\\$PARAM", tm$code)))
    expect_true(simtte:::.validate_ode_model_contract(tm$mod))
})

# ---------------------------------------------------------------------
# 7. sim_tte_ode() acceptance [fast: reuses the cached fixture].
# ---------------------------------------------------------------------
test_that("sim_tte_ode() rejects a bare (unconverted) compiled mrgmod", {
    skip_if_not_installed("mrgsolve")
    mod <- .modlib_model("pk2cmt", compile = FALSE)
    expect_error(sim_tte_ode(model = mod, n = 5, end = 10, delta = 2),
        "convert it first with tte_model")
})

test_that("sim_tte_ode() runs a converted simtte_model and returns well-formed $events", {
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    sim <- sim_tte_ode(model = tm, n = 20, end = 20, delta = 2,
        data = .pkpd_dose_data(20), seed = 1)
    expect_s3_class(sim, "simtte_ode_sim")
    expect_named(sim$events, c("ID", "sim_time", "sim_status", "sim_reason"))
    expect_equal(nrow(sim$events), 20)
})

test_that("sim_tte_ode() rejects knots/boundary_knots/coefs for a converted model", {
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    expect_error(sim_tte_ode(model = tm, knots = c(1, 2), coefs = rep(1, 5),
        n = 5, end = 10, delta = 2), "only used when model = \"mspline\"")
})

# ---------------------------------------------------------------------
# 8. Six-backbone self-consistency [slow: 12 real compiles].
# ---------------------------------------------------------------------
test_that("tte_model() reproduces every shipped *_hazard.cpp model's $events exactly [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    dose <- .pkpd_dose_data(20)
    for (nm in names(.PKPD_CONVERTER_SPECS)) {
        sp <- .PKPD_CONVERTER_SPECS[[nm]]
        tm <- tte_model(.modlib_model(sp$backbone, compile = FALSE),
            hazard = sp$hazard, params = .PKPD_TEST_DEFAULT_PARAM[[nm]])
        sim_conv <- sim_tte_ode(model = tm, n = 20, end = 20, delta = 2,
            data = dose, seed = 1)
        sim_ship <- sim_tte_ode(model = nm,
            param = .PKPD_TEST_DEFAULT_PARAM[[nm]], n = 20, end = 20,
            delta = 2, data = dose, seed = 1)
        expect_identical(sim_conv$events, sim_ship$events, info = nm)
    }
})

# ---------------------------------------------------------------------
# 9. A converted model's own declared $OMEGA dispatches through omat(),
#    end to end via the public API [slow: closes
#    reports/12_bsv_implementation_report.md open risk 1].
# ---------------------------------------------------------------------
test_that("a converted model with its own declared $OMEGA goes through omat() via sim_tte_ode() [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    base <- .modlib_model("pk2cmt", compile = FALSE)
    code <- base@code
    code <- sub("^CL   :  1  : Clearance \\(volume/time\\)$",
        "TVCL :  1  : Clearance (volume/time)", code)
    code <- sub("^V2   : 20  : Central volume \\(volume\\)$",
        "TVV2 : 20  : Central volume (volume)", code)
    cmt_idx <- which(code == "$CMT  @annotated")[1]
    code <- append(code, c(
        "$OMEGA @labels ETA_CL ETA_V2", "  0 0", "",
        "$MAIN",
        "  double CL = TVCL*exp(ETA_CL);",
        "  double V2 = TVV2*exp(ETA_V2);", ""), after = cmt_idx - 1L)

    tm <- tte_model(code, hazard = "H0 * exp(lp + beta_cp * CP)",
        params = list(H0 = 0.01, beta_cp = 0.3))
    expect_gt(nrow(mrgsolve::omat(tm$mod, make = TRUE)), 0L)

    n <- 1500
    dose <- .pkpd_dose_data(n, amt = 100)
    om <- diag(c(2, 2))
    dimnames(om) <- list(c("ETA_CL", "ETA_V2"), c("ETA_CL", "ETA_V2"))
    sim_pop <- sim_tte_ode(model = tm, n = n, end = 20, delta = 2,
        data = dose, seed = 1)
    sim_bsv <- sim_tte_ode(model = tm, omega = om, n = n, end = 20,
        delta = 2, data = dose, seed = 1)
    expect_gt(sd(sim_bsv$events$sim_time), sd(sim_pop$events$sim_time))
})

# ---------------------------------------------------------------------
# 10. A converted model's own bsv_targets goes through idata [slow].
#     tmdd_hazard's own conversion is used (author decision, "After the
#     BSV implementation": tmdd's spread effect is the robust one, not
#     pk_hazard's).
# ---------------------------------------------------------------------
test_that("a converted model's bsv_targets go through the idata route [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    # Same param/omega/dose/seed as the already-validated shipped-model
    # spread test (test-sim-tte-ode-bsv.R "event-time spread is
    # dramatically larger ... [slow, tmdd_hazard]") -- the point here
    # is that a *converted* model reproduces that same robust effect
    # via its own bsv_targets, not a fresh parameter search.
    sp <- .PKPD_CONVERTER_SPECS$tmdd_hazard
    tm <- tte_model(.modlib_model(sp$backbone, compile = FALSE),
        hazard = sp$hazard, params = list(H0 = 0.02, beta_rc = 0.3),
        bsv_targets = .ODE_BSV_TARGETS$tmdd_hazard)
    expect_equal(nrow(mrgsolve::omat(tm$mod, make = TRUE)), 0L)

    n <- 2000
    dose <- .pkpd_dose_data(n, amt = 500)
    om <- matrix(1, dimnames = list("V2", "V2"))
    sim_pop <- sim_tte_ode(model = tm, n = n, end = 30, delta = 1,
        data = dose, seed = 42)
    sim_bsv <- sim_tte_ode(model = tm, omega = om, n = n, end = 30,
        delta = 1, data = dose, seed = 42)
    sd_pop <- sd(sim_pop$events$sim_time[sim_pop$events$sim_status == 1])
    sd_bsv <- sd(sim_bsv$events$sim_time[sim_bsv$events$sim_status == 1])
    expect_gt(sd_bsv, 10 * sd_pop)
})

# ---------------------------------------------------------------------
# 11. A converted model composes with time-varying covariates and
#     dosing data [slow: exercises the full argument surface together].
# ---------------------------------------------------------------------
test_that("a converted model composes with covariates/beta and dosing data [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    tm <- .tte_model_fixture()
    n <- 40
    dose <- .pkpd_dose_data(n)
    cov <- data.frame(time = c(0, 5), sex = c(0, 1))
    # Used to warn "Parameter column lp must not contain missing
    # values" here (the dosing rows from 'data' carried no 'lp' column
    # while the covariate-update rows did) -- fixed by
    # .merge_ode_covariate_rows() (simttepower feedback 1: lp/data
    # merge, reports/29_simttepower_feedback_lp_merge.md), which fills
    # 'lp' on every dosing row by LOCF before the two are combined, for
    # a tte_model()-converted model exactly as for a shipped
    # *_hazard.cpp one.
    sim <- expect_no_warning(sim_tte_ode(model = tm, n = n, end = 20,
        delta = 2, data = dose, covariates = cov, beta = c(sex = 0.5),
        seed = 1))
    expect_equal(nrow(sim$events), n)
    expect_true(all(sim$events$sim_status %in% c(0L, 1L)))
})
