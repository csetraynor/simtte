# BSV/omega implementation (reports/11_bsv_review.md, reports/12_bsv_implementation_report.md):
# per-subject PK/PD parameters via idata, no model file changes. See
# R/helpers.R (.ODE_BSV_TARGETS, .build_ode_bsv_idata()) and
# sim_tte_ode()'s omega dispatch. Fast/slow split per the runbook
# convention: cheap/deterministic checks (exactness, reproducibility,
# zero-omega, input validation) run fast; the two checks that need a
# large n to be non-flaky (spread, the distributional check) are slow.

.bsv_om <- function(vars, values) {
    k <- length(vars)
    m <- diag(values, k, k)
    dimnames(m) <- list(vars, vars)
    m
}

# ---------------------------------------------------------------------
# 1. .ODE_BSV_TARGETS itself: sanity that every entry matches what the
#    compiled model actually declares (defensive regression -- these
#    are typed by hand from mread() output, reports/12_bsv_implementation_report.md).
# ---------------------------------------------------------------------
test_that(".ODE_BSV_TARGETS matches each model's own declared PARAM set (minus the hazard scaffold)", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    scaffold <- c("H0", "lp", "beta_cp", "beta_r", "beta_rc", "U", "END")
    for (model in names(simtte:::.ODE_BSV_TARGETS)) {
        mod <- simtte:::.read_ode_library_model(model)
        actual <- setdiff(names(mrgsolve::param(mod)), scaffold)
        expect_setequal(simtte:::.ODE_BSV_TARGETS[[model]], actual)
    }
})

# ---------------------------------------------------------------------
# 2. Exactness: .build_ode_bsv_idata() is the function that draws and
#    carries the values through; tested directly (not via
#    sim_tte_ode()'s public return value) because none of the six
#    models captures CL/V2/etc. in $CAPTURE and carry_out is fixed
#    inside sim_tte_ode() -- there is no public way to observe the
#    per-subject value exactly outside this function, so this is the
#    right level to test the "carried value equals the drawn value"
#    property at [CRAN].
# ---------------------------------------------------------------------
test_that(".build_ode_bsv_idata() carries the drawn value through exactly [CRAN]", {
    om <- .bsv_om(c("CL", "V2"), c(0.09, 0.16))
    idata <- data.frame(ID = 1:8, U = 0.5, END = 30)
    set.seed(1)
    eta <- mrgsolve::mvgauss(om, n = 8)
    set.seed(1)
    out <- simtte:::.build_ode_bsv_idata(om, "pk_hazard", n = 8, idata = idata,
        param = list(CL = 1, V2 = 20))
    expect_equal(out$CL, 1 * exp(eta[, 1]), tolerance = 0)
    expect_equal(out$V2, 20 * exp(eta[, 2]), tolerance = 0)
})

test_that("BSV reaches the ODE: sd(CP) is nonzero with omega, zero without [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    om <- .bsv_om(c("CL", "V2"), c(0.09, 0.16))
    dose <- .pkpd_dose_data(60)
    sim_bsv <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om, n = 60, end = 30, delta = 1,
        data = dose, seed = 1, keep_trajectory = TRUE)
    sim_pop <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), n = 60, end = 30, delta = 1, data = dose,
        seed = 1, keep_trajectory = TRUE)
    sd_bsv <- stats::sd(sim_bsv$trajectory$CP[sim_bsv$trajectory$time == 10])
    sd_pop <- stats::sd(sim_pop$trajectory$CP[sim_pop$trajectory$time == 10])
    expect_gt(sd_bsv, 0)
    expect_equal(sd_pop, 0)
})

# ---------------------------------------------------------------------
# 3. Reproducibility.
# ---------------------------------------------------------------------
test_that("same seed reproduces $events exactly with omega; different seed differs [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    om <- .bsv_om(c("CL", "V2"), c(0.09, 0.16))
    dose <- .pkpd_dose_data(30)
    s1 <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om, n = 30, end = 30, delta = 1,
        data = dose, seed = 99)
    s2 <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om, n = 30, end = 30, delta = 1,
        data = dose, seed = 99)
    s3 <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om, n = 30, end = 30, delta = 1,
        data = dose, seed = 100)
    expect_identical(s1$events, s2$events)
    expect_false(isTRUE(all.equal(s1$events, s3$events)))
})

# ---------------------------------------------------------------------
# 4. Additivity: omega = NULL is untouched. Not a new assertion here so
#    much as a pointer -- the real evidence is the full
#    test-sim-tte-ode-pkpd.R suite (90 tests) passing unchanged; this
#    is a lightweight in-file smoke check that omega = NULL still skips
#    the whole BSV/omat() branch (no idata columns injected, no error).
# ---------------------------------------------------------------------
test_that("omega = NULL is unaffected [CRAN]", {
    skip_if_not_installed("mrgsolve")
    sim <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), n = 10, end = 20, delta = 2,
        data = .pkpd_dose_data(10), seed = 1)
    expect_identical(names(sim$events), c("ID", "sim_time", "sim_status"))
})

# ---------------------------------------------------------------------
# 5. Zero-variance omega reproduces the population run exactly.
# ---------------------------------------------------------------------
test_that("omega of all zeros reproduces the population (no-omega) run exactly [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    om0 <- .bsv_om(c("CL", "V2"), c(0, 0))
    dose <- .pkpd_dose_data(40)
    sim_zero <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om0, n = 40, end = 30, delta = 1,
        data = dose, seed = 42)
    sim_pop <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), n = 40, end = 30, delta = 1, data = dose,
        seed = 42)
    expect_identical(sim_zero$events, sim_pop$events)
})

# ---------------------------------------------------------------------
# 6. Named subset, unnamed positional, and input-validation errors.
# ---------------------------------------------------------------------
test_that("a named omega subset (CL only) is accepted and applied [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    om <- .bsv_om("CL", 0.09)
    sim <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.02), omega = om, n = 30, end = 20, delta = 2,
        data = .pkpd_dose_data(30), seed = 1, keep_trajectory = TRUE)
    expect_gt(stats::sd(sim$trajectory$CP[sim$trajectory$time == 10]), 0)
})

test_that("an unnamed omega must match the full target-list dimension [CRAN]", {
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01),
        omega = matrix(c(0.09, 0, 0, 0.16), 2, 2), n = 5, end = 10,
        delta = 2, data = .pkpd_dose_data(5), seed = 1),
        "must have dimension 8 x 8")
})

test_that("an omega naming an unknown parameter errors, listing the valid targets [CRAN]", {
    skip_if_not_installed("mrgsolve")
    om <- .bsv_om("NOT_A_PARAM", 0.09)
    expect_error(sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01),
        omega = om, n = 5, end = 10, delta = 2, data = .pkpd_dose_data(5),
        seed = 1), "not in model.*between-subject-variability targets.*Valid targets")
})

test_that("omega on a model without a registry entry keeps the existing informative error [CRAN]", {
    skip_if_not_installed("mrgsolve")
    expect_error(sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        omega = matrix(0.1, 1, 1), n = 5, end = 10, delta = 2, seed = 1),
        "Could not apply 'omega'.*does not declare a matching")
})

test_that("an idata column clashing with an omega target errors [CRAN]", {
    skip_if_not_installed("mrgsolve")
    om <- .bsv_om("CL", 0.09)
    idata <- data.frame(ID = 1:5, CL = 1)
    expect_error(sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01),
        omega = om, idata = idata, end = 10, delta = 2,
        data = .pkpd_dose_data(5), seed = 1),
        "already has column.*CL")
})

test_that("mismatched omega row/column names error [CRAN]", {
    skip_if_not_installed("mrgsolve")
    om <- matrix(c(0.09, 0, 0, 0.16), 2, 2,
        dimnames = list(c("CL", "V2"), c("CL", "WRONG")))
    expect_error(sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01),
        omega = om, n = 5, end = 10, delta = 2, data = .pkpd_dose_data(5),
        seed = 1), "row and column names must be identical")
})

# ---------------------------------------------------------------------
# 7. Coexistence: a user-supplied model with its own declared $OMEGA
#    still goes through omat(), not the idata registry. sim_tte_ode()'s
#    own `model` argument only accepts the registered public names
#    (Phase 5's general/mrgmod dispatch is not built yet), so this
#    tests the exact boolean condition sim_tte_ode()'s own dispatch
#    branches on (nrow(omat(mod, make = TRUE)) > 0) directly, plus
#    .apply_ode_matlist()'s mechanics on such a model -- the same
#    experiment reports/experiments/pk_hazard_with_omega.cpp already
#    demonstrated (that file lives under the gitignored reports/, so
#    the model here is rebuilt inline from the identical minimal edit --
#    pk_hazard.cpp verbatim plus one $OMEGA block -- rather than
#    depending on a path this test suite does not ship).
# ---------------------------------------------------------------------
test_that("a model with its own declared $OMEGA block dispatches through omat(), not the idata registry [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    src <- readLines(system.file("models", "library", "pk_hazard.cpp",
        package = "simtte"))
    main_idx <- which(src == "[MAIN]")[1]
    src <- append(src, c("[OMEGA] @labels ETA_CL ETA_V2", "  0 0", ""),
        after = main_idx - 1L)
    tmp_dir <- tempfile("bsv_coexist_")
    dir.create(tmp_dir)
    writeLines(src, file.path(tmp_dir, "pk_hazard_coexist.cpp"))
    mod_bare <- simtte:::.read_ode_library_model("pk_hazard")
    mod_with_omega <- mrgsolve::mread(model = "pk_hazard_coexist",
        project = tmp_dir)

    expect_equal(nrow(mrgsolve::omat(mod_bare, make = TRUE)), 0L)
    expect_equal(nrow(mrgsolve::omat(mod_with_omega, make = TRUE)), 2L)

    updated <- simtte:::.apply_ode_matlist(mod_with_omega,
        matrix(c(0.09, 0, 0, 0.16), 2, 2), "omega", mrgsolve::omat)
    expect_equal(diag(mrgsolve::omat(updated, make = TRUE)), c(0.09, 0.16))
})

# ---------------------------------------------------------------------
# 8. Spread and the distributional check: both need a large n to be
#    non-flaky (verified directly during implementation -- at n = 500
#    on pk_hazard, the naive with-vs-without sd(sim_time) comparison
#    flipped direction on 3/5 tried seeds; at n = 2000 it was TRUE on
#    5/5 -- so these are the two checks explicitly called out as slow,
#    not left fast "where possible").
# ---------------------------------------------------------------------
test_that("event-time spread is larger with omega than without, same seed [slow, pk_hazard]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    om <- .bsv_om(c("CL", "V2"), c(2, 2))
    dose <- .pkpd_dose_data(2000)
    sim_bsv <- suppressMessages(sim_tte_ode(model = "pk_hazard",
        param = list(H0 = 0.01, beta_cp = 0.3), omega = om, n = 2000,
        end = 30, delta = 1, data = dose, seed = 42))
    sim_pop <- sim_tte_ode(model = "pk_hazard", param = list(H0 = 0.01,
        beta_cp = 0.3), n = 2000, end = 30, delta = 1, data = dose,
        seed = 42)
    et_bsv <- sim_bsv$events$sim_time[sim_bsv$events$sim_status == 1]
    et_pop <- sim_pop$events$sim_time[sim_pop$events$sim_status == 1]
    expect_gt(stats::sd(et_bsv), stats::sd(et_pop))
})

test_that("event-time spread is dramatically larger with omega than without [slow, tmdd_hazard]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    om <- .bsv_om("V2", 1)
    dose <- .pkpd_dose_data(2000, amt = 500)
    sim_bsv <- sim_tte_ode(model = "tmdd_hazard", param = list(H0 = 0.02,
        beta_rc = 0.3), omega = om, n = 2000, end = 30, delta = 1,
        data = dose, seed = 42)
    sim_pop <- sim_tte_ode(model = "tmdd_hazard", param = list(H0 = 0.02,
        beta_rc = 0.3), n = 2000, end = 30, delta = 1, data = dose,
        seed = 42)
    sd_bsv <- stats::sd(sim_bsv$events$sim_time[sim_bsv$events$sim_status == 1])
    sd_pop <- stats::sd(sim_pop$events$sim_time[sim_pop$events$sim_status == 1])
    expect_gt(sd_bsv, 10 * sd_pop)
})

test_that("sd(log(CL_i)) recovers the requested omega variance [slow, distributional]", {
    skip_on_cran()
    skip_if_not_slow()
    om <- .bsv_om("CL", 0.09)
    idata <- data.frame(ID = 1:2000, U = 0.5, END = 30)
    out <- simtte:::.build_ode_bsv_idata(om, "pk_hazard", n = 2000,
        idata = idata, param = list(CL = 1))
    expect_equal(stats::sd(log(out$CL)), sqrt(0.09), tolerance = 0.05)
})
