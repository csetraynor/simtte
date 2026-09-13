# sim_tte_ode_models(): public API freeze (reports/27_public_api_freeze.md).
# Pure data lookup -- no ODE solve, so no fast/slow split needed. Matches
# `^test-sim-tte-ode-` so it joins the `ode`/`all` groups automatically,
# with no dev/run-tests.R change required.

test_that("sim_tte_ode_models() lists every built-in model name, matching sim_tte_ode()'s own match.arg() choices", {
    models <- sim_tte_ode_models()
    expect_s3_class(models, "data.frame")
    expect_identical(names(models), c("model", "bsv_targets"))
    expect_setequal(models$model,
        c("exponential", "weibull", "gompertz", "mspline", "pk_hazard",
            "irm1_hazard", "irm2_hazard", "irm3_hazard", "irm4_hazard",
            "tmdd_hazard"))
})

test_that("every listed model name is actually accepted by sim_tte_ode()'s match.arg()", {
    models <- sim_tte_ode_models()$model
    for (m in models) {
        expect_no_error(match.arg(m, choices = c(names(.ODE_LIBRARY_FILES), "mspline")))
    }
})

test_that("bsv_targets is NA for models with no BSV route, comma-separated names otherwise", {
    models <- sim_tte_ode_models()
    no_bsv <- models$bsv_targets[models$model %in% c("exponential", "weibull",
        "gompertz", "mspline")]
    expect_true(all(is.na(no_bsv)))

    pkpd <- models$bsv_targets[models$model == "irm1_hazard"]
    expect_false(is.na(pkpd))
    expect_true(grepl("CL", pkpd))
})

test_that("sim_tte_ode_models()'s bsv_targets agrees with .ODE_BSV_TARGETS exactly", {
    models <- sim_tte_ode_models()
    for (m in names(.ODE_BSV_TARGETS)) {
        expected <- paste(.ODE_BSV_TARGETS[[m]], collapse = ", ")
        expect_identical(models$bsv_targets[models$model == m], expected)
    }
})
