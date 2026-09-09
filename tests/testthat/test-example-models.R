# Phase H: shipped example mrgsolve models (inst/models/examples/) and
# their discovery/loading helpers, simtte_example_models() and
# simtte_example_model(). These are user-facing templates for
# sim_tte_df(), entirely separate from the internal
# weibull/weibull_tv/ms engine loaded by .read_model_static_cache()
# (untouched by this phase; see test-protected-dots.R /
# test-weibull-analytical.R / test-ms-hazard-carry.R for that engine's
# own coverage).

test_that("simtte_example_models() lists exactly the two shipped examples", {
    expect_setequal(simtte_example_models(),
        c("pkpd_idr_hazard", "pkpd_linear_hazard"))
})

test_that("simtte_example_model() loads each shipped example without error", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    for (nm in simtte_example_models()) {
        mod <- simtte_example_model(nm)
        expect_s4_class(mod, "mrgmod")
    }
})

test_that("simtte_example_model() rejects an unknown model name informatively", {
    expect_error(simtte_example_model("not_a_real_model"),
        "Unknown example model 'not_a_real_model'")
    expect_error(simtte_example_model("not_a_real_model"),
        "pkpd_idr_hazard")
    expect_error(simtte_example_model("not_a_real_model"),
        "pkpd_linear_hazard")
})

test_that("simtte_example_model() validates 'name' before touching the filesystem", {
    expect_error(simtte_example_model(c("pkpd_idr_hazard", "pkpd_linear_hazard")),
        "single character string")
    expect_error(simtte_example_model(1), "single character string")
    expect_error(simtte_example_model(NA_character_), "single character string")
})

test_that("mrgsim() smoke test: pkpd_idr_hazard produces a valid p11 trajectory", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mod <- simtte_example_model("pkpd_idr_hazard")
    data <- mrgsolve::expand.ev(ID = 1:5, amt = 100, cmt = 1, ii = 1,
        addl = 23, time = 0)
    # obsonly = TRUE drops the internal dosing/bookkeeping rows mrgsim()
    # otherwise reports (which would duplicate time = 0: one row before
    # and one after the dose is applied) -- required for output destined
    # for sim_tte_df(), whose trajectory contract rejects duplicated
    # times within a subject; see the "M-spline hazard carry convention"
    # section of ?sim_tte for the same obsonly = TRUE convention used
    # internally throughout the package.
    out <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = 24,
        delta = 1, obsonly = TRUE))
    expect_true(all(c("ID", "time", "p11") %in% names(out)))
    expect_equal(anyDuplicated(out$time[out$ID == 1]), 0L)
    expect_true(all(is.finite(out$p11)))
    expect_true(all(out$p11 >= 0 & out$p11 <= 1))
    for (id in unique(out$ID)) {
        p <- out$p11[out$ID == id]
        expect_true(all(diff(p) <= 1e-8)) # non-increasing (solver tolerance)
    }
})

test_that("mrgsim() smoke test: pkpd_linear_hazard produces a valid p11 trajectory", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    mod <- simtte_example_model("pkpd_linear_hazard")
    data <- mrgsolve::ev(amt = 100, cmt = 1, time = 0)
    out <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = -1,
        add = seq(0, 10, by = 1), obsonly = TRUE))
    expect_true(all(c("ID", "time", "p11") %in% names(out)))
    expect_true(all(is.finite(out$p11)))
    expect_true(all(out$p11 >= 0 & out$p11 <= 1))
    expect_true(all(diff(out$p11) <= 1e-8))
})

test_that("both shipped examples run end to end through sim_tte_df()", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")

    mod_idr <- simtte_example_model("pkpd_idr_hazard")
    data_idr <- mrgsolve::expand.ev(ID = 1:8, amt = 100, cmt = 1, ii = 1,
        addl = 23, time = 0)
    out_idr <- as.data.frame(mrgsolve::mrgsim(mod_idr, data = data_idr,
        end = 24, delta = 1, obsonly = TRUE))
    set.seed(1)
    res_idr <- sim_tte_df(out_idr[, c("ID", "time", "p11")])
    expect_equal(nrow(res_idr), 8)
    expect_true(all(res_idr$sim_status %in% c(0, 1)))
    expect_true(all(is.finite(res_idr$sim_time)))

    mod_lin <- simtte_example_model("pkpd_linear_hazard")
    data_lin <- mrgsolve::ev(amt = 100, cmt = 1, time = 0)
    out_lin <- as.data.frame(mrgsolve::mrgsim(mod_lin, data = data_lin,
        end = -1, add = seq(0, 60, by = 0.5), obsonly = TRUE))
    traj_lin <- out_lin[rep(seq_len(nrow(out_lin)), 5),
        c("ID", "time", "p11")]
    traj_lin$ID <- rep(1:5, each = nrow(out_lin))
    set.seed(1)
    res_lin <- sim_tte_df(traj_lin)
    expect_equal(nrow(res_lin), 5)
    expect_true(all(res_lin$sim_status %in% c(0, 1)))
})

test_that("simtte_example_model() is not usable as a sim_tte() 'type'", {
    # Documents/enforces the integration contract: example models are
    # for sim_tte_df(), never for sim_tte(), which only ever simulates
    # its own built-in weibull/ms models via .read_model_static_cache().
    expect_error(
        sim_tte(pi = 0, mu = -1, coefs = 1, time = c(0, 1),
            type = "pkpd_idr_hazard"),
        "'arg' should be" # match.arg() failure on an unsupported 'type'
    )
})
