# Validate 'n' (or 'idata') against a dosing 'data' frame's own unique
# 'ID' set ("simttepower feedback 3" Task B, finding 5.3.2,
# reports/31_simttepower_feedback_process_loading.md): before this
# fix, 'n' larger than length(unique(data$ID)) silently simulated and
# returned only length(unique(data$ID)) subjects (no warning), and 'n'
# smaller failed with mrgsolve's own opaque "ID found in the data set,
# but not in idata." Both directions now error immediately from
# sim_tte_ode() itself, naming both numbers.

test_that("'n' larger than 'data's unique ID count errors immediately, naming both numbers", {
    dosing <- .pkpd_dose_data(5)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 10, end = 20,
            delta = 2, data = dosing, seed = 1),
        "'n' \\(10\\) does not match.*'data'.*\\(5\\)"
    )
})

test_that("'n' smaller than 'data's unique ID count errors immediately, naming both numbers (not mrgsolve's opaque error)", {
    dosing <- .pkpd_dose_data(5)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 3, end = 20,
            delta = 2, data = dosing, seed = 1),
        "'n' \\(3\\) does not match.*'data'.*\\(5\\)"
    )
})

test_that("'n' matching 'data's unique ID count still works, all subjects present", {
    dosing <- .pkpd_dose_data(5)
    sim <- sim_tte_ode(model = "pk_hazard",
        param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 5, end = 20,
        delta = 2, data = dosing, seed = 1)
    expect_setequal(sim$events$ID, 1:5)
})

test_that("'data' without an explicit 'n' infers 'n' from data's unique ID count", {
    dosing <- .pkpd_dose_data(5)
    sim <- sim_tte_ode(model = "pk_hazard",
        param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, end = 20, delta = 2,
        data = dosing, seed = 1)
    expect_setequal(sim$events$ID, 1:5)
    # matches the explicit n = 5 call exactly, same seed
    sim_explicit <- sim_tte_ode(model = "pk_hazard",
        param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 5, end = 20,
        delta = 2, data = dosing, seed = 1)
    expect_identical(sim$events, sim_explicit$events)
})

test_that("an explicit 'idata' with a row count mismatched against 'data' errors immediately", {
    dosing <- .pkpd_dose_data(5)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard,
            idata = data.frame(ID = 1:3), end = 20, delta = 2,
            data = dosing, seed = 1),
        "sim_tte_ode\\(\\) resolved 3 subject"
    )
})

test_that("'data' IDs and idata IDs of matching count but different values errors, naming the unknown ID(s)", {
    dosing <- data.frame(ID = 3:7, time = 0, cmt = 1, amt = 100, evid = 1)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 5, end = 20,
            delta = 2, data = dosing, seed = 1),
        "not present in the resolved 'idata'.*6, 7"
    )
})

test_that("'data' with no 'ID' column and more than one resolved subject errors clearly (previously silently dropped subjects)", {
    dosing <- data.frame(time = 0, cmt = 1, amt = 100, evid = 1)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 2, end = 20,
            delta = 2, data = dosing, seed = 1),
        "'data' has no 'ID' column"
    )
})

test_that("'data' with no 'ID' column and exactly one resolved subject still works (unambiguous)", {
    dosing <- data.frame(time = 0, cmt = 1, amt = 100, evid = 1)
    sim <- sim_tte_ode(model = "pk_hazard",
        param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 1, end = 20,
        delta = 2, data = dosing, seed = 1)
    expect_equal(nrow(sim$events), 1)
})

test_that("'n' mismatch against 'data' is checked even without covariates/beta", {
    dosing <- .pkpd_dose_data(4)
    expect_error(
        sim_tte_ode(model = "exponential", param = list(H0 = 0.1), n = 6,
            end = 20, delta = 2, data = dosing, seed = 1),
        "'n' \\(6\\) does not match.*'data'.*\\(4\\)"
    )
})
