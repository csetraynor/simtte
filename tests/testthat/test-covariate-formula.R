# "Before Phase 7: covariate interface" (reports/04_author_decisions.md;
# reports/25_covariate_interface_review.md; reports/26_covariate_interface_report.md):
# the formula-based generalization of sim_tte_ode()'s covariates/beta,
# additive to the Phase 2 mechanism tested in
# test-sim-tte-ode-covariates.R (unchanged, still passing, not
# duplicated here).
#
# Categories:
#   1. formula = NULL is the unchanged legacy default
#   2. the intercept rule (drop + message, vs. ~ 0 + ... silence)
#   3. baseline vs. time-varying data shapes
#   4. coefficient-name validation (two-directional)
#   5. missing-value handling (na.action = na.fail, not silent na.omit)
#   6. statistical recovery: baseline covariate, factor arm, transformed term
#   7. the fixed-EC50 Emax boundary: formula vs. tte_model(hazard = ...)
#   8. formula requires covariates/beta

# ---- 1. formula = NULL is the unchanged legacy default ----

test_that("formula = NULL and formula omitted give identical() output", {
    cov <- data.frame(time = c(0, 10), x = c(0, 1))
    a <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05), n = 10,
        end = 20, delta = 2, covariates = cov, beta = c(x = 0.8), seed = 5)
    b <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05), n = 10,
        end = 20, delta = 2, covariates = cov, beta = c(x = 0.8),
        formula = NULL, seed = 5)
    expect_identical(a$events, b$events)
})

test_that("the identity formula (~ 0 + x) reproduces the legacy covariates/beta path exactly", {
    cov <- data.frame(time = c(0, 5, 15), x = c(0, 1, 0.5))
    legacy <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 25, end = 20, delta = 1, covariates = cov, beta = c(x = 0.8),
        seed = 11)
    formula_form <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
        n = 25, end = 20, delta = 1, covariates = cov, formula = ~ 0 + x,
        beta = c(x = 0.8), seed = 11)
    expect_identical(legacy$events, formula_form$events)
})

test_that("a multi-covariate identity formula reproduces the additive-combination test's result", {
    covdat <- data.frame(time = c(0, 1), X1 = c(1, 2), X2 = c(3, 4))
    beta <- c(X1 = 0.5, X2 = -0.2)
    legacy_rows <- .build_ode_covariate_rows(covdat, beta, n_subjects = 2, end = 5)
    formula_rows <- .build_ode_covariate_rows_formula(~ 0 + X1 + X2, covdat,
        beta, n_subjects = 2, end = 5)
    expect_equal(legacy_rows, formula_rows)
})

# ---- 2. The intercept rule ----

test_that("a formula with an implied intercept drops it and messages", {
    d <- data.frame(age = c(40, 50, 60, 70))
    expect_message(X <- .build_lp_design_matrix(~ age, d),
        "dropped the intercept column")
    expect_identical(colnames(X), "age")
})

test_that("~ 0 + ... suppresses the intercept message and keeps the column set", {
    d <- data.frame(age = c(40, 50, 60, 70))
    expect_no_message(X <- .build_lp_design_matrix(~ 0 + age, d))
    expect_identical(colnames(X), "age")
})

test_that("a factor gets k-1 treatment-contrast columns with an intercept present, k with ~ 0 +", {
    d <- data.frame(arm = factor(c("A", "B", "A", "B")))
    expect_message(X1 <- .build_lp_design_matrix(~ arm, d),
        "dropped the intercept column")
    expect_identical(colnames(X1), "armB")
    # ~ 0 + arm gives k = 2 columns instead of the usual k - 1 = 1 --
    # its own, different message naming the factor (decision, "After
    # the covariate interface"), not silence.
    expect_message(X2 <- .build_lp_design_matrix(~ 0 + arm, d),
        "'arm'.*expanded to one design-matrix column per level")
    expect_identical(colnames(X2), c("armA", "armB"))
})

test_that("with two factors and no intercept, only the first gets full-level coding, and only it messages", {
    d <- data.frame(arm = factor(c("A", "B", "A", "B")),
        site = factor(c("X", "Y", "X", "Y")))
    msgs <- character(0)
    X <- withCallingHandlers(
        .build_lp_design_matrix(~ 0 + arm + site, d),
        message = function(m) {
            msgs[[length(msgs) + 1L]] <<- conditionMessage(m)
            invokeRestart("muffleMessage")
        })
    expect_identical(colnames(X), c("armA", "armB", "siteY"))
    expect_length(msgs, 1L)
    expect_match(msgs[[1]], "'arm'")
})

test_that("~ 0 + age (no factor at all) is silent", {
    d <- data.frame(age = c(1, 2, 3))
    expect_no_message(X <- .build_lp_design_matrix(~ 0 + age, d))
    expect_identical(colnames(X), "age")
})

test_that("a formula reducing to an intercept-only design matrix errors", {
    d <- data.frame(age = c(1, 2, 3))
    expect_error(.build_lp_design_matrix(~ 1, d), "empty design matrix")
})

# ---- 3. Baseline vs. time-varying data shapes ----

test_that("baseline data (no time column) is one row per subject, not a shared trajectory", {
    baseline <- data.frame(age = c(30, 40, 50, 60))
    rows <- suppressMessages(.build_ode_covariate_rows_formula(~ age, baseline,
        c(age = 0.1), n_subjects = 4, end = 10))
    expect_identical(sort(rows$ID), 1:4)
    expect_true(all(rows$time == 0))
    # one row per subject, each carrying that subject's own age -- not
    # every row of `baseline` recycled onto every subject.
    expect_equal(sort(rows$lp), sort(0.1 * baseline$age))
})

test_that("baseline data with an explicit ID column is honored", {
    baseline <- data.frame(ID = 4:1, age = c(60, 50, 40, 30))
    rows <- suppressMessages(.build_ode_covariate_rows_formula(~ age, baseline,
        c(age = 0.1), n_subjects = 4, end = 10))
    expect_equal(rows$lp[rows$ID == 4], 0.1 * 60)
    expect_equal(rows$lp[rows$ID == 1], 0.1 * 30)
})

test_that("baseline data with the wrong row count errors via the existing ID-coverage check", {
    baseline <- data.frame(age = c(30, 40, 50))
    expect_error(
        suppressMessages(.build_ode_covariate_rows_formula(~ age, baseline,
            c(age = 0.1), n_subjects = 4, end = 10)),
        "exactly one trajectory per")
})

test_that("time-varying data (has a time column) keeps the legacy no-ID recycle-to-everyone rule", {
    cov <- data.frame(time = c(0, 10), x = c(0, 1))
    rows <- .build_ode_covariate_rows_formula(~ 0 + x, cov, c(x = 1),
        n_subjects = 3, end = 20)
    expect_identical(sort(unique(rows$ID)), 1:3)
    expect_equal(nrow(rows), 6) # 2 time points x 3 subjects
})

# ---- 3b. Baseline single-row recycling and n inference (decision,
# reports/04_author_decisions.md "After the covariate interface") ----

test_that("a 1-row baseline frame is a population-level constant, recycled to n subjects", {
    baseline <- data.frame(age = 55)
    rows <- suppressMessages(.build_ode_covariate_rows_formula(~ age,
        baseline, c(age = 0.1), n_subjects = 5, end = 10))
    expect_identical(sort(rows$ID), 1:5)
    expect_true(all(rows$lp == 0.1 * 55))
})

test_that("a >1-row baseline frame infers n when n is left at its default", {
    baseline <- data.frame(age = c(20, 30, 40))
    sim <- suppressMessages(sim_tte_ode(model = "weibull",
        param = list(mu = -1, shape = 1.5), end = 10, delta = 1,
        covariates = baseline, formula = ~ age, beta = c(age = 0.01),
        seed = 1))
    expect_equal(sort(unique(sim$events$ID)), 1:3)
})

test_that("a >1-row baseline frame accepts an explicit, matching n", {
    baseline <- data.frame(age = c(20, 30, 40))
    sim <- suppressMessages(sim_tte_ode(model = "weibull",
        param = list(mu = -1, shape = 1.5), n = 3, end = 10, delta = 1,
        covariates = baseline, formula = ~ age, beta = c(age = 0.01),
        seed = 1))
    expect_equal(sort(unique(sim$events$ID)), 1:3)
})

test_that("a >1-row baseline frame rejects a mismatched, explicit n", {
    baseline <- data.frame(age = c(20, 30, 40))
    expect_error(
        suppressMessages(sim_tte_ode(model = "weibull",
            param = list(mu = -1, shape = 1.5), n = 5, end = 10, delta = 1,
            covariates = baseline, formula = ~ age, beta = c(age = 0.01),
            seed = 1)),
        "does not match the number of rows")
})

test_that("n alone, with no covariates at all, is untouched by the n-resolution logic", {
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05), n = 4,
        end = 10, delta = 1, seed = 1)
    expect_equal(sort(unique(sim$events$ID)), 1:4)
})

# ---- 4. Coefficient-name validation ----

test_that("a beta missing a required name errors, naming the expected columns", {
    d <- data.frame(age = c(1, 2, 3), arm = factor(c("A", "B", "A")))
    X <- suppressMessages(.build_lp_design_matrix(~ age + arm, d))
    expect_error(.validate_beta_matches_design(c(age = 1), X),
        "Missing from 'beta': armB")
})

test_that("a beta with an extra, unmatched name errors, naming it", {
    d <- data.frame(age = c(1, 2, 3))
    X <- suppressMessages(.build_lp_design_matrix(~ age, d))
    expect_error(
        .validate_beta_matches_design(c(age = 1, bogus = 2), X),
        "Not a column of the design matrix: bogus")
})

test_that("an unnamed beta errors", {
    d <- data.frame(age = c(1, 2, 3))
    X <- suppressMessages(.build_lp_design_matrix(~ age, d))
    expect_error(.validate_beta_matches_design(c(1), X),
        "must be a named numeric vector")
})

# ---- 5. Missing-value handling ----

test_that("NA in a covariate errors instead of silently dropping the row", {
    d <- data.frame(age = c(1, NA, 3))
    expect_error(.build_lp_design_matrix(~ age, d), "Could not build")
})

test_that("'formula' must actually be a formula", {
    d <- data.frame(age = c(1, 2, 3))
    expect_error(.build_lp_design_matrix("age", d), "must be a formula")
})

# ---- 6. Statistical recovery ----

check_baseline_cox_recovery <- function(n, beta_true, seed) {
    set.seed(seed)
    baseline <- data.frame(x = rnorm(n))
    sim <- suppressMessages(sim_tte_ode(model = "exponential",
        param = list(H0 = 0.05), n = n, end = 50, delta = 5,
        covariates = baseline, formula = ~ 0 + x, beta = c(x = beta_true),
        seed = seed))
    # $events never echoes covariate values back (?sim_tte_ode "Between-
    # subject variability" makes the same point for omega) -- rejoin by
    # ID, the way the vignette's own "irm1-idata" section does.
    dat <- sim$events
    dat$x <- baseline$x[dat$ID]
    fit <- survival::coxph(survival::Surv(sim_time, sim_status) ~ x, data = dat)
    co <- summary(fit)$coefficients
    expect_lt(abs(co["x", "coef"] - beta_true), 4 * co["x", "se(coef)"])
}

test_that("a baseline covariate's Cox coefficient recovers beta_x [fast]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    check_baseline_cox_recovery(n = 1000, beta_true = 0.6, seed = 101)
})
test_that("a baseline covariate's Cox coefficient recovers beta_x [slow]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    skip_if_not_slow()
    check_baseline_cox_recovery(n = 6000, beta_true = 0.6, seed = 102)
})

check_factor_arm_recovery <- function(n, beta_true, seed) {
    set.seed(seed)
    baseline <- data.frame(arm = factor(rep(c("control", "treated"), n / 2)))
    sim <- suppressMessages(sim_tte_ode(model = "exponential",
        param = list(H0 = 0.05), n = n, end = 50, delta = 5,
        covariates = baseline, formula = ~ arm,
        beta = c(armtreated = beta_true), seed = seed))
    dat <- sim$events
    dat$arm <- baseline$arm[dat$ID]
    fit <- survival::coxph(survival::Surv(sim_time, sim_status) ~ arm, data = dat)
    co <- summary(fit)$coefficients
    expect_lt(abs(co["armtreated", "coef"] - beta_true),
        4 * co["armtreated", "se(coef)"])
}

test_that("a factor arm's Cox coefficient recovers the between-arm log-HR [fast]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    check_factor_arm_recovery(n = 1000, beta_true = -0.7, seed = 201)
})
test_that("a factor arm's Cox coefficient recovers the between-arm log-HR [slow]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    skip_if_not_slow()
    check_factor_arm_recovery(n = 6000, beta_true = -0.7, seed = 202)
})

check_transformed_term_recovery <- function(n, beta_true, seed) {
    set.seed(seed)
    baseline <- data.frame(x = rnorm(n))
    sim <- suppressMessages(sim_tte_ode(model = "exponential",
        param = list(H0 = 0.05), n = n, end = 50, delta = 5,
        covariates = baseline, formula = ~ 0 + I(x^2),
        beta = c("I(x^2)" = beta_true), seed = seed))
    dat <- sim$events
    dat$x2 <- baseline$x[dat$ID]^2
    fit <- survival::coxph(survival::Surv(sim_time, sim_status) ~ x2, data = dat)
    co <- summary(fit)$coefficients
    expect_lt(abs(co["x2", "coef"] - beta_true), 4 * co["x2", "se(coef)"])
}

test_that("a transformed term (I(x^2)) recovers its coefficient [fast]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    check_transformed_term_recovery(n = 1000, beta_true = 0.5, seed = 301)
})
test_that("a transformed term (I(x^2)) recovers its coefficient [slow]", {
    skip_on_cran()
    skip_if_not_installed("survival")
    skip_if_not_slow()
    check_transformed_term_recovery(n = 6000, beta_true = 0.5, seed = 302)
})

# ---- 7. The fixed-EC50 Emax boundary: formula vs. tte_model(hazard = ...) ----

test_that("a fixed-EC50 Emax term expressed via formula matches the same link written in tte_model(hazard = ...)", {
    skip_if_not_installed("mrgsolve")
    n <- 20
    set.seed(41)
    dose <- runif(n, 10, 200)
    EC50 <- 50
    U <- runif(n, 0, 0.999)

    # Route A: a fixed, known EC50 folded into a covariate column,
    # supplied via formula/beta on an ordinary library model.
    logm <- log(dose / (EC50 + dose))
    idata_a <- data.frame(ID = seq_len(n), U = U)
    sim_a <- suppressMessages(sim_tte_ode(model = "exponential",
        param = list(H0 = 0.05), idata = idata_a, end = 20, delta = 1,
        covariates = data.frame(logm = logm), formula = ~ 0 + logm,
        beta = c(logm = 1)))

    # Route B: the same fixed EC50 written directly into a model's own
    # HAZ expression via tte_model() -- the "out of scope for formula"
    # side of the documented boundary, here with EC50 held fixed so the
    # two routes can be compared directly.
    model_code <- paste("$PARAM H0 = 0.05, EC50 = 50, DOSE = 0",
        "$CMT DUMMY", "$ODE", "dxdt_DUMMY = 0;", sep = "\n")
    mymod <- tte_model(model = model_code,
        hazard = "H0 * (DOSE / (EC50 + DOSE))")
    idata_b <- data.frame(ID = seq_len(n), U = U, DOSE = dose)
    sim_b <- sim_tte_ode(model = mymod, idata = idata_b, end = 20, delta = 1)

    expect_equal(sim_a$events$sim_time, sim_b$events$sim_time,
        tolerance = 1e-6)
    expect_identical(sim_a$events$sim_status, sim_b$events$sim_status)
})

# ---- 8. formula requires covariates/beta ----

test_that("'formula' without 'covariates'/'beta' errors", {
    expect_error(
        sim_tte_ode(model = "exponential", param = list(H0 = 0.05), n = 5,
            end = 10, formula = ~ age),
        "requires 'covariates' and 'beta'")
})
