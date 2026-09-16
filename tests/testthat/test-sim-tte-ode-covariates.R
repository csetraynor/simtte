# Phase 2: sim_tte_ode() time-varying named covariates (`covariates`/
# `beta`), generalizing sim_tte()'s single-covariate `lp_data` mechanism
# (Phase F/G) via the shared helpers in R/helpers.R
# (.canonicalize_lp_data()/.validate_lp_data_trajectories()/
# .check_lp_data_coverage(), now parameterized by `value_cols` instead
# of hardcoding "lp"; .build_ode_covariate_rows() is the new,
# sim_tte_ode()-specific combination step). See
# reports/06_phase2_report.md for the full account.
#
# Categories mirrored from test-time-varying-lp.R (Phase G), adapted to
# sim_tte_ode()'s shape (a named covariate set + beta, not a single lp
# column):
#   1. covariates = NULL is the unchanged default
#   2. population-level covariates (no ID column)
#   3. subject-specific covariates (ID column)
#   4. analytical agreement for a piecewise-constant multi-covariate lp(t)
#   5. multiple covariates combine additively, as documented
#   6. validation of malformed covariates/beta
#   7. sim_tte()'s lp_data behavior is unchanged by the shared-helper
#      refactor (explicit old-vs-new comparison, not only re-running the
#      pre-existing Phase G suite)
#   8. covariates/beta cannot leak through '...'

# ---- 1. covariates = NULL is the unchanged default ----

test_that("omitting covariates and passing covariates = NULL, beta = NULL give identical() output", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    s1 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 50, end = 20, delta = 2, seed = 11)
    s2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 50, end = 20, delta = 2, seed = 11, covariates = NULL,
        beta = NULL)
    expect_identical(s1$events, s2$events)
})

# ---- 2. Population-level covariates (no ID column) ----

test_that("population-level covariates give identical trajectories across subjects with equal idata", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    covdat <- data.frame(time = c(0, 5, 10), X1 = c(0, 1, -1))
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 3, end = 15, delta = 1, covariates = covdat, beta = c(X1 = 1),
        idata = data.frame(ID = 1:3, U = c(0.3, 0.3, 0.3)),
        keep_trajectory = TRUE, seed = 1)
    traj <- sim$trajectory
    p1 <- traj$p11[traj$ID == 1]
    p2 <- traj$p11[traj$ID == 2]
    p3 <- traj$p11[traj$ID == 3]
    expect_equal(p1, p2, tolerance = 1e-10)
    expect_equal(p1, p3, tolerance = 1e-10)
})

# ---- 3. Subject-specific covariates (ID column) ----

test_that("subject-specific covariates produce distinct trajectories per subject", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    covdat <- data.frame(ID = c(1, 1, 2, 2), time = c(0, 5, 0, 5),
        X1 = c(0, 2, 0, -2))
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 2, end = 15, delta = 1, covariates = covdat, beta = c(X1 = 1),
        idata = data.frame(ID = 1:2, U = c(0.3, 0.3)),
        keep_trajectory = TRUE, seed = 1)
    traj <- sim$trajectory
    p1 <- traj$p11[traj$ID == 1 & traj$time == 10]
    p2 <- traj$p11[traj$ID == 2 & traj$time == 10]
    expect_false(isTRUE(all.equal(p1, p2)))
    # Subject 1's X1 jumps up (higher hazard, beta = 1) -> lower survival
    # than subject 2's, whose X1 jumps down.
    expect_true(p1 < p2)
})

# ---- 4. Analytical agreement: piecewise-constant multi-covariate lp(t) ----
# Fast/slow split (interlude session, reports/07_test_runbook.md).

check_multi_covariate_agreement <- function(n) {
    H0 <- 0.2
    # lp(t) = 1*X1(t) - 2*X2(t): 0 on [0, 3), 1 on [3, 6), -1.5 on [6, ...)
    covdat <- data.frame(time = c(0, 3, 6), X1 = c(0, 1, 1), X2 = c(0, 0, 1))
    beta <- c(X1 = 1, X2 = -2)
    end <- 10
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = H0), n = n,
        end = end, delta = 1, covariates = covdat, beta = beta, seed = 42)

    # H(t) = integral of H0 * exp(lp(s)) ds over piecewise-constant lp
    H_at <- function(t) {
        segs <- c(0, 3, 6, Inf)
        lps <- c(0, 1, -1.5)
        total <- 0
        for (j in seq_along(lps)) {
            lo <- segs[j]; hi <- min(segs[j + 1], t)
            if (hi > lo) total <- total + H0 * exp(lps[j]) * (hi - lo)
        }
        total
    }
    for (t_j in c(2, 5, 9)) {
        p_a <- 1 - exp(-H_at(t_j))
        p_e <- mean(sim$events$sim_time <= t_j & sim$events$sim_status == 1)
        expect_lt(abs(p_a - p_e), binom_tol(p_a, n), label = paste0("t=", t_j))
    }
}

test_that("multi-covariate lp(t) matches a hand-derived piecewise-constant cumulative hazard (exponential model) [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_multi_covariate_agreement(n = 600)
})
test_that("multi-covariate lp(t) matches a hand-derived piecewise-constant cumulative hazard (exponential model) [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_multi_covariate_agreement(n = 4000)
})

# ---- 5. Multiple covariates combine additively ----

test_that(".build_ode_covariate_rows() combines multiple covariates as sum(beta * X)", {
    covdat <- data.frame(time = c(0, 1), X1 = c(1, 2), X2 = c(3, 4))
    rows <- simtte:::.build_ode_covariate_rows(covdat, beta = c(X1 = 2, X2 = -1),
        n_subjects = 1, end = 1)
    expect_equal(rows$lp, c(2 * 1 - 1 * 3, 2 * 2 - 1 * 4))
    expect_equal(rows$cmt, c(1, 1))
    expect_equal(rows$evid, c(1, 1))
    expect_equal(rows$amt, c(0, 0))
})

test_that(".build_ode_covariate_rows()'s cmt is a parameter, not hardcoded", {
    covdat <- data.frame(time = 0, X1 = 1)
    rows <- simtte:::.build_ode_covariate_rows(covdat, beta = c(X1 = 1),
        n_subjects = 1, end = 1, cmt = 3L)
    expect_equal(rows$cmt, 3L)
})

# ---- 6. Validation of malformed covariates/beta ----

test_that("covariates and beta must be supplied together", {
    expect_error(
        sim_tte_ode(model = "exponential", n = 2, end = 10,
            covariates = data.frame(time = 0, X1 = 0)),
        "supplied together"
    )
    expect_error(
        sim_tte_ode(model = "exponential", n = 2, end = 10,
            beta = c(X1 = 1)),
        "supplied together"
    )
})

test_that("beta names not present in covariates are rejected", {
    covdat <- data.frame(time = 0, X1 = 0)
    expect_error(
        simtte:::.build_ode_covariate_rows(covdat, beta = c(X2 = 1),
            n_subjects = 1, end = 1),
        "missing column"
    )
})

test_that("covariates$ID not equal to 1:n_subjects is rejected", {
    covdat <- data.frame(ID = c(1, 3), time = c(0, 0), X1 = c(0, 0))
    expect_error(
        simtte:::.build_ode_covariate_rows(covdat, beta = c(X1 = 1),
            n_subjects = 2, end = 1),
        "1:2"
    )
})

test_that("malformed covariates propagate the underlying .canonicalize_lp_data()/.validate_lp_data_trajectories() errors", {
    unsorted <- data.frame(time = c(1, 0), X1 = c(0, 1))
    expect_error(
        simtte:::.build_ode_covariate_rows(unsorted, beta = c(X1 = 1),
            n_subjects = 1, end = 1),
        "sorted in ascending order"
    )
    no_zero <- data.frame(time = c(0.5, 1), X1 = c(0, 1))
    expect_error(
        simtte:::.build_ode_covariate_rows(no_zero, beta = c(X1 = 1),
            n_subjects = 1, end = 1),
        "time = 0"
    )
})

# ---- 7. sim_tte()'s lp_data behavior is unchanged by the shared-helper refactor ----

test_that("the generalized .canonicalize_lp_data()/.validate_lp_data_trajectories() reproduce the pre-Phase-2 (single-lp) implementations exactly", {
    # The pre-Phase-2 bodies, reproduced verbatim (see PHASE_G_REPORT.md
    # and the source history of R/helpers.R before this session), kept
    # here only as a comparison oracle -- not used by any package code.
    old_canonicalize <- function(lp_data, n_subjects) {
        lp_data <- as.data.frame(lp_data)
        out <- if ("ID" %in% names(lp_data)) {
            data.frame(ID = lp_data$ID, time = lp_data$time, lp = lp_data$lp)
        } else {
            do.call(rbind, lapply(seq_len(n_subjects), function(i) {
                data.frame(ID = i, time = lp_data$time, lp = lp_data$lp)
            }))
        }
        out
    }
    old_validate <- function(lp_canonical) {
        by_id <- split(lp_canonical, lp_canonical$ID)
        for (id in names(by_id)) {
            times_i <- by_id[[id]]$time
            lp_i <- by_id[[id]]$lp
            dup <- anyDuplicated(times_i)
            if (dup) {
                dup_time <- times_i[dup]
                vals <- lp_i[times_i == dup_time]
                if (length(unique(vals)) > 1L) return(FALSE)
            }
        }
        TRUE
    }

    scenarios <- list(
        list(df = data.frame(time = c(0, 1, 2), lp = c(0, 0.5, -0.5)), n = 3),
        list(df = data.frame(ID = c(1, 1, 2, 2), time = c(0, 1, 0, 1),
            lp = c(0, 1, 0, -1)), n = 2),
        list(df = data.frame(time = c(0, 1, 1), lp = c(0, 1, 1)), n = 1)
    )
    for (s in scenarios) {
        old_out <- old_canonicalize(s$df, s$n)
        new_out <- simtte:::.canonicalize_lp_data(s$df, s$n)
        expect_identical(old_out, new_out)
        expect_identical(old_validate(old_out),
            isTRUE(tryCatch(simtte:::.validate_lp_data_trajectories(new_out),
                error = function(e) FALSE)))
    }
})

# ---- 8. covariates/beta cannot leak through '...' ----

test_that("covariates/beta are formal arguments and cannot be supplied a second time via '...'", {
    covdat <- data.frame(time = 0, X1 = 0)
    expect_error(
        do.call(sim_tte_ode, list(model = "exponential", n = 1, end = 10,
            covariates = covdat, beta = c(X1 = 1),
            covariates = covdat)),
        "formal argument"
    )
})

# ---- 9. covariates/beta combined with a caller-supplied dosing `data`
# (simttepower feedback 1: lp/data merge,
# reports/04_author_decisions.md; reports/29_simttepower_feedback_lp_merge.md) ----
#
# Before the fix, plain dplyr::bind_rows(data, cov_rows) left every
# dosing row's `lp` NA (only cov_rows has that column), which is not
# merely an mrgsolve::valid_data_set() warning: mrgsolve does not carry
# a NA $PARAM value forward the way it carries a real one, so once `lp`
# went NA at a data-set record, the p11 ODE state itself became (and
# permanently stayed) NaN for the rest of that subject's trajectory --
# silently turning a real later event into administrative censoring.
# See reports/experiments/29_lp_na_mechanism.R for the direct
# mrgsolve-level demonstration. .merge_ode_covariate_rows() now fills
# every dosing row's `lp` by LOCF against that subject's own covariate
# trajectory, and sorts the merged frame by ID/time (also fixing a
# separate failure mode: mrgsolve erroring "the data set is not sorted
# by time", or silently splitting a subject's rows into two disjoint
# blocks, when the unsorted merge put a dosing row after a
# time-earlier covariate row in the raw bind_rows() order).

# Time-varying covariates (age constant, trt_on switches on at t = 8)
# combined with a 4-dose repeat regimen (t = 0, 4, 12, 16) -- dosing
# rows straddle the covariate update, exercising the LOCF fill on both
# sides of it.
.lp_merge_fixture <- function(n, seed = 1) {
    set.seed(seed)
    age <- stats::runif(n, 40, 80)
    covdat <- data.frame(ID = rep(seq_len(n), each = 2),
        time = rep(c(0, 8), n), age = rep(age, each = 2),
        trt_on = rep(c(0, 1), n))
    dosing <- data.frame(ID = rep(seq_len(n), each = 4),
        time = rep(c(0, 4, 12, 16), n), evid = 1L, amt = 100, cmt = 1L)
    list(covdat = covdat, dosing = dosing, beta = c(age = 0.01, trt_on = 0.3))
}

check_lp_merge_matches_hand_built_oracle <- function(n) {
    fx <- .lp_merge_fixture(n)
    param <- .PKPD_TEST_DEFAULT_PARAM$pk_hazard

    actual <- expect_no_warning(
        sim_tte_ode(model = "pk_hazard", param = param, n = n, end = 20,
            delta = 1, data = fx$dosing, covariates = fx$covdat,
            beta = fx$beta, keep_trajectory = TRUE, seed = 11))

    # Independent oracle: pre-fill `lp` on the dosing rows by hand
    # (LOCF against the covariate-update times/lp values, not via
    # .merge_ode_covariate_rows()) and pass the merged, sorted data
    # straight through `data` with `covariates = NULL` -- a code path
    # the merge fix never touches.
    cov_rows <- simtte:::.build_ode_covariate_rows(fx$covdat, fx$beta,
        n_subjects = n, end = 20)
    dosing_lp <- fx$dosing
    dosing_lp$lp <- vapply(seq_len(nrow(dosing_lp)), function(i) {
        traj <- cov_rows[cov_rows$ID == dosing_lp$ID[i], ]
        traj$lp[max(which(traj$time <= dosing_lp$time[i]))]
    }, numeric(1))
    oracle_data <- rbind(dosing_lp, cov_rows)
    oracle_data <- oracle_data[order(oracle_data$ID, oracle_data$time), ]
    oracle <- sim_tte_ode(model = "pk_hazard", param = param, n = n,
        end = 20, delta = 1, data = oracle_data, keep_trajectory = TRUE,
        seed = 11)

    expect_equal(actual$trajectory$HAZ, oracle$trajectory$HAZ,
        tolerance = 1e-10)
    expect_equal(actual$events, oracle$events)
}

test_that("covariates + beta + custom dosing data: trajectory matches a hand-built LOCF oracle, zero warnings [fast]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_lp_merge_matches_hand_built_oracle(n = 5)
})
test_that("covariates + beta + custom dosing data: trajectory matches a hand-built LOCF oracle, zero warnings [slow]", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_slow()
    check_lp_merge_matches_hand_built_oracle(n = 60)
})

test_that("baseline (single-row-per-subject) covariates + custom dosing data also produce zero warnings", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    n <- 4
    covdat <- data.frame(ID = 1:n, time = 0, age = c(40, 50, 60, 70))
    dosing <- .pkpd_dose_data(n)
    expect_no_warning(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = n, end = 10,
            delta = 1, data = dosing, covariates = covdat,
            beta = c(age = 0.01), seed = 3))
})

test_that("a 'data' column named 'lp' conflicting with covariates/beta errors clearly", {
    dosing <- .pkpd_dose_data(2)
    dosing$lp <- 0
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 2, end = 10,
            delta = 1, data = dosing,
            covariates = data.frame(ID = 1:2, time = 0, age = c(50, 60)),
            beta = c(age = 0.01), seed = 1),
        "already has an 'lp' column"
    )
})

test_that("'data' with an unlabeled (no-ID) dosing row errors, regardless of whether covariates vary by subject", {
    dosing <- data.frame(time = 0, cmt = 1, amt = 100, evid = 1) # no ID
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 2, end = 10,
            delta = 1, data = dosing,
            covariates = data.frame(ID = 1:2, time = 0, age = c(50, 90)),
            beta = c(age = 0.01), seed = 1),
        "must have an 'ID' column"
    )
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 3, end = 10,
            delta = 1, data = dosing,
            covariates = data.frame(time = 0, age = 60),
            beta = c(age = 0.01), seed = 1),
        "must have an 'ID' column"
    )
})

test_that("'data' with an ID not present in covariates errors", {
    dosing <- data.frame(ID = 3, time = 0, cmt = 1, amt = 100, evid = 1)
    expect_error(
        sim_tte_ode(model = "pk_hazard",
            param = .PKPD_TEST_DEFAULT_PARAM$pk_hazard, n = 2, end = 10,
            delta = 1, data = dosing,
            covariates = data.frame(ID = 1:2, time = 0, age = c(50, 60)),
            beta = c(age = 0.01), seed = 1),
        "not present in 'covariates'"
    )
})
