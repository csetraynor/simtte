# Pre-Phase-B audit, Phase 9A: Weibull public-pipeline distribution
# validation. Uses the actual sim_tte() public workflow (not internal
# helpers) with a fixed seed and a cohort large enough that the
# empirical grid-crossing proportions and the analytical survival
# function agree within a statistically defensible tolerance (4 binomial
# standard errors, which is stable across machines/R versions and
# extremely unlikely to spuriously fail for a correct implementation).

binom_tol <- function(p, n, z = 4) {
    z * sqrt(p * (1 - p) / n)
}

check_distribution <- function(shape, mu = -0.5, n = 3000, seed = 20260826) {
    grid <- seq(0.1, 5, by = 0.1)
    set.seed(seed)
    lp <- matrix(rep(0, n), nrow = n)
    result <- sim_tte(pi = lp, mu = mu, coefs = shape, time = grid,
        type = "weibull", end_time = 5)

    eta <- exp(mu)
    analytic_S <- function(t) exp(-eta * t^shape)

    # 1. Empirical censoring proportion at end_time vs. analytical
    #    S(end_time).
    p_cens_analytic <- analytic_S(5)
    p_cens_empirical <- mean(result$sim_status == 0)
    expect_lt(abs(p_cens_empirical - p_cens_analytic),
        binom_tol(p_cens_analytic, n))

    # 2. Empirical P(T_grid <= t_j) vs. analytical 1 - S(t_j), at a few
    #    grid points spanning the follow-up.
    for (t_j in c(0.5, 1, 2, 3, 4)) {
        p_event_analytic <- 1 - analytic_S(t_j)
        p_event_empirical <- mean(result$sim_time <= t_j &
            result$sim_status == 1)
        expect_lt(abs(p_event_empirical - p_event_analytic),
            binom_tol(p_event_analytic, n),
            label = paste0("shape=", shape, " t_j=", t_j))
    }
}

test_that("Weibull public pipeline distribution matches theory: shape < 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_distribution(shape = 0.7)
})

test_that("Weibull public pipeline distribution matches theory: shape = 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_distribution(shape = 1)
})

test_that("Weibull public pipeline distribution matches theory: shape > 1", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    check_distribution(shape = 1.6)
})
