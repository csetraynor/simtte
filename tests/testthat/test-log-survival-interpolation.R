# Phase B: event_time_method = "log_survival" -- core algorithm,
# invariants, and regression tests. See PHASE_B_DESIGN_AUDIT.md for the
# full derivation of every edge case checked here.

# ---- Unit tests -----------------------------------------------------

test_that("simple two-point interpolation matches a hand-calculated answer", {
    # H_i = -log(exp(-1)) = 1, H_ip1 = -log(exp(-2)) = 2,
    # H_u = -log(exp(-1.5)) = 1.5 -> fraction = 0.5 -> t = 0.5
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = exp(-1), s_ip1 = exp(-2), u = exp(-1.5))
    expect_equal(got, 0.5, tolerance = 1e-12)
})

test_that("regular grid: interior crossing interpolates strictly between knots", {
    # Force a crossing strictly inside [1, 2]: S(1) = exp(-1) > u >
    # exp(-2) = S(2).
    got <- simtte:::.interpolate_log_survival(t_i = 1, t_ip1 = 2,
        s_i = exp(-1), s_ip1 = exp(-2), u = exp(-1.5))
    expect_equal(got, 1.5, tolerance = 1e-12)
})

test_that("irregular grid uses the actual spacing, not a uniform assumption", {
    got <- simtte:::.interpolate_log_survival(t_i = 0.2, t_ip1 = 5,
        s_i = exp(-1), s_ip1 = exp(-3), u = exp(-2))
    # H_i=1, H_ip1=3, H_u=2 -> fraction=0.5 -> t = 0.2 + 0.5*(5-0.2) = 2.6
    expect_equal(got, 2.6, tolerance = 1e-12)
})

test_that("exact crossing at a reported survival value is a direct grid hit, not an interpolated fraction", {
    # u exactly equal to a reported survival value: p11 = c(1, 0.5, 0.5, 0)
    # deterministic construction using flat-then-drop avoided; instead use
    # .get_tte() + .interpolate_log_survival() directly to confirm the
    # boundary fraction is exactly 0 or 1.
    got_left <- simtte:::.interpolate_log_survival(t_i = 1, t_ip1 = 2,
        s_i = 0.5, s_ip1 = 0.2, u = 0.5)
    expect_equal(got_left, 1, tolerance = 1e-9) # fraction == 0 -> t_i
    got_right <- simtte:::.interpolate_log_survival(t_i = 1, t_ip1 = 2,
        s_i = 0.5, s_ip1 = 0.2, u = 0.2)
    expect_equal(got_right, 2, tolerance = 1e-9) # fraction == 1 -> t_ip1
})

test_that("strictly interior crossing gives a time strictly between the two knots", {
    dat <- data.frame(ID = 1, time = c(1, 2), p11 = c(1, 0))
    # p11 = c(1, 0): deterministic event, crossing at index 2, i = 1
    # (S_i = 1 > 0, so H_i = 0); S_ip1 = 0 triggers the zero-survival rule
    # (covered separately below); use a strictly-interior example instead.
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 10,
        s_i = 0.9, s_ip1 = 0.1, u = 0.5)
    expect_true(got > 0 && got < 10)
})

test_that("no crossing censors at the subject's final reported time, both methods agree", {
    dat <- data.frame(ID = 1, time = seq(0.1, 5, by = 0.1), p11 = 1)
    r_grid <- sim_tte_df(dat, event_time_method = "grid")
    r_log <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_equal(r_grid$sim_status, 0)
    expect_equal(r_log$sim_status, 0)
    expect_equal(r_grid$sim_time, 5)
    expect_equal(r_log$sim_time, 5)
})

test_that("first-observation crossing is identical between grid and log_survival", {
    dat <- data.frame(ID = 1, time = c(5, 6, 7), p11 = c(0, 0, 0))
    r_grid <- sim_tte_df(dat, event_time_method = "grid")
    r_log <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_equal(r_grid$sim_time, 5)
    expect_equal(r_log$sim_time, 5)
    expect_equal(r_grid$sim_status, 1)
    expect_equal(r_log$sim_status, 1)
})

test_that("a flat survival segment elsewhere in the trajectory is never selected for interpolation", {
    # p11 = c(1, 1, 1, 0): the flat run of 1s (indices 1-3) can never
    # straddle any achievable U (1 <= U requires U >= 1, never true), so
    # the crossing is always found at index 4 (S = 0), with i = 3
    # pointing at the *last* 1 -- deterministic and free of any
    # zero-denominator from the flat run itself.
    dat <- data.frame(ID = 1, time = c(1, 2, 3, 4), p11 = c(1, 1, 1, 0))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_equal(r$sim_time, 4) # S[etime] == 0 rule
    expect_equal(r$sim_status, 1)
})

test_that("S[etime] == 0 returns t[etime] exactly, not the naive finite/Inf trap", {
    dat <- data.frame(ID = 1, time = c(1, 2, 3), p11 = c(1, 1, 0))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    # Deterministic: p11 = 1 never crosses (U < 1 always), so the
    # crossing is always found at index 3 (S = 0), i = 2 (S = 1).
    # A naive formula would give finite/Inf == 0 -> t_event = t[2] = 2.
    expect_equal(r$sim_time, 3)
    expect_false(isTRUE(all.equal(r$sim_time, 2)))
    expect_equal(r$sim_status, 1)
})

test_that("extremely small positive survival at the crossing endpoint is finite, no NaN/Inf", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = 0.5, s_ip1 = 1e-12, u = 1e-13)
    expect_true(is.finite(got))
    expect_true(got >= 0 && got <= 1)
})

test_that("U near zero gives a finite, in-range interpolated time", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = 0.5, s_ip1 = 1e-10, u = 1e-15)
    expect_true(is.finite(got))
    expect_true(got >= 0 && got <= 1)
})

test_that("U near one gives a finite, in-range interpolated time", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = 1 - 1e-10, s_ip1 = 0.1, u = 1 - 1e-9)
    expect_true(is.finite(got))
    expect_true(got >= 0 && got <= 1)
})

test_that("one-row trajectory resolves as a direct hit for both methods, never interpolated", {
    dat_event <- data.frame(ID = 1, time = 5, p11 = 0)
    r_grid <- sim_tte_df(dat_event, event_time_method = "grid")
    r_log <- sim_tte_df(dat_event, event_time_method = "log_survival")
    expect_identical(r_grid, r_log)

    dat_cens <- data.frame(ID = 1, time = 5, p11 = 1)
    r_grid2 <- sim_tte_df(dat_cens, event_time_method = "grid")
    r_log2 <- sim_tte_df(dat_cens, event_time_method = "log_survival")
    expect_identical(r_grid2, r_log2)
})

test_that("multiple subjects each interpolate independently", {
    dat <- data.frame(ID = rep(1:5, each = 20),
        time = rep(seq(0.5, 10, by = 0.5), 5),
        p11 = rep(exp(-0.4 * seq(0.5, 10, by = 0.5)), 5))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_equal(nrow(r), 5)
    expect_true(all(r$sim_time > 0))
})

test_that("character IDs work with log_survival", {
    dat <- data.frame(ID = rep(c("ctrl", "trt"), each = 20),
        time = rep(seq(0.5, 10, by = 0.5), 2),
        p11 = rep(exp(-0.4 * seq(0.5, 10, by = 0.5)), 2))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_setequal(r$ID, c("ctrl", "trt"))
})

test_that("non-contiguous numeric IDs work with log_survival", {
    dat <- data.frame(ID = rep(c(3, 17, 42), each = 20),
        time = rep(seq(0.5, 10, by = 0.5), 3),
        p11 = rep(exp(-0.4 * seq(0.5, 10, by = 0.5)), 3))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_setequal(r$ID, c(3, 17, 42))
})

# ---- Fraction/invariant checks ---------------------------------------

test_that("interpolation fraction is clamped to [0, 1] (defensive boundary check)", {
    # u exactly at s_i and exactly at s_ip1: fraction must be exactly 0/1,
    # not a hair outside due to log/exp round-trip.
    t0 <- simtte:::.interpolate_log_survival(t_i = 3, t_ip1 = 9,
        s_i = 0.37, s_ip1 = 0.02, u = 0.37)
    expect_equal(t0, 3)
    t1 <- simtte:::.interpolate_log_survival(t_i = 3, t_ip1 = 9,
        s_i = 0.37, s_ip1 = 0.02, u = 0.02)
    expect_equal(t1, 9)
})

test_that("interpolated event times never exceed the subject's final reported time", {
    set.seed(4242)
    dat <- data.frame(ID = rep(1:200, each = 30),
        time = rep(seq(0.1, 3, length.out = 30), 200),
        p11 = rep(exp(-1.2 * seq(0.1, 3, length.out = 30)), 200))
    r <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_true(all(r$sim_time <= 3 + 1e-9))
})

# ---- RNG invariant ----------------------------------------------------

test_that("exactly one runif(1) draw is consumed per subject, same as grid", {
    dat <- data.frame(ID = rep(1:10, each = 15),
        time = rep(seq(0.2, 3, length.out = 15), 10),
        p11 = rep(exp(-0.5 * seq(0.2, 3, length.out = 15)), 10))

    set.seed(99)
    u_before <- runif(1)
    set.seed(99)
    sim_tte_df(dat, event_time_method = "grid")
    u_after_grid <- runif(1)

    set.seed(99)
    sim_tte_df(dat, event_time_method = "log_survival")
    u_after_log <- runif(1)

    # If both methods consume exactly 10 draws (one per subject), the
    # *next* draw after each call must be identical.
    expect_equal(u_after_grid, u_after_log)
})

# ---- Regression tests --------------------------------------------------

test_that("default event_time_method is 'grid'", {
    dat <- data.frame(ID = rep(1:5, each = 20),
        time = rep(seq(0.5, 10, by = 0.5), 5),
        p11 = rep(exp(-0.3 * seq(0.5, 10, by = 0.5)), 5))
    set.seed(1)
    r_default <- sim_tte_df(dat)
    set.seed(1)
    r_grid <- sim_tte_df(dat, event_time_method = "grid")
    expect_identical(r_default, r_grid)
})

test_that("explicit 'grid' does not change any existing grid-method test outcome", {
    dat <- data.frame(ID = 1, time = 1:5, p11 = c(1, 1, 0, 0, 0))
    result <- sim_tte_df(dat, event_time_method = "grid")
    expect_equal(result$sim_time, 3)
    expect_equal(result$sim_status, 1)
})

test_that("event classification is identical between grid and log_survival for the same seed", {
    dat <- data.frame(ID = rep(1:100, each = 25),
        time = rep(seq(0.1, 5, length.out = 25), 100),
        p11 = rep(exp(-0.6 * seq(0.1, 5, length.out = 25)), 100))
    grid_times <- unique(dat$time)

    set.seed(777)
    r_grid <- sim_tte_df(dat, event_time_method = "grid")
    set.seed(777)
    r_log <- sim_tte_df(dat, event_time_method = "log_survival")

    # Same seed, same trajectories -> identical classification and ID
    # ordering for every subject (RNG invariant, §5/§6 of the audit).
    expect_identical(r_grid$sim_status, r_log$sim_status)
    expect_identical(r_grid$ID, r_log$ID)

    # sim_time may differ only where the grid time is *not* the first
    # grid point (a first-observation crossing) and status == 1 (a
    # strictly interior crossing); censoring (status == 0) always
    # matches exactly, since both methods censor at the same final time.
    is_censored <- r_grid$sim_status == 0
    expect_true(all(r_grid$sim_time[is_censored] == r_log$sim_time[is_censored]))

    # Wherever the two methods *do* differ, log_survival's time must be
    # strictly earlier than grid's (grid always reports the later,
    # right-hand endpoint of the crossing interval; interpolation gives
    # a point within that interval), and never before the earliest
    # reported time.
    differs <- r_grid$sim_time != r_log$sim_time
    if (any(differs)) {
        expect_true(all(r_log$sim_time[differs] < r_grid$sim_time[differs]))
        expect_true(all(r_log$sim_time[differs] >= min(grid_times)))
    }
})

test_that("event_time_method is validated with match.arg() semantics", {
    dat <- data.frame(ID = 1, time = 1, p11 = 1)
    expect_error(sim_tte_df(dat, event_time_method = "bogus"),
        "'arg' should be one of")
    # partial matching allowed
    expect_error(sim_tte_df(dat, event_time_method = "log"), NA)
})

# ---- Numerical tests ---------------------------------------------------

test_that("S = 1 at the crossing endpoints introduces no numerical issue", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1, s_i = 1,
        s_ip1 = 0.5, u = 0.7)
    expect_true(is.finite(got))
})

test_that("S close to 1 (post-normalization) is numerically safe", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1,
        s_i = 1 - 1e-9, s_ip1 = 0.3, u = 0.5)
    expect_true(is.finite(got))
    expect_true(got >= 0 && got <= 1)
})

test_that("S = 0 as the trajectory's final reported point is always an event, never censoring", {
    dat <- data.frame(ID = 1, time = c(1, 2, 3), p11 = c(1, 1, 0))
    r_grid <- sim_tte_df(dat, event_time_method = "grid")
    r_log <- sim_tte_df(dat, event_time_method = "log_survival")
    expect_equal(r_grid$sim_status, 1)
    expect_equal(r_log$sim_status, 1)
    expect_equal(r_grid$sim_time, 3)
    expect_equal(r_log$sim_time, 3)
})

test_that("very small time intervals do not introduce NaN/Inf", {
    got <- simtte:::.interpolate_log_survival(t_i = 1, t_ip1 = 1 + 1e-10,
        s_i = 0.6, s_ip1 = 0.5999999, u = 0.59999995)
    expect_true(is.finite(got))
    expect_true(got >= 1 && got <= 1 + 1e-10)
})

test_that("very large time intervals do not introduce NaN/Inf", {
    got <- simtte:::.interpolate_log_survival(t_i = 0, t_ip1 = 1e8,
        s_i = 0.9, s_ip1 = 0.1, u = 0.5)
    expect_true(is.finite(got))
    expect_true(got >= 0 && got <= 1e8)
})

test_that("extreme Weibull parameters through the public sim_tte() pipeline produce no NaN/Inf event times", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    # mu = 710: previously caused a NaN survival value at t = 0 before
    # the pre-Phase-B Weibull numerical fix; confirms interpolation
    # inherits that fix cleanly.
    result <- sim_tte(pi = 0, mu = 710, coefs = 1, time = c(0, 1),
        event_time_method = "log_survival")
    expect_true(is.finite(result$sim_time))

    result2 <- sim_tte(pi = 700, mu = 0, coefs = 1, time = c(0, 1),
        event_time_method = "log_survival")
    expect_true(is.finite(result2$sim_time))

    result3 <- sim_tte(pi = -700, mu = 0, coefs = 1,
        time = seq(0, 1000, by = 10), event_time_method = "log_survival")
    expect_true(is.finite(result3$sim_time))

    result4 <- sim_tte(pi = 0, mu = -1, coefs = 1e-4,
        time = seq(0, 100, by = 1), event_time_method = "log_survival")
    expect_true(is.finite(result4$sim_time))
})
