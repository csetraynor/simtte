# Phase A4: event-selection tests. Grid-based crossing detection only;
# interpolation between grid points is out of scope for Phase A (see
# NEWS/Phase A implementation report).

test_that("exact crossing at a reported grid point is returned deterministically", {
    # p11 drops to 0 exactly at time = 3; 0 <= U for any U in [0, 1),
    # so the event is always resolved at time = 3, never earlier.
    dat <- data.frame(ID = 1, time = 1:5, p11 = c(1, 1, 0, 0, 0))
    result <- sim_tte_df(dat)
    expect_equal(result$sim_time, 3)
    expect_equal(result$sim_status, 1)
})

test_that("crossing between reported points resolves to the grid point, not an interpolated time", {
    set.seed(99)
    grid <- seq(0.5, 10, by = 0.5)
    dat <- data.frame(ID = rep(1:20, each = length(grid)),
        time = rep(grid, 20),
        p11 = rep(exp(-0.4 * grid), 20))
    result <- sim_tte_df(dat)
    events <- result[result$sim_status == 1, ]
    expect_true(nrow(events) > 0)
    expect_true(all(events$sim_time %in% grid))
})

test_that("no crossing results in censoring at the subject's last reported time", {
    dat <- data.frame(ID = 1, time = seq(0.1, 5, by = 0.1), p11 = 1)
    result <- sim_tte_df(dat)
    expect_equal(result$sim_status, 0)
    expect_equal(result$sim_time, 5)
})

test_that("subject-specific final observation times are respected independently", {
    dat <- rbind(
        data.frame(ID = "a", time = seq(1, 3, by = 1), p11 = 1),
        data.frame(ID = "b", time = seq(1, 9, by = 1), p11 = 1)
    )
    result <- sim_tte_df(dat)
    expect_equal(result$sim_time[result$ID == "a"], 3)
    expect_equal(result$sim_time[result$ID == "b"], 9)
})

test_that("multiple subjects are all present exactly once in the output", {
    dat <- data.frame(ID = rep(1:8, each = 6),
        time = rep(1:6, 8),
        p11 = rep(seq(1, 0, length.out = 6), 8))
    result <- sim_tte_df(dat)
    expect_equal(sort(result$ID), 1:8)
    expect_equal(nrow(result), 8)
})

test_that("nonconsecutive numeric IDs are preserved", {
    dat <- data.frame(ID = rep(c(3, 17, 42), each = 5),
        time = rep(1:5, 3),
        p11 = rep(c(1, 0.8, 0.6, 0.4, 0.2), 3))
    result <- sim_tte_df(dat)
    expect_setequal(result$ID, c(3, 17, 42))
})

test_that("character IDs work end to end", {
    dat <- data.frame(ID = rep(c("ctrl", "trt"), each = 5),
        time = rep(1:5, 2),
        p11 = rep(c(1, 0.8, 0.6, 0.4, 0.2), 2))
    result <- sim_tte_df(dat)
    expect_setequal(result$ID, c("ctrl", "trt"))
})

test_that("sim_tte_df is reproducible with a fixed seed", {
    dat <- data.frame(ID = rep(1:10, each = 50),
        time = rep(seq(0.1, 5, length.out = 50), 10),
        p11 = rep(exp(-0.5 * seq(0.1, 5, length.out = 50)), 10))
    set.seed(2026)
    r1 <- sim_tte_df(dat)
    set.seed(2026)
    r2 <- sim_tte_df(dat)
    expect_identical(r1, r2)
})

test_that("sim_tte is reproducible with a fixed seed across full pipeline", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    lp <- matrix(c(-0.3, 0, 0.3, 0.6), nrow = 4)
    set.seed(7)
    r1 <- sim_tte(pi = lp, mu = -1, coefs = 1.3,
        time = seq(0.1, 20, by = 0.2), type = "weibull", end_time = 20)
    set.seed(7)
    r2 <- sim_tte(pi = lp, mu = -1, coefs = 1.3,
        time = seq(0.1, 20, by = 0.2), type = "weibull", end_time = 20)
    expect_identical(r1, r2)
})
