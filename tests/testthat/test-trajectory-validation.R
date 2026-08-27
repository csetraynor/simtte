# Phase A2/A3: input and trajectory validation for sim_tte_df(), and the
# custom-trajectory censoring contract.

valid_dat <- function() {
    data.frame(
        ID = rep(1:2, each = 5),
        time = rep(1:5, 2),
        p11 = rep(c(1, 0.8, 0.6, 0.4, 0.2), 2)
    )
}

test_that("empty input is rejected", {
    dat <- data.frame(ID = integer(0), time = numeric(0), p11 = numeric(0))
    expect_error(sim_tte_df(dat), "no rows")
})

test_that("NA in time or survival is rejected", {
    dat <- valid_dat()
    dat_t <- dat; dat_t$time[3] <- NA_real_
    expect_error(sim_tte_df(dat_t), "must not contain NA")

    dat_s <- dat; dat_s$p11[3] <- NA_real_
    expect_error(sim_tte_df(dat_s), "must not contain NA")
})

test_that("NaN in time or survival is rejected", {
    dat <- valid_dat()
    dat_t <- dat; dat_t$time[3] <- NaN
    expect_error(sim_tte_df(dat_t), "must not contain NA")

    dat_s <- dat; dat_s$p11[3] <- NaN
    expect_error(sim_tte_df(dat_s), "must not contain NA")
})

test_that("Inf and -Inf in time or survival are rejected", {
    dat <- valid_dat()
    dat_t <- dat; dat_t$time[3] <- Inf
    expect_error(sim_tte_df(dat_t), "Inf")

    dat_t2 <- dat; dat_t2$time[3] <- -Inf
    expect_error(sim_tte_df(dat_t2), "Inf")

    dat_s <- dat; dat_s$p11[3] <- Inf
    expect_error(sim_tte_df(dat_s), "Inf")
})

test_that("survival outside [0, 1] is rejected", {
    dat <- valid_dat()
    dat_neg <- dat; dat_neg$p11[3] <- -0.1
    expect_error(sim_tte_df(dat_neg), "\\[0, 1\\]")

    dat_big <- dat; dat_big$p11[3] <- 1.1
    expect_error(sim_tte_df(dat_big), "\\[0, 1\\]")
})

test_that("survival probabilities within floating-point tolerance of [0, 1] are accepted", {
    dat <- valid_dat()
    dat$p11[1] <- 1 + 5e-9
    dat$p11[nrow(dat)] <- -5e-9
    expect_error(sim_tte_df(dat), NA)
})

test_that("negative time is rejected", {
    dat <- valid_dat()
    dat$time[1] <- -1
    expect_error(sim_tte_df(dat), "negative")
})

test_that("unsorted time within a subject is rejected, not silently sorted", {
    dat <- data.frame(
        ID = rep(1, 4),
        time = c(1, 3, 2, 4),
        p11 = c(1, 0.8, 0.9, 0.4)
    )
    expect_error(sim_tte_df(dat), "ascending order")
})

test_that("duplicate time within a subject is rejected", {
    dat <- data.frame(
        ID = rep(1, 4),
        time = c(1, 2, 2, 3),
        p11 = c(1, 0.8, 0.7, 0.4)
    )
    expect_error(sim_tte_df(dat), "duplicated value")
})

test_that("non-monotone (increasing) survival is rejected, not silently repaired", {
    dat <- data.frame(
        ID = rep(1, 4),
        time = c(1, 2, 3, 4),
        p11 = c(1, 0.5, 0.7, 0.3)
    )
    expect_error(sim_tte_df(dat), "non-increasing")
})

test_that("flat (constant) survival is accepted as valid non-monotone-free input", {
    dat <- data.frame(ID = rep(1, 5), time = 1:5, p11 = rep(0.999, 5))
    result <- sim_tte_df(dat)
    expect_equal(nrow(result), 1)
})

test_that("survival noise within tolerance does not trigger a monotonicity error", {
    dat <- data.frame(
        ID = rep(1, 4),
        time = c(1, 2, 3, 4),
        # tiny non-monotone bump well inside the numerical tolerance
        p11 = c(1, 0.6000000005, 0.6, 0.3)
    )
    expect_error(sim_tte_df(dat), NA)
})

test_that("missing (NA) IDs are rejected", {
    dat <- valid_dat()
    dat$ID[1] <- NA
    expect_error(sim_tte_df(dat), "missing values")
})

test_that("non-finite numeric IDs are rejected", {
    dat <- valid_dat()
    dat$ID[1] <- Inf
    expect_error(sim_tte_df(dat), "Inf")
})

test_that("a single-observation trajectory is a valid, explicitly accepted input", {
    # Minimum trajectory length is one row; this is not an error.
    dat <- data.frame(ID = 1, time = 5, p11 = 0.3)
    result <- sim_tte_df(dat)
    expect_equal(nrow(result), 1)
    expect_equal(result$sim_time, 5)
})

test_that("a trajectory that does not start at survival = 1 is explicitly accepted", {
    # sim_tte_df() never invents an implicit S(0) = 1 state.
    dat <- data.frame(ID = 1, time = c(5, 6, 7), p11 = c(0, 0, 0))
    result <- sim_tte_df(dat)
    # p11 = 0 at the very first observation is <= any U in [0, 1), so
    # the event is deterministically returned at the first reported
    # time, regardless of the random draw.
    expect_equal(result$sim_time, 5)
    expect_equal(result$sim_status, 1)
})

test_that("a trajectory whose first observation is already <= U returns that first time", {
    dat <- data.frame(ID = c(1, 1, 1), time = c(2, 4, 6), p11 = c(0, 0, 0))
    result <- sim_tte_df(dat)
    expect_equal(result$sim_time, 2)
    expect_equal(result$sim_status, 1)
})

test_that("subject-specific final observation times define subject-specific censoring", {
    # p11 == 1 never crosses U (support of runif() is [0, 1)), so
    # censoring is deterministic here.
    dat <- rbind(
        data.frame(ID = 1, time = 1:5, p11 = 1),
        data.frame(ID = 2, time = 1:10, p11 = 1)
    )
    result <- sim_tte_df(dat)
    expect_equal(result$sim_status, c(0, 0))
    expect_equal(result$sim_time[result$ID == 1], 5)
    expect_equal(result$sim_time[result$ID == 2], 10)
})

test_that("valid multi-subject trajectories with different ID types all work", {
    dat_numeric <- valid_dat()
    expect_error(sim_tte_df(dat_numeric), NA)

    dat_noncontig <- valid_dat()
    dat_noncontig$ID <- rep(c(7, 42), each = 5)
    result <- sim_tte_df(dat_noncontig)
    expect_setequal(result$ID, c(7, 42))

    dat_char <- valid_dat()
    dat_char$ID <- rep(c("subjA", "subjB"), each = 5)
    result_char <- sim_tte_df(dat_char)
    expect_setequal(result_char$ID, c("subjA", "subjB"))
})
