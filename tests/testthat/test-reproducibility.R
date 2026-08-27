# Pre-Phase-B audit, Phase 7: reproducibility semantics.
#
# sim_tte_df() draws one stats::runif(1) per subject, in row-encounter
# order. This is documented in the "Reproducibility" section of
# ?sim_tte_df as reproducible for a *fixed row order*, not invariant to
# reordering subject blocks. This is standard R RNG behavior (a single
# shared stream consumed sequentially) and is not changed here -- the
# test below documents/demonstrates it rather than "fixing" it.

test_that("a fixed seed reproduces results for a fixed row order", {
    dat <- data.frame(
        ID = rep(1:5, each = 20),
        time = rep(seq(0.5, 10, by = 0.5), 5),
        p11 = rep(exp(-0.4 * seq(0.5, 10, by = 0.5)), 5)
    )
    set.seed(11)
    r1 <- sim_tte_df(dat)
    set.seed(11)
    r2 <- sim_tte_df(dat)
    expect_identical(r1, r2)
})

test_that("reversing subject block order changes draw assignment under the same seed", {
    # Two subjects, identical (non-degenerate) trajectories, so any
    # difference in outcome must come from which runif() draw each
    # subject receives -- not from the trajectories themselves.
    traj <- data.frame(time = seq(0.5, 10, by = 0.5),
        p11 = exp(-0.4 * seq(0.5, 10, by = 0.5)))
    forward <- rbind(
        cbind(ID = 1, traj),
        cbind(ID = 2, traj)
    )
    reversed <- rbind(
        cbind(ID = 2, traj),
        cbind(ID = 1, traj)
    )

    set.seed(2026)
    r_forward <- sim_tte_df(forward)
    set.seed(2026)
    r_reversed <- sim_tte_df(reversed)

    id1_forward <- r_forward[r_forward$ID == 1, ]
    id1_reversed <- r_reversed[r_reversed$ID == 1, ]

    # This is the documented, expected behavior (not a bug): subject 1
    # gets the *first* runif() draw when it appears first in the data,
    # and the *second* draw when subject 2 appears first instead. The
    # two draws differ, so (in general) the resulting sim_time/status
    # for subject 1 differs between the two row orders.
    expect_false(isTRUE(all.equal(id1_forward$sim_time,
        id1_reversed$sim_time)))
})

test_that("reordering rows does not change the *set* of subjects or the trajectory contract", {
    dat <- data.frame(
        ID = rep(1:3, each = 10),
        time = rep(1:10, 3),
        p11 = rep(seq(1, 0.1, length.out = 10), 3)
    )
    dat_reordered <- dat[order(-dat$ID), ]

    r1 <- sim_tte_df(dat)
    r2 <- sim_tte_df(dat_reordered)
    expect_setequal(r1$ID, r2$ID)
    expect_equal(nrow(r1), nrow(r2))
})
