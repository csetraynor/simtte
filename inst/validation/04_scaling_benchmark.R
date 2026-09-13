## Phase D validation, component 4: internal computational scaling
## benchmark.
##
## Purely descriptive: measures wall-clock time for sim_tte() (Weibull
## model, the fastest/simplest built-in model) as (a) cohort size and
## (b) output-grid density increase. No package code is modified or
## tuned as a result of this benchmark; it is a snapshot for future
## comparison, not a performance-optimization exercise.
##
## Internal benchmark only: this script uses base R timing
## (base::system.time(), repeated and summarized with the median to
## reduce noise) rather than an external package, so it introduces no
## new dependency and remains reproducible with only what simtte already
## Suggests/Imports. No fair, readily reproducible external comparator
## (e.g. an equivalent simulation in another package with a directly
## comparable model) was identified as being clearly in scope for this
## phase, so none is attempted here.
##
## Reproducibility: timings are inherently machine-dependent (CPU,
## load, R build); the *relative* scaling pattern reported here (how
## time changes as n or grid density changes on this machine, at this
## R/mrgsolve version) is the reproducible artifact, not the absolute
## timings. Each configuration is timed 5 times and the median is
## reported to reduce sensitivity to transient system load. A fixed
## seed is used throughout, though timing does not depend on the RNG
## stream.

suppressMessages(library(simtte))

set.seed(20260827)
n_reps <- 5

time_call <- function(n, delta) {
    lp <- matrix(stats::rnorm(n, 0, 0.5), nrow = n)
    grid <- seq(0.1, 20, by = delta)
    times <- replicate(n_reps, {
        t <- system.time({
            sim_tte(pi = lp, mu = -1, coefs = 1.2, time = grid,
                type = "weibull", end_time = 20)
        })
        t[["elapsed"]]
    })
    stats::median(times)
}

## ---- Scaling with cohort size (n), fixed grid -------------------------

cat("Benchmarking cohort-size scaling (grid fixed at delta = 0.1)...\n")
n_values <- c(10, 50, 100, 500, 1000, 2000, 5000)
n_scaling <- data.frame(
    n = n_values,
    median_seconds = vapply(n_values, time_call, numeric(1), delta = 0.1)
)
n_scaling$seconds_per_1000_subjects <-
    n_scaling$median_seconds / n_scaling$n * 1000
print(n_scaling, row.names = FALSE, digits = 4)

## ---- Scaling with output-grid density, fixed cohort size --------------

cat("\nBenchmarking grid-density scaling (n fixed at 500 subjects)...\n")
delta_values <- c(1, 0.5, 0.2, 0.1, 0.05, 0.02)
grid_scaling <- data.frame(
    delta = delta_values,
    n_grid_points = vapply(delta_values,
        function(d) length(seq(0.1, 20, by = d)), integer(1)),
    median_seconds = vapply(delta_values, time_call, numeric(1), n = 500)
)
print(grid_scaling, row.names = FALSE, digits = 4)

## ---- Empirical scaling exponent (log-log slope) ------------------------

slope <- function(x, y) {
    fit <- stats::lm(log(y) ~ log(x))
    unname(stats::coef(fit)[2])
}
cat("\nEmpirical scaling exponent (time ~ n^k), k =",
    round(slope(n_scaling$n, n_scaling$median_seconds), 2), "\n")
cat("Empirical scaling exponent (time ~ n_grid_points^k), k =",
    round(slope(grid_scaling$n_grid_points, grid_scaling$median_seconds), 2),
    "\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
cat("Platform:", R.version$platform, "\n")
cat("Replicates per configuration:", n_reps, "(median reported)\n")
