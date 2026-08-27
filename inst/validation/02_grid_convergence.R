## Phase D validation, component 2: grid convergence.
##
## Quantifies how event-time error (relative to an independent
## closed-form analytical reference) behaves as the reported `time`
## grid is refined, for both event_time_method = "grid" (the package
## default) and "log_survival" (opt-in interpolation).
##
## This script does NOT modify or reimplement interpolation/root
## finding: it calls the package's own internal building blocks exactly
## as sim_tte_df()/.simulate_survival_id() do --
##   simtte:::.sim_surv_df()          (generates S(t) on a grid)
##   simtte:::.get_tte()              (locates the crossing index)
##   simtte:::.interpolate_log_survival()  (the log_survival refinement)
## -- driven by uniform draws (U) generated directly in this script, so
## that the true continuous-time analytical event time
##   T_true = (-log(U) / eta)^(1 / shape)
## is known exactly for every subject and does not depend on any
## behavior internal to sim_tte()/sim_tte_df() (in particular, this
## sidesteps needing to recover individual U draws from the public
## sim_tte() API).
##
## Reproducibility: a single seed controls the U draws; grid definitions
## are deterministic given that seed. All other quantities are
## deterministic (the Weibull S(t) trajectory is closed-form, validated
## in 01_weibull_validation.R).

suppressMessages(library(simtte))

set.seed(20260827)

## ---- Scenario definitions -----------------------------------------------
## Two shapes bracketing shape = 1 (the case where log_survival is
## exact, per PHASE_B_DESIGN_AUDIT.md Section 7-8): a convex-hazard case
## (shape > 1) and a concave-hazard case (shape < 1).
mu <- -1
lp <- 0.3
eta <- exp(mu + lp)
t_max <- 20
n_subjects <- 800

deltas <- c(2, 1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005)

U <- stats::runif(n_subjects)

analytical_T <- function(U, eta, shape) (-log(U) / eta)^(1 / shape)

evaluate_grid <- function(shape, delta) {
    grid <- seq(0, t_max, by = delta)
    out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
        type = "weibull", times = grid)
    S <- out$p11
    t <- out$time

    T_true <- analytical_T(U, eta, shape)
    censored_by_range <- T_true > t_max # true event beyond simulated range

    T_grid <- rep(NA_real_, length(U))
    T_log <- rep(NA_real_, length(U))
    for (i in seq_along(U)) {
        etime <- simtte:::.get_tte(U[i], S)
        if (etime == -99L) next # censored on this grid
        if (etime > 1L) {
            T_log[i] <- simtte:::.interpolate_log_survival(
                t_i = t[etime - 1L], t_ip1 = t[etime],
                s_i = S[etime - 1L], s_ip1 = S[etime], u = U[i])
        } else {
            T_log[i] <- t[etime]
        }
        T_grid[i] <- t[etime]
    }

    keep <- !censored_by_range & !is.na(T_grid)
    data.frame(
        shape = shape, delta = delta, n_points = length(grid),
        n_compared = sum(keep),
        mae_grid = mean(abs(T_grid[keep] - T_true[keep])),
        mae_log_survival = mean(abs(T_log[keep] - T_true[keep])),
        max_ae_grid = max(abs(T_grid[keep] - T_true[keep])),
        max_ae_log_survival = max(abs(T_log[keep] - T_true[keep]))
    )
}

cat("Evaluating grid convergence for shape = 0.5 (concave hazard, H)...\n")
res_concave <- dplyr::bind_rows(lapply(deltas, evaluate_grid, shape = 0.5))

cat("Evaluating grid convergence for shape = 2 (convex hazard, H)...\n")
res_convex <- dplyr::bind_rows(lapply(deltas, evaluate_grid, shape = 2))

cat("Evaluating grid convergence for shape = 1 (exponential; log_survival exact)...\n")
res_exp <- dplyr::bind_rows(lapply(deltas, evaluate_grid, shape = 1))

cat("\n==== shape = 0.5 ====\n")
print(res_concave, row.names = FALSE, digits = 5)

cat("\n==== shape = 2 ====\n")
print(res_convex, row.names = FALSE, digits = 5)

cat("\n==== shape = 1 ====\n")
print(res_exp, row.names = FALSE, digits = 5)

## ---- Empirical convergence order ----------------------------------------
## log(error) ~ log(delta): the slope estimates the observed order of
## convergence (grid method expected O(delta); log_survival expected
## O(delta^2) locally, per PHASE_B_DESIGN_AUDIT.md Section 7).
convergence_slope <- function(res, col) {
    ok <- res[[col]] > 0
    fit <- stats::lm(log(res[[col]][ok]) ~ log(res$delta[ok]))
    unname(stats::coef(fit)[2])
}

cat("\n==== Empirical convergence order (log-log slope of MAE vs delta) ====\n")
for (nm in c("concave (shape=0.5)", "convex (shape=2)")) {
    res <- if (grepl("concave", nm)) res_concave else res_convex
    cat(sprintf("%-22s grid: %5.2f   log_survival: %5.2f\n", nm,
        convergence_slope(res, "mae_grid"),
        convergence_slope(res, "mae_log_survival")))
}

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
cat("Seed: 20260827; n_subjects:", n_subjects, "; mu:", mu, "; lp:", lp,
    "; t_max:", t_max, "\n")
