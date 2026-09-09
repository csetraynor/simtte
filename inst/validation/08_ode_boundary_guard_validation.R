## Phase 1 validation: sim_tte_ode() in-solver event detection
## (exponential/constant-hazard library model only).
##
## Package-scale reproduction of reports/experiments/
## 01_etime_boundary_experiment.R, run against the shipped
## `sim_tte_ode(model = "exponential")` (not an ad-hoc mcode() string),
## recording:
##   (a) quantization -- how many distinct internal SOLVERTIME values
##       actually resolve events, vs. the requested output grid;
##   (b) raw (unrefined TEVT) vs. refined (log_survival-interpolated)
##       accuracy against the closed-form exponential quantile;
##   (c) rtol/atol sensitivity -- one and two orders of magnitude
##       tighter and looser than the mrgsolve default (1e-8), per
##       reports/03_implementation_plan.md risk R2.
##
## Reproducibility: fixed seeds throughout; session info printed at the
## end. This script is excluded from the package build/tarball by the
## existing `^inst/validation$` .Rbuildignore rule (unchanged from every
## prior validation script).

suppressMessages(library(simtte))
suppressMessages(library(dplyr))

H0 <- 0.12
n <- 3000
end <- 40
delta <- 4   # deliberately coarse output grid, to make (a)/(b) informative
seed <- 20260910

cat("==== (a) Quantization: distinct internal SOLVERTIME values used ====\n")
sim <- sim_tte_ode(model = "exponential", param = list(H0 = H0), n = n,
    end = end, delta = delta, keep_trajectory = TRUE, seed = seed)
traj <- sim$trajectory
last_rows <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
n_events <- sum(last_rows$event_found == 1)
n_distinct_tevt <- length(unique(last_rows$TEVT[last_rows$event_found == 1]))
n_grid_points <- length(unique(traj$time))
cat("Requested output grid points:", n_grid_points, "\n")
cat("Subjects with a latched event:", n_events, "of", n, "\n")
cat("Distinct internal TEVT values actually used:", n_distinct_tevt, "\n")

cat("\n==== (b) Raw vs. refined accuracy vs. closed-form exponential quantile ====\n")
idata_u <- data.frame(ID = seq_len(n))
# Recover each subject's own U from the trajectory (carried via
# carry_out) so the analytical reference uses the *same* draw.
u_by_id <- last_rows$U[match(seq_len(n), last_rows$ID)]
analytic_T <- -log(u_by_id) / H0
is_event <- sim$events$sim_status == 1L & analytic_T < end

raw_tevt <- last_rows$TEVT[match(sim$events$ID, last_rows$ID)]
err_raw <- abs(raw_tevt[is_event] - analytic_T[is_event])
err_refined <- abs(sim$events$sim_time[is_event] - analytic_T[is_event])

cat("n compared:", sum(is_event), "\n")
cat("Raw TEVT      -- mean abs err:", round(mean(err_raw), 5),
    " max abs err:", round(max(err_raw), 5), "\n")
cat("Refined       -- mean abs err:", round(mean(err_refined), 5),
    " max abs err:", round(max(err_refined), 5), "\n")
cat("Refinement reduces mean abs error by:",
    round(100 * (1 - mean(err_refined) / mean(err_raw)), 1), "%\n")
stopifnot(mean(err_refined) < mean(err_raw))

cat("\n==== (c) rtol/atol sensitivity (risk R2) ====\n")
run_with_tol <- function(tol) {
    s <- sim_tte_ode(model = "exponential", param = list(H0 = H0), n = n,
        end = end, delta = delta, keep_trajectory = TRUE, seed = seed,
        rtol = tol, atol = tol)
    last <- s$trajectory[!duplicated(s$trajectory$ID, fromLast = TRUE), ]
    u_i <- last$U[match(seq_len(n), last$ID)]
    t_true <- -log(u_i) / H0
    ok <- s$events$sim_status == 1L & t_true < end
    raw_t <- last$TEVT[match(s$events$ID, last$ID)]
    c(mean_abs_err_raw = mean(abs(raw_t[ok] - t_true[ok])),
        mean_abs_err_refined = mean(abs(s$events$sim_time[ok] - t_true[ok])))
}

tol_default <- 1e-8
tols <- c(loose_2x = tol_default * 100, loose_1x = tol_default * 10,
    default = tol_default, tight_1x = tol_default / 10,
    tight_2x = tol_default / 100)
sens <- t(sapply(tols, run_with_tol))
sens <- data.frame(rtol_atol = tols, sens, row.names = names(tols))
print(sens)

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
cat("Platform:", R.version$platform, "\n")
cat("Seed:", seed, "; n:", n, "; H0:", H0, "; end:", end, "; delta:", delta,
    "\n")
