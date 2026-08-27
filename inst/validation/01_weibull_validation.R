## Phase D validation, component 1: Weibull survival function validation
## against an independent analytical reference.
##
## simtte's built-in Weibull model computes
##   S(t) = exp(-eta * t^shape),  eta = exp(mu + lp)
## as a closed-form expression inside a compiled mrgsolve model
## (inst/models/weibull.cpp). This script re-derives S(t) two ways that
## do NOT call that compiled model file, and compares both against
## simtte's actual output:
##
##   (a) a direct, hand-written R re-implementation of the closed-form
##       expression (independent of the C++ source, though the same
##       underlying mathematics);
##   (b) R's own base-distribution function stats::pweibull(), via the
##       reparameterization scale = eta^(-1/shape). This is a genuinely
##       independent reference implementation: a different codebase
##       (R's C-level distribution routines) maintained outside this
##       package, not derived from or aware of simtte's source.
##
## Reproducibility: no randomness is used anywhere in this script (the
## Weibull survival function is deterministic given its parameters), so
## no seed is required. R version and mrgsolve version used to produce
## the results in PHASE_D_REPORT.md are recorded at the end.

suppressMessages(library(simtte))

## ---- Independent references -------------------------------------------

## (a) Direct hand-written re-implementation (not sourced from the
## package; mirrors the mathematical definition only).
analytical_S_direct <- function(t, mu, lp, shape) {
    eta <- exp(mu + lp)
    exp(-eta * t^shape)
}

## (b) R base stats::pweibull(), an independent codebase.
##   S(t) = exp(-eta * t^shape) = exp(-(t/scale)^shape),
##   scale = eta^(-1/shape)
analytical_S_pweibull <- function(t, mu, lp, shape) {
    eta <- exp(mu + lp)
    scale <- eta^(-1 / shape)
    stats::pweibull(t, shape = shape, scale = scale, lower.tail = FALSE)
}

## ---- Parameter grid: representative settings + edge cases -------------

shapes <- c(0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3, 5, 10)
mus <- c(-5, -1, 0, 1, 5)
lps <- c(-2, -0.5, 0, 0.5, 2)
times <- c(0, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 50, 100,
    1000)

## Edge cases identified in prior phases as numerically hazardous for a
## naive implementation: extreme mu/lp causing eta to overflow, and
## t = 0 exactly for shape < 1 (divergent instantaneous hazard).
edge_cases <- data.frame(
    mu    = c(710, 710, -700, -700, 0, 0),
    lp    = c(0, 0, 0, 0, 700, -700),
    shape = c(1, 0.3, 1, 0.3, 0.5, 2)
)

run_grid <- function() {
    results <- list()
    for (shape in shapes) {
        for (mu in mus) {
            for (lp in lps) {
                out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu,
                    shape = shape, type = "weibull", times = times)
                ref_direct <- analytical_S_direct(out$time, mu, lp, shape)
                ref_pweib <- analytical_S_pweibull(out$time, mu, lp, shape)
                results[[length(results) + 1]] <- data.frame(
                    shape = shape, mu = mu, lp = lp, time = out$time,
                    p11 = out$p11, ref_direct = ref_direct,
                    ref_pweibull = ref_pweib,
                    abs_err_direct = abs(out$p11 - ref_direct),
                    abs_err_pweibull = abs(out$p11 - ref_pweib)
                )
            }
        }
    }
    dplyr::bind_rows(results)
}

run_edge_cases <- function() {
    out_list <- lapply(seq_len(nrow(edge_cases)), function(i) {
        mu <- edge_cases$mu[i]
        lp <- edge_cases$lp[i]
        shape <- edge_cases$shape[i]
        out <- simtte:::.sim_surv_df(log_hr = lp, mu = mu, shape = shape,
            type = "weibull", times = c(0, 1e-6, 1e-3, 1, 100))
        ref <- analytical_S_direct(out$time, mu, lp, shape)
        data.frame(mu = mu, lp = lp, shape = shape, time = out$time,
            p11 = out$p11, ref = ref, ok_finite = is.finite(out$p11),
            abs_err = ifelse(is.finite(ref), abs(out$p11 - ref), NA))
    })
    dplyr::bind_rows(out_list)
}

cat("Running main parameter grid (", length(shapes) * length(mus) *
    length(lps), " parameter combinations x ", length(times),
    " time points)...\n", sep = "")
main_results <- run_grid()

cat("Running numerical edge cases...\n")
edge_results <- run_edge_cases()

## ---- Summaries ----------------------------------------------------------

cat("\n==== Main grid: error summary (vs. direct R re-implementation) ====\n")
print(summary(main_results$abs_err_direct))
cat("Max absolute error:", max(main_results$abs_err_direct), "\n")
cat("Max absolute error location:\n")
print(main_results[which.max(main_results$abs_err_direct),
    c("shape", "mu", "lp", "time", "p11", "ref_direct")])

cat("\n==== Main grid: error summary (vs. stats::pweibull) ====\n")
print(summary(main_results$abs_err_pweibull))
cat("Max absolute error:", max(main_results$abs_err_pweibull), "\n")

cat("\n==== Near-zero times only (t <= 0.01), by shape ====\n")
near_zero <- main_results[main_results$time <= 0.01, ]
print(stats::aggregate(abs_err_pweibull ~ shape, data = near_zero, FUN = max))

cat("\n==== Edge cases ====\n")
print(edge_results)
cat("All edge-case p11 values finite:", all(edge_results$ok_finite), "\n")
finite_ref <- edge_results[is.finite(edge_results$ref), ]
cat("Max abs error where analytical reference itself is finite:",
    max(finite_ref$abs_err), "\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
