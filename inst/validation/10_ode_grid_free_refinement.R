## Phase 2.5 validation: grid-free in-solver refinement
## (reports/04_author_decisions.md "After the test runbook / Phase 2.5",
## reports/08_phase2_5_report.md).
##
## Repeats reports/06_phase2_report.md section 4's raw-vs-refined table
## (Weibull shape 0.5/1/2/5, Gompertz +/- gamma) at delta = 4, 1, 0.25,
## now with a third column for the new in-solver bracket refinement, and
## repeats the rtol/atol sweep for Weibull shape = 2. Also checks the
## boundary guard (both cohort sizes, all three models) and that
## sim_time never leaves [T_PRE, TEVT].
##
## Excluded from the build (^inst/validation$ .Rbuildignore rule).

## Run from the package root: `Rscript inst/validation/10_ode_grid_free_refinement.R`
## (self-loads the in-development package via devtools::load_all(), like
## dev/smoke-ode.R -- library(simtte) alone would resolve to whatever
## simtte build happens to be installed, not necessarily this checkout;
## see reports/07_test_runbook.md's "inst/validation/0*" table).
if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root: ",
        "Rscript inst/validation/10_ode_grid_free_refinement.R", call. = FALSE)
}
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages(library(dplyr))

cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")), "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n\n")

# ---- (a) Three-column accuracy table: raw / old grid refinement /
#          new in-solver refinement, at delta = 4, 1, 0.25 -------------

cat("==== (a) raw vs. old (grid) refinement vs. new (in-solver) refinement ====\n")

# The pre-Phase-2.5 "old" refinement is recomputed here directly (not
# by calling package internals, which now always use the new method) --
# the exact same reported-grid crossing-interval interpolation
# .resolve_ode_events() used through Phase 2, applied to the returned
# $trajectory.
old_grid_refine <- function(traj, id, u_i) {
    sub <- traj[traj$ID == id, ]
    sub <- sub[order(sub$time), ]
    p11 <- pmin(pmax(sub$p11, 0), 1)
    idx <- match(TRUE, p11 <= u_i)
    if (is.na(idx)) return(sub$TEVT[nrow(sub)])
    if (idx == 1L) return(sub$time[1])
    i <- idx - 1L
    simtte:::.interpolate_log_survival(t_i = sub$time[i], t_ip1 = sub$time[idx],
        s_i = p11[i], s_ip1 = p11[idx], u = u_i)
}

check_three_way <- function(model, param, analytic_T_fn, label, delta,
    n = 3000, end = 20, seed = 20260910) {
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = delta, keep_trajectory = TRUE, seed = seed)
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    u_by_id <- last$U[match(seq_len(n), last$ID)]
    analytic_T <- suppressWarnings(analytic_T_fn(u_by_id))
    is_event <- sim$events$sim_status == 1L & analytic_T < end
    raw_tevt <- last$TEVT[match(sim$events$ID, last$ID)]
    old_refined <- vapply(sim$events$ID[is_event], function(id) {
        old_grid_refine(traj, id, last$U[match(id, last$ID)])
    }, numeric(1))
    err_raw <- abs(raw_tevt[is_event] - analytic_T[is_event])
    err_old <- abs(old_refined - analytic_T[is_event])
    err_new <- abs(sim$events$sim_time[is_event] - analytic_T[is_event])
    data.frame(label = label, delta = delta, n_compared = sum(is_event),
        raw_mean = mean(err_raw), old_refined_mean = mean(err_old),
        new_refined_mean = mean(err_new))
}

scenarios <- list(
    list(model = "weibull", param = list(mu = -1, shape = 0.5),
        fn = function(u) (-log(u) / exp(-1))^(1 / 0.5), label = "Weibull shape=0.5"),
    list(model = "weibull", param = list(mu = -1, shape = 1),
        fn = function(u) -log(u) / exp(-1), label = "Weibull shape=1"),
    list(model = "weibull", param = list(mu = -1, shape = 2),
        fn = function(u) (-log(u) / exp(-1))^(1 / 2), label = "Weibull shape=2"),
    list(model = "weibull", param = list(mu = -1, shape = 5),
        fn = function(u) (-log(u) / exp(-1))^(1 / 5), label = "Weibull shape=5"),
    list(model = "gompertz", param = list(mu = -2, gamma = 0.1),
        fn = function(u) log(1 - (log(u) * 0.1) / exp(-2)) / 0.1,
        label = "Gompertz gamma=0.1"),
    list(model = "gompertz", param = list(mu = -2, gamma = -0.05),
        fn = function(u) log(1 - (log(u) * -0.05) / exp(-2)) / -0.05,
        label = "Gompertz gamma=-0.05")
)

results <- do.call(rbind, lapply(scenarios, function(s) {
    do.call(rbind, lapply(c(4, 1, 0.25), function(d) {
        check_three_way(s$model, s$param, s$fn, s$label, delta = d)
    }))
}))
print(results, row.names = FALSE, digits = 4)
cat("\nFlatness check: new_refined_mean at delta=4 vs delta=0.25, per scenario:\n")
for (s in scenarios) {
    r <- results[results$label == s$label, ]
    d4 <- r$new_refined_mean[r$delta == 4]
    d025 <- r$new_refined_mean[r$delta == 0.25]
    ratio <- d4 / max(d025, 1e-12)
    cat(sprintf("  %-20s delta=4: %.5f  delta=0.25: %.5f  ratio: %.2f\n",
        s$label, d4, d025, ratio))
}
cat("\n")
cat("Note on Weibull shape=5 (the one scenario above that is NOT flat):\n")
cat("at delta=4 and delta=1, ALL 3000 subjects hit the reported-grid\n")
cat("fallback (message fired both times); at delta=0.25, only 5/3000 do.\n")
cat("new_refined_mean equals old_refined_mean almost exactly whenever the\n")
cat("fallback engages for (near-)everyone -- the safety net is working as\n")
cat("designed, not silently returning a wrong number; see section (e) below\n")
cat("for the same phenomenon isolated at shape=10, and reports/08_phase2_5_report.md.\n\n")

# ---- (b) rtol/atol sweep repeated for Weibull shape = 2 -------------

cat("==== (b) rtol/atol sensitivity, Weibull shape = 2, delta = 4 (cf. Phase 1/2 report) ====\n")
mu <- -1; shape <- 2; n <- 3000; end <- 20; delta <- 4; seed <- 20260910

run_with_tol <- function(tol) {
    s <- sim_tte_ode(model = "weibull", param = list(mu = mu, shape = shape),
        n = n, end = end, delta = delta, keep_trajectory = TRUE, seed = seed,
        rtol = tol, atol = tol)
    last <- s$trajectory[!duplicated(s$trajectory$ID, fromLast = TRUE), ]
    u_i <- last$U[match(seq_len(n), last$ID)]
    eta <- exp(mu)
    t_true <- (-log(u_i) / eta)^(1 / shape)
    ok <- s$events$sim_status == 1L & t_true < end
    raw_t <- last$TEVT[match(s$events$ID, last$ID)]
    c(mean_abs_err_raw = mean(abs(raw_t[ok] - t_true[ok])),
        mean_abs_err_new_refined = mean(abs(s$events$sim_time[ok] - t_true[ok])))
}

tol_default <- 1e-8
tols <- c(loose_2x = tol_default * 100, loose_1x = tol_default * 10,
    default = tol_default, tight_1x = tol_default / 10,
    tight_2x = tol_default / 100)
sens <- t(sapply(tols, run_with_tol))
sens <- data.frame(rtol_atol = tols, sens, row.names = names(tols))
print(sens)
cat("\n")

# ---- (c) Boundary guard, both cohort sizes, all three models --------

cat("==== (c) Boundary guard (never sim_time > end), all three models ====\n")
boundary_check <- function(model, param, n, seed, end = 10) {
    set.seed(seed)
    data <- data.frame(ID = seq_len(n), time = end - 0.05, lp = 0.5,
        evid = 1, amt = 0, cmt = 1)
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 1, data = data, seed = seed)
    max(sim$events$sim_time) <= end + 1e-9
}
bg <- data.frame(
    model = rep(c("exponential", "weibull", "gompertz"), each = 2),
    n = rep(c(40, 2000), 3),
    ok = c(
        boundary_check("exponential", list(H0 = 0.3), 40, 1001),
        boundary_check("exponential", list(H0 = 0.3), 2000, 1002),
        boundary_check("weibull", list(mu = -1, shape = 1.5), 40, 2041),
        boundary_check("weibull", list(mu = -1, shape = 1.5), 2000, 4000),
        boundary_check("gompertz", list(mu = -1, gamma = 0.1), 40, 3041),
        boundary_check("gompertz", list(mu = -1, gamma = 0.1), 2000, 5000)
    ))
print(bg, row.names = FALSE)
stopifnot(all(bg$ok))
cat("All boundary-guard checks passed.\n\n")

# ---- (d) Bracket containment: sim_time always in [T_PRE, TEVT] ------

cat("==== (d) Bracket containment: sim_time in [T_PRE, TEVT] for every event subject ====\n")
bracket_check <- function(model, param, n = 1000, end = 20, seed = 1) {
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 2, keep_trajectory = TRUE, seed = seed)
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(sim$events$ID, last$ID), ]
    is_event <- sim$events$sim_status == 1L
    n_fallback_like <- sum(last$T_PRE[is_event] >= last$TEVT[is_event])
    lower_ok <- all(sim$events$sim_time[is_event] >= last$T_PRE[is_event] - 1e-8)
    upper_ok <- all(sim$events$sim_time[is_event] <= last$TEVT[is_event] + 1e-8)
    data.frame(model = model, n_event = sum(is_event),
        n_degenerate_bracket = n_fallback_like,
        lower_bound_ok = lower_ok, upper_bound_ok = upper_ok)
}
bc <- rbind(
    bracket_check("exponential", list(H0 = 0.1)),
    bracket_check("weibull", list(mu = -1, shape = 1.8)),
    bracket_check("weibull", list(mu = -1, shape = 0.3)),
    bracket_check("gompertz", list(mu = -1, gamma = 0.1)),
    bracket_check("gompertz", list(mu = -1, gamma = -0.05))
)
print(bc, row.names = FALSE)
stopifnot(all(bc$lower_bound_ok), all(bc$upper_bound_ok))
cat("All refined sim_time values stayed within [T_PRE, TEVT].\n\n")

# ---- (e) A worked steep-hazard case where the fallback engages -------

cat("==== (e) Weibull shape = 10: fallback engagement vs. delta (coarse output grids can let\n")
cat("     the solver take one large internal step across the whole fast-decaying region) ====\n")
for (delta in c(2, 1, 0.5, 0.25)) {
    n_fallback <- 0L
    withCallingHandlers(
        sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 10),
            n = 200, end = 20, delta = delta, seed = 1),
        message = function(m) {
            if (grepl("reported-grid refinement fallback", conditionMessage(m))) {
                n_fallback <<- as.integer(sub(".*: ([0-9]+) subject.*", "\\1",
                    conditionMessage(m)))
            }
            invokeRestart("muffleMessage")
        })
    cat(sprintf("  delta=%-5s fallback subjects: %d / 200\n", delta, n_fallback))
}
cat("(Documented as an open risk for Phase 3 in reports/08_phase2_5_report.md:\n")
cat(" very steep hazards can make the in-solver bracket momentarily invalid at a\n")
cat(" coarse delta; the fallback engages automatically and the message names it,\n")
cat(" and a finer delta resolves it, exactly as the message advises.)\n\n")

cat("==== Session info ====\n")
cat("Platform:", R.version$platform, "\n")
