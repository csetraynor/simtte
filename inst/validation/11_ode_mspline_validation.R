## Phase 3 validation: sim_tte_ode(model = "mspline").
## reports/09_phase3_report.md.
##
## Records, against the shipped package (not ad-hoc mcode() strings):
##   (a) analytical agreement per shipped knot-count variant (3, 5, 7),
##       a smooth and a sharply peaked coefficient vector each;
##   (b) the Phase 2.5 three-column table (raw / old grid refinement /
##       new in-solver refinement) at delta = 4, 1, 0.25;
##   (c) the steepness stress test (open risk 1 from Phase 2.5): fallback
##       rate across delta AND across rtol/atol, for a sharply peaked
##       coefficient vector;
##   (d) cross-method comparison vs sim_tte(type = "ms",
##       event_time_method = "log_survival") at three grid resolutions;
##   (e) boundary guard, bracket containment;
##   (f) rtol/atol sweep for one configuration.
##
## Excluded from the build (^inst/validation$ .Rbuildignore rule). Run
## from the package root: Rscript inst/validation/11_ode_mspline_validation.R

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root: ",
        "Rscript inst/validation/11_ode_mspline_validation.R", call. = FALSE)
}
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages(library(dplyr))
suppressMessages(library(splines2))

cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")), "\n")
cat("splines2 version:", as.character(utils::packageVersion("splines2")), "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n\n")

mu <- -0.5
DEGREE <- simtte:::.MSPLINE_DEGREE

variants <- list(
    list(K = 3, knots = c(5, 10, 15)),
    list(K = 5, knots = c(3.33, 6.67, 10, 13.33, 16.67)),
    list(K = 7, knots = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5))
)
bk <- c(0, 20)

analytic_S <- function(t, knots, coefs) {
    eta <- exp(mu)
    I <- iSpline(t, knots = knots, Boundary.knots = bk, degree = DEGREE,
        intercept = TRUE)
    exp(-eta * as.numeric(I %*% coefs))
}

# ---- (a) Analytical agreement, smooth and peaked coefficients -------

cat("==== (a) Analytical agreement per knot-count variant ====\n")
n <- 3000
probs <- c(0.2, 0.4, 0.6, 0.8)

check_analytical <- function(K, knots, coefs, label) {
    M <- length(coefs)
    sim <- sim_tte_ode(model = "mspline", knots = knots, coefs = coefs,
        boundary_knots = bk, param = list(mu = mu), n = n, end = bk[2],
        delta = 1, seed = 20260910)
    # Quantile-based checkpoints (t_j solved so 1 - S(t_j) = prob
    # exactly), same rationale as the Weibull tests: avoids saturated,
    # uninformative checkpoints for a peaked hazard.
    t_grid <- seq(0, bk[2], length.out = 20001)
    S_grid <- analytic_S(t_grid, knots, coefs)
    rows <- lapply(probs, function(p) {
        idx <- which.min(abs((1 - S_grid) - p))
        t_j <- t_grid[idx]
        p_emp <- mean(sim$events$sim_time <= t_j & sim$events$sim_status == 1)
        data.frame(label = label, prob_target = p, t = round(t_j, 3),
            empirical = round(p_emp, 4), analytic = round(p, 4),
            abs_diff = round(abs(p_emp - p), 4),
            tol_4se = round(4 * sqrt(p * (1 - p) / n), 4))
    })
    do.call(rbind, rows)
}

a_results <- do.call(rbind, lapply(variants, function(v) {
    M <- v$K + DEGREE + 1
    smooth <- rep(1, M)
    peaked <- rep(0.05, M)
    peaked[ceiling(M / 2)] <- 8
    rbind(
        check_analytical(v$K, v$knots, smooth, paste0("K=", v$K, " smooth")),
        check_analytical(v$K, v$knots, peaked, paste0("K=", v$K, " peaked"))
    )
}))
print(a_results, row.names = FALSE)
cat("\n")

# ---- (b) Three-column table: raw / old grid refinement / new in-solver

cat("==== (b) raw vs. old (grid) refinement vs. new (in-solver) refinement ====\n")

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

check_three_way <- function(K, knots, coefs, label, delta, n = 3000) {
    sim <- suppressMessages(sim_tte_ode(model = "mspline", knots = knots,
        coefs = coefs, boundary_knots = bk, param = list(mu = mu), n = n,
        end = bk[2], delta = delta, keep_trajectory = TRUE, seed = 20260910))
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    u_by_id <- last$U[match(seq_len(n), last$ID)]
    t_grid <- seq(0, bk[2], length.out = 200001)
    S_grid <- analytic_S(t_grid, knots, coefs)
    analytic_T <- approx(x = rev(S_grid), y = rev(t_grid), xout = u_by_id,
        ties = "ordered")$y
    is_event <- sim$events$sim_status == 1L & !is.na(analytic_T) &
        analytic_T < bk[2]
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

b_results <- do.call(rbind, lapply(variants, function(v) {
    M <- v$K + DEGREE + 1
    smooth <- rep(1, M)
    do.call(rbind, lapply(c(4, 1, 0.25), function(d) {
        check_three_way(v$K, v$knots, smooth, paste0("K=", v$K, " smooth"), d)
    }))
}))
print(b_results, row.names = FALSE, digits = 4)
cat("\n")

# ---- (c) Steepness stress test: fallback rate vs. delta AND rtol/atol

cat("==== (c) Steepness stress test (open risk 1, reports/08_phase2_5_report.md section 4) ====\n")
knots7 <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5)
peak_coefs <- rep(0.01, 10)
peak_coefs[6] <- 5

fallback_count <- function(...) {
    n_fb <- 0L
    withCallingHandlers({
        sim_tte_ode(model = "mspline", knots = knots7, coefs = peak_coefs,
            boundary_knots = bk, param = list(mu = -1), n = 200, end = bk[2],
            seed = 1, ...)
    }, message = function(m) {
        txt <- conditionMessage(m)
        if (grepl("fallback", txt)) {
            n_fb <<- as.integer(regmatches(txt,
                regexpr("[0-9]+(?= subject)", txt, perl = TRUE)))
        }
        invokeRestart("muffleMessage")
    })
    n_fb
}

cat("-- fallback vs. delta (rtol/atol at mrgsolve default) --\n")
for (d in c(2, 1, 0.5, 0.25, 0.1, 0.05)) {
    cat(sprintf("  delta=%-5s fallback: %d / 200\n", d, fallback_count(delta = d)))
}
cat("-- fallback vs. rtol/atol (delta = 0.25 fixed) --\n")
for (tol in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12)) {
    cat(sprintf("  rtol=atol=%-8s fallback: %d / 200\n", tol,
        fallback_count(delta = 0.25, rtol = tol, atol = tol)))
}
cat("\nInterpretation: unlike Weibull shape >= 5 (Phase 2.5 section 4, where a\n")
cat("finer delta alone resolved the fallback), this M-spline spike's fallback\n")
cat("does not resolve with delta at all (flat across delta = 2 down to 0.05),\n")
cat("but responds (non-monotonically) to rtol/atol -- a same-timestamp,\n")
cat("multiple-corrector-iteration degeneracy, not a step-size/delta issue. See\n")
cat("reports/09_phase3_report.md 'Open risks' for the full explanation.\n\n")

# ---- (d) Cross-method comparison vs sim_tte(type = "ms") ------------

cat("==== (d) Cross-method comparison vs sim_tte(type = 'ms', event_time_method = 'log_survival') ====\n")
knots3 <- c(5, 10, 15)
coefs3 <- rep(1, 6)
n <- 3000
seed <- 5151

km_at <- function(events, t) {
    ord <- order(events$sim_time)
    et <- events$sim_time[ord]
    st <- events$sim_status[ord]
    at_risk <- length(et):1
    surv <- cumprod(ifelse(st == 1, 1 - 1 / at_risk, 1))
    idx <- findInterval(t, et)
    ifelse(idx == 0, 1, surv[idx])
}

checkpoints <- seq(1, 19, by = 2)
analytic_chk <- analytic_S(checkpoints, knots3, coefs3)

lp <- matrix(rep(0, n), nrow = n)
grid_resolutions <- c(coarse = 2, medium = 0.5, fine = 0.1)
ms_tab <- data.frame(t = checkpoints, analytic = round(analytic_chk, 4))
for (nm in names(grid_resolutions)) {
    time_grid <- seq(0, bk[2], by = grid_resolutions[[nm]])
    basis <- mSpline(time_grid, knots = knots3, Boundary.knots = bk,
        degree = DEGREE, intercept = TRUE)
    set.seed(seed)
    ref <- sim_tte(pi = lp, mu = mu, basis = basis, coefs = coefs3,
        time = time_grid, type = "ms", end_time = bk[2],
        event_time_method = "log_survival")
    ms_tab[[nm]] <- round(km_at(ref, checkpoints), 4)
}
set.seed(seed)
sim_ode <- sim_tte_ode(model = "mspline", knots = knots3, coefs = coefs3,
    boundary_knots = bk, param = list(mu = mu), n = n, end = bk[2],
    delta = 1, seed = seed)
ms_tab$ode <- round(km_at(sim_ode$events, checkpoints), 4)
print(ms_tab, row.names = FALSE)

cat("\nMax abs diff vs analytic:\n")
for (nm in c(names(grid_resolutions), "ode")) {
    cat(sprintf("  %-8s %.4f\n", nm,
        max(abs(ms_tab[[nm]] - ms_tab$analytic))))
}
cat("\n")

# ---- (e) Boundary guard and bracket containment ----------------------

cat("==== (e) Boundary guard and bracket containment ====\n")
boundary_check <- function(K, knots, n, seed, end = 18, boundary_knots = c(0, 20)) {
    M <- K + DEGREE + 1
    coefs <- rep(1, M)
    set.seed(seed)
    data <- data.frame(ID = seq_len(n), time = end - 0.05, lp = 0.5,
        evid = 1, amt = 0, cmt = 1)
    sim <- sim_tte_ode(model = "mspline", knots = knots, coefs = coefs,
        boundary_knots = boundary_knots, param = list(mu = mu), n = n,
        end = end, delta = 1, data = data, seed = seed)
    max(sim$events$sim_time) <= end + 1e-9
}
bg <- data.frame(K = rep(c(3, 5, 7), each = 2), n = rep(c(40, 2000), 3),
    ok = c(
        boundary_check(3, c(5, 10, 15), 40, 8001),
        boundary_check(3, c(5, 10, 15), 2000, 8002),
        boundary_check(5, c(3.33, 6.67, 10, 13.33, 16.67), 40, 8003),
        boundary_check(5, c(3.33, 6.67, 10, 13.33, 16.67), 2000, 8004),
        boundary_check(7, c(2.5, 5, 7.5, 10, 12.5, 15, 17.5), 40, 8005),
        boundary_check(7, c(2.5, 5, 7.5, 10, 12.5, 15, 17.5), 2000, 8006)
    ))
print(bg, row.names = FALSE)
stopifnot(all(bg$ok))
cat("All boundary-guard checks passed.\n\n")

bracket_check <- function(K, knots, n = 1000, seed = 1) {
    M <- K + DEGREE + 1
    coefs <- rep(1, M)
    sim <- sim_tte_ode(model = "mspline", knots = knots, coefs = coefs,
        boundary_knots = bk, param = list(mu = mu), n = n, end = bk[2],
        delta = 2, keep_trajectory = TRUE, seed = seed)
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(sim$events$ID, last$ID), ]
    is_event <- sim$events$sim_status == 1L
    lower_ok <- all(sim$events$sim_time[is_event] >= last$T_PRE[is_event] - 1e-8)
    upper_ok <- all(sim$events$sim_time[is_event] <= last$TEVT[is_event] + 1e-8)
    data.frame(K = K, n_event = sum(is_event), lower_bound_ok = lower_ok,
        upper_bound_ok = upper_ok)
}
bc <- do.call(rbind, lapply(variants, function(v) bracket_check(v$K, v$knots)))
print(bc, row.names = FALSE)
stopifnot(all(bc$lower_bound_ok), all(bc$upper_bound_ok))
cat("All refined sim_time values stayed within [T_PRE, TEVT].\n\n")

# ---- (f) rtol/atol sweep for one configuration ------------------------

cat("==== (f) rtol/atol sensitivity, K=3 smooth, delta = 4 ====\n")
n <- 3000
run_with_tol <- function(tol) {
    s <- sim_tte_ode(model = "mspline", knots = knots3, coefs = coefs3,
        boundary_knots = bk, param = list(mu = mu), n = n, end = bk[2],
        delta = 4, keep_trajectory = TRUE, seed = 20260910, rtol = tol,
        atol = tol)
    last <- s$trajectory[!duplicated(s$trajectory$ID, fromLast = TRUE), ]
    u_i <- last$U[match(seq_len(n), last$ID)]
    t_grid <- seq(0, bk[2], length.out = 200001)
    S_grid <- analytic_S(t_grid, knots3, coefs3)
    t_true <- approx(x = rev(S_grid), y = rev(t_grid), xout = u_i,
        ties = "ordered")$y
    ok <- s$events$sim_status == 1L & !is.na(t_true) & t_true < bk[2]
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

cat("==== Session info ====\n")
cat("Platform:", R.version$platform, "\n")
