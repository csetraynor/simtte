## Phase 2 validation: sim_tte_ode() Weibull and Gompertz library models.
##
## Records, against the shipped package (not ad-hoc mcode() strings):
##   (a) risk R1 shape-sweep results (weibull_ode.cpp, shape 0.01-10);
##   (b) raw-vs-refined accuracy for Weibull at shapes 0.5, 1, 2, 5, and
##       Gompertz at two gamma values;
##   (c) Kaplan-Meier overlays vs. sim_tte(type = "weibull");
##   (d) the rtol/atol sweep (Phase 1 report section 4) repeated for
##       Weibull shape = 2, to check whether that finding still holds
##       when the hazard is not constant.
##
## Excluded from the build (^inst/validation$ .Rbuildignore rule).

suppressMessages(library(simtte))
suppressMessages(library(dplyr))

cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")), "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n\n")

# ---- (a) Risk R1: Weibull-ODE shape sweep ---------------------------

cat("==== (a) Risk R1: weibull_ode.cpp shape sweep (0.01-10) ====\n")
shapes <- c(0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 5, 10)
mu <- -1
n <- 500
end <- 20

r1 <- lapply(shapes, function(shape) {
    sim <- sim_tte_ode(model = "weibull", param = list(mu = mu, shape = shape),
        n = n, end = end, delta = 2, keep_trajectory = TRUE, seed = 1)
    traj <- sim$trajectory
    n_nan <- sum(!is.finite(traj$p11))
    by_id <- split(traj$p11, traj$ID)
    n_nonmono <- sum(vapply(by_id, function(p) any(diff(p) > 1e-8), logical(1)))
    eta <- exp(mu)
    analytic <- exp(-eta * traj$time^shape)
    max_err <- max(abs(traj$p11 - analytic))
    data.frame(shape = shape, n_nan = n_nan, n_nonmonotone = n_nonmono,
        max_abs_err_p11 = max_err,
        status = if (n_nan > 0 || n_nonmono > 0) "NUMERICAL_ISSUE" else "OK")
})
r1 <- do.call(rbind, r1)
print(r1, row.names = FALSE)
stopifnot(all(r1$status == "OK"))
cat("All shapes 0.01-10: no error, no NaN, no non-monotone p11.\n\n")

# ---- (b) Raw vs. refined accuracy: Weibull and Gompertz -------------

cat("==== (b) Raw vs. refined accuracy ====\n")
n <- 3000
end <- 20

check_raw_vs_refined <- function(model, param, analytic_T_fn, label) {
    sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
        delta = 4, keep_trajectory = TRUE, seed = 20260910)
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    u_by_id <- last$U[match(seq_len(n), last$ID)]
    # For gamma < 0 (a decreasing Gompertz hazard), cumulative hazard is
    # bounded as t -> Inf, so the analytic inverse is NaN for any U
    # below that bound (mathematically: no finite event time exists --
    # censored at infinity). Benign; `sim_status == 1L & NA` short-
    # circuits to FALSE in R (never selects those rows), so this never
    # affects the comparison below.
    analytic_T <- suppressWarnings(analytic_T_fn(u_by_id))
    is_event <- sim$events$sim_status == 1L & analytic_T < end
    raw_tevt <- last$TEVT[match(sim$events$ID, last$ID)]
    err_raw <- abs(raw_tevt[is_event] - analytic_T[is_event])
    err_refined <- abs(sim$events$sim_time[is_event] - analytic_T[is_event])
    cat(label, "-- n compared:", sum(is_event), "\n")
    cat(label, "-- raw      mean/max abs err:", round(mean(err_raw), 5), "/",
        round(max(err_raw), 5), "\n")
    cat(label, "-- refined  mean/max abs err:", round(mean(err_refined), 5), "/",
        round(max(err_refined), 5), "\n\n")
}

check_raw_vs_refined("weibull", list(mu = -1, shape = 0.5),
    function(u) (-log(u) / exp(-1))^(1 / 0.5), "Weibull shape=0.5")
check_raw_vs_refined("weibull", list(mu = -1, shape = 1),
    function(u) -log(u) / exp(-1), "Weibull shape=1")
check_raw_vs_refined("weibull", list(mu = -1, shape = 2),
    function(u) (-log(u) / exp(-1))^(1 / 2), "Weibull shape=2")
check_raw_vs_refined("weibull", list(mu = -1, shape = 5),
    function(u) (-log(u) / exp(-1))^(1 / 5), "Weibull shape=5")
check_raw_vs_refined("gompertz", list(mu = -2, gamma = 0.1),
    function(u) log(1 - (log(u) * 0.1) / exp(-2)) / 0.1, "Gompertz gamma=0.1")
check_raw_vs_refined("gompertz", list(mu = -2, gamma = -0.05),
    function(u) log(1 - (log(u) * -0.05) / exp(-2)) / -0.05, "Gompertz gamma=-0.05")

# ---- (c) Kaplan-Meier overlay vs. sim_tte(type = "weibull") ---------

cat("==== (c) Kaplan-Meier overlay: sim_tte_ode() vs. sim_tte() ====\n")
mu <- -1; shape <- 1.8; n <- 3000; end <- 15; seed <- 5050

set.seed(seed)
lp <- matrix(rep(0, n), nrow = n)
ref_grid <- sim_tte(pi = lp, mu = mu, coefs = shape,
    time = seq(0.05, end, by = 0.05), type = "weibull", end_time = end,
    event_time_method = "grid")
set.seed(seed)
ref_log <- sim_tte(pi = lp, mu = mu, coefs = shape,
    time = seq(0.05, end, by = 0.05), type = "weibull", end_time = end,
    event_time_method = "log_survival")
sim_ode <- sim_tte_ode(model = "weibull", param = list(mu = mu, shape = shape),
    n = n, end = end, delta = 0.25, seed = seed)

km_at <- function(events, t) {
    # Simple KM survival estimate via product-limit over observed
    # distinct event times <= t (no ties adjustment needed at this n).
    ord <- order(events$sim_time)
    et <- events$sim_time[ord]
    st <- events$sim_status[ord]
    at_risk <- length(et):1
    surv <- cumprod(ifelse(st == 1, 1 - 1 / at_risk, 1))
    idx <- findInterval(t, et)
    ifelse(idx == 0, 1, surv[idx])
}

eta <- exp(mu)
analytic_S <- function(t) exp(-eta * t^shape)
checkpoints <- seq(1, 14, by = 1)
km_tab <- data.frame(t = checkpoints,
    analytic = analytic_S(checkpoints),
    grid = km_at(ref_grid, checkpoints),
    log_survival = km_at(ref_log, checkpoints),
    ode = km_at(sim_ode$events, checkpoints))
print(km_tab, row.names = FALSE)
cat("Max abs diff vs analytic -- grid:",
    round(max(abs(km_tab$grid - km_tab$analytic)), 4),
    " log_survival:", round(max(abs(km_tab$log_survival - km_tab$analytic)), 4),
    " sim_tte_ode:", round(max(abs(km_tab$ode - km_tab$analytic)), 4), "\n\n")

# ---- (d) rtol/atol sweep repeated for Weibull shape = 2 -------------

cat("==== (d) rtol/atol sensitivity, Weibull shape = 2 (cf. Phase 1 report section 4) ====\n")
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
cat("Platform:", R.version$platform, "\n")
