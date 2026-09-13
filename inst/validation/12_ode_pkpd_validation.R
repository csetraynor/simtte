## Phase 4 validation: sim_tte_ode() PK/PD-linked hazard library
## (pk_hazard, irm1_hazard-irm4_hazard, tmdd_hazard). See
## reports/10_phase4_report.md.
##
## No closed form exists for a PK/PD-linked hazard in general, so the
## reference throughout is numerical: a very fine (delta = 0.05)
## mrgsim() trajectory of the SAME model/dose, read off at the SAME
## per-subject U draw (design report section 6 item 4) -- not the raw
## sim_tte_df() quantile comparison alone, which uses its own
## independently-drawn U and so cannot give a per-subject "true" time.
##
## Records:
##   (a) KM-style quantile comparison (median, IQR, 90th) vs the
##       fine-grid reference, plus mean-p11-vs-empirical-event-rate
##       checkpoints, per model, with a realistic dosing regimen
##       (repeated oral dosing; one IV bolus case for pk_hazard);
##   (b) mechanistic direction checks (dose/beta) and the beta=0/beta_r=0
##       -> exponential analytical check;
##   (c) BSV propagation with omega, and the disclosed limitation found
##       while building this (none of the six backbones wires ETA(n)
##       into any parameter);
##   (d) three-column table (raw / old grid-fallback / new in-solver
##       refinement) at delta = 4, 1, 0.25 vs the fine-grid reference,
##       plus fallback-classification counts per model at delta = 4 --
##       the steepness characterization Phase 3's open risk 1 asked for;
##   (e) boundary guard at 40/2000 subjects, bracket containment, one
##       rtol/atol sweep (tmdd_hazard, the intentionally stiff model);
##   (f) timing: cold compile time per model (already reported in
##       reports/10_phase4_report.md from isolated fresh-process runs,
##       not repeated here to keep this script's own runtime down).
##
## Excluded from the build (^inst/validation$ .Rbuildignore rule). Run
## from the package root: Rscript inst/validation/12_ode_pkpd_validation.R

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root: ",
        "Rscript inst/validation/12_ode_pkpd_validation.R", call. = FALSE)
}
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages(library(testthat))
source("tests/testthat/helper-ode-models.R")

cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")), "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n\n")

MODELS <- .PKPD_MODELS <- c("pk_hazard", "irm1_hazard", "irm2_hazard",
    "irm3_hazard", "irm4_hazard", "tmdd_hazard")
PARAM <- list(
    pk_hazard = list(H0 = 0.01, beta_cp = 0.02),
    irm1_hazard = list(H0 = 0.01, beta_r = 1),
    irm2_hazard = list(H0 = 0.01, beta_r = 1),
    irm3_hazard = list(H0 = 0.01, beta_r = 1),
    irm4_hazard = list(H0 = 0.01, beta_r = 1),
    tmdd_hazard = list(H0 = 0.02, beta_rc = 0.3)
)
END <- 30

oral_dose <- function(n, amt = 100, ii = 12, addl = 9) {
    data.frame(ID = seq_len(n), time = 0, cmt = 1, amt = amt, evid = 1,
        ii = ii, addl = addl)
}
iv_dose <- function(n, amt = 100) {
    data.frame(ID = seq_len(n), time = 0, cmt = 2, amt = amt, evid = 1)
}

# Generic "first crossing" lookup against ANY trajectory (fine or
# coarse) -- the same logic .resolve_ode_events()'s own reported-grid
# fallback uses, factored out here so it can be applied to a very fine
# reference trajectory too (giving a near-exact per-subject "true"
# event time without a closed form).
grid_event_time <- function(traj, id, u_i) {
    sub <- traj[traj$ID == id, ]
    sub <- sub[order(sub$time), ]
    p11 <- pmin(pmax(sub$p11, 0), 1)
    idx <- match(TRUE, p11 <= u_i)
    if (is.na(idx)) return(NA_real_)
    if (idx == 1L) return(sub$time[1])
    i <- idx - 1L
    simtte:::.interpolate_log_survival(t_i = sub$time[i], t_ip1 = sub$time[idx],
        s_i = p11[i], s_ip1 = p11[idx], u = u_i)
}

FINE_DELTA <- 0.05  # 20x finer than the coarsest delta used below; keeps
                     # this script's own runtime reasonable while still
                     # giving a near-exact per-subject reference (an
                     # initial delta = 0.01 attempt was computationally
                     # far too heavy at n = 2000-3000 -- see report).
fine_reference <- function(model, param, dose, n, seed) {
    suppressMessages(sim_tte_ode(model = model, param = param, n = n,
        end = END, delta = FINE_DELTA, data = dose, keep_trajectory = TRUE,
        seed = seed))
}

cat("==== (a) Quantile comparison vs. fine-grid (delta =", FINE_DELTA, ") reference ====\n")
N_A <- 1200
a_results <- do.call(rbind, lapply(MODELS, function(m) {
    dose <- oral_dose(N_A)
    fine <- fine_reference(m, PARAM[[m]], dose, n = N_A, seed = 909)
    sim <- suppressMessages(sim_tte_ode(model = m, param = PARAM[[m]],
        n = N_A, end = END, delta = 1, data = dose, seed = 909))
    ev <- sim$events
    fe <- fine$events[match(ev$ID, fine$events$ID), ]
    qs <- c(0.25, 0.5, 0.75, 0.9)
    q_sim <- stats::quantile(ev$sim_time[ev$sim_status == 1], qs)
    q_fine <- stats::quantile(fe$sim_time[fe$sim_status == 1], qs)
    data.frame(model = m, metric = c(paste0("q", qs * 100), "event_rate"),
        sim_tte_ode = c(q_sim, mean(ev$sim_status == 1)),
        fine_grid_ref = c(q_fine, mean(fe$sim_status == 1)))
}))
a_results$abs_diff <- abs(a_results$sim_tte_ode - a_results$fine_grid_ref)
print(a_results, row.names = FALSE, digits = 4)
cat("\n")

cat("---- mean-p11-vs-empirical-event-rate checkpoints (pk_hazard) ----\n")
{
    m <- "pk_hazard"
    sim <- suppressMessages(sim_tte_ode(model = m, param = PARAM[[m]],
        n = N_A, end = END, delta = 1, data = oral_dose(N_A),
        keep_trajectory = TRUE, seed = 909))
    traj <- sim$trajectory
    checkpoints <- c(5, 10, 15, 20, 25)
    for (t_j in checkpoints) {
        rows <- traj[traj$time == t_j, ]
        mean_p11 <- mean(rows$p11)
        p_event_empirical <- mean(sim$events$sim_time <= t_j &
            sim$events$sim_status == 1)
        cat(sprintf("t=%2d  mean(p11)=%.4f  1-mean(p11)=%.4f  emp.event.rate=%.4f  diff=%.4f\n",
            t_j, mean_p11, 1 - mean_p11, p_event_empirical,
            abs((1 - mean_p11) - p_event_empirical)))
    }
}
cat("\n")

cat("==== (a-IV) pk_hazard IV bolus case ====\n")
{
    m <- "pk_hazard"
    dose_iv <- iv_dose(N_A)
    fine <- fine_reference(m, PARAM[[m]], dose_iv, n = N_A, seed = 707)
    sim <- suppressMessages(sim_tte_ode(model = m, param = PARAM[[m]],
        n = N_A, end = END, delta = 1, data = dose_iv, seed = 707))
    p_event_sim <- mean(sim$events$sim_status == 1)
    p_event_fine <- mean(fine$events$sim_status == 1)
    cat(sprintf("event rate sim=%.4f fine=%.4f abs.diff=%.4f\n", p_event_sim,
        p_event_fine, abs(p_event_sim - p_event_fine)))
}
cat("\n")

cat("==== (b) Mechanistic direction checks ====\n")
{
    n <- 1200
    dose_lo <- oral_dose(n, amt = 20)
    dose_hi <- oral_dose(n, amt = 200)
    sim_lo <- suppressMessages(sim_tte_ode(model = "pk_hazard",
        param = PARAM$pk_hazard, n = n, end = END, delta = 1,
        data = dose_lo, seed = 1))
    sim_hi <- suppressMessages(sim_tte_ode(model = "pk_hazard",
        param = PARAM$pk_hazard, n = n, end = END, delta = 1,
        data = dose_hi, seed = 1))
    cat(sprintf("pk_hazard: p(event) low dose=%.4f  high dose=%.4f  (higher dose -> higher event rate: %s)\n",
        mean(sim_lo$events$sim_status == 1), mean(sim_hi$events$sim_status == 1),
        mean(sim_hi$events$sim_status == 1) > mean(sim_lo$events$sim_status == 1)))

    dose <- oral_dose(n)
    sim_harm <- suppressMessages(sim_tte_ode(model = "irm3_hazard",
        param = list(H0 = 0.01, beta_r = 3), n = n, end = END, delta = 1,
        data = dose, seed = 1))
    sim_prot <- suppressMessages(sim_tte_ode(model = "irm3_hazard",
        param = list(H0 = 0.01, beta_r = -3), n = n, end = END, delta = 1,
        data = dose, seed = 1))
    cat(sprintf("irm3_hazard: p(event) beta_r=+3: %.4f  beta_r=-3: %.4f  (sign flips direction: %s)\n",
        mean(sim_harm$events$sim_status == 1), mean(sim_prot$events$sim_status == 1),
        mean(sim_harm$events$sim_status == 1) > mean(sim_prot$events$sim_status == 1)))

    # beta = 0 -> exponential with hazard H0 (analytical: p(event by end)
    # = 1 - exp(-H0*end)).
    sim0 <- suppressMessages(sim_tte_ode(model = "pk_hazard",
        param = list(H0 = 0.02, beta_cp = 0), n = n, end = END, delta = 1,
        data = dose, seed = 1))
    p_analytic <- 1 - exp(-0.02 * END)
    p_emp <- mean(sim0$events$sim_status == 1)
    cat(sprintf("pk_hazard beta_cp=0: p(event) empirical=%.4f analytical(H0=0.02)=%.4f abs.diff=%.4f\n",
        p_emp, p_analytic, abs(p_emp - p_analytic)))
}
cat("\n")

cat("==== (c) BSV propagation (omega) ====\n")
{
    om <- matrix(c(0.3, 0, 0, 0.3), 2, 2)
    res <- tryCatch({
        sim_tte_ode(model = "pk_hazard", param = PARAM$pk_hazard,
            omega = om, n = 500, end = END, delta = 1,
            data = oral_dose(500), seed = 1)
        "no error (unexpected)"
    }, error = function(e) conditionMessage(e))
    cat("sim_tte_ode(model = 'pk_hazard', omega = <2x2 matrix>) ->\n  ", res, "\n")
    cat("FINDING: none of the six modlib() backbones wires ETA(n) into\n",
        "any parameter or declares an OMEGA block, so 'omega=' cannot add\n",
        "BSV to any of the six as shipped -- confirmed above via the\n",
        "informative error sim_tte_ode() now raises instead of mrgsolve's\n",
        "own cryptic one (R/helpers.R .apply_ode_matlist()). See\n",
        "reports/10_phase4_report.md open risks.\n", sep = "")
}
cat("\n")

cat("==== (d) Three-column table (raw / old grid fallback / new in-solver) + fallback classification ====\n")
N_D <- 800
# One fine-grid reference per model, reused across all three delta
# values below (the reference does not itself depend on the coarse
# run's delta) -- computing it per (model, delta) combination instead
# was needlessly 3x the cost for no extra information.
check_three_way <- function(model, param, dose, delta, ref_T, n = N_D,
    seed = 20260911) {
    sim <- suppressWarnings(suppressMessages(sim_tte_ode(model = model,
        param = param, n = n, end = END, delta = delta, data = dose,
        keep_trajectory = TRUE, seed = seed)))
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(seq_len(n), last$ID), ]

    is_event <- sim$events$sim_status == 1L & !is.na(ref_T) & ref_T < END
    raw_tevt <- last$TEVT
    old_refined <- vapply(seq_len(n), function(id)
        grid_event_time(traj, id, last$U[id]), numeric(1))

    err_raw <- abs(raw_tevt[is_event] - ref_T[is_event])
    err_old <- abs(old_refined[is_event] - ref_T[is_event])
    err_new <- abs(sim$events$sim_time[is_event] - ref_T[is_event])

    sig <- rep(NA_character_, n)
    for (id in which(is_event)) {
        refined <- simtte:::.refine_ode_event_time_insolver(
            t_pre = last$T_PRE[id], p_pre = last$P_PRE[id],
            tevt = last$TEVT[id], p_post = last$P_POST[id], u = last$U[id])
        if (is.na(refined)) {
            sig[id] <- simtte:::.classify_ode_fallback(t_pre = last$T_PRE[id],
                p_pre = last$P_PRE[id], tevt = last$TEVT[id],
                p_post = last$P_POST[id])
        }
    }
    n_fallback <- sum(!is.na(sig))
    list(row = data.frame(model = model, delta = delta,
            n_compared = sum(is_event), n_fallback = n_fallback,
            raw_mean_err = mean(err_raw), old_refined_mean_err = mean(err_old),
            new_refined_mean_err = mean(err_new)),
        sig_table = table(sig[which(is_event)], useNA = "no"))
}

d_rows <- list()
d_sig <- list()
for (m in MODELS) {
    dose_d <- oral_dose(N_D)
    fine <- fine_reference(m, PARAM[[m]], dose_d, n = N_D, seed = 20260911)
    fine_traj <- fine$trajectory
    u_by_id <- fine_traj$U[!duplicated(fine_traj$ID, fromLast = TRUE)]
    ref_T <- vapply(seq_len(N_D), function(id)
        grid_event_time(fine_traj, id, u_by_id[id]), numeric(1))
    for (delta in c(4, 1, 0.25)) {
        out <- check_three_way(m, PARAM[[m]], dose_d, delta, ref_T)
        d_rows[[length(d_rows) + 1L]] <- out$row
        if (delta == 4) d_sig[[m]] <- out$sig_table
    }
}
d_results <- do.call(rbind, d_rows)
print(d_results, row.names = FALSE, digits = 4)
cat("\nFallback signature counts at delta = 4 (steepness characterization):\n")
for (m in MODELS) {
    cat(m, ": ", if (length(d_sig[[m]])) paste(names(d_sig[[m]]), d_sig[[m]],
        sep = "=", collapse = ", ") else "(none)", "\n", sep = "")
}
cat("\n")

cat("==== (e) Boundary guard, bracket containment, rtol/atol sweep (tmdd_hazard) ====\n")
for (n_bg in c(40, 2000)) {
    events <- .run_ode_boundary_guard_check("tmdd_hazard", PARAM$tmdd_hazard,
        n = n_bg, seed = 5555, end = 10)
    cat(sprintf("boundary guard n=%d: max(sim_time)=%.6f (<= 10? %s)\n",
        n_bg, max(events$sim_time), all(events$sim_time <= 10 + 1e-9)))
}
{
    sim <- suppressMessages(sim_tte_ode(model = "tmdd_hazard",
        param = PARAM$tmdd_hazard, n = 500, end = END, delta = 2,
        data = oral_dose(500), keep_trajectory = TRUE, seed = 1))
    traj <- sim$trajectory
    last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
    last <- last[match(sim$events$ID, last$ID), ]
    is_event <- sim$events$sim_status == 1L
    contained <- all(sim$events$sim_time[is_event] >= last$T_PRE[is_event] - 1e-8) &&
        all(sim$events$sim_time[is_event] <= last$TEVT[is_event] + 1e-8)
    cat("bracket containment (tmdd_hazard, n=500):", contained, "\n")
}
cat("\n---- rtol/atol sweep (tmdd_hazard, delta = 4, n = 500) ----\n")
for (tol in c(1e-6, 1e-8, 1e-10)) {
    msgs <- testthat::capture_messages(sim <- suppressWarnings(sim_tte_ode(
        model = "tmdd_hazard", param = PARAM$tmdd_hazard, n = 500, end = END,
        delta = 4, data = oral_dose(500), seed = 1, rtol = tol, atol = tol)))
    # Filter for the fallback message specifically -- capture_messages()
    # also picks up mrgsolve's own "Loading model from cache."/build
    # progress message() calls, which are not what this sweep measures.
    fb_msg <- grep("used the reported-grid", msgs, value = TRUE)
    n_fb <- if (length(fb_msg)) as.integer(sub("^sim_tte_ode\\(\\): ([0-9]+).*",
        "\\1", fb_msg[1])) else 0L
    cat(sprintf("rtol=atol=%-8g -> n_fallback = %d / 500\n", tol, n_fb))
}
cat("\n")

cat("==== Done ====\n")
