## Phase H validation: pkpd_idr_hazard, the shipped indirect-response
## PK/PD mechanistic example model
## (inst/models/examples/pkpd_idr_hazard.cpp).
##
## Unlike pkpd_linear_hazard (validated in
## 03_pkpd_mechanism_validation.R against an independent closed-form
## S(t)), this model's hazard is driven by a nonlinear indirect-response
## (turnover) mediator and has no closed-form cumulative hazard. Its
## validation therefore follows the mechanistic-direction and
## self-consistency approach the package manuscript itself proposes for
## this exact model ("ODE-coupled PK/PD hazard" section):
##
##   (a) p11 trajectory validity (the same acceptance criterion
##       sim_tte_df() itself already applies to any input trajectory);
##   (b) a Kaplan-Meier vs. mean-p11 self-consistency check, reusing the
##       existing 4-SE binomial tolerance already established in
##       03_pkpd_mechanism_validation.R's Check 2 (see the note at
##       Check (b) below for why that reuse is valid here -- no new
##       numerical criterion is introduced);
##   (c) the manuscript's own IC50 potency sweep, reproduced as a
##       programmatic monotonicity assertion rather than the prose
##       numbers currently only asserted informally in the manuscript
##       text;
##   (d) a between-subject-variability sanity check ($OMEGA is active
##       and actually produces heterogeneous individual trajectories,
##       not silently inert);
##   (e) fixed-seed reproducibility.
##
## Reproducibility: mrgsolve draws $OMEGA random effects from R's own
## global RNG stream, so a single seed controls both those draws and the
## uniform draws consumed by sim_tte_df() (verified directly: two
## mrgsim() calls with the same seed give identical() output).
##
## obsonly = TRUE is required in every mrgsim() call whose output feeds
## sim_tte_df(): without it, mrgsim() also reports the internal
## dosing/bookkeeping row mrgsolve inserts at each dose time, duplicating
## time = 0 (one row before, one after the dose is applied) and tripping
## sim_tte_df()'s "no duplicated times within a subject" trajectory
## contract. (This is the same obsonly = TRUE convention
## .sim_surv_df() already fixes internally throughout the package; see
## R/helpers.R's .RESERVED_MRGSIM_ARGS.)

suppressMessages(library(simtte))
suppressMessages(library(mrgsolve))
suppressMessages(library(dplyr))

mod <- simtte_example_model("pkpd_idr_hazard")

n_subjects <- 2000
dose_regimen <- function(n) {
    mrgsolve::expand.ev(ID = seq_len(n), amt = 100, cmt = 1, ii = 1,
        addl = 23, time = 0)
}
data_cohort <- dose_regimen(n_subjects)

set.seed(20260830)
out <- as.data.frame(mrgsim(mod, data = data_cohort, end = 24, delta = 0.5,
    obsonly = TRUE))

## ---- Check (a): p11 trajectory validity ---------------------------------

cat("==== Check (a): p11 trajectory validity ====\n")
cat("All finite:", all(is.finite(out$p11)), "\n")
cat("All in [0, 1]:", all(out$p11 >= 0 & out$p11 <= 1), "\n")
starts_at_one <- tapply(out$p11, out$ID, function(p) p[1] == 1)
cat("All subjects start at p11 = 1:", all(starts_at_one), "\n")
nonincreasing <- tapply(out$p11, out$ID, function(p) all(diff(p) <= 1e-8))
cat("All subjects non-increasing (1e-8 solver tolerance):",
    all(nonincreasing), "\n")
no_dup_time <- tapply(out$time, out$ID, function(t) anyDuplicated(t) == 0L)
cat("No duplicated times within any subject:", all(no_dup_time), "\n")

## ---- sim_tte_df() smoke run ----------------------------------------------

sim_events <- sim_tte_df(out[, c("ID", "time", "p11")])
cat("\n==== sim_tte_df() smoke run ====\n")
cat("Rows == n_subjects:", nrow(sim_events) == n_subjects, "\n")
cat("sim_status only 0/1:", all(sim_events$sim_status %in% c(0, 1)), "\n")
cat("sim_time all finite:", all(is.finite(sim_events$sim_time)), "\n")
cat("Event rate:", round(mean(sim_events$sim_status), 4), "\n")

## ---- Check (b): Kaplan-Meier vs. mean-p11 self-consistency --------------
## Every subject shares the same administrative horizon (end = 24, no
## additional random censoring), so the risk set is exactly n_subjects
## at every checkpoint strictly before 24 -- the Kaplan-Meier estimator
## therefore has the same sampling behaviour as a simple empirical
## proportion at each checkpoint, which is exactly the setting the
## existing 4-SE binomial tolerance was designed for
## (03_pkpd_mechanism_validation.R, Check 2). Reused here unchanged.

binom_tol <- function(p, n, z = 4) z * sqrt(p * (1 - p) / n)
mean_p11 <- out %>% dplyr::group_by(time) %>%
    dplyr::summarise(mean_p11 = mean(p11), .groups = "drop")

check_points <- c(2, 6, 12, 18, 23.5)
cat("\n==== Check (b): Kaplan-Meier vs. mean-p11 ====\n")
if (requireNamespace("survival", quietly = TRUE)) {
    fit_km <- survival::survfit(
        survival::Surv(sim_events$sim_time, sim_events$sim_status) ~ 1)
    km_at <- function(tj) {
        idx <- findInterval(tj, fit_km$time)
        if (idx == 0L) 1 else fit_km$surv[idx]
    }
    comparison_b <- lapply(check_points, function(tj) {
        s_km <- km_at(tj)
        s_mean_p11 <- mean_p11$mean_p11[which.min(abs(mean_p11$time - tj))]
        tol <- binom_tol(1 - s_mean_p11, n_subjects)
        data.frame(t = tj, km = s_km, mean_p11 = s_mean_p11,
            abs_diff = abs(s_km - s_mean_p11), tol_4se = tol,
            within_tol = abs(s_km - s_mean_p11) <= tol)
    })
    comparison_b <- dplyr::bind_rows(comparison_b)
    print(comparison_b, row.names = FALSE, digits = 4)
    cat("All checkpoints within 4-SE binomial tolerance:",
        all(comparison_b$within_tol), "\n")
} else {
    cat("'survival' package not available; skipping Check (b) ",
        "(the p11-only checks above and Check (c)/(d)/(e) below do not ",
        "require it).\n", sep = "")
}

## ---- Check (c): IC50 potency sweep (manuscript direction) ---------------
## Manuscript text: "IC50 = 0.3 -> event rate = 0.503; IC50 = 1.0 ->
## event rate = 0.560; IC50 = 5.0 -> event rate = 0.727" (asserted only
## as prose there). Reproduced here as a programmatic monotonicity
## check: lower IC50 = greater potency = greater mediator suppression =
## lower event rate, i.e. event rate must be non-decreasing in IC50.

cat("\n==== Check (c): event rate vs. IC50 (lower potency -> more events) ====\n")
ic50_values <- c(0.3, 1.0, 5.0)
event_rates_ic50 <- vapply(ic50_values, function(ic50) {
    set.seed(20260830) # identical BSV draws and U draws across IC50 values
    out_ic50 <- as.data.frame(mrgsim(mod, data = data_cohort, end = 24,
        delta = 0.5, obsonly = TRUE, param = list(IC50 = ic50)))
    sim_ic50 <- sim_tte_df(out_ic50[, c("ID", "time", "p11")])
    mean(sim_ic50$sim_status)
}, numeric(1))
ic50_sweep <- data.frame(IC50 = ic50_values, event_rate = event_rates_ic50)
print(ic50_sweep, row.names = FALSE)
cat("Event rate monotonically non-decreasing in IC50:",
    all(diff(ic50_sweep$event_rate) >= -1e-8), "\n")

## ---- Check (d): between-subject variability is active -------------------

cat("\n==== Check (d): between-subject variability sanity check ====\n")
cp_at_checkpoint <- out$CP[out$time == 6]
cat("Variance of CP across subjects at t = 6:", var(cp_at_checkpoint), "\n")
cat("BSV propagating into concentration (variance > 0):",
    var(cp_at_checkpoint) > 0, "\n")
event_time_var <- var(sim_events$sim_time)
cat("Variance of simulated event/censoring time:", event_time_var, "\n")
cat("BSV propagating into event times (variance > 0):",
    event_time_var > 0, "\n")

## ---- Check (e): fixed-seed reproducibility -------------------------------

cat("\n==== Check (e): fixed-seed reproducibility ====\n")
set.seed(20260830)
out_rep1 <- as.data.frame(mrgsim(mod, data = data_cohort, end = 24,
    delta = 0.5, obsonly = TRUE))
sim_rep1 <- sim_tte_df(out_rep1[, c("ID", "time", "p11")])

set.seed(20260830)
out_rep2 <- as.data.frame(mrgsim(mod, data = data_cohort, end = 24,
    delta = 0.5, obsonly = TRUE))
sim_rep2 <- sim_tte_df(out_rep2[, c("ID", "time", "p11")])

cat("mrgsim() output identical across repeated runs (same seed):",
    identical(out_rep1, out_rep2), "\n")
cat("sim_tte_df() output identical across repeated runs (same seed):",
    identical(sim_rep1, sim_rep2), "\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
cat("n_subjects:", n_subjects, "; seed: 20260830\n")
