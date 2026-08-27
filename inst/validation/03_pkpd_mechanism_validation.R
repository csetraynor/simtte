## Phase D validation, component 3: compact PK/PD mechanistic validation.
##
## Builds one small custom mrgsolve model (a one-compartment IV-bolus PK
## model linearly coupled to a hazard, dxdt_p11 = -p11 * HAZ, the same
## survival co-integration pattern documented in the manuscript) and
## checks that simtte's event-time simulation (via sim_tte_df(), the
## function intended for exactly this custom-model workflow) responds
## appropriately to the underlying mechanistic trajectory.
##
## Independent reference: the PK and hazard are chosen so that the
## cumulative hazard has a closed form, giving an S(t) that does not
## depend on simtte or on the mrgsolve ODE solve at all:
##
##   C(t)    = C0 * exp(-k * t),          C0 = Dose / V,  k = CL / V
##   h(t)    = H0 - SLOPE * C(t),         (requires H0 >= SLOPE * C0)
##   H(t)    = integral_0^t h(s) ds = H0 * t - SLOPE * C0 * (1 - exp(-k*t)) / k
##   S(t)    = exp(-H(t))
##
## This is a genuinely independent check: the mrgsolve model computes
## S(t) by numerically integrating the ODE system; the reference below
## computes it from closed-form calculus, using no simtte or mrgsolve
## code.
##
## Reproducibility: a single seed controls the uniform draws consumed by
## sim_tte_df(); the PK/PD model itself has no between-subject
## variability (kept out deliberately to make the S(t) reference
## unambiguous -- BSV is a modeling choice orthogonal to what Phase D is
## validating here).

suppressMessages(library(simtte))
suppressMessages(library(mrgsolve))
suppressMessages(library(dplyr))

## ---- Compact PK/PD model --------------------------------------------

code <- '
$PARAM CL = 1, V = 10, H0 = 0.3, SLOPE = 0.02

$CMT CENT

$INIT p11 = 1

$MAIN
double K = CL / V;

$ODE
double C   = CENT / V;
double HAZ = H0 - SLOPE * C;
dxdt_CENT = -K * CENT;
dxdt_p11  = -p11 * HAZ;

$CAPTURE C HAZ
'
mod <- mcode("pkpd_validation", code, quiet = TRUE)

analytical_S <- function(t, dose, V, CL, H0, SLOPE) {
    C0 <- dose / V
    k <- CL / V
    Ht <- H0 * t - SLOPE * C0 * (1 - exp(-k * t)) / k
    exp(-Ht)
}

dose <- 100
V <- 10
CL <- 1
H0 <- 0.3
SLOPE <- 0.02
stopifnot(H0 >= SLOPE * dose / V) # hazard must stay non-negative

grid <- seq(0, 60, by = 0.5)

## ---- Check 1: mrgsolve-simulated S(t) vs. the closed-form reference --

data_dose <- ev(amt = dose, cmt = 1, time = 0)
out <- as.data.frame(mrgsim(mod, data = data_dose, end = -1, add = grid,
    obsonly = TRUE, output = "df"))
ref_S <- analytical_S(out$time, dose, V, CL, H0, SLOPE)

cat("==== Check 1: mrgsolve ODE trajectory vs. closed-form S(t) ====\n")
cat("Max absolute error:", max(abs(out$p11 - ref_S)), "\n")
cat("Mean absolute error:", mean(abs(out$p11 - ref_S)), "\n")

## ---- Check 2: simulated event-time distribution vs. closed-form S(t) -

n_subjects <- 3000
set.seed(20260827)

traj <- data.frame(ID = rep(1, length(grid)), time = grid, p11 = out$p11)
traj_rep <- traj[rep(seq_len(nrow(traj)), n_subjects), ]
traj_rep$ID <- rep(seq_len(n_subjects), each = length(grid))

sim_events <- sim_tte_df(traj_rep)

binom_tol <- function(p, n, z = 4) z * sqrt(p * (1 - p) / n)

check_points <- c(5, 10, 20, 30, 45, 60)
cat("\n==== Check 2: empirical event probability vs. analytical 1 - S(t) ====\n")
comparison <- lapply(check_points, function(tj) {
    p_analytic <- 1 - analytical_S(tj, dose, V, CL, H0, SLOPE)
    p_empirical <- mean(sim_events$sim_time <= tj & sim_events$sim_status == 1)
    tol <- binom_tol(p_analytic, n_subjects)
    data.frame(t = tj, p_analytic = p_analytic, p_empirical = p_empirical,
        abs_diff = abs(p_analytic - p_empirical), tol_4se = tol,
        within_tol = abs(p_analytic - p_empirical) <= tol)
})
comparison <- dplyr::bind_rows(comparison)
print(comparison, row.names = FALSE, digits = 4)
cat("All comparisons within 4-SE binomial tolerance:",
    all(comparison$within_tol), "\n")

p_cens_analytic <- analytical_S(max(grid), dose, V, CL, H0, SLOPE)
p_cens_empirical <- mean(sim_events$sim_status == 0)
cat("\nCensoring proportion at end of follow-up: analytical =",
    round(p_cens_analytic, 4), " empirical =",
    round(p_cens_empirical, 4), "\n")

## ---- Check 3: monotone dose-response sanity check ---------------------
## Increasing SLOPE (a stronger protective drug effect) should
## monotonically decrease the simulated event rate.

cat("\n==== Check 3: event rate as a function of drug effect (SLOPE) ====\n")
## A shorter follow-up horizon than Checks 1-2 is used here deliberately:
## by t = 60 the baseline hazard alone (H0 = 0.3) makes an event nearly
## certain regardless of SLOPE (see the t = 45/60 rows in Check 2), which
## would saturate every SLOPE value at an event rate of 1 and hide any
## drug effect. At the shorter horizon below, baseline survival is not
## yet saturated, so a genuine SLOPE gradient is observable.
grid_short <- seq(0, 8, by = 0.5)
slope_values <- c(0, 0.005, 0.01, 0.02, 0.03)
event_rates <- vapply(slope_values, function(s) {
    out_s <- as.data.frame(mrgsim(mod, data = data_dose, end = -1,
        add = grid_short, obsonly = TRUE, output = "df",
        param = list(SLOPE = s)))
    traj_s <- data.frame(ID = rep(1, length(grid_short)), time = grid_short,
        p11 = out_s$p11)
    traj_s_rep <- traj_s[rep(seq_len(nrow(traj_s)), n_subjects), ]
    traj_s_rep$ID <- rep(seq_len(n_subjects), each = length(grid_short))
    set.seed(20260827) # identical U draws across SLOPE values
    ev <- sim_tte_df(traj_s_rep)
    mean(ev$sim_status)
}, numeric(1))
dose_response <- data.frame(SLOPE = slope_values, event_rate = event_rates)
print(dose_response, row.names = FALSE)
cat("Event rate monotonically non-increasing in SLOPE:",
    all(diff(dose_response$event_rate) <= 1e-8), "\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
cat("Model parameters: dose =", dose, ", V =", V, ", CL =", CL, ", H0 =",
    H0, ", SLOPE =", SLOPE, "\n")
cat("Seed: 20260827; n_subjects:", n_subjects, "\n")
