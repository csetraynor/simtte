## Phase E (Internal Solver Threshold Tracking) investigation.
##
## REJECTED APPROACH -- NOT SHIPPED. This script documents, and
## reproducibly demonstrates, why $DES/$ODE internal-solver
## threshold-tracking was investigated and NOT adopted as a third
## event_time_method for the built-in Weibull/M-spline models. It does
## not modify, and is never sourced by, any package file; it exists
## purely as a reproducible record of the investigation described in
## PHASE_E_THRESHOLD_TRACKING_REPORT.md.
##
## The draft mechanism (as specified): sample U before simulation, pass
## it into the model as a parameter, and inside $DES/$ODE record
## SOLVERTIME the first time an internal solver evaluation satisfies
## p11 <= U.
##
## Five findings, each reproduced below:
##   1. BLOCKER (Weibull): inst/models/weibull.cpp has no $ODE/$DES
##      block at all (Phase A replaced it with a closed-form $TABLE
##      expression specifically because ODE integration is unstable for
##      shape < 1) -- there is nothing to attach threshold tracking to
##      without reverting that fix.
##   2. Where structurally possible (M-spline, which still integrates a
##      real ODE), tracked event times are quantized onto the SOLVER'S
##      OWN internal adaptive step grid, not continuous -- and that grid
##      depends on solver tolerance and, to a lesser extent, on the
##      requested output grid itself.
##   3. For the one case where it is exact (piecewise-constant hazard,
##      i.e. the M-spline model), the existing event_time_method =
##      "log_survival" is ALREADY exact (machine precision), making
##      threshold tracking strictly worse there.
##   4. CORRECTNESS DEFECT: the naive mechanism as specified can report
##      an event PAST the administrative end_time boundary, because the
##      solver evaluates at least one internal step beyond the final
##      requested output time before finalizing, and a one-shot flag
##      cannot "un-latch" if that step is later superseded.
##   5. Reproducibility itself is not a concern (deterministic given
##      fixed data), but the mechanism only works because mrgsolve
##      currently simulates individuals strictly sequentially, single
##      threaded, with per-individual C++ `static` globals reset at
##      NEWIND <= 1 -- an implementation detail, not part of mrgsolve's
##      documented public contract.
##
## Reproducibility: mrgsolve models are built here with mcode() (not
## from any package file); a fixed seed controls all U draws.

suppressMessages(library(mrgsolve))

cat("==== Finding 1: Weibull has no $ODE/$DES block to attach to ====\n")
## Resolves against an installed/loaded simtte first (system.file()),
## falling back to the conventional source-tree-relative path used by
## the other inst/validation/ scripts (run from the package root).
candidate_paths <- c(
    system.file("models", "weibull.cpp", package = "simtte"),
    "inst/models/weibull.cpp"
)
candidate_paths <- candidate_paths[nzchar(candidate_paths)]
weibull_path <- candidate_paths[file.exists(candidate_paths)][1]
if (!is.na(weibull_path)) {
    weibull_src <- readLines(weibull_path)
    has_ode <- any(grepl("^\\s*\\[ODE\\]|^\\s*\\[DES\\]", weibull_src))
    has_table <- any(grepl("^\\s*\\[TABLE\\]", weibull_src))
    cat("weibull.cpp contains an [ODE]/[DES] block:", has_ode, "\n")
    cat("weibull.cpp contains a [TABLE] block (closed-form, no solver):",
        has_table, "\n")
    stopifnot(!has_ode, has_table)
    cat("CONFIRMED: no internal solver evaluations exist to track for",
        "the built-in Weibull model as currently implemented.\n")
} else {
    cat("(weibull.cpp not found from this working directory; see",
        "PHASE_E_THRESHOLD_TRACKING_REPORT.md for the confirmed result:",
        "no [ODE]/[DES] block, only [TABLE].)\n")
}

## ---- Findings 2-4: M-spline-structure model with $DES tracking ------

des_code <- paste(
    "$PARAM lp = 0, mu = 0, basehaz = 1, U = 0.5",
    "$INIT p11 = 1",
    "$GLOBAL",
    "static int event_found = 0;",
    "static double event_time = -1.0;",
    "$MAIN",
    "double eta = exp(mu + lp);",
    "if (NEWIND <= 1) { event_found = 0; event_time = -1.0; }",
    "$DES",
    "if (p11 <= U && event_found == 0) {",
    "    event_found = 1;",
    "    event_time = SOLVERTIME;",
    "}",
    "dxdt_p11 = -p11 * basehaz * eta;",
    "$TABLE",
    "double tracked_time = event_time;",
    "double tracked_found = event_found;",
    "$CAPTURE tracked_time tracked_found",
    sep = "\n"
)
mod <- mcode("phase_e_des_investigation", des_code, quiet = TRUE)

## The c(1, 2, 4) on c(0, 1, 2) hazard example used throughout prior
## phases (test-ms-hazard-carry.R, PHASE_B_DESIGN_AUDIT.md): h = 1 on
## [0, 1), h = 2 on [1, 2). Analytical piecewise-constant-hazard event
## time given U:
T_true_fn <- function(U) {
    ifelse(U >= exp(-1), -log(U),
        ifelse(U >= exp(-3), 1 + (-log(U) - 1) / 2, NA_real_))
}

n <- 1000
set.seed(20260828)
U <- runif(n)
T_true <- T_true_fn(U)

data <- data.frame(
    ID = rep(seq_len(n), each = 3),
    time = rep(c(0, 1, 2), n),
    amt = 0, evid = 1, cmt = 1,
    basehaz = rep(c(1, 2, 4), n),
    mu = 0, lp = 0,
    U = rep(U, each = 3)
)

out <- as.data.frame(mrgsim(data_set(mod, data), tgrid = c(0, 2),
    obsonly = TRUE, output = "df", nocb = FALSE))
last_row <- out[out$time == max(out$time), ]

known <- !is.na(T_true)
err <- abs(last_row$tracked_time[known] - T_true[known])

cat("\n==== Finding 2/3: accuracy vs. the M-spline log_survival exactness ====\n")
cat("Subjects with finite T_true:", sum(known), "/", n, "\n")
cat("Mean absolute error (DES tracking):", mean(err), "\n")
cat("Max absolute error (DES tracking):", max(err), "\n")
cat("log_survival's error on this exact scenario is 0 (machine precision;",
    "see test-log-survival-ms.R and PHASE_B_DESIGN_AUDIT.md Section 9).\n")
stopifnot(mean(err) > 1e-6) # DES tracking is measurably NOT exact here

cat("\n==== Finding 4: administrative end_time boundary can be violated ====\n")
censored <- !known
n_censored <- sum(censored)
n_misclassified <- sum(last_row$tracked_found[censored] == 1)
cat("Subjects that should be censored (true event time > end_time = 2):",
    n_censored, "\n")
cat("Of those, incorrectly reported as events by DES tracking:",
    n_misclassified, "\n")
if (n_misclassified > 0) {
    bad_times <- last_row$tracked_time[censored][
        last_row$tracked_found[censored] == 1]
    cat("Their reported event time(s), which exceed end_time = 2:\n")
    print(unique(bad_times))
    stopifnot(any(bad_times > 2))
    cat("CONFIRMED: this is a genuine correctness defect in the naive",
        "single-flag mechanism as specified in the draft.\n")
} else {
    cat("(Not reproduced in this run; see report for the run in which",
        "it was observed -- the underlying cause, a solver evaluation",
        "one step beyond the requested end, is deterministic given the",
        "model/data/tolerances, but its manifestation as a reportable",
        "misclassification depends on exactly which U values are drawn.)\n")
}

## ---- Finding 5: reproducibility and tolerance/grid dependence -------

const_code <- paste(
    "$PARAM H = 0.5, U = 0.5",
    "$INIT p11 = 1",
    "$GLOBAL",
    "static int event_found = 0;",
    "static double event_time = -1.0;",
    "$MAIN",
    "if (NEWIND <= 1) { event_found = 0; event_time = -1.0; }",
    "$DES",
    "if (p11 <= U && event_found == 0) {",
    "    event_found = 1;",
    "    event_time = SOLVERTIME;",
    "}",
    "dxdt_p11 = -p11 * H;",
    "$TABLE",
    "double tracked_time = event_time;",
    "$CAPTURE tracked_time",
    sep = "\n"
)
mod2 <- mcode("phase_e_des_investigation_const", const_code, quiet = TRUE)

n2 <- 300
set.seed(1)
U2 <- runif(n2)
H <- 0.5
data2 <- data.frame(ID = seq_len(n2), time = 0, amt = 0, evid = 1, cmt = 1,
    H = H, U = U2)

run_once <- function(rtol = 1e-8, atol = 1e-8, tgrid = c(0, 30)) {
    out <- as.data.frame(mrgsim(data_set(mod2, data2), tgrid = tgrid,
        obsonly = TRUE, output = "df", rtol = rtol, atol = atol))
    out[out$time == max(out$time), "tracked_time"]
}

r1 <- run_once()
r2 <- run_once()
cat("\n==== Finding 5: reproducibility (fixed data, repeated calls) ====\n")
cat("Identical across repeated calls:", identical(r1, r2), "\n")
stopifnot(identical(r1, r2))

r_loose <- run_once(rtol = 1e-4, atol = 1e-4)
r_tight <- run_once(rtol = 1e-12, atol = 1e-12)
cat("\n==== Solver-tolerance dependence ====\n")
cat("Mean |loose - default|:", mean(abs(r_loose - r1)), "\n")
cat("Mean |tight - default|:", mean(abs(r_tight - r1)), "\n")
stopifnot(mean(abs(r_loose - r1)) > 1e-3) # tolerance materially changes results

r_fine_grid <- run_once(tgrid = seq(0, 30, by = 0.01))
cat("\n==== Requested-output-grid dependence (should ideally be none) ====\n")
cat("Mean |fine_tgrid - coarse_tgrid|:", mean(abs(r_fine_grid - r1)), "\n")
cat("(Non-zero: tracked_time is not fully decoupled from the requested",
    "output grid, contrary to the hoped-for grid-independent property.)\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
