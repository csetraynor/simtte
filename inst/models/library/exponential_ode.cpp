Model file:  exponential_ode.cpp
[PROB]
# Model: `Exponential (constant) hazard, in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Survival function: S(t) = exp(-eta * t), eta = H0 * exp(lp)
  Hazard function:   h(t) = eta (constant)

  This is the Phase 1 library entry for `sim_tte_ode()`'s in-solver
  event-detection mechanism (see reports/02_technical_design.md section
  2 and reports/05_phase1_report.md). Unlike the closed-form
  `weibull.cpp`/`weibull_tv.cpp` engine models, `p11` here is a genuine
  ODE state, and the event/censoring time is located *during*
  integration (at the solver's own internal evaluations), not resolved
  afterwards from a materialized trajectory the way `sim_tte_df()` does.

  NAMING RULE (do not remove this comment; copy it into every future
  library model): never name a $GLOBAL/tracking variable `ETIME` --
  `<sys/errno.h>` defines `ETIME` as a macro (STREAMS "ioctl timeout",
  value 101) on macOS/BSD-derived toolchains, which is pulled in
  transitively by mrgsolve's generated model source and makes the model
  fail to compile. Verified directly (reports/experiments/
  01_etime_boundary_experiment.R); use `TEVT` (event time) and
  `event_found` (one-shot latch flag) instead, as done below. Before
  introducing any new all-caps $GLOBAL name, check it isn't a libc/POSIX
  errno macro.

[PARAM] @annotated
  H0  :  0.1 : Baseline (constant) hazard (1/time)
  lp  :  0   : Log hazard ratio (linear predictor); eta = H0 * exp(lp)
  U   :  0   : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END :  1e9 : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[INIT]
  // `p11` is not declared in a separate $CMT block: mrgsolve infers the
  // compartment from this $INIT assignment alone, exactly like
  // inst/models/ms.cpp and inst/models/weibull_tv.cpp already do
  // (declaring it in both places raises "Duplicated model names").
  p11 = 1

[GLOBAL]
  // -- BEGIN simtte survival scaffolding -----------------------------
  // General-purpose in-solver event-detection scaffolding, identical
  // in every library model this package ships. Written as a single,
  // self-contained, copy-pasteable unit (this $GLOBAL declaration, the
  // $MAIN reset below, and the $ODE latch below, excluding only the
  // model-specific `double HAZ = ...;` line) so a future model-
  // converter utility (reports/04_author_decisions.md section 6) can
  // lift it verbatim into a user-supplied mrgsolve model.
  static int    event_found = 0;
  static double TEVT        = 0.0;
  // Grid-free refinement bracket (Phase 2.5, reports/04_author_decisions.md
  // "After the test runbook / Phase 2.5"): the solver's own last
  // pre-crossing evaluation (T_PRE, P_PRE) and the crossing evaluation
  // itself (TEVT, P_POST), one internal solver step apart -- narrow
  // enough that interpolating between them makes the constant-hazard-
  // within-bracket assumption close to exact for any smooth hazard,
  // independent of the reported output grid (`delta`).
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;

[MAIN]
  if (NEWIND <= 1) {
    event_found = 0;
    TEVT = 0.0;
    T_PRE = 0.0;
    P_PRE = 1.0;
    P_POST = 0.0;
  }

[ODE]
  double eta = H0 * exp(lp);
  double HAZ = eta;                 // <- the only model-specific line
  dxdt_p11 = -p11 * HAZ;
  // Record the last pre-crossing evaluation. The SOLVERTIME >= T_PRE
  // guard is a monotone update: a rejected/retried step that revisits
  // an earlier time must not move the bracket backwards.
  if (!event_found && p11 > U && SOLVERTIME <= END && SOLVERTIME >= T_PRE) {
    T_PRE = SOLVERTIME;
    P_PRE = p11;
  }
  if (!event_found && p11 <= U && SOLVERTIME <= END) {
    event_found = 1;
    TEVT = SOLVERTIME;
    P_POST = p11;
  }
  // -- END simtte survival scaffolding (ODE half) --------------------

[CAPTURE] @annotated
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
