Model file:  weibull_ode.cpp
[PROB]
# Model: `Weibull hazard, in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Survival function: S(t) = exp(-eta * t^shape), eta = exp(mu + lp)
  Hazard function:   h(t) = shape * eta * t^(shape - 1)

  Same `mu`/`lp`/`shape` parameterization as the closed-form
  `inst/models/weibull.cpp`, so results are directly comparable
  one-to-one; unlike that model, `p11` here is a genuine `$ODE` state
  and the event time is located in-solver (see
  reports/02_technical_design.md section 2 and
  reports/06_phase2_report.md).

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro on macOS/BSD-derived toolchains and the model
  fails to compile.
  Use `TEVT`/`event_found` instead (see exponential_ode.cpp).

  WEIBULL-ODE SHAPE RANGE (risk R1, reports/03_implementation_plan.md):
  h(t) diverges as t -> 0+ for shape < 1 -- a genuine mathematical
  property of the instantaneous Weibull hazard, not a floating-point
  artifact (PHASE_A_REPORT.md section 4: this is exactly why the
  closed-form weibull.cpp evaluates S(t) directly instead of
  integrating an ODE). Integrating the raw hazard as an ODE derivative
  therefore fails outright (lsoda error, confirmed directly for this
  file's naive form across shape in {0.05, ..., 0.7}) unless the time
  used *inside the hazard expression only* is floored away from exactly
  0. `T_FLOOR` below implements that floor, evidenced (see
  reports/06_phase2_report.md, risk R1 section) to:
    - eliminate solver failure/NaN/non-monotone p11 for every shape in
      {0.01, ..., 10} tested, at every T_FLOOR from 1e-4 down to 1e-100
      (only T_FLOOR = 0, i.e. no floor, fails);
    - accuracy vs. the closed-form S(t) improves monotonically as
      T_FLOOR shrinks, with no stability cost observed down to 1e-100;
    - at T_FLOOR = 1e-100 specifically: max abs error is at or below
      mrgsolve's own default solver tolerance (~1e-8) for shape >= 0.05,
      growing to a still-finite, still-monotone but larger bias for
      shape < 0.05 (e.g. ~0.03 at shape = 0.01).
  This model therefore RUNS (no error) for the entire tested range
  (0.01-10); accuracy is solver-tolerance-scale for shape >= 0.05 and
  degrades gracefully (not catastrophically) below that. See
  ?sim_tte_ode's "Weibull shape support" section for the documented
  policy this implements.

[PARAM] @annotated
  mu    :  0.1 : Intercept
  lp    :  0   : Log hazard ratio (linear predictor); eta = exp(mu + lp)
  shape :  1   : Weibull shape parameter (must be > 0)
  U     :  0   : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END   :  1e9 : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  static int    event_found = 0;
  static double TEVT        = 0.0;
  // Evidenced floor on the time value used *only* inside the hazard
  // expression below (never on SOLVERTIME/END/TEVT themselves, which
  // must stay exact); see the shape-range note in [PROB] above.
  #define T_FLOOR 1e-100
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
  double eta = exp(mu + lp);
  double t_eff = SOLVERTIME > T_FLOOR ? SOLVERTIME : T_FLOOR;
  double HAZ = shape * eta * pow(t_eff, shape - 1);   // <- model-specific line
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

[CAPTURE] @annotated
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
