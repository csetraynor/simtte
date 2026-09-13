Model file:  gompertz_ode.cpp
[PROB]
# Model: `Gompertz hazard, in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Survival function: S(t) = exp[-(eta / gamma) * (exp(gamma * t) - 1)],
                      eta = exp(mu + lp), for gamma != 0
                      (gamma = 0 reduces to the exponential model)
  Hazard function:   h(t) = eta * exp(gamma * t)

  Unlike the Weibull hazard, h(t) is finite everywhere, including
  t = 0 (h(0) = eta) and for gamma < 0 (a decreasing hazard) -- there
  is no t -> 0+ singularity of the kind risk R1 addresses for
  weibull_ode.cpp, so no time-floor clamp is needed here.

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro on macOS/BSD-derived toolchains and the model
  fails to compile.
  Use `TEVT`/`event_found` instead (see exponential_ode.cpp).

[PARAM] @annotated
  mu    :  0.1 : Intercept
  lp    :  0   : Log hazard ratio (linear predictor); eta = exp(mu + lp)
  gamma :  0.1 : Gompertz shape parameter (rate of exponential change in the hazard; may be negative)
  U     :  0   : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END   :  1e9 : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  // -- BEGIN simtte survival scaffolding -----------------------------
  // General-purpose in-solver event-detection scaffolding, identical
  // in every library model this package ships; see exponential_ode.cpp
  // for the full rationale comment (kept there only, to avoid repeating
  // it verbatim in every file).
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
  double eta = exp(mu + lp);
  double HAZ = eta * exp(gamma * SOLVERTIME);   // <- model-specific line
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
