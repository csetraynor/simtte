Model file:  tmdd_hazard.cpp
[PROB]
# Model: `Target-engagement-driven hazard (TMDD), in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Backbone: `mrgsolve::modlib("tmdd")` (target-mediated drug disposition:
  free drug CP, free target REC, and drug-target complex RC, with
  target turnover KSYN/KDEG and binding KON/KOFF) -- every
  PARAM/CMT/GLOBAL(PK part)/MAIN(PK part)/ODE(PK part)/TABLE/
  CAPTURE(PK part) block below is that backbone's own code, unedited,
  INCLUDING its plain (non-annotated, comma-form) PARAM/CMT/CAPTURE
  style, which is how mrgsolve ships tmdd itself; "verbatim" is taken
  literally here, not normalized to the annotated style. See
  irm1_hazard.cpp for the full mechanical-edit list (identical here,
  minus the RESP0 baseline step -- see "Hazard link" below for why)
  and for why this doc block avoids dollar-prefixed block names and
  the literal "at-annotated" tag text.

  Hazard link -- chosen deliberately different in kind from the other
  five models, not just a different backbone:
    HAZ = H0 * exp(beta_rc * RC)
  RC (the bound drug-target complex) is the quantity unique to
  target-mediated disposition: linking the hazard to CP (free drug)
  instead would make this model differ from pk_hazard.cpp only in its
  PK backbone, not in the *kind* of link, and would not exercise
  anything TMDD-specific. RC is an absolute level (not relative to a
  baseline, unlike the RESP-based IRM link): REC_0 = R0 but RC_0 = 0
  by construction (no pre-dose complex exists), so there is no nonzero
  baseline to divide by the way RESP0 = KIN/KOUT serves the IRM models
  -- an absolute log-linear link (mirroring pk_hazard.cpp's own CP
  link) is therefore the natural, well-defined form, and is why this
  model's survival scaffold has no RESP0-equivalent baseline variable.
  `beta_rc = 0` reduces exactly to the exponential model with hazard
  `H0`, independent of RC. TMDD is stiff by nature (fast binding
  kinetics next to slow turnover), so this model is also this
  library's designated steepness stress case (reports/10_phase4_report.md).

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro (and ERANGE/EDOM-style names are libc macros
  too) on macOS/BSD-derived toolchains and the model fails to compile.
  Use `TEVT`/`event_found` instead. None of
  H0/beta_rc/KPT/KTP/V2/KA/KA2/KEL/R0/KDEG/KINT/KON/KOFF/KSYN collide
  with an errno.h name. The backbone's own `tmdd::` namespace (used to
  avoid a name clash between the `CP` macro and the PK derivative) is
  untouched.

  BETWEEN-SUBJECT VARIABILITY (reports/11_bsv_review.md,
  reports/12_bsv_implementation_report.md): `sim_tte_ode(model =
  "tmdd_hazard", omega = ...)` applies BSV via per-subject `idata`
  columns, not a declared OMEGA block -- see `.ODE_BSV_TARGETS` in
  R/helpers.R. Valid targets (every structural parameter of this
  backbone, taken from this file's own $PARAM block below, in this
  order): KPT, KTP, V2, KA, KA2, KEL, R0, KDEG, KINT, KON, KOFF.
  H0/lp/beta_rc/U/END (the hazard-scaffold parameters) are not BSV
  targets.

[PARAM]
  KPT = 0.064, KTP = 0.123, V2=0.032, KA = 0.142, KA2 = 0.6
  KEL = 0.106
  R0 = 64.31, KDEG = 0.079, KINT = 2, KON=0.101, KOFF = 10.1

[PARAM] @annotated
  H0      : 0.01 : Baseline (RC == 0, lp = 0) hazard (1/time)
  lp      : 0    : Log hazard ratio (linear predictor); composes with beta_rc * RC for sim_tte_ode()'s covariates/beta mechanism (see ?sim_tte_ode "Time-varying covariates")
  beta_rc : 0    : Log hazard ratio per unit drug-target complex concentration; HAZ(t) = H0 * exp(lp + beta_rc * RC(t))
  U       : 0    : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END     : 1e9  : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[CMT]
  EV CENT TISS REC RC EV2

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  namespace tmdd {
    double _dxdt_CP=0;
    double TMDDR0 = 0;
  }
  #define KSYN (R0*KDEG)
  #define CP (CENT/V2)
  // -- BEGIN simtte survival scaffolding -----------------------------
  static int    event_found = 0;
  static double TEVT        = 0.0;
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;
  // -- END simtte survival scaffolding (GLOBAL half) ------------------

[MAIN]
  REC_0 = R0;
  tmdd::TMDDR0 = _R(3);
  // -- BEGIN simtte survival scaffolding -----------------------------
  if (NEWIND <= 1) {
    event_found = 0;
    TEVT = 0.0;
    T_PRE = 0.0;
    P_PRE = 1.0;
    P_POST = 0.0;
  }
  // -- END simtte survival scaffolding (MAIN half) --------------------

[ODE]
  dxdt_EV = -KA*EV;
  dxdt_EV2 = -KA2*EV2;
  tmdd::_dxdt_CP = (tmdd::TMDDR0+KA*EV + KA2*EV2)/V2 - (KEL+KPT)*CP - KON*CP*REC + KOFF*RC + KTP*TISS/V2;
  dxdt_CENT = tmdd::_dxdt_CP * V2;
  dxdt_TISS = KPT*CP*V2 - KTP*TISS;
  dxdt_REC = KSYN - KDEG*REC - KON*CP*REC + KOFF*RC;
  dxdt_RC = KON*CP*REC - (KINT+KOFF)*RC;
  // -- BEGIN simtte survival scaffolding -----------------------------
  double HAZ = H0 * exp(lp + beta_rc * RC);   // <- the only model-specific line
  dxdt_p11 = -p11 * HAZ;
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

[TABLE]
  double TOTAL = REC+RC;

[CAPTURE]
  CP TOTAL

[CAPTURE] @annotated
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
  HAZ         : Instantaneous hazard, H0 * exp(lp + beta_rc * RC)
