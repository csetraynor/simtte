Model file:  irm4_hazard.cpp
[PROB]
# Model: `Response-driven hazard (IRM4), in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Backbone: `mrgsolve::modlib("irm4")` (indirect response model, type 4
  -- stimulation of response elimination; two-compartment PK, optional
  nonlinear clearance) -- every PARAM/CMT/GLOBAL(PK part)/MAIN(PK
  part)/ODE(PK part)/CAPTURE(PK part) block below is that backbone's
  own code, unedited, INCLUDING its plain (non-annotated, comma-form)
  PARAM/CMT/CAPTURE style, which -- unlike irm1/irm2/irm3's own
  backbones -- is how mrgsolve ships irm4 itself; "verbatim" is taken
  literally here, not normalized to the annotated style. See
  irm1_hazard.cpp for the full mechanical-edit list (identical here,
  applied to a different backbone dynamic) and for why this doc block
  avoids dollar-prefixed block names and the literal "at-annotated" tag
  text.

  Hazard link (log-linear in the response, relative to its own
  pre-dose baseline) -- identical form to irm1_hazard.cpp/
  irm2_hazard.cpp/irm3_hazard.cpp, only the backbone RESP dynamics
  differ:
    HAZ = H0 * exp(beta_r * (RESP / RESP0 - 1)), RESP0 = KIN / KOUT.
  The sign of `beta_r` gives the direction: positive makes a response
  *above* baseline harmful (raises the hazard); negative makes it
  protective. `beta_r = 0` reduces exactly to the exponential model
  with hazard `H0`, independent of RESP.

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro (and ERANGE/EDOM-style names are libc macros
  too) on macOS/BSD-derived toolchains and the model fails to compile.
  Use `TEVT`/`event_found` instead. None of
  H0/beta_r/RESP0/CL/V2/KA/KA2/Q/V3/KIN/KOUT/EC50/EMAX/n/VMAX/KM
  collide with an errno.h name.

  BETWEEN-SUBJECT VARIABILITY (reports/11_bsv_review.md,
  reports/12_bsv_implementation_report.md): `sim_tte_ode(model =
  "irm4_hazard", omega = ...)` applies BSV via per-subject `idata`
  columns, not a declared OMEGA block -- see `.ODE_BSV_TARGETS` in
  R/helpers.R. Valid targets (every structural PK/PD parameter of
  this backbone, taken from this file's own $PARAM block below, in
  this order -- irm4's own plain-comma-form parameter order, unlike
  irm1-3): CL, V2, KA, KA2, Q, V3, KIN, KOUT, EC50, EMAX, VMAX, KM, n.
  H0/lp/beta_r/U/END (the hazard-scaffold parameters) are not BSV
  targets.

[PARAM]
  CL=1, V2=10, KA=0.5, KA2=0.5
  Q = 0, V3=10
  KIN = 10, KOUT=2, EC50 = 2, EMAX=1
  VMAX = 0, KM=2, n=1

[PARAM] @annotated
  H0     : 0.01 : Baseline (RESP == RESP0, lp = 0) hazard (1/time)
  lp     : 0    : Log hazard ratio (linear predictor); composes with beta_r * (RESP/RESP0 - 1) for sim_tte_ode()'s covariates/beta mechanism (see ?sim_tte_ode "Time-varying covariates")
  beta_r : 0    : Log hazard ratio per unit relative response deviation; HAZ(t) = H0 * exp(lp + beta_r * (RESP(t)/RESP0 - 1))
  U      : 0    : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END    : 1e9  : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[CMT]
  EV CENT PERIPH RESP EV2

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  #define CP (CENT/V2)
  #define CT (PERIPH/V3)
  #define CLNL (VMAX/(KM+CP))
  #define STIM (EMAX*pow(CP,n)/(pow(EC50,n)+pow(CP,n)))
  // -- BEGIN simtte survival scaffolding -----------------------------
  static int    event_found = 0;
  static double TEVT        = 0.0;
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;
  static double RESP0       = 0.0;
  // -- END simtte survival scaffolding (GLOBAL half) ------------------

[MAIN]
  RESP_0 = KIN/KOUT;
  // -- BEGIN simtte survival scaffolding -----------------------------
  RESP0 = KIN/KOUT;
  if (NEWIND <= 1) {
    event_found = 0;
    TEVT = 0.0;
    T_PRE = 0.0;
    P_PRE = 1.0;
    P_POST = 0.0;
  }
  // -- END simtte survival scaffolding (MAIN half) --------------------

[ODE]
  dxdt_EV     = -KA*EV;
  dxdt_EV2    = -KA2*EV2;
  dxdt_CENT   =  KA*EV + KA2*EV2 - (CL+CLNL+Q)*CP  + Q*CT;
  dxdt_PERIPH =  Q*CP - Q*CT;
  dxdt_RESP   =  KIN - KOUT*(1+STIM)*RESP;
  // -- BEGIN simtte survival scaffolding -----------------------------
  double HAZ = H0 * exp(lp + beta_r * (RESP/RESP0 - 1.0));   // <- the only model-specific line
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

[CAPTURE]
  CP

[CAPTURE] @annotated
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
  HAZ         : Instantaneous hazard, H0 * exp(lp + beta_r * (RESP/RESP0 - 1))
