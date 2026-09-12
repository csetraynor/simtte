Model file:  irm1_hazard.cpp
[PROB]
# Model: `Response-driven hazard (IRM1), in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Backbone: `mrgsolve::modlib("irm1")` (indirect response model, type 1
  -- inhibition of response input; two-compartment PK, optional
  nonlinear clearance) -- every PARAM/CMT/GLOBAL(PK part)/MAIN(PK
  part)/ODE(PK part)/CAPTURE(PK part) block below is that backbone's
  own code, unedited. (NOTE: this doc block intentionally spells block
  names without a leading dollar sign and avoids the literal
  "at-annotated" tag text -- mrgsolve's spec scraper scans the whole
  file for those tokens, not just line starts, and a prose mention of
  them elsewhere in the file breaks parsing; see pk_hazard.cpp, where
  this was found.)

  Hazard link (log-linear in the response, relative to its own
  pre-dose baseline):
    HAZ = H0 * exp(beta_r * (RESP / RESP0 - 1)), RESP0 = KIN / KOUT
    (the backbone's own steady-state baseline, copied once per subject
    from the RESP_0 initial-condition assignment already present in
    the backbone's MAIN block).
  The sign of `beta_r` gives the direction: positive makes a
  response *above* baseline harmful (raises the hazard); negative
  makes it protective. `beta_r = 0` reduces exactly to the exponential
  model with hazard `H0`, independent of RESP.

  MECHANICAL EDITS applied to the `irm1` backbone to reach this file
  (identical list and order used for irm2_hazard.cpp/irm3_hazard.cpp/
  irm4_hazard.cpp; see pk_hazard.cpp for the same list applied to a
  backbone without a response compartment -- this list is the working
  specification for the future model-converter utility,
  reports/04_author_decisions.md section 6):
    1. Copy the backbone's PARAM block(s) verbatim.
    2. Copy the backbone's CMT block verbatim.
    3. Copy the backbone's GLOBAL block verbatim, then append the
       survival scaffold (event_found/TEVT/T_PRE/P_PRE/P_POST, plus
       RESP0 for the response-relative link) to the end of that same
       block (GLOBAL may only appear once per model).
    4. Copy the backbone's MAIN block verbatim (it already assigns
       RESP_0 = KIN/KOUT), then append `RESP0 = KIN/KOUT;` (the plain
       baseline copy the ODE link reads) and the survival scaffold's
       reset logic, to the end of that same block (MAIN may only
       appear once per model).
    5. Copy the backbone's ODE block verbatim, then append: one
       model-specific `double HAZ = ...;` line (the only line that
       differs across the six PK/PD hazard models -- the response-
       relative form is shared by all four IRM variants), `dxdt_p11 =
       -p11 * HAZ;`, and the T_PRE/P_PRE/TEVT/P_POST latch (identical
       in every model in this package).
    6. Copy the backbone's CAPTURE block verbatim, then add a second
       annotated CAPTURE block for TEVT/event_found/T_PRE/P_PRE/
       P_POST/HAZ (CAPTURE, unlike GLOBAL/MAIN, may appear more than
       once and is concatenated). RESP/RESP0 are not captured
       separately: RESP is already a CMT compartment and therefore
       already a reported column.
    7. Add one new annotated PARAM block: H0, beta_r, U, END.
    8. Add one new INIT block: p11 = 1 (not CMT -- see the naming
       rule below).

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro (and ERANGE/EDOM-style names are libc macros
  too) on macOS/BSD-derived toolchains and the model fails to compile.
  Use `TEVT`/`event_found` instead (see exponential_ode.cpp). None of
  H0/beta_r/RESP0/CL/V2/Q/V3/KA/KA2/KIN/KOUT/IC50/IMAX/n/VMAX/KM
  collide with an errno.h name.

  BETWEEN-SUBJECT VARIABILITY (reports/11_bsv_review.md,
  reports/12_bsv_implementation_report.md): `sim_tte_ode(model =
  "irm1_hazard", omega = ...)` applies BSV via per-subject `idata`
  columns, not a declared OMEGA block -- see `.ODE_BSV_TARGETS` in
  R/helpers.R. Valid targets (every structural PK/PD parameter of
  this backbone, taken from this file's own $PARAM block below, in
  this order): CL, V2, Q, V3, KA, KA2, KIN, KOUT, IC50, IMAX, n, VMAX,
  KM. H0/lp/beta_r/U/END (the hazard-scaffold parameters) are not BSV
  targets.

[PARAM] @annotated
  CL   :  1  : Clearance (volume/time)
  V2   : 20  : Central volume (volume)
  Q    :  2  : Inter-compartmental clearance (volume/time)
  V3   : 10  : Peripheral volume of distribution (volume)
  KA   :  1  : Absorption rate constant 1 (1/time)
  KA2  :  1  : Absorption rate constant 2 (1/time)
  KIN  : 10  : Response in rate constant (1/time)
  KOUT :  2  : Response out rate constant (1/time)
  IC50 :  2  : Concentration for 50% of max inhibition (mass/volume)
  IMAX :  1  : Maximum inhibition
  n    :  1  : Emax model sigmoidicity
  VMAX :  0  : Maximum reaction velocity (mass/time)
  KM   :  2  : Michaelis constant (mass/volume)

[PARAM] @annotated
  H0     : 0.01 : Baseline (RESP == RESP0, lp = 0) hazard (1/time)
  lp     : 0    : Log hazard ratio (linear predictor); composes with beta_r * (RESP/RESP0 - 1) for sim_tte_ode()'s covariates/beta mechanism (see ?sim_tte_ode "Time-varying covariates")
  beta_r : 0    : Log hazard ratio per unit relative response deviation; HAZ(t) = H0 * exp(lp + beta_r * (RESP(t)/RESP0 - 1))
  U      : 0    : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END    : 1e9  : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[CMT] @annotated
  EV     : First extravascular compartment (mass)
  CENT   : Central compartment (mass)
  PERIPH : Peripheral compartment (mass)
  RESP   : Response compartment
  EV2    : Second extravascular compartment (mass)

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  #define CP (CENT/V2)
  #define CT (PERIPH/V3)
  #define CLNL (VMAX/(KM+CP))
  #define INH (IMAX*pow(CP,n)/(pow(IC50,n)+pow(CP,n)))
  // -- BEGIN simtte survival scaffolding -----------------------------
  static int    event_found = 0;
  static double TEVT        = 0.0;
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;
  // Pre-dose response baseline, copied once per subject from RESP_0
  // (see MAIN below); the hazard link is relative to this value.
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
  dxdt_RESP   =  KIN*(1-INH) - KOUT*RESP;
  // -- BEGIN simtte survival scaffolding -----------------------------
  double HAZ = H0 * exp(lp + beta_r * (RESP/RESP0 - 1.0));   // <- the only model-specific line
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
  CP : Plasma concentration (mass/volume)

[CAPTURE] @annotated
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
  HAZ         : Instantaneous hazard, H0 * exp(lp + beta_r * (RESP/RESP0 - 1))
