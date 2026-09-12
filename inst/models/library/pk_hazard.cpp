Model file:  pk_hazard.cpp
[PROB]
# Model: `Concentration-driven hazard, in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Backbone: `mrgsolve::modlib("pk2cmt")` (two-compartment PK, dual
  first-order absorption, optional nonlinear clearance from `CENT`) --
  every PARAM/CMT/GLOBAL(PK part)/ODE(PK part)/CAPTURE(PK part) block
  below is that backbone's own code, unedited. (NOTE: this doc block
  intentionally spells block names without a leading dollar sign and
  avoids the literal "at-annotated" tag text -- mrgsolve's spec scraper
  scans the whole file for those tokens, not just line starts, and a
  prose mention of them elsewhere in this file breaks parsing; verified
  directly while drafting this file.)

  Hazard link (log-linear in plasma concentration):
    HAZ = H0 * exp(beta_cp * CP), CP = CENT / V2 (already defined by the
    backbone's own GLOBAL block).
  `beta_cp = 0` reduces exactly to the exponential model
  (`exponential_ode.cpp`) with hazard `H0`.

  MECHANICAL EDITS applied to the `pk2cmt` backbone to reach this file
  (the same edit list, in the same order, is applied to every
  `*_hazard.cpp` model in this library -- this list is the working
  specification for the future model-converter utility,
  reports/04_author_decisions.md section 6):
    1. Copy the backbone's PARAM block(s) verbatim.
    2. Copy the backbone's CMT block verbatim.
    3. Copy the backbone's GLOBAL block verbatim, then append the
       survival scaffold (event_found/TEVT/T_PRE/P_PRE/P_POST) to the
       end of that same block (GLOBAL may only appear once per model).
    4. Copy the backbone's MAIN block verbatim (if any), then append
       the survival scaffold's reset logic to the end of that same
       block (MAIN may only appear once per model); if the backbone has
       no MAIN block, add one containing only the reset logic.
    5. Copy the backbone's ODE block verbatim, then append: one
       model-specific `double HAZ = ...;` line (the only line that
       differs across the six PK/PD hazard models), `dxdt_p11 = -p11 *
       HAZ;`, and the T_PRE/P_PRE/TEVT/P_POST latch (identical in every
       model in this package).
    6. Copy the backbone's CAPTURE block verbatim, then add a second
       annotated CAPTURE block for TEVT/event_found/T_PRE/P_PRE/
       P_POST/HAZ (CAPTURE, unlike GLOBAL/MAIN, may appear more than
       once and is concatenated -- verified directly before building
       this file).
    7. Add one new annotated PARAM block: H0, beta_cp, U, END.
    8. Add one new INIT block: p11 = 1 (not CMT -- see the naming
       rule below).

  KNOWN LIMITATION (disclosed, not fixed here -- reports/10_phase4_report.md
  open risks): none of the six modlib() backbones wire ETA(n) into any
  parameter themselves, and none declares an OMEGA block -- so
  `sim_tte_ode(model = "pk_hazard", omega = <matrix>)` cannot add
  between-subject variability to this model as shipped (verified
  directly, including against an unmodified mrgsolve::mread() of
  pk2cmt itself: mrgsolve's `omat(mod, matrix)` only updates an
  *already-declared* OMEGA block, it does not create one from nothing).
  `sim_tte_ode()` now raises an informative error for this case instead
  of mrgsolve's own cryptic "improper signature: omat" (see
  R/sim_tte_ode.R). Wiring in genuine IIV would need a non-mechanical,
  per-model edit (renaming the backbone's own CL/V2 references, since
  mrgsolve declares PARAM-sourced names as const references that
  cannot be reassigned in MAIN under the same name) -- out of scope for
  the mechanical, verbatim-backbone edit list above.

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- <sys/errno.h>
  defines it as a macro (and ERANGE/EDOM-style names are libc macros
  too) on macOS/BSD-derived toolchains and the model fails to compile.
  Use `TEVT`/`event_found` instead (see exponential_ode.cpp). None of
  H0/beta_cp/CP/CL/V2/Q/V3/KA/KA2/VMAX/KM collide with an errno.h name.

  BETWEEN-SUBJECT VARIABILITY (reports/11_bsv_review.md,
  reports/12_bsv_implementation_report.md): `sim_tte_ode(model =
  "pk_hazard", omega = ...)` applies BSV via per-subject `idata`
  columns, not a declared OMEGA block (this file declares none,
  unchanged from Phase 4) -- see `.ODE_BSV_TARGETS` in R/helpers.R.
  Valid targets (every structural PK parameter of this backbone,
  taken from this file's own $PARAM blocks below, in this order):
  CL, V2, Q, V3, KA, KA2, VMAX, KM. H0/lp/beta_cp/U/END (the hazard-
  scaffold parameters) are not BSV targets.

[PARAM] @annotated @input
  CL   :  1  : Clearance (volume/time)
  V2   : 20  : Central volume (volume)
  Q    :  2  : Inter-compartmental clearance (volume/time)
  V3   : 10  : Peripheral volume of distribution (volume)
  KA   :  1  : Absorption rate constant 1 (1/time)

[PARAM] @annotated
  KA2  :  1  : Absorption rate constant 2 (1/time)
  VMAX :  0  : Maximum velocity (mass/time)
  KM   :  2  : Michaelis Constant (mass/volume)

[PARAM] @annotated
  H0      : 0.01 : Baseline (concentration = 0, lp = 0) hazard (1/time)
  lp      : 0    : Log hazard ratio (linear predictor); composes with beta_cp * CP for sim_tte_ode()'s covariates/beta mechanism (see ?sim_tte_ode "Time-varying covariates")
  beta_cp : 0    : Log hazard ratio per unit plasma concentration; HAZ(t) = H0 * exp(lp + beta_cp * CP(t))
  U       : 0    : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END     : 1e9  : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[CMT] @annotated
  EV     : First extravascular compartment (mass)
  CENT   : Central compartment (mass)
  PERIPH : Peripheral compartment (mass)
  EV2    : Second extravascular compartment (mass)

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  #define CP (CENT/V2)
  #define CT (PERIPH/V3)
  #define CLNL (VMAX/(KM+CP))
  // -- BEGIN simtte survival scaffolding -----------------------------
  static int    event_found = 0;
  static double TEVT        = 0.0;
  // Grid-free refinement bracket (Phase 2.5, reports/04_author_decisions.md
  // "After the test runbook / Phase 2.5"): the solver's own last
  // pre-crossing evaluation (T_PRE, P_PRE) and the crossing evaluation
  // itself (TEVT, P_POST), one internal solver step apart.
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;
  // -- END simtte survival scaffolding (GLOBAL half) ------------------

[MAIN]
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
  dxdt_EV     = -KA*EV;
  dxdt_EV2    = -KA2*EV2;
  dxdt_CENT   =  KA*EV + KA2*EV2 - (CL+CLNL+Q)*CP  + Q*CT;
  dxdt_PERIPH =  Q*CP - Q*CT;
  // -- BEGIN simtte survival scaffolding -----------------------------
  double HAZ = H0 * exp(lp + beta_cp * CP);   // <- the only model-specific line
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
  HAZ         : Instantaneous hazard, H0 * exp(lp + beta_cp * CP)
