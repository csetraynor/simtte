Model file:  pkpd_idr_hazard.txt
[PROB]
# Model: `Oral PK / indirect-response mediator hazard (shipped example)`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  This is simtte's flagship shipped example of the general
  "survival-compartment co-integration" pattern
  (dxdt_p11 = -p11 * HAZ) that a user's own custom mrgsolve model is
  expected to follow when using sim_tte_df() directly (see the package
  vignettes and PHASE_H_IMPLEMENTATION_PLAN.md). It is the same model
  described in the manuscript's "ODE-coupled PK/PD hazard" section,
  reproduced here verbatim (same compartments, parameters, and
  equations) as a standalone, loadable, user-editable file instead of
  an inline code string, so it can actually be discovered and used
  without retyping it from the paper.

  Structure: a one-compartment oral-absorption PK model (GUT -> CENT)
  feeds an indirect-response (turnover) model for a mediator R, whose
  production is inhibited by drug concentration through an Emax-form
  term INH(t). The event hazard is proportional to the mediator level
  relative to its untreated steady state (HAZ = H0 * R / R0), so
  suppressing R lowers the hazard below its baseline value H0.
  Between-subject variability (log-normal) is placed on clearance,
  volume, and IC50 via $OMEGA.

  This file is loaded via simtte_example_model("pkpd_idr_hazard"), not
  through the package's internal .read_model_static_cache() (which
  only ever loads the three built-in weibull/weibull_tv/ms models used
  by sim_tte()). It is a user-facing template/example, not part of the
  core survival-simulation engine, and is never used internally by
  sim_tte() or sim_tte_df().

  Intended workflow (see ?simtte_example_model):
    mod  <- simtte_example_model("pkpd_idr_hazard")
    data <- mrgsolve::expand.ev(ID = 1:n, amt = 100, cmt = 1, ii = 1,
               addl = 23, time = 0)
    out  <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = 24,
               delta = 0.1))
    sim_tte_df(out[, c("ID", "time", "p11")])

[PARAM] @annotated
  CL   :  2    : Clearance (volume/time)
  V    :  20   : Central volume of distribution (volume)
  KA   :  1    : First-order absorption rate constant (1/time)
  KIN  :  1    : Mediator (R) zero-order production rate (1/time)
  KOUT :  0.5  : Mediator (R) first-order elimination rate constant (1/time)
  IC50 :  1    : Concentration for 50% inhibition of R production (mass/volume)
  IMAX :  0.9  : Maximum fractional inhibition of R production
  H0   :  0.15 : Baseline event hazard (1/time), attained when R = R0

[CMT] @annotated
  GUT  : First-order absorption depot (mass)
  CENT : Central compartment (mass)
  R    : Mediator level driving the event hazard

[INIT]
  p11 = 1

[OMEGA] @labels ECL EV EIC50
  0.09 0.09 0.16

[MAIN]
  // Log-normal between-subject variability on clearance, volume, and
  // potency (IC50): CV% approx sqrt(exp(omega^2) - 1), i.e. ~30% for
  // ECL/EV (omega^2 = 0.09) and ~42% for EIC50 (omega^2 = 0.16).
  double CLi   = CL * exp(ECL);
  double Vi    = V  * exp(EV);
  double IC50i = IC50 * exp(EIC50);
  R_0 = KIN / KOUT; // mediator steady state at t = 0, untreated

[ODE]
  // HAZ is recomputed at every solver evaluation from the current
  // mediator level R, so the survival compartment p11 stays
  // numerically consistent with the pharmacological state throughout
  // the integration (the general pattern documented for user-defined
  // mechanistic hazard models in the package vignettes).
  double CP  = CENT / Vi;
  double INH = IMAX * CP / (IC50i + CP);
  double HAZ = H0 * (R / (KIN / KOUT));
  dxdt_GUT  = -KA * GUT;
  dxdt_CENT =  KA * GUT - (CLi / Vi) * CENT;
  dxdt_R    =  KIN * (1 - INH) - KOUT * R;
  dxdt_p11  = -p11 * HAZ;

[CAPTURE] @annotated
  CP  : Plasma drug concentration (mass/volume)
  HAZ : Instantaneous event hazard (1/time)
