Model file:  pkpd_linear_hazard.txt
[PROB]
# Model: `IV-bolus PK / linear-drug-effect hazard (shipped example)`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  A minimal, single-compartment instance of simtte's general
  "survival-compartment co-integration" pattern
  (dxdt_p11 = -p11 * HAZ): one IV-bolus PK compartment, and a hazard
  that decreases linearly in the current drug concentration from a
  constant baseline H0. No between-subject variability and no
  turnover/indirect-response dynamics -- deliberately the simplest
  possible working instance of the pattern, intended as a short,
  readable starting point to copy and modify (see
  pkpd_idr_hazard.cpp for a more realistic, mechanistic example).

  Because both the PK and the hazard are linear, this model has an
  independent closed-form cumulative hazard and survival function,
  which is what makes it the package's validation anchor for this
  example library (see inst/validation/03_pkpd_mechanism_validation.R,
  which validates exactly this model -- loaded from this file, not a
  duplicated inline definition -- against that closed form):

    C(t) = (Dose / V) * exp(-k * t),      k = CL / V
    HAZ(t) = H0 - SLOPE * C(t)            (requires H0 >= SLOPE * C(0))
    H(t) = integral_0^t HAZ(s) ds
         = H0 * t - SLOPE * (Dose / V) * (1 - exp(-k * t)) / k
    S(t) = exp(-H(t))

  This file is loaded via simtte_example_model("pkpd_linear_hazard"),
  not through the package's internal .read_model_static_cache() (which
  only ever loads the three built-in weibull/weibull_tv/ms models used
  by sim_tte()). It is a user-facing template/example, not part of the
  core survival-simulation engine, and is never used internally by
  sim_tte() or sim_tte_df().

  Intended workflow (see ?simtte_example_model):
    mod  <- simtte_example_model("pkpd_linear_hazard")
    data <- mrgsolve::ev(amt = 100, cmt = 1, time = 0)
    out  <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = -1,
               add = seq(0, 60, by = 0.5), obsonly = TRUE))
    sim_tte_df(out[, c("ID", "time", "p11")])

[PARAM] @annotated
  CL    :  1    : Clearance (volume/time)
  V     :  10   : Central volume of distribution (volume)
  H0    :  0.3  : Baseline event hazard (1/time), attained at zero concentration
  SLOPE :  0.02 : Reduction in hazard per unit drug concentration ((1/time)/(mass/volume))

[CMT] @annotated
  CENT : Central compartment (mass)

[INIT]
  p11 = 1

[MAIN]
  double K = CL / V;

[ODE]
  double C   = CENT / V;
  double HAZ = H0 - SLOPE * C;
  dxdt_CENT = -K * CENT;
  dxdt_p11  = -p11 * HAZ;

[CAPTURE] @annotated
  C   : Plasma drug concentration (mass/volume)
  HAZ : Instantaneous event hazard (1/time)
