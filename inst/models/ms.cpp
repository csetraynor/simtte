Model file:  ms_prop_haz.cpp
[PROB]
# Model: `M-spline proportional hazard model`
- Author: Carlos S Traynor
- Date: `r Sys.Date()`
- Version: `r packageVersion("mrgsolve")`

[PARAM] @annotated
  lp     :  0.1 : linear predictor
  mu      :  -9  : intercept
  basehaz :  1   : M-splines baseline hazard

[INIT]
  p11 = 1

[ODE]
  // eta = exp(mu + lp) is recomputed at every solver evaluation (not
  // cached once in $MAIN, which fires only once per individual at
  // simulation start) so that a time-varying `lp` covariate -- carried
  // forward under the same nocb = FALSE convention as `basehaz` -- is
  // actually reflected during integration. For the constant-lp case
  // (lp supplied as a single, unchanging value, as in every existing
  // use of this model), eta is identical at every evaluation either
  // way, so this is a pure relocation with no numerical effect on
  // existing behavior (verified in test-ms-hazard-carry.R and
  // PHASE_G_REPORT.md).
  double eta = exp(mu + lp);

  if(SOLVERTIME > 10E-10) {
    dxdt_p11 =  - p11 * basehaz * eta;
  } else {
    dxdt_p11 = 0;
  }
