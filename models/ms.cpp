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

[MAIN]
  double eta = exp(mu + lp);

[ODE]

  if(SOLVERTIME > 10E-10) {
    dxdt_p11 =  - p11 * basehaz * eta;
  } else {
    dxdt_p11 = 0;
  }
