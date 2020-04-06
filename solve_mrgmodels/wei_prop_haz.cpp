Model file:  wei_prop_haz.cpp
[PROB]
# Model: `Simulate Weibull parametric proportional hazard model`
  - Forward Kolmogorov differential equation
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

[PARAM] @annotated
  lp   :  0.15  : linear predictor
  mu   :  0.1   : intercept
  shape :  1     : shape parameter

[INIT]
  p11 = 1

[MAIN]
  double eta = exp(mu + lp);

[ODE]
  if(SOLVERTIME > 10E-10) {
    dxdt_p11 =  - p11 * shape * pow (SOLVERTIME * eta, shape - 1) * eta; // note af * time = caf
  } else {
    dxdt_p11 = 0;
  }
