Model file:  weibull.txt
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
if(SOLVERTIME >= 1E-1) {
  dxdt_p11 =  - p11 * shape * pow (SOLVERTIME, shape - 1) * eta;
} else {
  dxdt_p11 = - p11 * eta; // approximate via exponential model
}
