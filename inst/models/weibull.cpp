Model file:  weibull.txt
[PROB]
# Model: `Simulate Weibull parametric proportional hazard model`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Survival function:  S(t) = exp[-eta * t^shape],  eta = exp(mu + lp)
  Hazard function:    h(t) = eta * shape * t^(shape - 1)

[PARAM] @annotated
  lp   :  0.15  : linear predictor
  mu   :  0.1   : intercept
  shape :  1     : shape parameter

[TABLE]
  // The Weibull survival function has a closed form that depends only
  // on (mu, lp, shape) and the current output TIME; it is evaluated
  // directly here rather than by integrating dS/dt = -h(t) * S(t) as an
  // ODE compartment.
  //
  // Why not integrate the ODE: the instantaneous hazard
  // h(t) = eta * shape * t^(shape - 1) diverges as t -> 0+ whenever
  // shape < 1. Although the *cumulative* hazard (and hence S(t)) stays
  // perfectly finite for every t > 0, a numerical ODE solver has to
  // evaluate h(t) at times arbitrarily close to 0 while probing step
  // sizes, and no finite step size keeps the local error below
  // tolerance there: the adaptive solver shrinks its step toward zero
  // and the integration fails outright (verified empirically for
  // shape < 1; this is not a hypothetical edge case). Any fix that
  // keeps the ODE formulation has to either (a) cap h(t) near t = 0,
  // trading a tunable, shape/eta-dependent amount of bias for
  // stability, or (b) special-case a neighbourhood of t = 0 with an
  // unrelated approximation, which is exactly the scientifically
  // incorrect behaviour this model previously had (a hard-coded
  // exponential-hazard branch for SOLVERTIME < 0.1, regardless of
  // shape).
  //
  // Because the exact solution is available in closed form and does
  // not depend on any other ODE state, computing it directly is exact
  // to floating-point precision at every reported time, for every
  // shape > 0, with no tuning parameter and no risk of solver failure.
  // This is the numerically safest option for this specific model.
  // (User-defined mrgsolve hazard models that lack a closed form
  // should keep using the dxdt_p11 = -p11 * HAZ ODE pattern described
  // in the package vignettes.)
  //
  // Numerical robustness for extreme but finite (mu + lp): the naive
  // expression eta = exp(mu + lp); p11 = exp(-eta * pow(TIME, shape))
  // computes eta as an intermediate value that can itself overflow to
  // Inf for large finite (mu + lp) (e.g. mu + lp > ~709.78, the double
  // overflow threshold), *before* it is multiplied by pow(TIME, shape).
  // At TIME == 0, pow(TIME, shape) == 0 exactly, so that product
  // becomes Inf * 0 = NaN, even though the true survival is exactly
  // S(0) = 1 for every valid (finite mu, lp, shape > 0). Two changes
  // fix this without approximation:
  //   1. TIME == 0 is handled directly: p11 = 1 always, without ever
  //      touching eta or pow().
  //   2. For TIME > 0, the calculation works in log-cumulative-hazard
  //      space: logH = (mu + lp) + shape * log(TIME), so the two
  //      finite terms are summed *before* exponentiating, rather than
  //      exponentiating (mu + lp) alone first. This defers the only
  //      exponentiation to the very end. If logH is large enough that
  //      exp(logH) overflows to Inf, p11 = exp(-Inf) = 0, which is the
  //      mathematically correct saturated answer (not NaN); if logH is
  //      very negative, exp(logH) underflows to 0 and p11 = exp(-0) =
  //      1, again correct.
  double p11;
  if (TIME <= 0.0) {
    p11 = 1.0;
  } else {
    double logH = (mu + lp) + shape * log(TIME);
    p11 = exp(-exp(logH));
  }

[CAPTURE]
  p11
