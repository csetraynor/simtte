Model file:  weibull_tv.txt
[PROB]
# Model: `Time-varying-lp Weibull parametric proportional hazard model`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Piecewise-constant lp(t) on its own supplied grid; cumulative hazard
  by segment:
    H(t_k) = sum_j exp(mu + lp_j) * (t_j^shape - t_{j-1}^shape)
    p11    = exp(-H(t))
  See the "Time-varying Weibull model" section of PHASE_G_REPORT.md for
  the full derivation and the empirical verification that a $TABLE-only
  model sees exactly the "closing" (nocb = FALSE / LOCF) covariate value
  needed for this formula, at every reported row -- no $ODE/$DES is
  used, deliberately: this preserves the closed-form, solver-free
  numerical robustness established in Phase A for the constant-lp
  baseline (inst/models/weibull.cpp, which this file does not modify or
  replace), and avoids every failure mode documented for internal
  solver evaluations in PHASE_E_THRESHOLD_TRACKING_REPORT.md.

[PARAM] @annotated
  lp    :  0    : linear predictor (time-varying covariate)
  mu    :  0.1  : intercept
  shape :  1    : shape parameter

[GLOBAL]
  static double cum_H = 0.0;
  static double prev_time = 0.0;

[MAIN]
  if (NEWIND <= 1) {
    cum_H = 0.0;
    prev_time = 0.0;
  }

[TABLE]
  // Called once per reported row, in increasing TIME order, for each
  // individual (never at internal/rejected solver evaluations -- there
  // is no ODE here at all). `lp` is supplied as a nocb = FALSE
  // (last-observation-carried-forward) time-varying covariate, exactly
  // like the M-spline model's `basehaz`. Verified empirically that at a
  // row whose TIME exactly coincides with a new lp knot, the value
  // visible here is the OLD (pre-update) value -- i.e. exactly the
  // value that was active over the segment ENDING at this row, which is
  // precisely what the summation below needs; the NEW value only
  // becomes visible starting at the next reported row after the knot.
  // This mirrors, and was verified against, the already-established
  // M-spline hazard-carry convention (nocb = FALSE): "the value at
  // time[i] applies from time[i] until time[i+1]" is the mirror-image
  // statement of "the value seen when computing the segment ending at
  // time[i] is the value active up to but not including time[i]".
  //
  // Because the internal mrgsolve output grid used by .sim_surv_df()
  // for this model is the union of the user's requested `time` grid and
  // every supplied lp(t) knot time (see R/simtte.R), no segment is ever
  // skipped: consecutive $TABLE calls are always adjacent points on
  // that merged grid, so at most one true segment lies between them.
  double p11;
  if (TIME <= 0.0) {
    p11 = 1.0;
    cum_H = 0.0;
    prev_time = 0.0;
  } else {
    double p_i = pow(TIME, shape);
    double p_prev = pow(prev_time, shape);
    double diff = p_i - p_prev;
    // Defensive clamp: TIME > prev_time always holds by construction of
    // the merged, deduplicated, sorted output grid, so this difference
    // is mathematically non-negative; the clamp only guards against
    // floating-point noise when TIME and prev_time are extremely close,
    // which would otherwise risk log() of a tiny negative number.
    if (diff < 0.0) diff = 0.0;
    // Log-space, exactly as in the constant-lp baseline
    // (inst/models/weibull.cpp): defers the only exponentiation to the
    // very end, so a large finite (mu + lp) cannot overflow the
    // intermediate eta = exp(mu + lp) before being combined with a
    // possibly small segment width. diff == 0 gives log(diff) = -Inf,
    // hence delta_H = exp(-Inf) = 0 -- a zero-width or first (t = 0)
    // segment contributes nothing, with no NaN anywhere.
    double log_delta_H = (mu + lp) + log(diff);
    double delta_H = exp(log_delta_H);
    cum_H += delta_H;
    p11 = exp(-cum_H);
    prev_time = TIME;
  }

[CAPTURE]
  p11
