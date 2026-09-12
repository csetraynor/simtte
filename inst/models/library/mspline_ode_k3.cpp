Model file:  mspline_ode_k3.cpp
[PROB]
# Model: `M-spline baseline hazard (3 interior knots), in-solver event detection`
  - Author: Carlos Traynor
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`

  Hazard function:   h(t) = eta * sum_m c_m * M_m(t), eta = exp(mu + lp)
  M_m(t): degree-2 (quadratic) normalized M-spline basis functions
  (Ramsay 1988), 3 interior knots (`k1`-`k3`) plus boundary knots
  (`bk_lo`, `bk_hi`), evaluated by the Cox-de Boor recursion below --
  6 basis functions / coefficients (`c1`-`c6`) total (interior knots +
  degree + 1).

  CONVENTION (reports/09_phase3_report.md section 1): matches
  `splines2::mSpline(x, knots = <interior>, Boundary.knots = c(bk_lo,
  bk_hi), degree = 2, intercept = TRUE)` exactly -- the recursion below
  was verified against that reference to floating-point precision
  (max abs diff ~3e-17) across the full domain before being shipped.
  This is the same convention `sim_tte(type = "ms")`'s own example data
  was generated with (see the M-spline demonstration script referenced
  from `?ms_data`), so knots/coefficients are directly portable between
  `sim_tte(type = "ms")` and `sim_tte_ode(model = "mspline")`.

  This is a FIXED-VARIANT library file (one of three shipped interior-
  knot counts: 3, 5, 7; `reports/02_technical_design.md` section 4
  option 1's "moderate, fixed set" compromise). Users needing an
  arbitrary knot count use the existing, piecewise-constant
  `sim_tte(type = "ms")` path instead; see `?sim_tte_ode` "M-spline
  knot-count variants".

  NO TIME FLOOR NEEDED (unlike `weibull_ode.cpp`'s risk R1): the
  M-spline basis is bounded and finite everywhere on `[bk_lo, bk_hi]`,
  including exactly at the boundary (verified directly: the first basis
  function at `t = bk_lo` evaluates to a finite, nonzero value, never a
  singularity) -- there is no `t -> 0+`-style divergence to guard
  against, so `SOLVERTIME` is used directly, unfloored, inside the
  hazard expression.

  DOMAIN: the basis (and therefore the hazard) is identically zero for
  `t` outside `[bk_lo, bk_hi]` (a consequence of the clamped-knot
  construction, not a special case coded here). `sim_tte_ode()` requires
  `end <= bk_hi` at the R level specifically so this boundary is never
  silently crossed during a real call -- see `?sim_tte_ode`.

  NAMING RULE (copy into every future library model): never name a
  tracking variable in the global-state block ETIME -- `<sys/errno.h>`
  defines it as a macro on macOS/BSD-derived toolchains and the model
  fails to compile. Use `TEVT`/`event_found` instead (see
  `exponential_ode.cpp`). The M-spline order macro below is named
  `MS_ORDER`, not `ERANGE`/`EDOM`/any other libc errno-family name --
  checked against `<cerrno>`'s macro list before use.

[PARAM] @annotated
  mu    :  -1  : Intercept
  lp    :  0   : Log hazard ratio (linear predictor); eta = exp(mu + lp)
  bk_lo :  0   : M-spline lower boundary knot
  bk_hi :  20  : M-spline upper boundary knot
  k1    :  5   : M-spline interior knot 1
  k2    :  10  : M-spline interior knot 2
  k3    :  15  : M-spline interior knot 3
  c1    :  1   : M-spline coefficient 1 (must be >= 0)
  c2    :  1   : M-spline coefficient 2 (must be >= 0)
  c3    :  1   : M-spline coefficient 3 (must be >= 0)
  c4    :  1   : M-spline coefficient 4 (must be >= 0)
  c5    :  1   : M-spline coefficient 5 (must be >= 0)
  c6    :  1   : M-spline coefficient 6 (must be >= 0)
  U     :  0   : Uniform(0,1) draw for in-solver event detection, per subject via idata (see ?sim_tte_ode)
  END   :  1e9 : Administrative censoring horizon for the boundary guard, per subject via idata (see ?sim_tte_ode)

[INIT]
  // `p11` is not declared in a separate $CMT block; see exponential_ode.cpp.
  p11 = 1

[GLOBAL]
  // -- BEGIN simtte survival scaffolding -----------------------------
  static int    event_found = 0;
  static double TEVT        = 0.0;
  // Grid-free refinement bracket (Phase 2.5, reports/04_author_decisions.md
  // "After the test runbook / Phase 2.5"): the solver's own last
  // pre-crossing evaluation (T_PRE, P_PRE) and the crossing evaluation
  // itself (TEVT, P_POST), one internal solver step apart -- narrow
  // enough that interpolating between them makes the constant-hazard-
  // within-bracket assumption close to exact for any smooth hazard,
  // independent of the reported output grid (`delta`).
  static double T_PRE       = 0.0;
  static double P_PRE       = 1.0;
  static double P_POST      = 0.0;

  // -- M-spline basis (Phase 3, reports/09_phase3_report.md) ---------
  // Degree-2 normalized M-spline basis via the Cox-de Boor recursion
  // (Ramsay 1988 scaling, not B-spline scaling), evaluated against a
  // clamped (open-uniform) extended knot vector: `MS_ORDER` copies of
  // `bk_lo`, the interior knots in order, `MS_ORDER` copies of `bk_hi`.
  // Identical in every mspline_ode_k*.cpp variant except the fixed
  // array sizes below (n_interior); expected duplication, not
  // consolidated (see reports/09_phase3_report.md section on the
  // ponytail review of this file set).
  #define MS_DEGREE 2
  #define MS_ORDER (MS_DEGREE + 1)

  inline void mspline_basis(double t, const double* interior,
      int n_interior, double bk_lo, double bk_hi, double* out) {
    int n_ext = n_interior + 2 * MS_ORDER;
    double ext[32];
    int idx = 0;
    for (int j = 0; j < MS_ORDER; j++) ext[idx++] = bk_lo;
    for (int j = 0; j < n_interior; j++) ext[idx++] = interior[j];
    for (int j = 0; j < MS_ORDER; j++) ext[idx++] = bk_hi;
    int n0 = n_ext - 1;
    double M[32][MS_ORDER];
    for (int i = 0; i < n0; i++) {
      double lo = ext[i];
      double hi = ext[i + 1];
      double v = 0.0;
      if (hi > lo) {
        // Right-closed only for the segment reaching the overall domain
        // maximum (ext[n_ext - 1]) -- NOT simply the last array index,
        // which (under a clamped/repeated-boundary-knot extended knot
        // vector) is a zero-width degenerate segment, not the last
        // segment with positive width. Verified against
        // splines2::mSpline() exactly at t == bk_hi before this fix was
        // applied (reports/09_phase3_report.md section 1).
        if (t >= lo && (t < hi || (hi == ext[n_ext - 1] && t <= hi))) {
          v = 1.0 / (hi - lo);
        }
      }
      M[i][0] = v;
    }
    for (int k = 2; k <= MS_ORDER; k++) {
      int n_k = n0 - k + 1;
      for (int i = 0; i < n_k; i++) {
        double denom = ext[i + k] - ext[i];
        double val = 0.0;
        if (denom > 0) {
          double left = M[i][k - 2] * (t - ext[i]);
          double right = (i + 1 < n0) ? M[i + 1][k - 2] * (ext[i + k] - t) : 0.0;
          val = k * (left + right) / ((k - 1) * denom);
        }
        M[i][k - 1] = val;
      }
    }
    int n_basis = n_interior + MS_ORDER;
    for (int m = 0; m < n_basis; m++) out[m] = M[m][MS_ORDER - 1];
  }

[MAIN]
  if (NEWIND <= 1) {
    event_found = 0;
    TEVT = 0.0;
    T_PRE = 0.0;
    P_PRE = 1.0;
    P_POST = 0.0;
  }

[ODE]
  double eta = exp(mu + lp);
  double interior[3] = {k1, k2, k3};
  double basis[6];
  mspline_basis(SOLVERTIME, interior, 3, bk_lo, bk_hi, basis);
  double coef[6] = {c1, c2, c3, c4, c5, c6};
  double basehaz = 0.0;
  for (int m = 0; m < 6; m++) basehaz += coef[m] * basis[m];
  double HAZ = eta * basehaz;          // <- the only model-specific line
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
  TEVT        : Latched in-solver event time (SOLVERTIME at first p11 <= U), or 0 if not yet found
  event_found : One-shot latch flag (1 once an event has been detected, 0 otherwise)
  T_PRE       : Last pre-crossing solver evaluation time (refinement bracket lower end)
  P_PRE       : p11 at T_PRE (refinement bracket lower end)
  P_POST      : p11 at TEVT (refinement bracket upper end), or 0 if not yet found
  HAZ         : Instantaneous hazard eta * sum_m c_m * M_m(SOLVERTIME), captured for direct convention-equivalence validation against the R-side splines2::mSpline() basis
