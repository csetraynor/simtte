# simtte 1.1.0

New features, fully backward compatible with CRAN 1.0.2 (verified
directly, `inst/validation/13_backcompat_v1_0_2.R`).

## New features

* **`sim_tte_ode()`: joint PK/PD and time-to-event simulation with
  in-solver event detection.** The survival probability is added to a
  PK/PD `mrgsolve` model as one more ODE compartment (`dxdt_p11 = -p11
  * HAZ`); each subject's event time is located *during* ODE
  integration, at the first internal solver evaluation where `p11`
  falls to or below a per-subject uniform draw `U`, rather than
  resolved afterward from a pre-simulated trajectory (the mechanism
  `sim_tte_df()` already used). Built-in models: `"exponential"`,
  `"weibull"`, `"gompertz"`, and `"mspline"` (a quadratic M-spline
  baseline hazard evaluated continuously inside the compiled model, in
  one of three shipped interior-knot-count variants: 3, 5, or 7).
* **Six built-in PK/PD-linked hazard models**: `"pk_hazard"`
  (concentration-driven), `"irm1_hazard"`-`"irm4_hazard"`
  (response-driven, the four indirect-response model types),
  `"tmdd_hazard"` (driven by drug-target complex, target-mediated drug
  disposition). Each reuses an unedited `mrgsolve::modlib()` PK/PD
  backbone plus a hazard link; dosing is ordinary `mrgsolve` event data
  via `sim_tte_ode()`'s `data` argument.
* **Time-varying covariates for `sim_tte_ode()`** (`covariates`/`beta`):
  an arbitrary named set of log-linear covariates, `lp(t) = sum_k
  beta_k * X_k(t)`, generalizing `sim_tte()`'s single-covariate
  `lp_data` mechanism.
* **Between-subject variability (`omega`) for the PK/PD hazard
  library**: a log-normal draw on any subset of a model's own
  structural PK/PD parameters, applied via per-subject `idata` columns
  -- no model file needs a declared `$OMEGA` block. A user-supplied
  model that already declares its own `$OMEGA` (the classic
  `TVCL`/`ETA()` idiom) goes through `mrgsolve::omat()` instead,
  automatically.
* **`tte_model()`: convert your own `mrgsolve` PK/PD model** for use
  with `sim_tte_ode()`, instead of being limited to the six built-in
  models. Accepts model code, a file path, or a compiled model; adds
  the survival compartment and in-solver event-detection scaffolding
  without editing the model's own `$ODE`/`$DES` content (only
  `$GLOBAL`/`$MAIN`, which `mrgsolve` allows only one of each, are
  edited in place). Verified to reproduce all six built-in hazard
  models exactly when converting their own `mrgsolve::modlib()`
  backbones.
* **`add_censoring()`: independent right censoring**, on top of any
  already-simulated events data frame (`sim_tte_ode()`, `sim_tte()`, or
  `sim_tte_df()` output alike). Draws a per-subject censoring time from
  an exponential, Weibull, uniform, lognormal, gamma, or user-supplied
  distribution and takes `min(event time, censoring time, administrative
  end)`; `sim_status`/`sim_time` are updated in place and a new
  `sim_reason` column (`"event"`/`"censored"`/`"administrative"`)
  records why. `censoring_rate_for()` solves for the distribution
  parameter (rate, Weibull scale at a fixed shape, lognormal meanlog at
  a fixed sdlog, or gamma rate at a fixed shape) giving a target
  censoring fraction, from the event times of an already-run uncensored
  simulation. `sim_tte_ode()` also accepts a `censoring = ` argument
  that applies this automatically, inside its own seeded draw, so one
  `seed` reproduces `U`, any between-subject draw, and censoring
  together. `$events` now always carries `sim_reason`, even with no
  `censoring` argument supplied (`"event"`/`"administrative"`) --
  `sim_status`/`sim_time` are unaffected either way. Dependent/informative
  censoring (a censoring hazard driven by a subject's own simulated
  PK/PD state) is out of scope for this release; see
  `reports/16_censoring_design.md` for why and what it would take.
* **`add_interval_censoring()`: interval censoring** from a visit/
  assessment-time schedule, on top of any already-simulated (and, if
  used, already right-censored) events data frame. An event is only
  known to have happened between the last event-free visit and the
  first visit at which it was detected, giving `(L, R]`; a censored
  subject is only known event-free through their last visit, giving
  `(L, Inf)`. Adds `sim_time_left`/`sim_time_right` columns --
  `sim_time`/`sim_status`/`sim_reason` are unchanged, so every existing
  consumer keeps working. `visit_schedule()` generates a per-subject
  schedule with fixed spacing and optional random jitter, uniform
  (default) or normal (`jitter_dist = `; `jitter` is the half-width for
  uniform, the standard deviation for normal) -- a `jitter` large
  enough relative to the spacing to cross visit order raises a
  `warning()` rather than silently reordering.
  `sim_tte_ode()` also accepts a `visits = ` argument (a schedule, or a
  jitter spec) that applies this automatically, always *after* any
  `censoring = `, inside its own seeded draw. Designed to feed
  `survival::Surv(time = sim_time_left, time2 = sim_time_right, type =
  "interval2")` directly (convert the `Inf` upper bound to `NA` first --
  see `?add_interval_censoring`). An event with `sim_time_left == 0`
  *is* a left-censored observation, with no separate mechanism needed
  (see `?add_interval_censoring` and the vignette's "Left censoring"
  subsection). A visit-process model reacting to a subject's own
  simulated state (`reports/18_interval_censoring_design.md`'s deferred
  option; "C3" below) is still out of scope.
* **`visit_schedule()`'s `jitter_dist = "normal"` is now a *truncated*
  normal**, truncated at `+/- jitter_trunc * jitter` (new argument,
  default `2`), drawn via inverse-CDF. Both `jitter_dist` options are
  now validated against `every` *before* any draw is made (uniform:
  `jitter < every / 2`; normal: `jitter_trunc * jitter < every / 2`),
  which guarantees jittered visits can never cross order or go
  negative -- a `jitter`/`jitter_trunc` combination violating this
  bound is now an error, not a `warning()` (the previous session's
  reordering `warning()` is gone; it is now an unreachable internal
  invariant for either built-in distribution). `"uniform"` output/seeds
  are unchanged from the previous session; `"normal"` output/seeds
  change (the draw mechanism itself changed to truncated inverse-CDF).
* **`thin_visits(visits, ...)`: missed visits and assessment dropout
  ("C1")** on top of a visit schedule (the `visits` argument -- the
  same schedule shape and name `add_interval_censoring()`/
  `sim_tte_ode()` already use) -- reads no simulated outcome, so it can
  never make an assessment schedule informative. Each non-baseline
  visit is independently missed with probability `p_miss`; each subject
  may independently start dropping out from a randomly chosen
  non-baseline visit onward with probability `p_dropout`. Composes with
  `add_censoring()` in either order; must run before
  `add_interval_censoring()`.
* **`visit_schedule_informative(visits, events, ...)`: an outcome-
  reactive visit schedule ("C2")** -- unlike `thin_visits()`, this reads
  `$events` and lets visit attendance/timing react to the simulated
  event time
  (`miss_near_event`: a visit shortly before an event is missed at an
  elevated rate; `extra_visit_after_event`: an unscheduled visit is
  added shortly after an event, using the same distribution spec shape
  as `add_censoring()`'s `censoring`). `miss_near_event$p` and
  `p_miss_base` are **additive/competing risks, not one overriding the
  other**: a near-window visit for an event subject is missed with
  probability `1 - (1 - p_miss_base) * (1 - p)` (collapses exactly to
  `p` when `p_miss_base = 0`) -- pre-release seeds for a call combining
  both changed accordingly. This makes the returned schedule
  informative by construction; a `message()` naming this is emitted on
  every call. See the vignette's "Informative assessment schedules"
  subsection for a worked bias demonstration against `thin_visits()`.
  A third generator reacting to a subject's own PK/PD trajectory ("C3")
  was evaluated and deferred; see
  `reports/20_visit_process_evaluation.md`.

## Improvements

* **New `event_time_method` argument** on `sim_tte_df()` and `sim_tte()`,
  with two values:
  * `"grid"` (the default, unchanged): the event time is the first
    reported trajectory time at which survival falls to or below the
    sampled uniform draw — exactly the existing behavior, byte-for-byte.
  * `"log_survival"` (opt-in): when a crossing occurs strictly after the
    first reported observation, the event time is refined by linear
    interpolation of cumulative hazard (`H = -log(S)`) between the two
    reported points surrounding the crossing, i.e. assuming the hazard
    is constant over that interval. For the M-spline model, whose hazard
    is genuinely piecewise-constant on the reported grid, this recovers
    the event time implied by that discretized hazard exactly. For
    trajectories from a continuously-varying hazard (e.g. Weibull with
    `shape != 1`, or a custom mechanistic model), it is an approximation
    that improves as the reported time grid is refined.
  * Censoring and first-observation crossings are identical between the
    two methods; exactly one `stats::runif(1)` draw is consumed per
    subject either way, so classification (event vs. censored) never
    differs between methods for the same seed and trajectory.
  * `sim_tte()` forwards `event_time_method` explicitly to its internal
    `sim_tte_df()` call; it is unrelated to and never passed to the
    `mrgsolve` simulation step.
* **Time-varying `lp(t)` for `sim_tte()`** (`lp_data` argument): a
  per-subject or population-level log hazard ratio trajectory,
  last-observation-carried-forward between supplied time points.
* **Example PK/PD-driven-hazard models** exported for use with
  `sim_tte_df()`: `simtte_example_model()`/`simtte_example_models()`
  load/list two bundled, user-editable `mrgsolve` `.cpp` files
  demonstrating the pattern by hand (superseded, for new work, by
  `sim_tte_ode()`'s in-solver mechanism and `tte_model()`, but kept as
  a documented, tested way to use `sim_tte_df()` directly).
* **Informative `omega`/`sigma` errors.** Calling `sim_tte_ode(...,
  omega = <matrix>)` (or `sigma =`) on a model with no matching
  declared block previously surfaced `mrgsolve`'s own cryptic
  `"improper signature: omat"`; it now explains that `omat()`/`smat()`
  only update an already-declared block and names what a fix looks
  like.

## Bug fixes / behavior changes versus CRAN 1.0.2

* **Weibull closed-form correctness for `shape < 1`.** CRAN 1.0.2's
  `inst/models/weibull.cpp` approximated the hazard as shape-independent
  (constant-hazard) for solver time below 0.1, which was substantially
  wrong for `shape != 1` -- **measured up to 0.16 absolute error** in
  survival probability for `shape < 1` near `t -> 0` in the tested
  range. The survival probability is now computed as the exact
  closed-form expression `S(t) = exp(-exp(mu + lp) * t^shape)`
  directly, matching the analytical formula to floating-point precision
  at every reported time, for every `shape > 0`.
* **`time` argument now actually controls the output grid.**
  `sim_tte()`'s `time` argument previously only set `end_time` when
  not supplied explicitly; the real output/event-time grid silently
  came from `mrgsolve`'s own default schedule (`delta = 1`) regardless
  of the spacing requested in `time`, for both Weibull and M-spline
  models. `time` now genuinely determines the simulation output grid
  via `mrgsolve`'s `tgrid` mechanism. This is a real, and often
  substantial, change in reported `sim_time` values whenever a caller
  supplied a `time` grid finer than `delta = 1` (the common case, e.g.
  the package's own README examples) -- **measured directly in this
  release's backward-compatibility audit**
  (`inst/validation/13_backcompat_v1_0_2.R`): 16.6% mean relative
  difference in `sim_time` for the README Weibull example (`shape =
  1.1`, itself unaffected by the shape `< 1` fix above), 24.4% for a
  `shape = 2` case, 3.6% for the README M-spline example. `sim_tte()`'s
  `time` default also changed from the scalar `100` to `seq(0, 100, by
  = 1)` (the same effective grid resolution as before when `time` is
  not supplied, consistent with this fix).
* **Protected `...`.** `tgrid`, `obsonly`, `nocb`,
  `carry_out`/`carry.out`, and `data` can no longer be overridden via
  `...` in `sim_tte()`/`explore_pi_tq_surv()`/`.sim_surv_df()`
  (or, new in this release, `sim_tte_ode()`): each was confirmed to
  silently defeat the package's output-grid or trajectory contract if
  supplied by a caller. Doing so now raises a clear error naming the
  argument.
* **M-spline hazard-carry convention fixed and documented.**
  `.sim_surv_df()` now calls `mrgsim(..., nocb = FALSE)`
  (last-observation-carried-forward): the hazard value at `time[i]`
  applies over `[time[i], time[i+1])`. Previously `mrgsolve`'s default
  (`nocb = TRUE`, next-observation-carried-backward) applied a hazard
  value to the *preceding* interval, producing an incorrect survival
  trajectory. See the "M-spline hazard carry convention" section of
  `?sim_tte`.
* **Survival tolerance is now clamped, not merely accepted.** Survival
  values within `1e-8` of `[0, 1]` are clamped to exactly `0`/`1` before
  monotonicity checking and event-time selection (values further outside
  `[0, 1]` remain a hard error). This makes every accepted, normalized
  trajectory safe for a `-log(S)` transform and is provably incapable
  of changing event/censoring classification, given `stats::runif()`'s
  `[0, 1)` support.
* **M-spline input validation hardened.** `basis`/`coefs`/`basehaz` are
  now validated (numeric, finite, matching dimensions, non-negative
  resulting hazard) before reaching `mrgsolve`. Duplicate `time` values
  are accepted only when they carry identical hazard values; conflicting
  duplicates are a hard error.
* **`end_time` boundary semantics resolved.** `end_time = 0` is now
  accepted for both model types (a degenerate zero-duration follow-up;
  every subject is censored at 0). For M-spline models, `end_time` must
  now satisfy `min(time) <= end_time <= max(time)`
  (`end_time < min(time)` was previously accepted and produced a
  censoring time before the supplied hazard trajectory begins).
* **Package hygiene.** Fixed a `.Rbuildignore` regex bug
  (`^\.manuscript$`, which matched nothing) so the `manuscript/`
  development directory is now correctly excluded from the built source
  tarball.
* **Reproducibility semantics documented.** `sim_tte_df()` now documents
  explicitly that a fixed seed reproduces results for a fixed row order
  of the input data, not invariant to reordering subject blocks
  (standard sequential-RNG behavior, unchanged); `sim_tte_ode()`'s own
  reproducibility contract (one `U` draw per `idata` row, insensitivity
  of `mrgsolve::mvgauss()`'s BSV draw to intervening ordinary RNG draws)
  is documented the same way.

## Internal

* Test runbook (`reports/07_test_runbook.md`), targeted test groups and
  a fast/slow split (`dev/run-tests.R`), and a `Makefile` for common
  development tasks -- none shipped in the built package
  (`.Rbuildignore`d).
* Extensive validation scripts under `inst/validation/` (also
  `.Rbuildignore`d, not part of the built package): closed-form
  correctness, grid convergence, PK/PD mechanism checks, in-solver
  boundary-guard/refinement accuracy, M-spline convention equivalence,
  and a CRAN-1.0.2 backward-compatibility audit.
* Added several hundred new tests across the legacy and new APIs (see
  `reports/07_test_runbook.md` for current counts and the fast/slow
  split).
* Pre-Phase-7 whole-package review (`reports/24_pre_phase7_review.md`):
  `sim_tte()`'s own input validation now uses `call. = FALSE`
  consistently with every other exported function; `?sim_tte_ode` is
  reordered into one coherent reading order (mechanism, models, inputs,
  BSV, censoring, reproducibility); `weibull_ode.cpp`/`gompertz_ode.cpp`
  gained the same scaffold marker comments every other library model
  already had, plus a new test asserting the survival scaffold is
  byte-identical across all 12 shipped ODE library models; a shared
  `.fake_events()` test fixture replaces three identical copies;
  `inst/WORDLIST` added for `devtools::spell_check()`; `.gitignore`
  lost a dead pattern and gained explicit exceptions for `NEWS.md`/
  `README.md`/`cran-comments.md`; `.Rbuildignore` dropped 15 entries for
  files that no longer exist (a pre-`reports/`-era filename list and
  unused CI/doc-tool boilerplate). No exported function's behavior
  changed as a result of this review (verified against
  `inst/validation/13_backcompat_v1_0_2.R` and the full test suite).

# simtte 1.0.2

* Fixed time-column position bug in `sim_tte_df()`: time is now resolved by
  column name, not position.
* Preserved original subject IDs in output (no longer renumbered).
* Corrected `.get_tte()` to use the conventional `<= U` inverse-transform
  boundary rule.
* Added optional `time_var` argument to `sim_tte_df()`.
* Changed package Title to "Simulate Bespoke Time-to-Event Models Using
  ODEs" and rewrote the Description to frame **simtte** around bespoke
  ODE-based time-to-event models (built-in Weibull/M-spline hazards plus
  fully custom `mrgsolve` ODE systems), rather than Weibull and spline
  models specifically.

# simtte 1.0.1

* Fixed issues from CRAN review:
  * Added `\value` tag to `pipe.Rd` documenting the return value of the pipe operator.
  * Reduced example run times and unwrapped short examples from `\donttest{}` so they run during automated checks; examples now execute in under 5 seconds.

# simtte 1.0.0

* Initial CRAN release.
* Simulate time-to-event datasets using Weibull and M-spline baseline hazard models.
* Inverse transform sampling from cumulative hazard functions via mrgsolve ODE solver.
* Exported functions: `sim_tte()`, `sim_tte_df()`, `explore_pi_tq_surv()`.
* Includes example M-spline dataset (`ms_data`).
* Two vignettes: "Getting Started" and "Advanced Simulation Scenarios".
