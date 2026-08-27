# simtte (development version)

Pre-Phase-B hardening pass, resolving findings from an audit of the
Phase A work (see `PRE_PHASE_B_REPORT.md` for full detail):

* **Protected `...`.** `tgrid`, `obsonly`, `nocb`,
  `carry_out`/`carry.out`, and `data` can no longer be overridden via
  `...` in `sim_tte()`/`explore_pi_tq_surv()`/`.sim_surv_df()`: each was
  confirmed to silently defeat the package's output-grid or trajectory
  contract if supplied by a caller. Doing so now raises a clear error
  naming the argument.
* **M-spline hazard carry convention fixed and documented.**
  `.sim_surv_df()` now calls `mrgsim(..., nocb = FALSE)`
  (last-observation-carried-forward): the hazard value at `time[i]`
  applies over `[time[i], time[i+1])`. Previously mrgsolve's default
  (`nocb = TRUE`, next-observation-carried-backward) applied a hazard
  value to the *preceding* interval, producing an incorrect survival
  trajectory. See the "M-spline hazard carry convention" section of
  `?sim_tte`.
* **Survival tolerance is now clamped, not merely accepted.** Survival
  values within `1e-8` of `[0, 1]` are clamped to exactly `0`/`1` before
  monotonicity checking and event-time selection (values further outside
  `[0, 1]` remain a hard error). This makes every accepted, normalized
  trajectory safe for a future `-log(S)` transform and is provably
  incapable of changing event/censoring classification, given
  `stats::runif()`'s `[0, 1)` support.
* **M-spline input validation hardened.** `basis`/`coefs`/`basehaz` are
  now validated (numeric, finite, matching dimensions, non-negative
  resulting hazard) before reaching mrgsolve. Duplicate `time` values
  are accepted only when they carry identical hazard values; conflicting
  duplicates are a hard error.
* **Weibull numerical robustness.** The closed-form Weibull calculation
  in `inst/models/weibull.cpp` now special-cases `t = 0` (`S(0) = 1`
  always) and computes in log-cumulative-hazard space for `t > 0`,
  eliminating a `NaN` that previously occurred for large but finite
  `mu + lp` (e.g. `mu = 710`) at `t = 0`. `sim_tte()` also now validates
  `pi`/`mu`/`coefs` are finite and non-empty before calling mrgsolve.
* **`end_time` boundary semantics resolved.** `end_time = 0` is now
  accepted for both model types (a degenerate zero-duration follow-up;
  every subject is censored at 0). For M-spline models, `end_time` must
  now satisfy `min(time) <= end_time <= max(time)`
  (`end_time < min(time)` was previously accepted and produced a
  censoring time before the supplied hazard trajectory begins).
* **Reproducibility semantics documented.** `sim_tte_df()` now documents
  explicitly that a fixed seed reproduces results for a fixed row order
  of the input data, not invariant to reordering subject blocks
  (standard sequential-RNG behavior, unchanged).
* **Package hygiene.** Fixed a `.Rbuildignore` regex bug
  (`^\.manuscript$`, which matched nothing) so the `manuscript/`
  development directory is now correctly excluded from the built source
  tarball.
* Added 139 new tests (258 total) covering all of the above.

Phase A correctness and robustness pass (see the implementation report
for full detail):

* **Weibull model correctness.** Removed the near-zero workaround in
  `inst/models/weibull.cpp` that approximated the hazard as
  shape-independent (constant-hazard) for solver time below 0.1,
  which was substantially wrong for `shape != 1`. The built-in
  Weibull survival probability is now computed as the closed-form
  expression `S(t) = exp(-exp(mu + lp) * t^shape)` directly (rather
  than by integrating a singular ODE for `shape < 1`), and matches
  the analytical formula to floating-point precision at every
  reported time, for every `shape > 0`.
* **`time` argument now controls the output grid.** `sim_tte()`'s
  `time` argument previously only set `end_time` when not supplied
  explicitly; the actual output/event-time grid silently came from
  mrgsolve's own default schedule (`delta = 1`), regardless of the
  spacing requested in `time`. `time` now genuinely determines the
  simulation output grid via mrgsolve's `tgrid` mechanism, for both
  Weibull and M-spline models.
* **`end_time` bug fix.** `sim_tte()`'s `end_time` argument was
  accepted but never forwarded to the internal simulation call, so it
  had no effect whenever it differed from `max(time)`. It is now
  applied correctly and documented.
* **Input and trajectory validation.** `sim_tte()` and `sim_tte_df()`
  now validate `time`/survival trajectories for missing values,
  non-finite values, out-of-range survival probabilities, negative
  times, unsorted or duplicated times within a subject, and
  non-monotone (increasing) survival, with informative errors rather
  than silent repair. See `?sim_tte_df` for the full contract.
* **Documented censoring contract.** `sim_tte_df()` now documents
  precisely that a subject whose trajectory never crosses the sampled
  threshold is censored at that subject's own last reported time
  (not necessarily a shared administrative cutoff), and that a
  trajectory is not required to start at survival = 1 or time = 0.
* Added extensive regression tests for the above, including tests
  that would have failed under the previous Weibull implementation.

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
