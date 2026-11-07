# Changelog

## simtte (development version)

- The blog post [Driving a hazard from a PK/PD ODE system in
  simtte](https://csetraynor.github.io/simtte/articles/pkpd-hazard-ode-simulation.html)
  is now an article on the package website.
- [`simtte_model_cache_clear()`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md)
  warns, and leaves the path out of its return value, when it cannot
  fully remove a directory. On Windows this happens for a model loaded
  in the current session, whose compiled file stays open until R exits;
  see “Windows” in
  [`?simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md).

## simtte 2.0.0

A major release. It adds in-solver time-to-event simulation, censoring
and visit schedules. It also corrects
[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
behaviours that change its results or checks compared with CRAN 1.0.2;
see “Breaking changes”.

### Breaking changes

- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  uses `time` as the output grid. Event times (`sim_time`) change
  whenever `time` is finer than 1, for example by a mean of 16.6% for
  the README Weibull example. Pass `time = seq(0, end, by = 1)` to
  reproduce the 1.0.2 grid. The default of `time` is now
  `seq(0, 100, by = 1)`, the same effective resolution as before.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  computes the Weibull survival probability in closed form,
  `S(t) = exp(-exp(mu + lp) * t^shape)`. Results for `shape != 1` change
  at small times; 1.0.2 was off by up to 0.16 for `shape < 1`.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  applies each M-spline hazard value over the interval that follows its
  time point (`nocb = FALSE`), so M-spline survival trajectories and
  event times change. See
  [`?sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
  “M-spline hazard carry convention”.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
  [`explore_pi_tq_surv()`](https://csetraynor.github.io/simtte/reference/explore_pi_tq_surv.md)
  and
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  raise an error naming the argument if `...` supplies `tgrid`,
  `obsonly`, `nocb`, `carry_out` or `data`, which would break the
  output-grid contract. Remove the argument.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  rejects an M-spline `end_time` outside
  `min(time) <= end_time <= max(time)`, and rejects non-finite or
  mismatched `basis`, `coefs` or `basehaz`, a negative resulting hazard,
  and duplicate `time` values with conflicting hazards.

### New features

- [`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
  applies independent right censoring to any events data frame from
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
  [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  or
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md).
  It draws censoring times from an exponential, Weibull, uniform,
  lognormal, gamma or user-supplied distribution, updates `sim_time` and
  `sim_status`, and records the outcome in a new `sim_reason` column
  (`"event"`, `"censored"` or `"administrative"`).
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  also takes `censoring`, so one `seed` reproduces the whole simulation.
  Informative (dependent) censoring is not available.
- [`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
  maps events onto a visit schedule, adding `sim_time_left` and
  `sim_time_right` columns for `(L, R]` intervals and `(L, Inf)` for
  subjects censored event-free. The columns feed
  `survival::Surv(type = "interval2")` after converting `Inf` to `NA`;
  an event with `sim_time_left == 0` is left-censored.
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  also takes `visits`, applied after any `censoring`.
- [`censoring_rate_for()`](https://csetraynor.github.io/simtte/reference/censoring_rate_for.md)
  solves for the censoring-distribution parameter that gives a target
  censoring fraction from the event times of an uncensored run.
- `draw_censoring_times(spec, n, seed = NULL)` draws directly from an
  [`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
  distribution spec, with no events data frame.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  and
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
  gain `event_time_method`. `"grid"` (the default) keeps the first
  reported time at which survival reaches the uniform draw.
  `"log_survival"` interpolates the cumulative hazard between the two
  surrounding reported times: exact for the M-spline model and an
  approximation that improves with a finer grid for continuously varying
  hazards. Event versus censored classification is identical for the
  same seed.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  takes a time-varying `lp_data`: a per-subject or population-level log
  hazard ratio trajectory, carried forward between supplied time points.
- [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  simulates joint PK/PD and time-to-event data with in-solver event
  detection. Survival is one more ODE compartment
  (`dxdt_p11 = -p11 * HAZ`), and each event time is located during
  integration, where `p11` first falls to or below a per-subject uniform
  draw. Built-in models are `"exponential"`, `"weibull"`, `"gompertz"`
  and `"mspline"` (quadratic M-spline hazard with 3, 5 or 7 interior
  knots).
- [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  includes six PK/PD-linked hazard models built on unedited
  [`mrgsolve::modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
  backbones: `"pk_hazard"` (concentration-driven), `"irm1_hazard"` to
  `"irm4_hazard"` (response-driven) and `"tmdd_hazard"` (target-mediated
  drug disposition). Dosing is ordinary `mrgsolve` event data passed as
  `data`. A `data` frame whose unique `ID`s do not match the resolved
  subject count errors, naming both numbers.
- [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  takes `covariates` and `beta` (optionally `formula`) for a log-linear
  linear predictor, `lp(t) = sum_k beta_k * X_k(t)`. `formula` (for
  example `~ age + arm`) builds the design matrix with
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  so factors, interactions and terms such as `I(age^2)` or
  `splines::ns(weight, 3)` work. `covariates` can be baseline (one row
  per subject) or time-varying. An implied intercept is dropped with a
  message, as in
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). A
  baseline frame with one row and no `ID` applies to every subject; with
  several rows it sets the subject count.
- [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  takes `omega` for log-normal between-subject variability on a model’s
  structural PK/PD parameters, applied through per-subject `idata`, so
  no `$OMEGA` block is needed. A model that declares its own `$OMEGA`
  uses
  [`mrgsolve::omat()`](https://mrgsolve.org/docs/reference/omega.html)
  instead; `omega` or `sigma` on a model with no matching block explains
  this rather than surfacing the `mrgsolve` error.
- [`sim_tte_ode_models()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
  lists the built-in model names accepted by `sim_tte_ode(model = )`
  with each model’s valid `omega` targets.
- [`simtte_example_model()`](https://csetraynor.github.io/simtte/reference/simtte_example_model.md)
  and
  [`simtte_example_models()`](https://csetraynor.github.io/simtte/reference/simtte_example_models.md)
  load and list two bundled, editable `mrgsolve` models for use with
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md).
- [`simtte_model_cache()`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
  reports where compiled models are stored (under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by default). Set
  `options(simtte.cache_dir = tools::R_user_dir("simtte", "cache"))` to
  keep them across sessions.
- `simtte_model_cache_clear(all = FALSE)` removes compiled models for
  the current R, simtte and `mrgsolve` versions, or for every version
  with `all = TRUE`. It only removes directories simtte created.
- `simtte_prepare_model(model, ...)` builds (or reuses) and loads a
  model once so
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  never rebuilds it. Call it once per
  [`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  worker, then pass the same `model` to later
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  calls. It rebuilds a
  [`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
  result from its stored code in a fresh process, so a saved object
  works in a separate R process. For `"mspline"` it takes `knots`;
  `coefs` and `boundary_knots` error. See
  [`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
  “Parallel simulation”.
- `thin_visits(visits, ...)` drops visits with probability `p_miss` and
  starts dropout with probability `p_dropout`. It reads no outcome, so
  it never makes a schedule informative. Run it before
  [`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md).
- [`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
  converts your own `mrgsolve` PK/PD model for
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
  from code, a file path or a compiled model. It adds the survival
  compartment without editing your `$ODE`/`$DES` and reproduces the six
  built-in hazard models from their
  [`modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
  backbones.
- [`visit_schedule()`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)
  builds a per-subject schedule with fixed spacing and jitter.
  `jitter_dist = "normal"` is truncated at `+/- jitter_trunc * jitter`
  (default `jitter_trunc = 2`), and a jitter large enough to reorder
  visits is an error.
- `visit_schedule_informative(visits, events, ...)` makes attendance
  react to the simulated event time. `miss_near_event` raises the chance
  of missing a visit shortly before an event, and
  `extra_visit_after_event` adds a visit after it. `miss_near_event$p`
  and `p_miss_base` combine as competing risks,
  `1 - (1 - p_miss_base) * (1 - p)`. The schedule is informative by
  construction, and every call says so with a
  [`message()`](https://rdrr.io/r/base/message.html).

### Minor improvements and fixes

- Help pages link related functions and say why they go together,
  [`?simtte`](https://csetraynor.github.io/simtte/reference/simtte-package.md)
  gains a “Which functions go together” section, and the `model`,
  `events` and `name` arguments say where their values come from.
- [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  accepts `end_time = 0`; every subject is censored at 0.
- [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
  clamps survival values within `1e-8` of `[0, 1]` to 0 or 1 before
  checking monotonicity. Values further outside still error.
- [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
  and
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  document that a fixed `seed` reproduces results for a fixed row order
  of the input.

## simtte 1.0.2

CRAN release: 2026-07-31

- Fixed time-column position bug in
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md):
  time is now resolved by column name, not position.
- Preserved original subject IDs in output (no longer renumbered).
- Corrected `.get_tte()` to use the conventional `<= U`
  inverse-transform boundary rule.
- Added optional `time_var` argument to
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md).
- Changed package Title to “Simulate Bespoke Time-to-Event Models Using
  ODEs” and rewrote the Description to frame **simtte** around bespoke
  ODE-based time-to-event models (built-in Weibull/M-spline hazards plus
  fully custom `mrgsolve` ODE systems), rather than Weibull and spline
  models specifically.

## simtte 1.0.1

CRAN release: 2026-07-11

- Fixed issues from CRAN review:
  - Added `\value` tag to `pipe.Rd` documenting the return value of the
    pipe operator.
  - Reduced example run times and unwrapped short examples from
    `\donttest{}` so they run during automated checks; examples now
    execute in under 5 seconds.

## simtte 1.0.0

- Initial CRAN release.
- Simulate time-to-event datasets using Weibull and M-spline baseline
  hazard models.
- Inverse transform sampling from cumulative hazard functions via
  mrgsolve ODE solver.
- Exported functions:
  [`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
  [`explore_pi_tq_surv()`](https://csetraynor.github.io/simtte/reference/explore_pi_tq_surv.md).
- Includes example M-spline dataset (`ms_data`).
- Two vignettes: “Getting Started” and “Advanced Simulation Scenarios”.
