# simtte: project status

One page, no phase-report reading required. For the detailed history,
see `reports/` (not shipped with the package -- a local, gitignored
directory of design/implementation reports, one or more per
development phase).

## What the package does

**simtte** simulates time-to-event (survival) data for clinical trial
design and analysis, two ways:

1. **Joint PK/PD and time-to-event simulation** (`sim_tte_ode()`): the
   hazard is driven by a pharmacokinetic/pharmacodynamic model. The
   survival probability is added to the model as one more ODE
   compartment, and each subject's event time is found *during* the
   ODE solve (not resolved afterward from a saved trajectory). Six
   built-in PK/PD hazard models ship out of the box (concentration-
   driven, four indirect-response variants, and a target-mediated-
   disposition model), plus a converter, `tte_model()`, that turns any
   user-supplied `mrgsolve` model into one of these.
2. **Bespoke parametric/flexible hazard simulation** (`sim_tte()`/
   `sim_tte_df()`, the original CRAN 1.0.2 API): closed-form Weibull
   and flexible M-spline baseline hazards, plus a fully model-agnostic
   entry point for any already-simulated survival trajectory from any
   source. Unchanged in this rework except for a small number of
   real, documented bug fixes (see `NEWS.md`).

Either way, `add_censoring()` can layer independent right censoring
(exponential, Weibull, uniform, lognormal, gamma, or a user-supplied
distribution) on top of the simulated event times; `censoring_rate_for()`
picks a distribution parameter for a target censoring fraction. On top
of that, `add_interval_censoring()` can map the (possibly
right-censored) outcome onto a visit/assessment-time interval `(L, R]`
-- `sim_time_left`/`sim_time_right` columns, ready for
`survival::Surv(..., type = "interval2")` (an event with
`sim_time_left == 0` is a left-censored observation, with no separate
mechanism needed); `visit_schedule()` generates a schedule with fixed
spacing and jitter, `thin_visits()` layers missed visits/dropout on top
(non-informative), and `visit_schedule_informative()` layers a visit
process that reacts to the simulated event time (informative by
construction, documented loudly). A visit process depending on a
subject's own simulated PK/PD state, and dependent/informative right
censoring, are the remaining out-of-scope pieces -- see
`reports/16_censoring_design.md`/`reports/18_interval_censoring_design.md`/
`reports/20_visit_process_evaluation.md`.

## Public API

| Function | What it's for |
|---|---|
| `sim_tte_ode()` | In-solver PK/PD + time-to-event simulation |
| `tte_model()` | Convert your own `mrgsolve` model for `sim_tte_ode()` |
| `add_censoring()` | Apply independent right censoring to any simulated events data frame |
| `censoring_rate_for()` | Solve for a censoring-distribution parameter hitting a target censoring fraction |
| `add_interval_censoring()` | Map a (possibly right-censored) outcome onto a visit-schedule interval |
| `visit_schedule()` | Generate a per-subject visit schedule (fixed spacing + jitter) |
| `thin_visits()` | Missed visits and assessment dropout on a schedule (non-informative) |
| `visit_schedule_informative()` | A visit schedule that reacts to the simulated event time (informative) |
| `sim_tte()` | Closed-form Weibull/M-spline simulation (original API) |
| `sim_tte_df()` | Inverse-transform sampling on any custom trajectory |
| `explore_pi_tq_surv()` | Survival-difference-at-a-quantile utility |
| `simtte_example_model()`/`simtte_example_models()` | Bundled example PK/PD-hazard models for `sim_tte_df()` |

Three vignettes: `vignette("pkpd-time-to-event")` (the new capability,
start here), `vignette("introduction")`, `vignette("advanced-usage")`
(the original API).

## How to run the tests

```
Rscript dev/run-tests.R              # fast: ~40s
Rscript dev/run-tests.R all --slow   # full: ~2m35s
Rscript dev/run-tests.R --check --fast  # R CMD check --as-cran, slow tests off: ~1m30s
Rscript dev/run-tests.R --check      # R CMD check --as-cran, slow tests on (pre-release): ~4min
```

See `reports/07_test_runbook.md` for targeted groups (one file/model at
a time), what each test file covers, and a "how do I check X" table.
Current counts: 1126 pass / 0 fail / 49 skip (fast), 1267 pass / 0 fail /
0 skip (slow). `R CMD check`: `Status: OK` (0 errors, 0 warnings, 0
notes).

## Where the validation evidence lives

`inst/validation/` (not shipped with the built package): one script
per major claim, each run directly and read for its printed output --
not a pass/fail suite, the citable numbers behind specific report
claims (closed-form correctness, grid convergence, PK/PD mechanism
checks, in-solver boundary-guard/refinement accuracy, M-spline
convention equivalence against `splines2`, and a full backward-
compatibility audit against the CRAN 1.0.2 tag). `reports/` has the
narrative for each: what was tested, why, and what the numbers mean.

## Known limitations

- **Weibull `shape < 0.05`** (in-solver ODE form): accuracy degrades
  near `t = 0`; a `message()` fires and points at `sim_tte(type =
  "weibull")` (exact there) as the alternative. Never an error.
- **A very steep hazard at a coarse `delta`** can trigger a safety-net
  fallback for some or all subjects; the fallback still runs and is
  disclosed via a `message()` naming the remedy (a finer `delta`, or
  tighter `rtol`/`atol`, depending on which of two known causes
  applies). Not silent, not wrong, just less precise until the remedy
  is applied.
- **`model = "mspline"` ships only 3 fixed interior-knot counts** (3,
  5, 7). `sim_tte(type = "ms")` accepts an arbitrary count, at the
  cost of a piecewise-constant (not continuous) hazard.
- **`tte_model()` compiles fresh each R session** by default (no
  persistent build cache); `options(mrgsolve.project = <a directory>)`
  makes it persist across sessions, at the cost of two disclosed
  caveats -- see `vignette("pkpd-time-to-event")` section 5.
- **An extreme, far-past-default parameter/dose combination can make
  the ODE solver itself fail** (not just trigger the safety-net
  fallback above) for the stiffer PK/PD models (`tmdd_hazard`
  specifically, documented). Not hit under any shipped default; no
  informative wrapper exists yet for this case, so the raw `mrgsolve`/
  `lsoda` error is what you would see.
- **Between-subject variability (`omega`) on `tmdd_hazard`** is
  verified thoroughly on one parameter (`V2`); its other ten
  registered targets are not individually stress-tested the same way.
- **Right censoring is independent only.** `add_censoring()` draws a
  censoring time with no dependence on a subject's own simulated
  PK/PD trajectory; informative/dependent censoring (e.g. dropout
  driven by toxicity) needs a model-based mechanism not built yet --
  see `reports/16_censoring_design.md` options B/C for the design and
  cost estimate.
- **`censoring_rate_for()`'s solved rate/scale is approximate**: it
  depends on the particular uncensored run supplied to it, not a
  closed-form property of the model. Re-check the realized fraction on
  a larger cohort after drawing with the solved parameter.
- **A visit process depending on a subject's own simulated PK/PD state
  ("C3") is not built.** `thin_visits()` (schedule-only) and
  `visit_schedule_informative()` (reacts to the simulated event time)
  cover the first two complexity levels evaluated; C3 needs
  `keep_trajectory = TRUE` plus interpolation between reported times --
  see `reports/20_visit_process_evaluation.md`.
- **Delayed entry / left truncation is not built.** A subject observed
  only from `entry > 0` onward, with those already past the event
  before `entry` excluded from the risk set, needs a different analysis
  dataset shape (`survival::Surv(start, stop, event)`) that nothing in
  this package currently produces -- see
  `reports/21_left_censoring_evaluation.md` for the design sketch
  (`add_delayed_entry(events, entry)`) and the open `end`-vs-`entry`
  alignment question a future implementation would need to resolve.

## Not started / open before a CRAN release

- See `reports/23_phase6c_report.md` "Questions for the author" for the
  specific open decisions and the release readiness assessment.
