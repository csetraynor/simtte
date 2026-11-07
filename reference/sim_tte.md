# Simulate Time-to-Event Data

Main function of simtte. Simulates a time-to-event dataset (survival
dataset) using either a Weibull parametric model or an M-spline flexible
baseline hazard model. Event times are generated via inverse transform
sampling from the cumulative hazard function computed by mrgsolve.

## Usage

``` r
sim_tte(
  pi,
  log_pi = TRUE,
  mu = -3,
  coefs = 0,
  basis = NULL,
  time = seq(0, 100, by = 1),
  end_time,
  type = "weibull",
  event_time_method = c("grid", "log_survival"),
  lp_data = NULL,
  ...
)
```

## Arguments

- pi:

  Numeric matrix or vector. Prognostic index (linear predictor) for each
  individual. Given a covariate matrix X and coefficient vector b,
  `pi = X %*% b`. Must be non-empty and finite (no
  `NA`/`NaN`/`Inf`/`-Inf`).

- log_pi:

  Logical; is the prognostic index already on the log scale? Default
  `TRUE`.

- mu:

  Numeric scalar. Intercept parameter of the model. Must be finite.
  Default `-3`.

- coefs:

  Numeric vector. For M-spline models, the coefficients of each spline
  basis function. For Weibull models, the shape parameter (scalar),
  which must be strictly positive. Must be finite in all cases.

- basis:

  Numeric matrix. Basis matrix for M-spline models, with one row per
  element of `time` and one column per coefficient in `coefs`. Must be
  finite. Ignored for Weibull models. Default `NULL`.

- time:

  Numeric vector. The output time grid: the exact set of times at which
  the survival trajectory is simulated and reported, and therefore the
  set of times at which an event can be detected (see "Time grid and
  event-time resolution" below). For M-spline models, `time` must
  additionally be sorted in non-decreasing order and have one element
  per row of `basis` (the baseline hazard is only defined at these
  points). Must be non-negative, finite, and contain no `NA`/`NaN`.
  Default `seq(0, 100, by = 1)`.

- end_time:

  Numeric scalar. Administrative censoring horizon. Defaults to
  `max(time)`. Any requested `time` points beyond `end_time` are
  dropped, and `end_time` itself is always added to the output grid so
  that censoring lands exactly on it. Must be non-negative;
  `end_time = 0` is accepted (a degenerate zero-duration follow-up: the
  output grid is just time 0, and since \\S(0) = 1\\ for every supported
  model, every subject is censored at 0). For Weibull models `end_time`
  may otherwise be any non-negative value, including below `min(time)`
  or above `max(time)` (the closed-form hazard is defined for all
  `t >= 0`). For M-spline models `end_time` must satisfy
  `min(time) <= end_time <= max(time)`: the baseline hazard is only
  defined on the supplied `time` grid and this package does not
  extrapolate it in either direction.

- type:

  Character string. Model type: `"weibull"` or `"ms"` (M-splines).
  Default `"weibull"`.

- event_time_method:

  Character string. `"grid"` (default) or `"log_survival"`. See
  "Event-time interpolation" below. Forwarded explicitly to
  [`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md);
  it is unrelated to the mrgsolve simulation step and is never passed
  via `...`.

- lp_data:

  `NULL` (default) or a data frame with columns `time`, `lp`, and
  optionally `ID`, giving a time-varying linear predictor `lp(t)` for
  the built-in Weibull/M-spline models. See "Time-varying `lp(t)`"
  below. `NULL` preserves every prior behavior of `sim_tte()`
  byte-for-byte.

- ...:

  Additional arguments passed to
  [`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html). `tgrid`,
  `obsonly`, `nocb`, `carry_out`/`carry.out`, and `data` are controlled
  internally to guarantee the contracts documented here and cannot be
  overridden; supplying them raises an error naming the argument.

## Value

A data frame with columns:

- sim_time:

  Simulated event or censoring time.

- sim_status:

  Event indicator (1 = event, 0 = censored).

- ID:

  Subject identifier.

- lp:

  Log prognostic index (linear predictor).

## Time grid and event-time resolution

The output grid used for simulation is `sort(unique(time))`, restricted
to values `<= end_time`, with `end_time` appended if it is not already
present. mrgsolve is instructed to report the survival trajectory *only*
at these times (via its `tgrid` mechanism), so the grid is deterministic
and is not silently replaced by mrgsolve's own default output schedule.

Event times are then resolved by
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
as the first grid point at which the survival probability falls to or
below a uniform random draw. For `event_time_method = "grid"` (the
default), no interpolation between grid points is performed: the
resolution of the simulated event time is exactly the spacing of `time`.
For `event_time_method = "log_survival"`, a crossing that falls strictly
between two grid points is refined by interpolation; see "Event-time
interpolation" below and the corresponding section of
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
for the full mathematical contract (this is the function `sim_tte()`
delegates to internally, forwarding `event_time_method` explicitly). A
subject whose survival never crosses the threshold on the grid is
censored at `end_time` under either method.

## Event-time interpolation

See the "Event-time interpolation" section of
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
for the full mathematical contract and edge-case rules; `sim_tte()`
simply forwards `event_time_method` to it. In summary: `"log_survival"`
assumes the hazard is constant between consecutive reported trajectory
points and linearly interpolates cumulative hazard accordingly. For the
M-spline model this recovers the event time implied by the discretized
(piecewise-constant, via the "M-spline hazard carry convention" below)
hazard exactly. For the Weibull model with `shape != 1`, whose true
hazard is not constant within an interval, it is an approximation whose
accuracy improves as `time` is refined (exact when `shape == 1`, the one
case where the Weibull hazard genuinely is constant).

## M-spline hazard carry convention

For `type = "ms"`, `basis %*% coefs` gives a baseline hazard value at
each element of `time`; this is supplied to mrgsolve as a time-varying
input, not evaluated as a continuous spline function at arbitrary times.
The package fixes the piecewise interpretation explicitly (mrgsolve's
`nocb = FALSE`, i.e. last-observation-carried-forward): the hazard value
reported at `time[i]` applies from `time[i]` (inclusive) until
`time[i + 1]` (exclusive) – a right-continuous step function in the
hazard. Concretely, for knots/hazards `time = c(0, 1, 2)`,
`basehaz = c(1, 2, 4)`: \\S(0.5) = \exp(-0.5)\\, \\S(1) = \exp(-1)\\,
\\S(1.5) = \exp(-2)\\, \\S(2) = \exp(-3)\\ (verified directly against
mrgsolve output; this is **not** the mrgsolve default, which is
next-observation-carried-backward and would instead give \\S(0.5) =
\exp(-1)\\).

A consequence of this convention: the hazard value associated with the
*last* element of `time` has no interval after it to apply to, so it
never affects the survival trajectory (it would only matter if
`end_time` extended past `max(time)`, which is not permitted; see
`end_time` above). `basis`/`coefs` therefore define `length(time) - 1`
active piecewise-constant hazard segments, not `length(time)`.

This is a grid-based piecewise-constant hazard, not continuous M-spline
evaluation: the package receives precomputed basis values at specific
times and never re-evaluates the spline basis at arbitrary times.

## Time-varying `lp(t)`

Supplying `lp_data` allows the linear predictor itself to vary over
time, for both built-in models. `lp_data` must contain `time` (finite,
non-negative) and `lp` (finite; unlike a hazard, negative values are
valid) columns, and, optionally, an `ID` column:

- no `ID`: a single population-level \\lp(t)\\ trajectory, applied
  identically to every subject;

- `ID` present: subject-specific trajectories.
  `sort(unique(lp_data$ID))` must equal `1:length(pi)` exactly (one full
  trajectory per subject; no partial coverage).

In both cases, `pi` keeps its existing role as each subject's baseline
linear-predictor offset: \\lp\_{\mathrm{effective}}(t) = pi_i + lp(t)\\.
`lp_data` rows are not sorted or deduplicated automatically (same
philosophy as
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)'s
trajectory contract): unsorted times are an error; duplicate times are
accepted only with an identical `lp` value, otherwise an error. `lp(t)`
is piecewise constant under last-observation-carried-forward (the same
convention as `basehaz`, above): the value at `lp_data$time[i]` applies
from that time until the next `lp_data` knot for that subject.
`lp_data`'s own time grid is independent of, and need not coincide with,
the main `time` argument.

Coverage requirement: every subject's `lp_data` must include an
observation at `time = 0` (there is no implicit \\lp(0)\\). For M-spline
models, which cannot extrapolate, it must additionally reach at least
`end_time`; for Weibull models this is not required, since the last
known `lp` value is carried forward for any remaining follow-up.

For `type = "weibull"`, a separate, distinct model (`weibull_tv`) is
used internally whenever `lp_data` is supplied; the constant-`lp`
closed-form Weibull model (`inst/models/weibull.cpp`) used whenever
`lp_data` is `NULL` is not modified or replaced. `weibull_tv` computes
\\S(t)\\ by an **exact** closed-form cumulative-hazard summation over
the piecewise-constant \\lp(t)\\ segments – it does not use ODE
integration or internal-solver event detection (the latter was
investigated and rejected; see `PHASE_E_THRESHOLD_TRACKING_REPORT.md`) –
so it retains the same solver-free numerical robustness established for
the baseline model. For `type = "ms"`, the existing `basehaz(t)`
mechanism and `lp(t)` combine as \\h(t) = basehaz(t) \exp\[\mu +
lp(t)\]\\, both under the same last-observation-carried-forward
convention.

## References

Bender R, Augustin T, Blettner M (2005). Generating survival times to
simulate Cox proportional hazards models. *Statistics in Medicine*,
24(11), 1713–1723.
[doi:10.1002/sim.2059](https://doi.org/10.1002/sim.2059)

Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
and proportional-odds models for censored survival data, with
application to prognostic modelling and estimation of treatment effects.
*Statistics in Medicine*, 21(15), 2175–2197.
[doi:10.1002/sim.1203](https://doi.org/10.1002/sim.1203)

## See also

[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
which turns any `(ID, time, p11)` trajectory, including one from your
own mrgsolve model, into event times.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
for joint PK/PD and time-to-event simulation with in-solver event
detection (also the route for Cox-model data; see its "Simulating from a
Cox model" section).
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
to add independent right censoring to this function's output (call it
directly on the returned data frame; there is no argument here for it),
then
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
to map the result onto an assessment-time interval.
[`explore_pi_tq_surv`](https://csetraynor.github.io/simtte/reference/explore_pi_tq_surv.md),
to see how the prognostic index moves survival before simulating.

## Examples

``` r
# Fast Weibull example with a small dataset
set.seed(1)
lp <- matrix(rnorm(5, 0, 0.5), nrow = 5)
result <- sim_tte(
  pi = lp, mu = -1, coefs = 1.1,
  time = seq(0.1, 10, by = 0.5),
  type = "weibull", end_time = 10
)
#> Loading model from cache.
head(result)
#> # A tibble: 5 × 4
#>   sim_time sim_status    ID      lp
#>      <dbl>      <dbl> <dbl>   <dbl>
#> 1      5.1          1     1 -0.313 
#> 2      4.1          1     2  0.0918
#> 3      1.6          1     3 -0.418 
#> 4      1.6          1     4  0.798 
#> 5      1.1          1     5  0.165 
# \donttest{
# Larger examples using bundled ms_data
data("ms_data")
mu <- ms_data$mu
basis <- ms_data$basis
coefs <- ms_data$coefs
time <- ms_data$time
lp <- matrix(runif(nrow(basis)), nrow = nrow(basis))
wei_sim <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
  time = time, type = "weibull", end_time = 100)
#> Loading model from cache.
ms_sim <- sim_tte(pi = lp, mu = mu, basis = basis,
  coefs = coefs, time = time, type = "ms")
#> Building ms ... 
#> done.
# }
```
