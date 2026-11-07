# Simulate Time-to-Event Data from a Custom mrgsolve Output

Applies inverse transform sampling to a data frame produced by an
mrgsolve simulation that contains a survival probability column. This
function is useful when you have a custom time-to-event model
implemented in mrgsolve and want to generate event times.

## Usage

``` r
sim_tte_df(
  dat,
  surv_var = "p11",
  id_var = "ID",
  time_var = "time",
  xdata = NULL,
  event_time_method = c("grid", "log_survival")
)
```

## Arguments

- dat:

  Data frame. Output from
  [`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html) containing
  at minimum columns for time, subject ID, and survival probability.
  Columns are resolved by name (not position).

- surv_var:

  Character string. Name of the column containing the survival
  probability (probability of remaining event-free). Default `"p11"`.

- id_var:

  Character string. Name of the subject ID column. Default `"ID"`.

- time_var:

  Character string. Name of the time column. Default `"time"`.

- xdata:

  Optional data frame of additional covariates to merge into the output.
  Must contain a column matching `id_var`.

- event_time_method:

  Character string. `"grid"` (default) or `"log_survival"`. See
  "Event-time interpolation" below. Partial matching is allowed
  (resolved via [`match.arg`](https://rdrr.io/r/base/match.arg.html)).

## Value

A data frame with columns:

- sim_time:

  Simulated event or censoring time.

- sim_status:

  Event indicator (1 = event, 0 = censored).

- ID:

  Subject identifier.

Plus any additional columns from `xdata`.

## Trajectory contract

`dat` must contain, for each subject, one ordered survival trajectory
satisfying:

- at least one row;

- `time` values that are non-negative, finite, and strictly increasing
  within the subject, with no duplicates. Rows are **not** sorted
  automatically: unsorted or duplicated times raise an error rather than
  being silently re-ordered or resolved by row order;

- survival probabilities in \\\[0, 1\]\\ (a small numerical tolerance of
  \\10^{-8}\\ is allowed for floating-point/ODE solver noise: values
  outside \\\[0, 1\]\\ but within this tolerance are accepted and then
  deterministically **clamped** to exactly \\0\\ or \\1\\ before
  monotonicity checking and event-time selection; values further outside
  \\\[0, 1\]\\ are a hard error, never silently clamped or dropped).
  This clamping cannot change whether a subject is classified as an
  event or censored, because
  [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html) draws from
  \\\[0, 1)\\: a raw value already below 0 and its clamped value of
  exactly 0 both satisfy \\value \le U\\ for every achievable \\U\\, and
  symmetrically a raw value already above 1 and its clamped value of
  exactly 1 both fail it;

- survival that is non-increasing within the subject (checked on the
  clamped values). A trajectory that increases by more than the same
  \\10^{-8}\\ tolerance at any step is rejected rather than silently
  repaired, because that pattern cannot arise from a valid survival
  function and usually indicates a data or model error; an increase
  within the tolerance is treated as numerical noise and accepted as-is
  (not further adjusted).

A trajectory is **not** required to start at survival = 1 or at time =
0: `sim_tte_df()` never invents an implicit \\S(0) = 1\\ state. Whatever
survival value is reported at the first observed time is used exactly as
given (see "Censoring and event-time resolution" below for what this
implies when that first value is already at or below the sampled
threshold). Subjects may have different time ranges; see below for what
this means for censoring.

## Censoring and event-time resolution

For each subject, a uniform random variate \\U \sim \mathrm{Uniform}(0,
1)\\ is drawn, and the **first reported time** at which survival is at
or below \\U\\ is located. If the first reported observation already
satisfies \\S(t_1) \le U\\, that first time is returned as the event
time; no earlier crossing is inferred, and this is identical for both
values of `event_time_method` (there is no earlier reported point to
interpolate from). What happens next depends on `event_time_method`: see
"Event-time interpolation" below. Exactly one `stats::runif(1)` is drawn
per subject regardless of `event_time_method`.

If a subject's trajectory never crosses \\U\\, that subject is censored
(`sim_status = 0`) at the **last reported time in their own
trajectory**, i.e. `max(time)` for that subject, for both values of
`event_time_method`. Because this censoring time is subject-specific,
`sim_tte_df()` does **not** inherently represent a common administrative
censoring cutoff unless every subject's trajectory happens to share the
same final time — each subject's final observation time is effectively
that subject's own censoring horizon. There is currently no separate
`end_time` argument to `sim_tte_df()` itself;
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
provides one for its built-in Weibull/M-spline models by ensuring every
subject's generated trajectory is extended to a common `end_time` (see
[`?sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)).

## Event-time interpolation

`event_time_method = "grid"` (the default) resolves event times purely
by grid-based lookup, not interpolation: the returned time is always one
of the times present in `dat` for that subject, never a value between
two reported times.

`event_time_method = "log_survival"` additionally refines a crossing
that occurs strictly after the first reported observation: it assumes
the hazard is constant between the two consecutive reported trajectory
points surrounding the crossing, and linearly interpolates cumulative
hazard (\\H = -\log(S)\\) between them accordingly. For trajectories
whose hazard is genuinely piecewise-constant on the reported grid (e.g.
this package's M-spline model; see
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)'s
"M-spline hazard carry convention"), this recovers the event time
implied by that discretized hazard exactly. For trajectories generated
from a continuously-varying hazard (e.g. the Weibull model with
`shape != 1`, or a custom mechanistic model), it is an approximation
whose accuracy improves as the reported time grid is refined.
`event_time_method` never changes which subjects are classified as
events vs. censored, and never changes the returned time for censoring
or a first-observation crossing — it can only change `sim_time` for a
crossing that falls strictly between two reported points, and the
interpolated time is always within that interval (never before the
earlier point, never after the subject's final reported time).

## Reproducibility

`sim_tte_df()` draws exactly one `stats::runif(1)` per subject, in the
order subjects are first encountered in `dat` (row order), not in
ID-sorted order. With a fixed
[`set.seed()`](https://rdrr.io/r/base/Random.html), results are
therefore reproducible **for a fixed row order** of `dat`, but
reordering the rows of `dat` so that subjects are first encountered in a
different sequence changes which random draw each subject receives, even
with the same seed and the same set of per-subject trajectories. This is
standard, unmodified R RNG behavior (a single shared random stream
consumed sequentially), not a bug; `sim_tte_df()` does not currently use
per-subject or per-ID RNG streams.

## See also

[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
which builds the Weibull or M-spline trajectory this function consumes.
[`simtte_example_model`](https://csetraynor.github.io/simtte/reference/simtte_example_model.md)
loads a bundled mrgsolve model whose simulated `p11` trajectory can be
passed as `dat`.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
for in-solver event detection with no separate trajectory step.
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
and
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
to post-process this function's output (call
[`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
first if using both).

## Examples

``` r
# Create a mock survival probability data frame (no mrgsolve required)
mock_dat <- data.frame(
  ID = rep(1:3, each = 50),
  time = rep(seq(0.1, 10, length.out = 50), 3),
  p11 = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 3)
)
result <- sim_tte_df(mock_dat)
head(result)
#> # A tibble: 3 × 3
#>   sim_time sim_status    ID
#>      <dbl>      <dbl> <int>
#> 1     1.31          1     1
#> 2     9.19          1     2
#> 3    10             0     3

# Custom column names
mock_dat2 <- data.frame(
  subj = rep(101:103, each = 50),
  tgrid = rep(seq(0.1, 10, length.out = 50), 3),
  surv = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 3)
)
result2 <- sim_tte_df(mock_dat2, surv_var = "surv",
  id_var = "subj", time_var = "tgrid")
head(result2)
#> # A tibble: 3 × 3
#>   sim_time sim_status    ID
#>      <dbl>      <dbl> <int>
#> 1     2.73          1   101
#> 2     6.16          1   102
#> 3     1.11          1   103
```
