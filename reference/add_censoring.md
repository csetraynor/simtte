# Apply independent right censoring to a simulated events data frame

Adds a second, independent censoring process to an already-simulated
time-to-event data frame: a censoring time \\C_i\\ is drawn per subject
from `censoring`, and the observed outcome becomes \\\min(T_i, C_i,
end)\\, with `sim_status = 1` only if \\T_i\\ (the original `sim_time`)
is the minimum. Works on the output of
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
or
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
alike – all three already share the same `ID`/`sim_time`/`sim_status`
shape.

## Usage

``` r
add_censoring(
  events,
  censoring,
  end,
  id_var = "ID",
  time_var = "sim_time",
  status_var = "sim_status",
  seed = NULL
)
```

## Arguments

- events:

  Data frame of simulated outcomes: typically the `$events` element
  returned by
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
  or the output of
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)/[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md).
  Needs `ID`, `sim_time`, `sim_status` columns by default – column names
  are configurable via `id_var`/`time_var`/`status_var` for a data frame
  using different names. Any other column is left untouched.

- censoring:

  A list describing the censoring-time distribution, with a `dist`
  element:

  `"exponential"`

  :   `rate` (positive scalar). `C ~ Exp(rate)`.

  `"weibull"`

  :   `shape`, `scale` (positive scalars). `C ~ Weibull(shape, scale)`
      ([`rweibull`](https://rdrr.io/r/stats/Weibull.html)
      parameterization).

  `"uniform"`

  :   `min`, `max` (`0 <= min < max`). `C ~ Uniform(min, max)` – the
      common case for a uniform-accrual administrative censoring time.

  `"lognormal"`

  :   `meanlog` (finite scalar), `sdlog` (positive scalar).
      `C ~ Lognormal(meanlog, sdlog)`
      ([`rlnorm`](https://rdrr.io/r/stats/Lognormal.html)
      parameterization).

  `"gamma"`

  :   `shape`, `rate` (positive scalars). `C ~ Gamma(shape, rate)`
      ([`rgamma`](https://rdrr.io/r/stats/GammaDist.html)
      parameterization) – `rate` only, matching `"exponential"`'s own
      rate-only convention (a `scale` element is an error, not silently
      accepted).

  `"user"`

  :   `fun`, a function `n -> numeric(n)` returning `n` non-negative,
      finite censoring times. The escape hatch for any distribution not
      covered above.

- end:

  Numeric scalar. Administrative censoring horizon – the third term of
  `min(T, C, end)`. Always wins: no observed time in the result exceeds
  `end`. Required (part of the censoring formula itself, not inferred
  from `events`).

- id_var, time_var, status_var:

  Character scalars naming the subject-ID, event/censoring-time, and
  status columns of `events`. Default `"ID"`/`"sim_time"`/`"sim_status"`
  (the column names every simtte simulation function already uses).

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before \\C\\
  is drawn. Leave `NULL` when calling this from inside a larger
  already-seeded workflow (this is how
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
  own `censoring =` argument uses it, so that one `seed` reproduces `U`,
  any between-subject draw, and the censoring draw together, in that
  fixed order – see
  [`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  "Reproducibility").

## Value

`events` with `time_var`/`status_var` updated in place, plus one new
column, `sim_reason`: `"event"`, `"censored"` (the independent draw, or
a pre-existing censoring time from `events` itself, was earlier than
both the event and `end`), or `"administrative"` (the observed time
equals `end` and the subject is not an event). Every other column of
`events` is returned unchanged.

## See also

[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
whose `$events` this usually post-processes (or pass its `censoring`
argument to do this in one call);
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
the next step if you also need assessment intervals.
[`censoring_rate_for`](https://csetraynor.github.io/simtte/reference/censoring_rate_for.md),
to choose `rate`/ `scale` for a target censoring fraction;
[`draw_censoring_times`](https://csetraynor.github.io/simtte/reference/draw_censoring_times.md),
to draw from the same `censoring`-shaped spec directly, with no `events`
data frame needed; the "Right censoring" section of
[`vignette("pkpd-time-to-event", package = "simtte")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
for a worked example, including the scale-helper workflow.

## Examples

``` r
# \donttest{
sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
  n = 100, end = 30, delta = 2, seed = 1)
#> Building exponential_ode ... 
#> done.
out <- add_censoring(sim$events, censoring = list(dist = "exponential",
  rate = 0.02), end = 30, seed = 2)
table(out$sim_reason)
#> 
#> administrative       censored          event 
#>              9             22             69 

# Also works directly on the legacy sim_tte() output.
lp <- matrix(rnorm(50, 0, 0.5), nrow = 50)
legacy <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
  time = seq(0, 30, by = 1), type = "weibull", end_time = 30)
#> Building weibull ... 
#> done.
out2 <- add_censoring(legacy, censoring = list(dist = "uniform",
  min = 5, max = 30), end = 30, seed = 3)
table(out2$sim_reason)
#> 
#> censored    event 
#>        1       49 
# }
```
