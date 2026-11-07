# Generate an outcome-informative visit schedule (C2)

Layers visit-process randomness that **reacts to the simulated event
time** on top of a fixed visit schedule: a scheduled visit shortly
before a subject's event may be missed at an elevated rate
(`miss_near_event`, e.g. "too unwell to attend"), and/or an unscheduled
extra visit may be added shortly after it (`extra_visit_after_event`,
e.g. "brought in for follow-up after a reported symptom"). Both are
optional and independent of each other; supplying neither is equivalent
to
[`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md)
with `p_miss = p_miss_base`.

## Usage

``` r
visit_schedule_informative(
  visits,
  events,
  end,
  miss_near_event = NULL,
  extra_visit_after_event = NULL,
  p_miss_base = 0,
  id_var = "ID",
  time_var = "sim_time",
  status_var = "sim_status",
  seed = NULL
)
```

## Arguments

- visits:

  A per-subject schedule (data frame with `ID`/ `time`), or a common
  numeric vector (expanded per subject using `events[[id_var]]` – no
  separate `ids` argument needed here, unlike
  [`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md),
  since `events` already supplies the ID list).

- events:

  Data frame of simulated outcomes (typically the `$events` element of
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md))
  with `ID`, `sim_time`, `sim_status` columns by default, **after** any
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
  call (see "Pipeline position" below) – column names configurable via
  `id_var`/`time_var`/`status_var`.

- end:

  Numeric scalar. An unscheduled extra visit (`extra_visit_after_event`)
  landing after `end` is dropped, not added – the same administrative
  horizon
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
  uses.

- miss_near_event:

  `NULL` (default, no effect), or
  `list(window = <positive scalar>, p = <scalar in [0, 1]>)`. For an
  event subject (`sim_status == 1`), a scheduled visit `v` with
  `t - window <= v < t` (`t` the event time) is missed with probability
  `1 - (1 - p_miss_base) * (1 - p)` – `p_miss_base`'s background miss
  chance and `p`'s near-event miss chance are independent, competing
  risks, so they add rather than one overriding the other (this
  collapses exactly to `p` when `p_miss_base = 0`, the common case). Has
  no effect on a censored subject, or on a visit outside the window –
  both use `p_miss_base` alone.

- extra_visit_after_event:

  `NULL` (default, no effect), or a distribution spec in
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)'s
  `censoring` shape (`list(dist = "exponential", rate = ...)`,
  `"weibull"`, `"uniform"`, `"lognormal"`, `"gamma"`, or `"user"`). For
  an event subject only, one delay is drawn and an unscheduled visit is
  added at `sim_time + delay`, unless that would exceed `end`.

- p_miss_base:

  Numeric scalar in `[0, 1]`, default `0`. Background miss probability
  applied to every non-baseline visit not otherwise covered by
  `miss_near_event` (i.e. every visit for a censored subject, and every
  out-of-window visit for an event subject) – the same mechanism as
  [`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md)'s
  `p_miss`, folded in here so a caller does not need to chain both
  functions just to add a uniform background miss rate on top of the
  near-event effect.

- id_var, time_var, status_var:

  Character scalars naming the subject-ID, event/censoring-time, and
  status columns of `events`. Default
  `"ID"`/`"sim_time"`/`"sim_status"`.

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before any
  draw.

## Value

A data frame with `ID`/`time` columns – directly usable as
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)'s
`visits` argument.

## Informative assessment

**This function's output is not a neutral assessment schedule.** By
construction, whether/when a subject is assessed now depends on
whether/when they had the simulated event – exactly the situation most
interval-censored analyses (Turnbull NPMLE, standard parametric
interval-censored likelihoods) assume does *not* hold when they treat
assessment times as ignorable. Using this function's schedule is
appropriate for *studying* that bias (e.g. comparing a fit against
[`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md)'s
non-informative schedule at the same nominal miss rate, as the
vignette's "Interval censoring" section does), not as a drop-in,
general-purpose visit generator. A
[`message()`](https://rdrr.io/r/base/message.html) naming this is
emitted on every call.

## Pipeline position

`events` should already reflect any independent right censoring: exact
simulation -\>
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
-\> `visit_schedule_informative()` -\>
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md).
Generating the schedule from pre-censoring event times would make the
near-event/extra-visit mechanics react to a `sim_time` that the final
`events` no longer has – the same ordering reasoning as
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)'s
own "Ordering with right censoring" section.

## See also

[`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md),
which builds the base schedule passed as `visits`;
[`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md),
the non-informative counterpart;
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
which takes the returned schedule as its `visits`;
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
which should run before this.

## Examples

``` r
# \donttest{
sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
  n = 50, end = 24, delta = 2, seed = 1)
#> Loading model from cache.
sched <- visit_schedule(n = 50, every = 4, end = 24, seed = 1)
realized <- visit_schedule_informative(sched, sim$events, end = 24,
  miss_near_event = list(window = 3, p = 0.5), seed = 2)
#> visit_schedule_informative(): the returned schedule depends on the simulated event times; interval-censored analyses assuming non-informative assessment times will be biased by design (see ?visit_schedule_informative "Informative assessment").
out <- add_interval_censoring(sim$events, visits = realized)
head(out)
#>   ID  sim_time sim_status     sim_reason sim_time_left sim_time_right
#> 1  1 24.000000          0 administrative            24            Inf
#> 2  2 19.770560          1          event            16             20
#> 3  3 11.142509          1          event             8             12
#> 4  4  1.925641          1          event             0              4
#> 5  5 24.000000          0 administrative            24            Inf
#> 6  6  2.143028          1          event             0              4
# }
```
