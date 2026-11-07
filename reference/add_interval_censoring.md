# Map simulated event/censoring times onto an assessment-time interval

Adds interval-censoring bounds to an already-simulated (and, if
applicable, already right-censored – see "Ordering" below) events data
frame: given a visit schedule, an event known to have occurred exactly
at `sim_time` is only *detectable* at the first visit at or after it,
and a censored subject is only known event-free through their last visit
at or before their own censoring time. Works on the output of
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
or
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
alike, same as
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md).

## Usage

``` r
add_interval_censoring(
  events,
  visits,
  id_var = "ID",
  time_var = "sim_time",
  status_var = "sim_status"
)
```

## Arguments

- events:

  Data frame of simulated outcomes: typically the `$events` element
  returned by
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
  or the output of
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)/[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
  optionally after
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md).
  Needs `ID`, `sim_time`, `sim_status` columns by default (rename via
  `id_var`/`time_var`/`status_var`). Every other column is left
  untouched, and so are these three: this function only adds
  `sim_time_left`/`sim_time_right` (see "Return" below).

- visits:

  A visit schedule, either:

  a numeric vector

  :   applied to every subject (e.g. `c(0, 4, 8, 12, 24)`). Must start
      at `0` and be strictly increasing (sorted, no duplicate visit
      times).

  a data frame

  :   with `ID`/`time` columns, for a per-subject schedule. Every ID in
      `events` must have a schedule (starting at `0`, strictly
      increasing), and every ID in `visits` must be one that appears in
      `events`.

  [`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)
  builds the data-frame form from a fixed spacing plus optional
  per-subject jitter, rather than this function accepting a jitter spec
  itself – see its own documentation.

- id_var, time_var, status_var:

  Character scalars naming the subject-ID, event/censoring-time, and
  status columns of `events`. Default
  `"ID"`/`"sim_time"`/`"sim_status"`.

## Value

`events` with two new columns:

- sim_time_left:

  Left (open) bound of the interval, `L`.

- sim_time_right:

  Right (closed) bound, `R`, or `Inf` for an interval with no known
  upper bound (a censored subject, or an event that occurred after the
  subject's last visit and so was never confirmed by the schedule –
  `sim_status` still reads `1` in that case: the exact columns are the
  simulation's ground truth, unaffected by this mapping. The two columns
  answer different questions: `sim_status` says what happened, the
  bounds say what the schedule could observe).

An event with `sim_time_left == 0` (an event before the first
post-baseline visit) **is** a left-censored observation – the event is
known to have occurred before `sim_time_right`, exact time unknown –
with no separate mechanism needed; see the vignette's "Left censoring"
subsection for a worked example. `sim_time`/`sim_status`/`sim_reason`
(if present) are returned unchanged – this function only adds columns,
so every existing consumer of those three keeps working. For
`survival::Surv(time = sim_time_left, time2 = sim_time_right, type = "interval2")`,
which expects `NA` rather than `Inf` for an open upper bound, convert
with `ifelse(is.infinite( sim_time_right), NA, sim_time_right)` – see
the vignette's "Interval censoring" section for a full worked example.

## Ordering with right censoring

Call this **after**
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
(or after the exact simulation, if no independent right censoring was
requested), not before: the interval mapping reads
`sim_time`/`sim_status` as given, so a dropout applied afterward would
leave the already-added interval columns stale and inconsistent with the
final exact columns.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
own `censoring`/ `visits` convenience arguments apply them in this order
internally.

## See also

[`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md),
to build a per-subject schedule with fixed spacing and optional jitter;
[`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md),
to drop visits at random (missed visits, dropout);
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md),
for visits that react to the simulated outcome.
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
which should run first.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
which can apply this step itself through its `visits` argument.

## Examples

``` r
# \donttest{
sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
  n = 50, end = 24, delta = 2, seed = 1)
#> Loading model from cache.
out <- add_interval_censoring(sim$events, visits = c(0, 4, 8, 12, 16, 20, 24))
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
