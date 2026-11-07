# Censoring and Realistic Assessment Schedules

``` r

library(simtte)
```

Every simulation in
[`vignette("pkpd-time-to-event")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
is administratively censored at `end` only. This vignette covers the
rest of the assessment toolkit, layered on top of any already-simulated
`$events`-shaped data frame
([`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
or
[`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
alike): independent right censoring, interval censoring from a visit
schedule, and realistic (jittered, missed-visit, outcome-reactive)
visit-schedule generators.

``` r

sim_uncensored <- sim_tte_ode(model = "exponential",
  param = list(H0 = 0.05), n = 200, end = 30, delta = 2, seed = 1)
#> Building exponential_ode ... done.
table(sim_uncensored$events$sim_reason)
#> 
#> administrative          event 
#>             37            163
```

## 1. Right censoring

[`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
layers a second, independent censoring process on top of any
`$events`-shaped output, drawing a per-subject censoring time and taking
`min(event time, censoring time, end)`:

``` r

censored <- add_censoring(sim_uncensored$events,
  censoring = list(dist = "exponential", rate = 0.02), end = 30, seed = 2)
table(censored$sim_reason)
#> 
#> administrative       censored          event 
#>             20             42            138
```

`sim_status`/`sim_time` are updated in place (so any existing code
reading those two columns keeps working); `sim_reason` records why –
`"event"`, `"censored"` (pulled earlier by the independent draw), or
`"administrative"` (hit `end`).
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
also accepts a `censoring =` argument directly, drawing inside the same
seeded call as `U`:

``` r

sim_censored <- sim_tte_ode(model = "exponential", param = list(H0 = 0.09),
  n = 200, end = 30, delta = 2,
  censoring = list(dist = "exponential", rate = 0.02), seed = 1)
#> Loading model from cache.
table(sim_censored$events$sim_reason)
#> 
#> administrative       censored          event 
#>              6             31            163
```

`"exponential"`, `"weibull"` (fixed `shape`, solves for `scale`),
`"uniform"` (a fixed interval – the common case for uniform-accrual
administrative censoring), `"lognormal"` (fixed `sdlog`), `"gamma"`
(fixed `shape`, `rate` only), and `"user"` (any function
`n -> numeric(n)`) are the six `dist` options; see
[`?add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md).
[`draw_censoring_times()`](https://csetraynor.github.io/simtte/reference/draw_censoring_times.md)
draws directly from any of these specs without needing an `$events` data
frame at all, e.g. to check a distribution’s own shape before applying
it. Censoring that depends on a subject’s own simulated PK/PD trajectory
(informative dropout, e.g. driven by toxicity) is out of scope.

### Choosing the censoring rate for a target censoring fraction

**[`censoring_rate_for()`](https://csetraynor.github.io/simtte/reference/censoring_rate_for.md)
solves for a `dist` parameter (rate, scale, meanlog, or gamma rate) that
would give a target independent-censoring fraction, estimated from the
event times of an already-run, uncensored simulation.** It is an
approximation, not an exact calculation – the estimate depends on the
particular uncensored run it is given, so always redraw on a larger
cohort afterward and check the realized fraction landed close to the
target, the way this example does:

``` r

# 1. Solve: what rate gives ~20% independently censored, given this
#    model's own event-time scale?
event_times <- sim_uncensored$events$sim_time[
  sim_uncensored$events$sim_status == 1]
rate_20pct <- censoring_rate_for(event_times, target = 0.2)
rate_20pct
#> [1] 0.0205978

# 2. Redraw: apply that rate to a larger cohort.
check <- sim_tte_ode(model = "exponential", param = list(H0 = 0.09),
  n = 2000, end = 30, delta = 2,
  censoring = list(dist = "exponential", rate = rate_20pct), seed = 1)
#> Loading model from cache.

# 3. Check: did it land close to 20%?
mean(check$events$sim_reason == "censored")
#> [1] 0.1895
```

If the realized fraction in step 3 is not close enough, re-solve with a
larger or different uncensored sample rather than trusting the first
answer – the check is the point, not an optional extra.

## 2. Interval censoring

Right censoring still gives an *exact* event time when an event happens.
In a real trial, many endpoints (progression, response) are only checked
at scheduled visits – an event is known to have happened somewhere
between the last visit where the subject was event-free and the first
visit where it was detected.
[`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
maps an already-simulated (and, if used, already right-censored) outcome
onto that interval:

``` r

visits <- c(0, 5, 10, 15, 20, 25, 30)
interval <- add_interval_censoring(sim_censored$events, visits = visits)
head(interval[, c("ID", "sim_time", "sim_status", "sim_reason",
  "sim_time_left", "sim_time_right")])
#>   ID  sim_time sim_status sim_reason sim_time_left sim_time_right
#> 1  1 14.734534          1      event            10             15
#> 2  2 10.983644          1      event            10             15
#> 3  3  6.190283          1      event             5             10
#> 4  4  1.069800          1      event             0              5
#> 5  5  6.278294          0   censored             5            Inf
#> 6  6  1.190571          1      event             0              5
```

`sim_time`/`sim_status`/`sim_reason` are untouched – `sim_time_left`/
`sim_time_right` are added on top, so anything already reading the first
three keeps working. `sim_time_right` is `Inf` whenever there is no
known upper bound: a censored subject (never detected as an event under
this schedule), or – more subtly – an event that genuinely happened but
only *after* the subject’s last scheduled visit, so the schedule itself
never confirmed it. `sim_status` still reads `1` for that second case;
the exact columns are the simulation’s ground truth, independent of what
the visit schedule happened to observe.

**Order matters**: apply
[`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
*after*
[`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
never before – a dropout changes `sim_time`/ `sim_status`, and the
interval mapping needs to see the final value to know there was no
assessment after the dropout time.
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
own `censoring =`/`visits =` arguments apply them in this order
automatically.

For a nonparametric interval-censored survival fit, convert
`sim_time_right`’s `Inf` to `NA` (the `survival` package’s own
convention for an open upper bound) and hand both bounds to
`survival::Surv(..., type = "interval2")` – this and the “Informative
assessment schedules” chunk below are the only places in this vignette
that need the optional `survival` package (`Suggests`, not `Imports`);
both are skipped if it is not installed:

``` r

right_for_surv <- ifelse(is.infinite(interval$sim_time_right), NA,
  interval$sim_time_right)
fit <- survival::survfit(survival::Surv(time = sim_time_left,
  time2 = right_for_surv, type = "interval2") ~ 1,
  data = data.frame(sim_time_left = interval$sim_time_left,
    right_for_surv = right_for_surv))
summary(fit, times = visits)
#> Call: survfit(formula = survival::Surv(time = sim_time_left, time2 = right_for_surv, 
#>     type = "interval2") ~ 1, data = data.frame(sim_time_left = interval$sim_time_left, 
#>     right_for_surv = right_for_surv))
#> 
#>  time n.risk n.event survival std.err lower 95% CI upper 95% CI
#>     0    200       0   1.0000  0.0000       1.0000       1.0000
#>     5    116      72   0.6170  0.0355       0.5513       0.6906
#>    10     66      43   0.3736  0.0360       0.3093       0.4512
#>    15     32      24   0.2135  0.0321       0.1589       0.2868
#>    20     17      13   0.1210  0.0265       0.0787       0.1860
#>    25     13       4   0.0925  0.0238       0.0559       0.1532
#>    30      6       7   0.0427  0.0169       0.0197       0.0926
```

### Left censoring

An event with `sim_time_left == 0` – one that happened before the first
post-baseline visit – is a **left-censored observation**: the event is
known to have occurred before `sim_time_right`, with the exact time
unknown.
[`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
already produces this with no separate mechanism:

``` r

ev <- data.frame(ID = 1:3, sim_time = c(5, 15, 25),
  sim_status = c(1L, 1L, 0L))
add_interval_censoring(ev, visits = c(0, 10))
#>   ID sim_time sim_status sim_time_left sim_time_right
#> 1  1        5          1             0             10
#> 2  2       15          1            10            Inf
#> 3  3       25          0            10            Inf
```

Subject 1’s row (`sim_time_left = 0`, `sim_time_right = 10`) *is* left
censoring – “no lower bound” is exactly what `sim_time_left == 0` means.
There is no dedicated `left_censor_at =` option: it would only ever
rebuild, less completely, what a plain
[`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
call already gives for the whole cohort at once.

### Delayed entry (not implemented)

A genuinely different concept, easy to confuse with left censoring:
**delayed entry / left truncation**, where a subject is only observed
from some entry time `entry > 0` onward (e.g. enrollment after treatment
already started), and a subject whose event already happened before
`entry` is excluded from the risk set entirely rather than bounded by an
interval. This changes the analysis dataset’s *shape* – a
counting-process `survival::Surv(start, stop, event)` rather than a
single time column – not just how one event is bounded, and nothing in
this package currently represents it.

## 3. Realistic visit schedules

A fixed schedule assumes every visit lands exactly on the nominal day
and every subject attends every one.
[`visit_schedule()`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)
builds a per-subject schedule with fixed spacing and, optionally, random
jitter (uniform or truncated-normal);
[`thin_visits()`](https://csetraynor.github.io/simtte/reference/thin_visits.md)
and
[`visit_schedule_informative()`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md)
layer more realistic visit-process randomness on top.

``` r

sched <- visit_schedule(n = 200, every = 4, end = 24, jitter = 1, seed = 1)
head(sched)
#>   ID      time
#> 1  1  0.000000
#> 2  1  3.531017
#> 3  1  7.744248
#> 4  1 12.145707
#> 5  1 16.816416
#> 6  1 19.403364

sim_jittered <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
  n = 200, end = 24, delta = 2, visits = list(every = 4, jitter = 1),
  seed = 1)
#> Loading model from cache.
head(sim_jittered$events[, c("sim_time", "sim_status", "sim_time_left",
  "sim_time_right")])
#>    sim_time sim_status sim_time_left sim_time_right
#> 1 24.000000          0     19.362337            Inf
#> 2 19.770560          1     16.435871      20.922820
#> 3 11.142509          1      8.895933      12.637269
#> 4  1.925641          1      0.000000       4.907465
#> 5 24.000000          0     19.835713            Inf
#> 6  2.143028          1      0.000000       3.694461
```

### Missed visits and assessment dropout

[`thin_visits()`](https://csetraynor.github.io/simtte/reference/thin_visits.md)
reads no simulated outcome at all – it only thins a schedule, so it can
never bias an interval-censored analysis (the visit process stays
independent of the outcome, the assumption most such analyses make).
Each non-baseline visit is independently missed with probability
`p_miss`; independently of that, a subject may also start dropping out
at a random visit and never return, with probability `p_dropout`:

``` r

realized <- thin_visits(sched, p_miss = 0.1, p_dropout = 0.05, seed = 2)
nrow(sched) - nrow(realized) # visits lost to missed/dropout
#> [1] 156
out <- add_interval_censoring(sim_censored$events, visits = realized)
head(out[, c("sim_time", "sim_status", "sim_time_left", "sim_time_right")])
#>    sim_time sim_status sim_time_left sim_time_right
#> 1 14.734534          1     12.145707      16.816416
#> 2 10.983644          1      8.321596      12.258228
#> 3  6.190283          1      4.374046       7.768207
#> 4  1.069800          1      0.000000       3.760070
#> 5  6.278294          0      3.534441            Inf
#> 6  1.190571          1      0.000000       3.964160
```

[`thin_visits()`](https://csetraynor.github.io/simtte/reference/thin_visits.md)
also accepts a plain vector schedule directly, given the subject IDs to
expand it across (`ids =`), and composes with
[`add_censoring()`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
in either relative order – it only ever touches the visit schedule,
never `events`.

### Informative assessment schedules

[`visit_schedule_informative()`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md)
is different in kind: it reads `$events` and lets the visit process
**react to the simulated outcome** – a visit shortly before a subject’s
event may be missed at an elevated rate (`miss_near_event`, e.g. “too
unwell to attend”), and an unscheduled extra visit may be added shortly
after it (`extra_visit_after_event`, e.g. “brought in for follow-up
after a reported symptom”). This makes the resulting schedule
**informative by construction** – whether/when a subject is assessed now
depends on whether/when they had the event – which every function in
this package says loudly, on every call, because it breaks the
assumption most interval-censored analyses (including the Turnbull fit
above) make about assessment timing. Use it to *study* that bias, not as
a general visit generator:

``` r

rate <- 0.08
sim_small <- sim_tte_ode(model = "exponential", param = list(H0 = rate),
  n = 800, end = 24, delta = 2, seed = 1)
base_sched <- visit_schedule(n = 800, every = 3, end = 24, seed = 1)
base_sched$ID <- sim_small$events$ID[base_sched$ID]
visit_times <- seq(0, 24, by = 3)
analytical_surv <- exp(-rate * visit_times)

fit_surv <- function(schedule) {
  out <- add_interval_censoring(sim_small$events, visits = schedule)
  right <- ifelse(is.infinite(out$sim_time_right), NA, out$sim_time_right)
  fit <- survival::survfit(survival::Surv(time = sim_time_left, time2 = right,
    type = "interval2") ~ 1,
    data = data.frame(sim_time_left = out$sim_time_left, right = right))
  summary(fit, times = visit_times, extend = TRUE)$surv
}

c1_sched <- thin_visits(base_sched, p_miss = 0.3, seed = 2)
c2_sched <- visit_schedule_informative(base_sched, sim_small$events,
  end = 24, p_miss_base = 0.1, miss_near_event = list(window = 3, p = 0.9),
  seed = 2)

data.frame(time = visit_times,
  truth = round(analytical_surv, 3),
  c1_non_informative = round(fit_surv(c1_sched), 3),
  c2_informative = round(fit_surv(c2_sched), 3))
#>   time truth c1_non_informative c2_informative
#> 1    0 1.000              1.000          1.000
#> 2    3 0.787              0.796          0.666
#> 3    6 0.619              0.623          0.576
#> 4    9 0.487              0.509          0.443
#> 5   12 0.383              0.372          0.332
#> 6   15 0.301              0.282          0.237
#> 7   18 0.237              0.220          0.209
#> 8   21 0.186              0.173          0.143
#> 9   24 0.147              0.132          0.140
```

At the same *nominal* miss rate (0.3 for C1’s uniform miss; 0.1
background + 0.9 near an event for C2 – additive, not overriding: a
near-window visit for an event subject is missed with probability
`1 - (1 - p_miss_base) * (1 - p)`, so the two mechanisms compete rather
than one replacing the other), C2’s estimate strays further from the
analytical curve than C1’s – missing visits specifically *because* a
subject is about to have the event pulls the estimated survival away
from the truth in a way a uniform miss rate does not. A visit-process
generator reacting to a subject’s own PK/PD trajectory, rather than the
exact event time, is out of scope for this release – it would need
`keep_trajectory = TRUE` and interpolation between reported times, the
same complexity class this package has consistently declined to add.
