# Thin a visit schedule with missed visits and dropout (C1)

Layers schedule-only randomness on top of a fixed visit schedule – no
simulated outcome is read, so this is the "safe", always-
non-informative member of the visit-process generators (contrast with
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md),
which reads the simulated outcome by design). Each non-baseline visit is
independently missed with probability `p_miss`; independently of that,
each subject may also start dropping out: with probability `p_dropout`,
one of their non-baseline visits is chosen uniformly at random as the
dropout onset, and it and every later visit is dropped. The baseline
visit (`time = 0`) is never dropped by either mechanism.

## Usage

``` r
thin_visits(visits, p_miss = 0, p_dropout = 0, ids = NULL, seed = NULL)
```

## Arguments

- visits:

  A per-subject schedule (a data frame with `ID`/ `time` columns, e.g.
  from
  [`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)),
  or a common numeric vector – in which case `ids` is required, and the
  vector is expanded to one copy per ID first (via the same internal
  validation
  [`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
  itself uses).

- p_miss, p_dropout:

  Numeric scalars in `[0, 1]`, default `0` (no thinning; `visits` is
  returned as given). The dropout-onset visit is drawn from the
  subject's *original* non-baseline visits, independently of any
  `p_miss` draws for that subject.

- ids:

  Required only when `visits` is a numeric vector; the subject IDs to
  expand it across (e.g. `events$ID`).

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before any
  draw. Leave `NULL` inside a larger already-seeded workflow.

## Value

A data frame with `ID`/`time` columns – a (possibly shorter, per
subject) schedule directly usable as
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)'s
`visits` argument. Must be applied before
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
and composes with
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
in either relative order (they touch disjoint objects: this function
only ever reads/writes a visit schedule, never `events`).

## See also

[`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md),
which builds the schedule to thin;
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md),
the outcome-reactive counterpart;
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
which takes the thinned schedule as `visits`.

## Examples

``` r
sched <- visit_schedule(n = 20, every = 4, end = 24, seed = 1)
thinned <- thin_visits(sched, p_miss = 0.1, p_dropout = 0.05, seed = 2)
nrow(thinned) <= nrow(sched)
#> [1] TRUE
```
