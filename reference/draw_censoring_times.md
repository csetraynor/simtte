# Draw independent times from a censoring/delay distribution spec

A thin, exported wrapper over the same distribution dispatch
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
uses internally for its own `censoring` argument, and
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md)
uses for `extra_visit_after_event`. Exported so a caller that only wants
`n` raw draws from a distribution spec – to plan a scenario's expected
censoring/dropout burden, or to pre-validate a spec before running many
replicates of it – does not need to build an `events`-shaped data frame
just to reach this logic.
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)/[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md)
are unaffected and still call the same, unexported dispatcher directly;
this function does not replace either.

## Usage

``` r
draw_censoring_times(spec, n, seed = NULL)
```

## Arguments

- spec:

  A list describing the distribution, in exactly the shape documented
  for
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)'s
  own `censoring` argument (`dist` one of `"exponential"`, `"weibull"`,
  `"uniform"`, `"lognormal"`, `"gamma"`, `"user"`, plus that
  distribution's own parameters).

- n:

  Integer. Number of draws.

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before the
  draw (same convention as
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)'s
  own `seed`).

## Value

Numeric vector of length `n`, all finite and non-negative.

## See also

[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
and
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md),
the two functions that use this same spec shape internally;
[`censoring_rate_for`](https://csetraynor.github.io/simtte/reference/censoring_rate_for.md),
to choose a spec's rate for a target censoring fraction.

## Examples

``` r
draw_censoring_times(list(dist = "exponential", rate = 0.05), n = 10,
  seed = 1)
#>  [1] 15.103637 23.632856  2.914135  2.795905  8.721373 57.899371 24.591241
#>  [8] 10.793657 19.131350  2.940920
```
