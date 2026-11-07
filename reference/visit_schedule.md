# Generate a per-subject visit schedule with fixed spacing and jitter

Builds the per-subject data-frame form of `visits` that
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
accepts: assessments every `every` time units from `0` to `end`,
optionally perturbed per subject by up to `jitter` (a common pattern in
trial simulation – scheduled visits rarely land exactly on the nominal
day). The baseline visit (`time = 0`) is never jittered, so the result
always satisfies
[`add_interval_censoring()`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)'s
"must start at 0" requirement.

## Usage

``` r
visit_schedule(
  n,
  every,
  end,
  jitter = 0,
  jitter_dist = c("uniform", "normal"),
  jitter_trunc = 2,
  seed = NULL
)
```

## Arguments

- n:

  Integer. Number of subjects; schedules are generated for `ID = 1:n`
  (remap the `ID` column afterward if your `events` uses different
  subject identifiers).

- every:

  Positive numeric scalar. Nominal spacing between visits.

- end:

  Non-negative numeric scalar. Last nominal visit time (the schedule is
  `seq(0, end, by = every)` before jitter).

- jitter:

  Non-negative numeric scalar (default `0`, no jitter). Each
  non-baseline visit is perturbed by an independent draw from
  `jitter_dist`: `jitter` is the draw's half-width for `"uniform"`
  (`stats::runif(1, -jitter, jitter)`) or its standard deviation for
  `"normal"`, which is truncated at `+/- jitter_trunc * jitter` (see
  `jitter_trunc`) and drawn via inverse-CDF
  ([`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html) on a
  [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html) restricted to
  the truncated range) – one uniform draw per visit either way, so both
  distributions consume the seeded stream the same shape of way.

- jitter_dist:

  `"uniform"` (default; unchanged behavior and seeds from before this
  argument existed) or `"normal"` (truncated, see `jitter_trunc`; its
  seeds changed when the truncation/bound was added – see `NEWS.md`).

- jitter_trunc:

  Positive numeric scalar, default `2`. Only used for
  `jitter_dist = "normal"`: the truncation half-width, in standard
  deviations, of the truncated normal jitter draw. Ignored for
  `"uniform"`.

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before any
  jitter is drawn. Leave `NULL` when calling this from inside a larger
  already-seeded workflow (this is how
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
  own `visits =` argument uses it when given a jitter spec, so one
  `seed` reproduces `U`, any between-subject draw, censoring, and the
  visit jitter together).

## Value

A data frame with `ID`/`time` columns, one row per subject per visit –
directly usable as
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)'s
`visits` argument.

## Details

`jitter` is validated against `every` *before* any draw is made: the
maximum possible deviation (`jitter` for `"uniform"`,
`jitter_trunc * jitter` for `"normal"`) must be less than `every / 2`,
which guarantees – by construction, not by a post-hoc check – that
jittered visits can never cross order or collide (adjacent nominal
visits are `every` apart; two deviations each smaller than half of that
cannot close the gap) or go negative (the first non-baseline visit is at
least `every / 2` above `0`). A `jitter` violating this bound is an
error, naming the bound, rather than a silently reordered or clamped
schedule.

## See also

[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md),
[`thin_visits`](https://csetraynor.github.io/simtte/reference/thin_visits.md)
and
[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md),
to layer missed visits/dropout on top of the schedule this function
returns.

## Examples

``` r
visit_schedule(n = 5, every = 4, end = 20, jitter = 1, seed = 1)
#>    ID      time
#> 1   1  0.000000
#> 2   1  3.531017
#> 3   1  7.744248
#> 4   1 12.145707
#> 5   1 16.816416
#> 6   1 19.403364
#> 7   2  0.000000
#> 8   2  4.796779
#> 9   2  8.889351
#> 10  2 12.321596
#> 11  2 16.258228
#> 12  2 19.123573
#> 13  3  0.000000
#> 14  3  3.411949
#> 15  3  7.353114
#> 16  3 12.374046
#> 17  3 15.768207
#> 18  3 20.539683
#> 19  4  0.000000
#> 20  4  3.995398
#> 21  4  8.435237
#> 22  4 12.983812
#> 23  4 15.760070
#> 24  4 20.554890
#> 25  5  0.000000
#> 26  5  4.869410
#> 27  5  7.424285
#> 28  5 12.303348
#> 29  5 15.251110
#> 30  5 19.534441
visit_schedule(n = 5, every = 4, end = 20, jitter = 0.4,
  jitter_dist = "normal", seed = 1)
#>    ID      time
#> 1   1  0.000000
#> 2   1  3.762307
#> 3   1  7.875644
#> 4   1 12.070080
#> 5   1 16.489834
#> 6   1 19.684673
#> 7   2  0.000000
#> 8   2  4.470521
#> 9   2  8.574241
#> 10  2 12.157895
#> 11  2 16.125600
#> 12  2 19.442576
#> 13  3  0.000000
#> 14  3  3.690248
#> 15  3  7.650715
#> 16  3 12.185418
#> 17  3 15.887623
#> 18  3 20.279395
#> 19  4  0.000000
#> 20  4  3.997798
#> 21  4  8.218691
#> 22  4 12.749539
#> 23  4 15.883566
#> 24  4 20.288759
#> 25  5  0.000000
#> 26  5  4.548691
#> 27  5  7.698154
#> 28  5 12.148498
#> 29  5 15.572502
#> 30  5 19.764259
```
