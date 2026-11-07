# Solve for a censoring-distribution parameter giving a target censoring fraction

Users think in terms of "about 20 `end`", not in the rate/scale units
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
needs, and the right units depend on the model's own time scale. This
solves for one parameter of a fixed-family censoring distribution – the
exponential/gamma `rate`, the Weibull `scale` at a fixed `shape`, or the
lognormal `meanlog` at a fixed `sdlog` – that would give a target
independent-censoring fraction, treating `times` – the event times from
an already-run, **uncensored** simulation – as an empirical sample of
the latent event time \\T\\: it solves `mean(F_C(times)) == target` via
[`uniroot`](https://rdrr.io/r/stats/uniroot.html), where `F_C` is the
censoring distribution's own CDF
([`pexp`](https://rdrr.io/r/stats/Exponential.html)/
[`pweibull`](https://rdrr.io/r/stats/Weibull.html)/[`plnorm`](https://rdrr.io/r/stats/Lognormal.html)/
[`pgamma`](https://rdrr.io/r/stats/GammaDist.html)). This is the Monte
Carlo estimate of \\P(C \< T)\\ under independence, using the supplied
sample in place of \\T\\'s true (generally not closed-form once
between-subject variability or covariates are in the model) distribution
– **approximate, and dependent on the run it is given**: a different
(even same-seed-different-`n`) uncensored run will give a slightly
different answer. Re-check the realized fraction on a large cohort after
drawing with the solved parameter; see
[`vignette("pkpd-time-to-event", package = "simtte")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
"Right censoring" for a worked before/after check.

## Usage

``` r
censoring_rate_for(
  times,
  target,
  dist = c("exponential", "weibull", "lognormal", "gamma"),
  shape = NULL,
  sdlog = NULL
)
```

## Arguments

- times:

  Numeric vector of event times (`sim_status == 1` only) from an
  uncensored simulation. Must be non-negative and finite.

- target:

  Numeric scalar in \\(0, 1)\\. The desired \\P(C \< T)\\.

- dist:

  `"exponential"` (default), `"weibull"`, `"lognormal"`, or `"gamma"`.

- shape:

  Required (positive numeric scalar) for `dist = "weibull"` (solves for
  `scale`) or `dist = "gamma"` (solves for `rate`); the shape is fixed
  by the caller in either case. Unused otherwise.

- sdlog:

  Required (positive numeric scalar) for `dist = "lognormal"` (solves
  for `meanlog`, with `sdlog` fixed by the caller). Unused otherwise.

## Value

Numeric scalar: `rate` for `dist = "exponential"` or `"gamma"`, `scale`
for `dist = "weibull"`, `meanlog` for `dist = "lognormal"`.

## Details

Only genuine event times are informative here: an
administratively-censored subject's true \\T_i\\ is unknown beyond
`end`, so it carries no information about where independent censoring
would actually bite. Pass
`times = events$sim_time[events$sim_status == 1]` from an uncensored
run, not the full `events$sim_time`.

## See also

[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
which takes the solved parameter in its `censoring` list;
[`draw_censoring_times`](https://csetraynor.github.io/simtte/reference/draw_censoring_times.md),
to check the draws directly;
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
whose uncensored `$events` supply `times`.

## Examples

``` r
# \donttest{
sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
  n = 500, end = 30, delta = 2, seed = 1)
#> Loading model from cache.
event_times <- sim$events$sim_time[sim$events$sim_status == 1]
rate <- censoring_rate_for(event_times, target = 0.2)
out <- add_censoring(sim$events, censoring = list(dist = "exponential",
  rate = rate), end = 30, seed = 2)
mean(out$sim_reason == "censored") # close to 0.2, on the events subset
#> [1] 0.24
# }
```
