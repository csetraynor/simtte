# Explore Prognostic Index and Quantile Survival Times

Calculates the difference in survival probability at a given quantile
time across a range of prognostic index values. Useful for understanding
how different log hazard ratios affect survival outcomes.

## Usage

``` r
explore_pi_tq_surv(q = 0.5, pi, mu, shape, type, times, basehaz, end_time, ...)
```

## Arguments

- q:

  Numeric scalar. Survival quantile of interest (default 0.5, i.e.
  median survival).

- pi:

  Numeric vector. Prognostic index values to evaluate. Must be finite.

- mu:

  Numeric scalar. Intercept parameter. Must be finite.

- shape:

  Numeric scalar or vector. Shape parameter(s) for the Weibull model.
  Must be finite and strictly positive.

- type:

  Character string. Model type: `"weibull"` or `"ms"`.

- times:

  Numeric vector. Output time grid; see the "Time grid and event-time
  resolution" and "M-spline hazard carry convention" sections of
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  for the exact contract. For M-spline models this must also be the time
  points corresponding to the rows of `basehaz`. If missing, defaults to
  `seq(0, end_time, by = 1)`.

- basehaz:

  Numeric matrix. Baseline hazard (M-spline models), one row per element
  of `times`. Must be finite and non-negative; see
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)'s
  "M-spline hazard carry convention" for how it is interpreted between
  supplied times.

- end_time:

  Numeric scalar. Administrative censoring time. See the `end_time`
  contract documented in
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  (including the M-spline `min(times) <= end_time <= max(times)`
  requirement).

- ...:

  Additional arguments passed to
  [`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html). `tgrid`,
  `obsonly`, `nocb`, `carry_out`/`carry.out`, and `data` cannot be
  overridden; see
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md).

## Value

A data frame including columns `lp`, `p11`, `survdiff_tq`, and model
parameters.

## See also

[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
the closed-form engine this function diagnoses, and the function to call
once you have picked `pi`, `mu` and `shape`;
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
for the underlying trajectory-to-event mechanism.

## Examples

``` r
# Small fast example
data_sim <- explore_pi_tq_surv(
  pi = seq(-1, 1, by = 0.5),
  mu = -1,
  shape = 1.1,
  end_time = 10,
  type = "weibull"
)
#> Loading model from cache.
head(data_sim)
#>   ID time basehaz_id   lp mu shape       p11 p11_baseline survdiff_tq
#> 1  1    2          1 -1.0 -1   1.1 0.7481906    0.4544957   0.2936949
#> 2  2    2          1 -0.5 -1   1.1 0.6198424    0.4544957   0.1653467
#> 3  3    2          1  0.0 -1   1.1 0.4544957    0.4544957   0.0000000
#> 4  4    2          1  0.5 -1   1.1 0.2724972    0.4544957  -0.1819985
#> 5  5    2          1  1.0 -1   1.1 0.1172383    0.4544957  -0.3372574
# \donttest{
# Larger range example
data_sim2 <- explore_pi_tq_surv(
  pi = seq(-3, 3, by = 0.1),
  mu = -1,
  shape = 1.1,
  end_time = 200,
  type = "weibull"
)
#> Loading model from cache.
plot(survdiff_tq ~ lp, data = data_sim2)

# }
```
