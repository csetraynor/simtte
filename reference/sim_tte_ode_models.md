# List sim_tte_ode()'s bundled library models

The built-in, bundled model names accepted by
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
`model` argument – the character-string case, not a
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
user model. Exported so a caller building a scenario grid over `model`
(e.g. a companion package comparing several hazard mechanisms) can
enumerate and validate model names programmatically, and look up each
model's valid `omega` targets (see
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
"Between-subject variability"), instead of hardcoding a copy of
`names(.ODE_LIBRARY_FILES)` (private, `@noRd`) or reaching in via `:::`.

## Usage

``` r
sim_tte_ode_models()
```

## Value

A data frame, one row per model: `model` (character, passable to
`sim_tte_ode(model = ...)`), `bsv_targets` (character, comma-separated
valid `omega` target names, or `NA` when the model takes no `omega` at
all).

## Details

`"mspline"` is included with `bsv_targets = NA`: it takes no
between-subject variability (its hazard is entirely determined by
`knots`/`boundary_knots`/`coefs`), which is a different situation from
the exponential/Weibull/Gompertz models' empty-but-valid target list
(also `NA` here, for the same reason: `omega` is rejected outright for
those three, exactly as for `"mspline"`, since none of the four has a
structural PK/PD parameter for BSV to attach to).

## See also

[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
where these names go as `model`.
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md),
to bring in a model beyond this list.
[`simtte_example_models`](https://csetraynor.github.io/simtte/reference/simtte_example_models.md)
lists a different set: editable templates for
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
not
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
models.

## Examples

``` r
sim_tte_ode_models()
#>          model                                                bsv_targets
#> 1  exponential                                                       <NA>
#> 2      weibull                                                       <NA>
#> 3     gompertz                                                       <NA>
#> 4    pk_hazard                           CL, V2, Q, V3, KA, KA2, VMAX, KM
#> 5  irm1_hazard CL, V2, Q, V3, KA, KA2, KIN, KOUT, IC50, IMAX, n, VMAX, KM
#> 6  irm2_hazard CL, V2, Q, V3, KA, KA2, KIN, KOUT, IC50, IMAX, n, VMAX, KM
#> 7  irm3_hazard CL, V2, Q, V3, KA, KA2, KIN, KOUT, EC50, EMAX, n, VMAX, KM
#> 8  irm4_hazard CL, V2, KA, KA2, Q, V3, KIN, KOUT, EC50, EMAX, VMAX, KM, n
#> 9  tmdd_hazard      KPT, KTP, V2, KA, KA2, KEL, R0, KDEG, KINT, KON, KOFF
#> 10     mspline                                                       <NA>
```
