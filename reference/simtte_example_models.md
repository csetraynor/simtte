# List simtte's bundled example mrgsolve models

Returns the names of the mechanistic mrgsolve models shipped as
user-facing examples/templates in `inst/models/examples/`. These are
plain, editable `.cpp` model source files, not part of the package's
internal Weibull/M-spline simulation engine (see
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md))
and not used by it. They demonstrate the `dxdt_p11 = -p11 * HAZ`
survival co-integration pattern that a user's own custom mrgsolve model
is expected to follow when working with
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
directly; see "Shipped example models" below.

## Usage

``` r
simtte_example_models()
```

## Value

Character vector of example model names, suitable for passing to
[`simtte_example_model`](https://csetraynor.github.io/simtte/reference/simtte_example_model.md).

## Details

Discovery is a directory listing of `*.cpp` files in
`inst/models/examples/`, not a maintained registry or catalog: adding a
new example model in a future package version requires only adding a new
`.cpp` file to that directory, with no change to this function.

## Shipped example models

- `pkpd_idr_hazard`:

  A one-compartment oral-absorption PK model coupled to an
  indirect-response (turnover) mediator model, in which drug exposure
  inhibits mediator production and the event hazard is proportional to
  the current mediator level. Includes between-subject variability
  (`$OMEGA`) on clearance, volume, and potency. This is the model
  described in the package's companion manuscript ("ODE-coupled PK/PD
  hazard"), shipped here as a standalone, loadable file.

- `pkpd_linear_hazard`:

  A minimal one-compartment IV-bolus PK model with a hazard that
  decreases linearly in drug concentration from a constant baseline. No
  between-subject variability or turnover dynamics. Both the PK and
  hazard are linear, so this model has an independent closed-form
  survival function, making it a short, readable starting point for a
  custom model as well as this package's validation anchor for
  mechanistic hazard models (see `inst/validation/` in the package's
  source repository – a developer-only directory, not installed with the
  package).

## See also

[`simtte_example_model`](https://csetraynor.github.io/simtte/reference/simtte_example_model.md),
to load one;
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
the function these models are for.
[`sim_tte_ode_models`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
lists the built-in models for
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
which are a separate set.

## Examples

``` r
simtte_example_models()
#> [1] "pkpd_idr_hazard"    "pkpd_linear_hazard"
```
