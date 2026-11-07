# Load a bundled simtte example mrgsolve model

Compiles (or loads from cache) one of simtte's shipped example
mechanistic mrgsolve models; see
[`simtte_example_models`](https://csetraynor.github.io/simtte/reference/simtte_example_models.md)
for the list of available names and what each model represents. This is
a thin wrapper around
[`mread_cache`](https://mrgsolve.org/docs/reference/mread.html) pointed
at the package's example-model directory; it adds nothing beyond
discoverability and an informative error for an unrecognized name.

## Usage

``` r
simtte_example_model(name, ...)
```

## Arguments

- name:

  Character scalar naming one bundled model; must be an element of
  [`simtte_example_models()`](https://csetraynor.github.io/simtte/reference/simtte_example_models.md).

- ...:

  Additional arguments passed to
  [`mread_cache`](https://mrgsolve.org/docs/reference/mread.html) (for
  example `quiet = TRUE`).

## Value

A compiled mrgsolve model object.

## Details

The returned model is a standard mrgsolve model object: use
[`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html) to simulate
it, exactly as with any other mrgsolve model. The underlying `.cpp`
source file can be found under
`system.file("models", "examples", package = "simtte")` and copied out
and modified freely with any mrgsolve tooling; loading it through this
function does not require continuing to use this function afterward.

## Intended workflow

These example models are for use with
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
the package's entry point for custom/mechanistic models with an
ODE-integrated survival compartment (see its "Trajectory contract"
section); they are not used by, and cannot be passed to,
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md),
which only ever simulates its own built-in Weibull and M-spline models.


    mod  <- simtte_example_model("pkpd_linear_hazard")
    data <- mrgsolve::ev(amt = 100, cmt = 1, time = 0)
    out  <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = -1,
              add = seq(0, 60, by = 0.5), obsonly = TRUE))
    sim_tte_df(out[, c("ID", "time", "p11")])

## See also

[`simtte_example_models`](https://csetraynor.github.io/simtte/reference/simtte_example_models.md),
to list available names;
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
which turns this model's simulated trajectory into event times. To
simulate a model of your own with
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
convert it with
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md).

## Examples

``` r
# \donttest{
mod <- simtte_example_model("pkpd_linear_hazard")
#> Building pkpd_linear_hazard ... 
#> done.
# }
```
