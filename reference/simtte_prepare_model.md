# Prepare a model for simulation on this process, once

Builds (or, if this process has already built the exact same model
before –
[`mread_cache`](https://mrgsolve.org/docs/reference/mread.html)'s/[`mcode_cache`](https://mrgsolve.org/docs/reference/mcode.html)'s
own in-session cache – reuses) and
[`loadso`](https://mrgsolve.org/docs/reference/loadso.html)s a model,
returning it wrapped in a class
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
recognizes and never rebuilds.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
calls this function internally for every `model` it is given (a
character library name, a
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
`simtte_model`, or an object this function already returned) – there is
exactly one model-loading path, whether or not a caller prepares
anything ahead of time.

## Usage

``` r
simtte_prepare_model(model, knots = NULL)
```

## Arguments

- model:

  A character library-model name (as
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
  own `model` argument), a `simtte_model` from
  [`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md),
  or an object this function already returned.

- knots:

  Only used, and required, when `model = "mspline"`: its *length*
  selects which of the three shipped knot-count variants
  ([`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  "M-spline knot-count variants") gets compiled – the knot *values*
  themselves, like `coefs`/`boundary_knots`, only ever affect `$PARAM`
  values, applied later by
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  itself, never which compiled variant is loaded (checked directly: this
  function's only use of `knots` is `length(knots)`). `coefs`/
  `boundary_knots` are therefore not parameters of this function at all
  – supplying either raises R's own "unused argument" error, naming it,
  rather than silently accepting a value that would not affect what gets
  built.

## Value

An object of class `"simtte_prepared_model"`. Pass it directly as
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
own `model` argument.

## Details

**The point of calling this yourself**: a parallel PSOCK worker is a
separate R process with no memory of what the parent process has already
built or loaded. Calling
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
directly, unprepared, on each worker is not wrong – every path through
it now calls [`loadso`](https://mrgsolve.org/docs/reference/loadso.html)
– but it does mean each worker independently pays its own full compile
cost the first time (harmless for a fast-compiling built-in model;
measured up to ~1.5s for a
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
one, "Parallel simulation" section of
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)).
Calling `simtte_prepare_model()` once per worker at cluster
initialization avoids that redundant work *and* is the pattern verified
race-free: a worker receiving an already-compiled model object (via
ordinary parallel serialization –
[`parallel::clusterCall()`](https://rdrr.io/r/parallel/clusterApply.html)/`clusterExport()`
both work) and only ever calling `loadso()` on it never races, unlike a
worker independently re-triggering `mcode_cache()`/`mread_cache()`
against a directory another process might be writing to at the same
time.

Idempotent: calling this on its own return value (exactly what happens
when the same `model =` argument is broadcast to every worker via
`clusterCall()`, each worker then calling this function on it again)
only re-`loadso()`s (cheap, safe to repeat) and returns the object
unchanged – never a second build.

**A genuinely separate process** (not a parallel cluster worker that
received the compiled object by serialization, but e.g. a
chunked/array-job task reading a
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
object back from a manifest file written by an earlier, already-exited
process): a `simtte_model`'s compiled shared object lives wherever the
*building* process's own cache put it, which a later process cannot see.
Calling this function there rebuilds it from `$code`/`$name`
automatically, into *this* process's own
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
– deterministic (the same code and name always produce the same compiled
model), and never a global `mrgsolve` option, so correctness never
depends on the two processes sharing a filesystem. A built-in library
model needs no such rebuild (it is always (re)built from the shipped
`.cpp` source into this process's own cache directly, every time).

## See also

[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
which takes the returned object as `model` (see "Parallel simulation"
there for the full cluster-init pattern).
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md),
which builds the `simtte_model` to pass in for your own model.
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md),
where the compiled model is stored, and
[`simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md),
to remove it.

## Examples

``` r
# \donttest{
prepared <- simtte_prepare_model("exponential")
#> Building exponential_ode ... 
#> done.
sim <- sim_tte_ode(model = prepared, param = list(H0 = 0.1), n = 10,
  end = 20, delta = 2, seed = 1)
head(sim$events)
#>   ID   sim_time sim_status sim_reason
#> 1  1 13.2610805          1      event
#> 2  2  9.8852799          1      event
#> 3  3  5.5712546          1      event
#> 4  4  0.9628203          1      event
#> 5  5 16.0106348          1      event
#> 6  6  1.0715139          1      event
# }
```
