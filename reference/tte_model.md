# Convert a PK/PD mrgsolve model into a \`sim_tte_ode()\`-ready model

Writes the survival scaffold (`p11`, `U`, `END`, the
`TEVT`/`event_found`/`T_PRE`/`P_PRE`/`P_POST` in-solver event-detection
latch –
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
own mechanism) onto a user's own mrgsolve PK/PD model, exactly the
mechanical edit list every shipped `*_hazard.cpp` library model already
documents by hand at the top of its `[PROB]` block – generalized here to
any model with an `$ODE`/`$DES` block, not just the six shipped
backbones. See
[`vignette("bring-your-own-model", package = "simtte")`](https://csetraynor.github.io/simtte/articles/bring-your-own-model.md)
for the full mechanics (block-locating/splicing, the name-collision
policy) and worked examples.

## Usage

``` r
tte_model(model, hazard, params = list(), name = NULL, bsv_targets = NULL)
```

## Arguments

- model:

  One of: a file path to a model file; a character vector of model code
  (either one string with embedded newlines, or already one line per
  element – both accepted); or a compiled `mrgmod` (its own source is
  read back from `mod@code`). Must declare an `$ODE` or `$DES` block – a
  `$PKMODEL`/`$PRED`-only model is rejected, since the survival
  compartment has to be co-integrated in the same ODE system; use the
  ODE-form equivalent from
  [`modlib`](https://mrgsolve.org/docs/reference/modlib.html) instead
  (e.g. `"pk2cmt"`, not `"pk2"`). A `simtte_model` (this function's own
  return class) is rejected – convert the original model again with a
  different `hazard`/`params` instead of converting an already-converted
  one.

- hazard:

  Character scalar: one C++ expression for the instantaneous hazard
  `HAZ`, written in terms of any name already in scope inside the
  model's own `$ODE`/`$DES` (states, `$GLOBAL` `#define`s, `$MAIN`
  locals, `$PARAM` values) plus whatever names `params` adds – e.g.
  `"H0 * exp(lp + beta_cp * CP)"`. Include `lp` in the expression if
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
  `covariates`/ `beta` mechanism should affect this model; `lp` is
  always declared (default 0) whether or not the expression uses it.

- params:

  Named list, the hazard's own parameters and their population-default
  values (e.g. `list(H0 = 0.01, beta_cp = 0)`). Must not name `lp`, `U`,
  or `END` – those three are always added automatically (defaults `0`,
  `0`, `1e9`, matching every shipped `*_hazard.cpp` model) and are
  rejected if present here.

- name:

  Character scalar or `NULL` (default): passed to
  [`mcode_cache`](https://mrgsolve.org/docs/reference/mcode.html) as the
  model name (caching and [`print()`](https://rdrr.io/r/base/print.html)
  only). Left at `NULL`, a name is derived from an MD5 hash of the
  converted code, so two different conversions never share a build
  directory/model name even if called back to back – `mcode_cache()` is
  normally content-hash aware under a repeated name (verified, verified
  directly), but two builds under the same name landing in the same
  filesystem-timestamp second can defeat its underlying `make`-based
  rebuild check (verified directly – exactly the pattern an automated
  test loop hits, converting many small variants back to back); hashing
  the default name sidesteps this rather than relying on
  `mcode_cache()`'s own detection. A caller-supplied `name` does not get
  this treatment (it is used verbatim, for predictable caching across
  sessions) – avoid reusing one `name` for genuinely different `hazard`/
  `params` content called in quick succession.

- bsv_targets:

  Character vector or `NULL` (default): the parameter names
  [`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
  `omega` argument may target via its per-subject-`idata` route. Ignored
  (with a [`message()`](https://rdrr.io/r/base/message.html)) if `model`
  already declares its own `$OMEGA` block – that model already has a
  working BSV route
  ([`omat`](https://mrgsolve.org/docs/reference/omega.html)), and the
  two never combine for the same model.

  The compiled model is built into
  [`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
  (never a global `mrgsolve` option). **Using the result in a separate
  process** (a cluster worker, a chunked/array job task – see
  [`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  "Parallel simulation"): pass this object to
  [`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
  there; it rebuilds from `$code` automatically if the compiled object
  this function returned does not resolve in that process (its shared
  object lives wherever *this* process's own cache put it), so
  correctness never depends on the two processes sharing a filesystem.

## Value

An object of class `"simtte_model"`, a list with `code` (the converted
model source, printable directly with `cat(x$code, sep = "\n")`), `mod`
(the compiled `mrgmod`), `hazard`, `hazard_params` (the fully resolved
parameter list, including the always-added `lp`/ `U`/`END`),
`bsv_targets`, `name`, and `source` (how `model` was supplied: `"code"`,
`"file"`, or `"mrgmod"`). Pass the whole object as
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)'s
`model` argument.

## Details

Call this first, then pass the result to
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
as its `model` argument. In a parallel or multi-process run, run the
result through
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
on each worker.

Nothing in the input model is edited in place except its
`$GLOBAL`/`$MAIN` blocks (both singleton in mrgsolve; every other
addition – a new `$PARAM`, `$INIT`, an `$ODE` block carrying only the
hazard link, a `$CAPTURE` block – is appended as a brand-new block at
the end of the file, never spliced into the model's own `$ODE`/`$DES`
content.

## See also

[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
which takes this function's return value as its `model` argument and
runs the simulation.
[`sim_tte_ode_models`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md),
if a built-in model already fits and no conversion is needed.
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md),
to make the converted model usable in another R process;
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
shows where it is compiled.

## Examples

``` r
# \donttest{
# mrgsolve::modlib() itself needs mrgsolve *attached* (library(mrgsolve)),
# not just installed; mrgsolve::mread() with an explicit 'project' reaches
# the same internal model library without that requirement.
irm1 <- mrgsolve::mread(model = "irm1",
  project = system.file("models", package = "mrgsolve"), compile = FALSE)
tm <- tte_model(irm1,
  hazard = "H0 * exp(lp + beta_r * (RESP / (KIN / KOUT) - 1.0))",
  params = list(H0 = 0.01, beta_r = 0),
  bsv_targets = c("CL", "V2", "KIN", "KOUT"))
#> Building simtte_tte_model_0f9e7be42471 ... 
#> done.
tm
#> <simtte_model> "simtte_tte_model_0f9e7be42471" (source: mrgmod)
#>  hazard: HAZ = H0 * exp(lp + beta_r * (RESP / (KIN / KOUT) - 1.0)) 
#>  parameters:
#>    H0 = 0.01 
#>    beta_r = 0 
#>    lp = 0 
#>    U = 0 
#>    END = 1e+09 
#>  between-subject variability: via idata, targets: CL, V2, KIN, KOUT 
#>  $code: 98 lines -- cat(x$code, sep = "\n") to view
sim <- sim_tte_ode(model = tm, n = 20, end = 30, delta = 5, seed = 1)
head(sim$events)
#>   ID  sim_time sim_status     sim_reason
#> 1  1 30.000000          0 administrative
#> 2  2 30.000000          0 administrative
#> 3  3 30.000000          0 administrative
#> 4  4  9.628203          1          event
#> 5  5 30.000000          0 administrative
#> 6  6 10.715139          1          event
# }
```
