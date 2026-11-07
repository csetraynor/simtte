# Bring Your Own PK/PD Model

``` r

library(simtte)
```

[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
six built-in PK/PD models
([`vignette("pkpd-time-to-event")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md))
cover common structures;
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
converts *any* `mrgsolve` model with an `$ODE`/`$DES` block into one
usable the same way. This vignette covers the conversion mechanics, how
it interacts with a model’s own between-subject variability, the
boundary between a covariate effect and a model parameter, and the
compile-cost/caching trade-off.

## 1. A minimal conversion

Here is a complete, short model – one-compartment PK with first-order
absorption – with no survival scaffolding of its own:

``` r

my_model <- "
$PARAM CL = 1, V = 20, KA = 1
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/V)*CENT;
$CAPTURE CL KA
"
```

Convert it, naming the hazard as a plain C++ expression in terms of the
model’s own quantities (`CENT`, `V`) and new hazard parameters (`H0`,
`beta`):

``` r

tm <- tte_model(my_model, hazard = "H0 * exp(lp + beta * (CENT/V))",
  params = list(H0 = 0.01, beta = 0.05))
#> Building simtte_tte_model_61f3815d2ede ... done.
tm
#> <simtte_model> "simtte_tte_model_61f3815d2ede" (source: code)
#>  hazard: HAZ = H0 * exp(lp + beta * (CENT/V)) 
#>  parameters:
#>    H0 = 0.01 
#>    beta = 0.05 
#>    lp = 0 
#>    U = 0 
#>    END = 1e+09 
#>  between-subject variability: none configured
#>  $code: 59 lines -- cat(x$code, sep = "\n") to view
```

`tm$code` is the full converted model, ready to inspect:

``` r

cat(tm$code, sep = "\n")
#> 
#> $PARAM CL = 1, V = 20, KA = 1
#> $CMT DEPOT CENT
#> $ODE
#> dxdt_DEPOT = -KA*DEPOT;
#> dxdt_CENT = KA*DEPOT - (CL/V)*CENT;
#> $CAPTURE CL KA
#> 
#> $GLOBAL
#>   // -- BEGIN simtte survival scaffolding (tte_model()) ----
#>   static int    event_found = 0;
#>   static double TEVT        = 0.0;
#>   static double T_PRE       = 0.0;
#>   static double P_PRE       = 1.0;
#>   static double P_POST      = 0.0;
#>   // HAZ is declared here (not as an $ODE-local) so that a
#>   // brand-new, separate $ODE block can still assign it and
#>   // have $CAPTURE see it -- mrgsolve's autodec only hoists
#>   // an $ODE-local into shared storage when it is declared
#>   // in the SAME $ODE occurrence as the rest of the model's
#>   // own equations, verified directly not to happen across
#>   // two separate $ODE (or $ODE + $DES) occurrences (see
#>   // reports/experiments/10_two_ode_blocks_capture_test.R).
#>   static double HAZ         = 0.0;
#>   // -- END simtte survival scaffolding --------------------
#> 
#> $MAIN
#>   // -- BEGIN simtte survival scaffolding (tte_model()) ----
#>   if (NEWIND <= 1) {
#>     event_found = 0;
#>     TEVT = 0.0;
#>     T_PRE = 0.0;
#>     P_PRE = 1.0;
#>     P_POST = 0.0;
#>   }
#>   // -- END simtte survival scaffolding --------------------
#> 
#> $ODE
#>   // -- BEGIN simtte survival scaffolding (tte_model()) ----
#>   HAZ = H0 * exp(lp + beta * (CENT/V));
#>   dxdt_p11 = -p11 * HAZ;
#>   if (!event_found && p11 > U && SOLVERTIME <= END && SOLVERTIME >= T_PRE) {
#>     T_PRE = SOLVERTIME;
#>     P_PRE = p11;
#>   }
#>   if (!event_found && p11 <= U && SOLVERTIME <= END) {
#>     event_found = 1;
#>     TEVT = SOLVERTIME;
#>     P_POST = p11;
#>   }
#>   // -- END simtte survival scaffolding --------------------
#> 
#> $PARAM H0 = 0.01, beta = 0.05, lp = 0, U = 0, END = 1000000000
#> 
#> $INIT
#>   p11 = 1
#> 
#> $CAPTURE
#>   TEVT event_found T_PRE P_PRE P_POST HAZ
```

Pass `tm` to
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
exactly like a built-in model name:

``` r

dose <- data.frame(ID = 1:30, time = 0, cmt = 1, amt = 100, evid = 1)
sim_tm <- sim_tte_ode(model = tm, n = 30, end = 20, delta = 2, data = dose,
  seed = 1)
head(sim_tm$events)
#>   ID  sim_time sim_status     sim_reason
#> 1  1 20.000000          0 administrative
#> 2  2 20.000000          0 administrative
#> 3  3 20.000000          0 administrative
#> 4  4  8.003627          1          event
#> 5  5 20.000000          0 administrative
#> 6  6  8.918469          1          event
```

A bare compiled `mrgsolve` model (not wrapped by
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md))
is not accepted directly by
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
– convert it first. This keeps there to exactly one way to bring a model
in.

## 2. What conversion actually does

[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
writes the same survival scaffold every built-in `*_hazard` model
already has – a `p11` compartment with `dxdt_p11 = -p11 * HAZ`, a
per-subject uniform draw `U`, an administrative horizon `END`, and the
in-solver latch that records the first crossing time – onto your model,
mechanically:

1.  Your model’s own `$ODE`/`$DES` content is required (a
    `$PKMODEL`/`$PRED`-only model is rejected, since the survival
    compartment must be co-integrated in the same ODE system; use an
    ODE-form
    [`mrgsolve::modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
    equivalent instead, e.g. `"pk2cmt"` not `"pk2"`).
2.  `lp`, `U`, and `END` are added as new `$PARAM` entries (defaults 0,
    0, 1e9) – always, whether or not your `hazard` expression uses `lp`.
    Include `lp` in `hazard` if
    [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
    `covariates`/`beta`/`formula` mechanism should affect this model.
3.  Your `hazard` expression becomes `HAZ`, computed inside a brand-new
    `$ODE` block appended after your model’s own – your model’s own
    `$ODE`/`$DES` content is never edited.
4.  `p11 = 1` is added as a new `$INIT` entry, and `p11` as a new
    compartment.
5.  The event-detection latch (`TEVT`, `event_found`, and the bracket
    columns `T_PRE`/`P_PRE`/`P_POST` that
    [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
    refinement step reads) is added to `$GLOBAL`/`$MAIN` – the only two
    blocks edited in place, since `mrgsolve` allows only one of each.
6.  New `$CAPTURE` entries expose the scaffold’s own columns.
7.  A reserved-name collision check runs before compiling: your model
    must not already declare `p11`, `U`, `END`, `HAZ`, `TEVT`,
    `event_found`, `T_PRE`, `P_PRE`, `P_POST`, or `lp` for anything
    else.
8.  The result compiles via
    [`mrgsolve::mcode_cache()`](https://mrgsolve.org/docs/reference/mcode.html)
    into a `"simtte_model"` object, ready for
    `sim_tte_ode(model = ...)`.

Verified to reproduce all six built-in hazard models exactly when
converting their own
[`mrgsolve::modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
backbones directly.

## 3. `$OMEGA` coexistence

If your own model already declares its own `$OMEGA` block (the classic
`TVCL`/`ETA()` idiom),
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
leaves it untouched and
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
`omega` argument goes through it directly instead of the `idata` route
[`vignette("pkpd-time-to-event")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
shows for the built-in models:

``` r

my_model_omega <- "
$PARAM TVCL = 1, TVV = 20, KA = 1
$OMEGA @labels ETA_CL ETA_V
0 0
$CMT DEPOT CENT
$MAIN
double CL = TVCL*exp(ETA_CL);
double V  = TVV*exp(ETA_V);
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/V)*CENT;
"
tm_omega <- tte_model(my_model_omega,
  hazard = "H0 * exp(lp + beta * (CENT/V))",
  params = list(H0 = 0.01, beta = 0.05))
#> Building simtte_tte_model_cc72360c2f85 ... done.
sim_tm_omega <- sim_tte_ode(model = tm_omega, omega = diag(c(0.09, 0.16)),
  n = 30, data = dose, end = 20, delta = 2, seed = 1)
head(sim_tm_omega$events)
#>   ID  sim_time sim_status     sim_reason
#> 1  1 20.000000          0 administrative
#> 2  2 20.000000          0 administrative
#> 3  3 20.000000          0 administrative
#> 4  4  7.101013          1          event
#> 5  5 20.000000          0 administrative
#> 6  6  9.650633          1          event
```

[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)’s
own `bsv_targets` argument is the alternative, for a model with *no*
declared `$OMEGA` block: it names which of the model’s structural
parameters
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
`omega` may target via the `idata` route instead, exactly like the
built-in PK/PD models. Supplying `bsv_targets` for a model that already
declares `$OMEGA` is ignored (with a
[`message()`](https://rdrr.io/r/base/message.html)) – the two routes
never combine for the same model.

## 4. The `lp`/hazard boundary

A covariate’s effect
([`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)’s
`covariates`/`beta`/`formula`) must be linear in `beta`; a link that is
nonlinear in a *parameter* you want to estimate or vary belongs directly
in `hazard` instead, composed with whatever PK/PD state, BSV draw, or
dosing your model already provides:

``` r

tm_emax <- tte_model(my_model, hazard = "H0 * (CENT / (EC50 + CENT))",
  params = list(H0 = 0.01, EC50 = 5))
#> Building simtte_tte_model_fc381db0bb37 ... done.
sim_emax <- sim_tte_ode(model = tm_emax, n = 30, end = 20, delta = 2,
  data = dose, seed = 1)
head(sim_emax$events)
#>   ID sim_time sim_status     sim_reason
#> 1  1 20.00000          0 administrative
#> 2  2 20.00000          0 administrative
#> 3  3 20.00000          0 administrative
#> 4  4 10.38735          1          event
#> 5  5 20.00000          0 administrative
#> 6  6 11.56369          1          event
```

Here `EC50` is a `param` you can override or, in a future variation,
sample per subject via `idata` – something a covariate column can never
do, since a covariate’s effect only ever multiplies a fixed column by a
coefficient. See
[`vignette("pkpd-time-to-event")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
section 3 for the same boundary from the covariate side, including a
demonstration that a *fixed* `EC50` gives identical results whichever
way it is written.

## 5. Compile cost and caching

Every
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
conversion is a real C++ compile (a few seconds).
[`mrgsolve::mcode_cache()`](https://mrgsolve.org/docs/reference/mcode.html)
(used internally) skips recompiling an unchanged, previously-seen model
within the same session. To make that survive across R sessions too
(worthwhile if you reconvert the same model many times, e.g. a
simulation study script re-run repeatedly), set
`options(simtte.cache_dir = tools::R_user_dir("simtte", "cache"))` – see
[`?simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
for the persistence/versioning details, and
[`?simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md)
to remove what accumulates there (e.g. after an R/simtte/mrgsolve
upgrade); avoid pointing `mrgsolve`’s own
`mrgsolve.project`/`mrgsolve.soloc` options at a shared directory
yourself (see the next paragraph for why). One caveat either way: a
caller-supplied `name` reused for genuinely different `hazard`/`params`
content, called back to back, can in principle collide within the same
build directory (leave `name` at its default, which avoids this by
construction, since it is derived from a hash of the converted code).

**Parallel workflow**: don’t have each worker convert/compile the model
independently, and don’t point every worker’s `mrgsolve.soloc` at one
shared directory hoping for a cross-worker cache hit – measured
directly: multiple workers compiling the *same* not-yet-built model into
a shared directory at the same moment fails intermittently, and even
after it is fully built, a worker *re-compiling* it against that shared
directory a second time still can. Convert the model once (in the
parent, or however many times you like sequentially), call
[`simtte_prepare_model()`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
on the result, and broadcast *that* to every worker – see
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
“Parallel simulation” for the full pattern. A worker that only ever
calls `loadso()` on an already-compiled object cannot race with
anything.
