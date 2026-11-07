# Joint PK/PD and Time-to-Event Simulation

``` r

library(simtte)
```

## 1. What this is for

**simtte** simulates time-to-event data where the hazard is driven by a
pharmacokinetic/pharmacodynamic (PK/PD) model rather than a fixed
parametric shape. The survival probability is added to the PK/PD model
as one more ODE compartment, `p11`, alongside the model’s own
concentration/response compartments:

``` math
\frac{d\,p_{11}}{dt} = -p_{11}(t) \cdot \mathrm{HAZ}(t)
```

where `HAZ(t)` is any expression built from the model’s own PK/PD
quantities (concentration, response, drug-target complex, …). Each
subject gets one draw $`U \sim \mathrm{Uniform}(0, 1)`$; the ODE solver
is watched *during* integration for the first time at which
$`p_{11}(t) \le U`$, and that is the subject’s event time (or the
subject is censored at the follow-up horizon `end` if no crossing
occurs). This happens inside the solver, not after the fact from a
pre-simulated trajectory – the mechanism
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
implements.

This vignette is the tour: the built-in models, covariates, and
between-subject variability, all through
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md).
Two more vignettes cover the rest of the toolkit in depth –
[`vignette("bring-your-own-model")`](https://csetraynor.github.io/simtte/articles/bring-your-own-model.md)
for converting your own PK/PD model with
[`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md),
and
[`vignette("censoring-and-assessment")`](https://csetraynor.github.io/simtte/articles/censoring-and-assessment.md)
for right/interval censoring and realistic visit schedules. Everything
below is real, executed code, not illustrative snippets.

## 2. Library models

The simplest case has no PK/PD at all: a built-in Weibull hazard.

``` r

sim <- sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 1.5),
  n = 30, end = 20, delta = 2, seed = 1)
#> Building weibull_ode ... done.
head(sim$events)
#>   ID  sim_time sim_status sim_reason
#> 1  1 2.3508265          1      event
#> 2  2 1.9326955          1      event
#> 3  3 1.3186634          1      event
#> 4  4 0.4090209          1      event
#> 5  5 2.6655244          1      event
#> 6  6 0.4393047          1      event
```

`sim$events` is always at least three columns: `ID`, `sim_time`,
`sim_status` (1 = event, 0 = censored at `end`). `"exponential"`,
`"weibull"`, `"gompertz"`, and `"mspline"` (a flexible, continuously
evaluated baseline hazard) are the closed-hazard-shape models;
[`sim_tte_ode_models()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
lists every built-in name, including the six PK/PD-linked ones below,
together with what between-subject variability each supports.

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

Six more models link the hazard to a genuine PK/PD mechanism instead of
a fixed shape. `irm1_hazard` links it to an indirect-response PD
compartment relative to its own pre-dose baseline; dosing is ordinary
`mrgsolve` event data, passed via `data`:

``` r

dose <- data.frame(ID = 1:30, time = 0, cmt = 1, amt = 100, evid = 1)
sim_pkpd <- sim_tte_ode(model = "irm1_hazard",
  param = list(H0 = 0.01, beta_r = 1), n = 30, end = 20, delta = 2,
  data = dose, seed = 1)
#> Building irm1_hazard ... done.
head(sim_pkpd$events)
#>   ID sim_time sim_status     sim_reason
#> 1  1 20.00000          0 administrative
#> 2  2 20.00000          0 administrative
#> 3  3 20.00000          0 administrative
#> 4  4 16.40492          1          event
#> 5  5 20.00000          0 administrative
#> 6  6 18.15148          1          event
```

The other five (`pk_hazard`, concentration-driven;
`irm2_hazard`-`irm4_hazard`, the other indirect-response types;
`tmdd_hazard`, target-mediated drug disposition) work the same way, each
reusing an unedited
[`mrgsolve::modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
PK/PD backbone plus a hazard link – see
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
“PK/PD-linked hazard models” for each one’s exact hazard expression.

## 3. Covariates and the linear predictor

Covariates compose with any model via a formula and a design matrix:
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
builds `X` with
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and
computes `lp = X %*% beta`.

``` r

baseline <- data.frame(age = rnorm(30, 50, 10),
  arm = factor(rep(c("control", "treated"), 15)))
sim_cov <- sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 1.5),
  n = 30, end = 20, delta = 2, covariates = baseline,
  formula = ~ age + arm, beta = c(age = 0.01, armtreated = -0.5),
  seed = 1)
#> Loading model from cache.
#> sim_tte_ode(): dropped the intercept column implied by 'formula' -- the model's own baseline hazard parameter (e.g. 'H0'/'mu') already plays that role, so an intercept in 'lp' would be additively unidentifiable against it. Write 'formula = ~ 0 + ...' to build the design matrix without this message. See ?sim_tte_ode "Covariates and the linear predictor".
head(sim_cov$events)
#>   ID  sim_time sim_status sim_reason
#> 1  1 1.6894386          1      event
#> 2  2 1.9347830          1      event
#> 3  3 0.8871804          1      event
#> 4  4 0.3872332          1      event
#> 5  5 1.8357952          1      event
#> 6  6 0.4132041          1      event
```

`baseline` has no `time` column, so it is a fixed, per-subject covariate
set (equivalent to a single observation at `time = 0`, in effect for the
whole follow-up – a multi-row frame like this one also determines `n`
when it is not supplied); `arm` is coded as an ordinary
treatment-contrast factor, and the intercept `formula` would otherwise
imply is dropped automatically (with a message) since the model’s own
`mu` already plays that role. `formula` accepts anything
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
does – `I(age^2)`, [`log()`](https://rdrr.io/r/base/Log.html),
interactions, [`splines::ns()`](https://rdrr.io/r/splines/ns.html) – as
long as the result is linear in `beta`.

The original interface, one raw covariate column per named coefficient
with no formula at all, still works unchanged:

``` r

cov <- data.frame(time = c(0, 10), sex = c(0, 1))
sim_legacy <- sim_tte_ode(model = "weibull", param = list(mu = -1, shape = 1.5),
  n = 30, end = 20, delta = 2, covariates = cov, beta = c(sex = 0.8),
  seed = 1)
#> Loading model from cache.
```

### Where a covariate ends and a model parameter begins

`formula`/`beta` only ever need to be linear in `beta`, not in the
covariates – a covariate column can already be an arbitrarily nonlinear,
*fixed* transform (`I(age^2)`, or a column computed as
`dose / (EC50 + dose)` for a literal numeric `EC50`). A link that is
nonlinear in a *parameter* you want to estimate or vary – an Emax effect
whose `EC50` is itself unknown – does not fit here: it belongs in the
model’s own hazard expression instead, via
`tte_model(hazard = "H0 * (CP / (EC50 + CP))", params = list(EC50 = 5))`
(see
[`vignette("bring-your-own-model")`](https://csetraynor.github.io/simtte/articles/bring-your-own-model.md)),
where `EC50` composes with PK/PD states, BSV, and dosing directly, and
can be varied like any other model parameter. A *fixed* `EC50` gives the
identical event-time distribution written either way; only the second
way lets `EC50` vary.

### Combining dosing data with covariates

A dosing regimen (`data`) and `covariates` on the same call is the
common case for a real trial simulation, and works exactly as either
does alone: `lp` is combined with the model’s own PK/PD-driven hazard,
carried onto every dosing record from the covariate value in force at
that time (last-observation-carried-forward, the same convention used
everywhere else in this vignette) – no `mrgsolve` warning is emitted.

``` r

dose_regimen <- mrgsolve::expand.ev(ID = 1:6, amt = 100, cmt = 1,
  ii = 12, addl = 3, time = 0)
age <- data.frame(age = c(40, 50, 60, 70, 80, 90))
sim_dose_cov <- sim_tte_ode(model = "pk_hazard",
  param = list(H0 = 0.01, beta_cp = 0.02), n = 6, end = 48, delta = 2,
  data = dose_regimen, covariates = age, formula = ~ age,
  beta = c(age = 0.01), seed = 1, keep_trajectory = TRUE)
#> Building pk_hazard ... done.
#> sim_tte_ode(): dropped the intercept column implied by 'formula' -- the model's own baseline hazard parameter (e.g. 'H0'/'mu') already plays that role, so an intercept in 'lp' would be additively unidentifiable against it. Write 'formula = ~ 0 + ...' to build the design matrix without this message. See ?sim_tte_ode "Covariates and the linear predictor".
head(sim_dose_cov$events)
#>   ID  sim_time sim_status     sim_reason
#> 1  1 48.000000          0 administrative
#> 2  2 48.000000          0 administrative
#> 3  3 28.170561          1          event
#> 4  4  4.492854          1          event
#> 5  5 48.000000          0 administrative
#> 6  6  4.095499          1          event
```

At `time = 0`, plasma concentration `CP` is still `0` for every subject
(the dose has just been given), so `HAZ = H0 * exp(lp + beta_cp * CP)`
reduces to `H0 * exp(lp)` there – directly showing each subject’s own
`age`-driven `lp` already in force at that very first, dosing record:

``` r

t0 <- subset(sim_dose_cov$trajectory, time == 0)
cbind(t0[, c("ID", "CP", "HAZ")], age = age$age)
#>     ID CP        HAZ age
#> 1    1  0 0.01491825  40
#> 26   2  0 0.01648721  50
#> 51   3  0 0.01822119  60
#> 76   4  0 0.02013753  70
#> 101  5  0 0.02225541  80
#> 126  6  0 0.02459603  90
```

## 4. Between-subject variability

Between-subject variability (BSV) on a structural PK/PD parameter (here,
clearance `CL`) is `omega`, a named matrix (any subset of the model’s
own targets – see
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
“Between-subject variability” for the full list per model):

``` r

om <- matrix(0.09, dimnames = list("CL", "CL"))
sim_bsv <- sim_tte_ode(model = "irm1_hazard", omega = om,
  param = list(H0 = 0.01, beta_r = 1), n = 30, end = 20, delta = 2,
  data = dose, seed = 1)
#> Loading model from cache.
head(sim_bsv$events)
#>   ID sim_time sim_status     sim_reason
#> 1  1 20.00000          0 administrative
#> 2  2 20.00000          0 administrative
#> 3  3 20.00000          0 administrative
#> 4  4 16.78841          1          event
#> 5  5 20.00000          0 administrative
#> 6  6 17.93413          1          event
```

`sim_bsv$events` does not echo back each subject’s own drawn `CL` – it
is consumed internally by the ODE solver, not one of the reported
columns. If you want to see the per-subject values, build `idata`
yourself instead of using `omega` (this is exactly the mechanism `omega`
uses internally, exposed as an ordinary R data frame):

``` r

set.seed(1)
n <- 30
CL_i <- 1 * exp(rnorm(n, 0, sqrt(0.09)))  # population CL = 1, log-normal BSV
idata <- data.frame(ID = seq_len(n), CL = CL_i)
head(idata)
#>   ID        CL
#> 1  1 0.8286676
#> 2  2 1.0566389
#> 3  3 0.7782647
#> 4  4 1.6137880
#> 5  5 1.1039033
#> 6  6 0.7818124

sim_idata <- sim_tte_ode(model = "irm1_hazard", idata = idata,
  param = list(H0 = 0.01, beta_r = 1), end = 20, delta = 2, data = dose,
  seed = 1)
#> Loading model from cache.
head(sim_idata$events)
#>   ID sim_time sim_status     sim_reason
#> 1  1 20.00000          0 administrative
#> 2  2 20.00000          0 administrative
#> 3  3 20.00000          0 administrative
#> 4  4 15.76692          1          event
#> 5  5 20.00000          0 administrative
#> 6  6 18.45315          1          event
```

Any `idata` column matching a declared model parameter overrides that
parameter per subject – `U`/`END` (drawn/filled automatically when
absent) work the same way. A user-supplied model with its own declared
`$OMEGA` block uses `omega` the ordinary `mrgsolve` way instead – see
[`vignette("bring-your-own-model")`](https://csetraynor.github.io/simtte/articles/bring-your-own-model.md).

## 5. Beyond this tour

Two more vignettes cover the rest in depth:

- **[`vignette("bring-your-own-model")`](https://csetraynor.github.io/simtte/articles/bring-your-own-model.md)**:
  [`tte_model()`](https://csetraynor.github.io/simtte/reference/tte_model.md)
  converts *any* `mrgsolve` model with an `$ODE`/`$DES` block into one
  usable here, for when the six built-in PK/PD models above do not cover
  your mechanism.
- **[`vignette("censoring-and-assessment")`](https://csetraynor.github.io/simtte/articles/censoring-and-assessment.md)**:
  right censoring (including choosing a rate for a target censoring
  fraction), interval censoring from a visit schedule (including left
  censoring), and realistic visit-schedule generators (jitter, missed
  visits/dropout, and outcome-reactive schedules).

## 6. Accuracy, messages, and when things get slow

`delta` (and `add`) only control how finely the returned event time is
*reported* – the event time itself is located inside the solver, at its
own internal step resolution, not snapped to the requested output grid.

Two situations trigger a
[`message()`](https://rdrr.io/r/base/message.html) (never a change in
what is simulated, just a heads-up):

- **A very steep hazard at a coarse `delta`** can force the solver
  across a region where the hazard changes by orders of magnitude in one
  internal step; when that happens for some subjects,
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  falls back to a coarser-grid event-time estimate for just those
  subjects and says so, naming which of two known causes applies and
  which of `delta` or `rtol`/`atol` fixes it.
- **Weibull `shape < 0.05`**: the ODE form of the Weibull hazard is
  mathematically singular as $`t \to 0`$ for small shapes; below 0.05,
  accuracy degrades near $`t = 0`$ specifically. The simulation still
  runs; `sim_tte(type = "weibull")` (the closed-form path) is exact
  there if that matters more than the in-solver mechanism’s other
  properties.

``` r

msgs <- capture.output(
  sim_small_shape <- sim_tte_ode(model = "weibull",
    param = list(mu = -1, shape = 0.01), n = 5, end = 10, delta = 1,
    seed = 1),
  type = "message"
)
#> Loading model from cache.
#> sim_tte_ode(): Weibull 'shape' = 0.01 is below 0.05; the in-solver mechanism's accuracy near t = 0 degrades below this threshold (see ?sim_tte_ode "Weibull shape support"). The simulation proceeds unchanged; if exactness near t = 0 matters more, use sim_tte(type = "weibull") instead.
msgs
#> character(0)
```

`model = "mspline"` evaluates a quadratic M-spline basis *inside* the
compiled model as a genuinely continuous hazard, at the cost of a fixed
set of shipped interior-knot counts (3, 5, or 7); `sim_tte(type = "ms")`
accepts an arbitrary knot count but treats the hazard as
piecewise-constant on its own reporting grid instead. Use whichever
property you need more.

The M-spline model is also how to **simulate from a Cox model**. A Cox
model leaves the baseline hazard unspecified, so simulating requires
choosing one; `"mspline"` gives a flexible baseline, and covariates
multiply it by `exp(lp)`, so the data are proportional-hazards by
construction and
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
recovers `beta`:

``` r

trt <- data.frame(ID = 1:400, time = 0, trt = rep(0:1, 200))
sim <- sim_tte_ode(model = "mspline", knots = c(5, 10, 15),
  coefs = c(0.2, 1, 0.5, 1.5, 0.8, 0.3), param = list(mu = -1.5),
  n = 400, end = 20, delta = 1, covariates = trt, beta = c(trt = -0.5),
  seed = 1)
#> Building mspline_ode_k3 ... done.
ev <- merge(sim$events, trt[c("ID", "trt")], by = "ID")
fit <- survival::coxph(survival::Surv(sim_time, sim_status) ~ trt, data = ev)
cbind(estimate = coef(fit), confint(fit))
#>       estimate      2.5 %     97.5 %
#> trt -0.6772759 -0.9634755 -0.3910762
```

The true `beta = -0.5` lies inside the 95% interval.

## 7. Backward compatibility

[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
and
[`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
– the original, closed-form/grid-based API – are unchanged:

``` r

set.seed(42)
lp <- matrix(rnorm(20, 0, 0.5), nrow = 20)
result <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
  time = seq(0.1, 20, by = 0.5), type = "weibull", end_time = 20)
#> Building weibull ... done.
head(result)
#> # A tibble: 6 × 4
#>   sim_time sim_status    ID      lp
#>      <dbl>      <dbl> <dbl>   <dbl>
#> 1      1.6          1     1  0.685 
#> 2      3.1          1     2 -0.282 
#> 3      6.6          1     3  0.182 
#> 4      0.1          1     4  0.316 
#> 5      2.1          1     5  0.202 
#> 6      0.6          1     6 -0.0531
```

Reasons to still reach for them instead of
[`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md):

- **An arbitrary M-spline knot count** (`sim_tte(type = "ms")`) –
  `sim_tte_ode(model = "mspline")` ships only 3 fixed variants.
- **Exact Weibull behavior for very small shapes** (`shape < 0.05`, see
  section 6) – the closed-form calculation has no singularity there.
- **An externally simulated `p11` trajectory** (from any source, not
  necessarily `mrgsolve`) –
  [`sim_tte_df()`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
  accepts any validated `(ID, time, p11)` data frame directly;
  [`sim_tte_ode()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
  always runs the ODE solver itself.
