# Simulate Time-to-Event Data via In-Solver Event Detection

The event/censoring time is located *during* ODE integration (at the
solver's own internal evaluations), by drawing a per-subject \\U \sim
\mathrm{Uniform}(0, 1)\\ and integrating the survival probability `p11`
as a state, latching the first solver time at which `p11 <= U`. This is
a different mechanism from
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)/[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
which resolve event times *after* the fact from an already-materialized
trajectory – the two are separate code paths, not two modes of one
function.

## Usage

``` r
sim_tte_ode(
  model,
  param = list(),
  omega = NULL,
  sigma = NULL,
  n = 1L,
  end = 100,
  delta = 1,
  add = NULL,
  idata = NULL,
  data = NULL,
  covariates = NULL,
  beta = NULL,
  formula = NULL,
  knots = NULL,
  boundary_knots = NULL,
  coefs = NULL,
  censoring = NULL,
  visits = NULL,
  seed = NULL,
  keep_trajectory = FALSE,
  ...
)
```

## Arguments

- model:

  What to simulate. One of: a character name from
  [`sim_tte_ode_models()`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
  (`"exponential"`, `"weibull"`, `"gompertz"`, `"mspline"`,
  `"pk_hazard"`, `"irm1_hazard"`–`"irm4_hazard"`, `"tmdd_hazard"`); a
  `"simtte_model"` returned by
  [`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
  (your own converted mrgsolve model); or the `"simtte_prepared_model"`
  returned by
  [`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md),
  for parallel or multi-process use. A bare compiled `mrgmod` is not
  accepted: convert it with
  [`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
  first.

- param:

  Named list or named numeric vector of population-level `$PARAM`
  overrides, forwarded to
  [`param`](https://mrgsolve.org/docs/reference/param.html). `U` and
  `END` are ordinarily supplied per subject via `idata` (see below), not
  here.

- omega:

  Optional square matrix. For a model that already declares a matching
  `$OMEGA` block (any user-supplied model with its own `ETA()` wiring),
  forwarded to [`omat`](https://mrgsolve.org/docs/reference/omega.html)
  as before. For a built-in PK/PD hazard-library model with no declared
  block, applied instead as per-subject between-subject variability via
  `idata` columns – see "Between-subject variability" below.

- sigma:

  Optional square matrix, forwarded to
  [`smat`](https://mrgsolve.org/docs/reference/sigma.html).
  User-model-only: no built-in library model's hazard depends on a
  residual-error-perturbed quantity, so there is no `idata`-based
  analogue for `sigma`.

- n:

  Integer. Number of subjects. Ignored (and inferred from `nrow(idata)`)
  when `idata` is supplied. Also ignored (and inferred from
  `nrow(covariates)` instead) when `idata` is `NULL` and a baseline,
  no-`ID` `formula` covariate frame (see "Covariates and the linear
  predictor" below) has more than one row – supplying `n` explicitly in
  that case is only valid if it matches that row count exactly (an error
  names both numbers otherwise). When a dosing `data` with an `ID`
  column is also supplied, `n` (or `idata`'s row count) must equal
  `length(unique(data$ID))` exactly – inferred from it when `n`/`idata`
  are both left unresolved, checked otherwise (an error names both
  numbers; see `data` below).

- end:

  Numeric scalar. Administrative censoring horizon. Written into the
  model's `END` parameter for every subject that does not already supply
  its own `END` column in `idata`, and used to guard in-solver event
  detection against latching a time past the requested follow-up (see
  "Boundary guard" below).

- delta:

  Numeric scalar. Reported-grid spacing, used to build the requested
  output grid `seq(0, end, by = delta)` when `add` is not supplied.
  Ignored if `add` is supplied.

- add:

  Numeric vector of explicit reported times. If `NULL` (default),
  `seq(0, end, by = delta)` is used. Either way, `end` itself is always
  included in the resolved grid (via `.resolve_output_grid()`), which
  the censoring rule and refinement step below both depend on.

- idata:

  `NULL` (default) or a data frame with an `ID` column and, optionally,
  `U` and/or `END` columns. Missing `U` is drawn as
  `stats::runif(nrow(idata))` (one independent draw per row/subject; see
  "Reproducibility" below). Missing `END` is filled with `end`. Any
  other column is passed through to mrgsolve as an ordinary per-subject
  covariate/parameter.

- data:

  Optional dosing event data, forwarded to
  [`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html)'s own
  `data` argument. If `covariates` is also supplied, the
  covariate-update rows built from it are appended to `data` (see
  "Covariates and the linear predictor" below); either or both may be
  `NULL`. When supplied together with more than one resolved subject,
  `data` must have an `ID` column whose unique values equal the resolved
  subject set (from `n`/`idata`) exactly – every subject dosed via
  `data` must have exactly one `idata`/ `n` row, and vice versa (an
  error, checked before simulating, names any mismatch and the fix).
  Before this check existed, `n` exceeding `data`'s own `ID` count
  silently dropped the extra subjects from the result with no warning,
  and `n` short of it failed with mrgsolve's own opaque
  `"ID found in the data set, but not in idata"` – both now error
  immediately from here instead.

- covariates:

  `NULL` (default) or a covariate data frame, either baseline (one row
  per subject, no `time` column) or time-varying (a `time` column, an
  optional `ID` column, and one numeric column per covariate) – see
  "Covariates and the linear predictor" below. Required together with
  `beta`; an error if only one of the two is supplied.

- beta:

  `NULL` (default) or a named numeric vector of log-linear coefficients.
  Without `formula`, one name per raw column of `covariates`
  (`lp(t) = sum(beta * covariates(t))`, the original interface,
  unchanged). With `formula`, one name per column of the design matrix
  `formula` implies (`lp(t) = X(t) %*% beta`) – see "Covariates and the
  linear predictor" below.

- formula:

  `NULL` (default, preserving the original `covariates`/`beta` behaviour
  exactly) or a one-sided formula, e.g. `~ age + arm`, built into a
  design matrix via
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) against
  `covariates`. See "Covariates and the linear predictor" below.

- knots, boundary_knots, coefs:

  Only used (and, together with `coefs`, required) when
  `model = "mspline"`; an error if any is supplied for another `model`.
  `knots`: numeric vector of interior knots, strictly increasing,
  strictly inside `boundary_knots` – its length selects which shipped
  knot-count variant is used (see "M-spline knot-count variants" below).
  `boundary_knots`: numeric `c(lower, upper)`, default `c(0, end)`.
  `coefs`: numeric vector of non-negative M-spline coefficients, length
  `length(knots) + 3` (2 = spline degree, fixed).

- censoring:

  `NULL` (default) or a censoring spec list (see
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)).
  When supplied, an independent censoring time is drawn per subject and
  applied to `$events` via
  [`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md),
  inside this call's own seeded block – see "Right censoring" below.

- visits:

  `NULL` (default), a visit schedule (a numeric vector or per-subject
  data frame, see
  [`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)),
  or a jitter spec
  `list(every = <spacing>, jitter = <optional>, jitter_dist = <optional, "uniform"/"normal">, jitter_trunc = <optional>)`.
  When supplied,
  [`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
  is applied to `$events` **after** any `censoring` above – see
  "Interval censoring" below. The jitter-spec form calls
  [`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)
  internally, inside this call's own seeded block
  ([`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
  itself does not accept this third form).

- seed:

  Optional integer. If supplied, `set.seed(seed)` is called before `U`
  is drawn.

- keep_trajectory:

  Logical. If `TRUE`, the full `mrgsim()` output is retained in the
  returned object's `$trajectory`. Default `FALSE` (the
  coarse-output/in-solver-detection design is partly motivated by *not*
  requiring a dense materialized trajectory for large cohorts).

- ...:

  Forwarded to
  [`mrgsim`](https://mrgsolve.org/docs/reference/mrgsim.html). `tgrid`,
  `obsonly`, `nocb`, and `carry_out`/`carry.out` are controlled
  internally to guarantee the contracts documented here and cannot be
  overridden; supplying them raises an error naming the argument (the
  same `.check_reserved_dots()` mechanism
  [`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
  already uses). This is where solver-tolerance arguments relevant to
  in-solver event-time accuracy (`rtol`, `atol`, `hmax`, ...) are set.

## Value

An object of class `"simtte_ode_sim"`, a list with:

- events:

  Data frame with columns `ID`, `sim_time`, `sim_status` (1 = event, 0 =
  censored), `sim_reason` (`"event"`, `"censored"`, or
  `"administrative"` – see "Right censoring" below); plus
  `sim_time_left`/ `sim_time_right` when `visits` is supplied (see
  "Interval censoring" below).

- trajectory:

  The full `mrgsim()` output data frame, or `NULL` unless
  `keep_trajectory = TRUE`.

- model:

  The compiled `mrgmod` object used.

- seed:

  The `seed` argument, as supplied (possibly `NULL`).

- call:

  The matched call, for reproducibility bookkeeping.

## Details

Typical order of use: pick a model (a built-in name from
[`sim_tte_ode_models`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md),
or convert your own mrgsolve model with
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
first), call this function, then optionally post-process the `$events`
it returns with
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
and
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
(or pass `censoring`/`visits` to do both in this call).

The sections below follow one reading order: the mechanism itself
("Boundary guard", "Event-time refinement"), then what varies by `model`
("Weibull shape support", "M-spline knot-count variants", "PK/PD-linked
hazard models"), then the cross-cutting inputs any model can take
("Covariates and the linear predictor", "Between-subject variability",
"Right censoring", "Interval censoring"), and finally "Reproducibility",
which ties every draw above together.

`model %in% c("exponential", "weibull", "gompertz")` are the three
closed-hazard-form library models; `model = "mspline"` is the flexible
baseline hazard (see "M-spline knot-count variants" below);
`"pk_hazard"`, `"irm1_hazard"`–`"irm4_hazard"`, and `"tmdd_hazard"` are
the six PK/PD-linked models (see "PK/PD-linked hazard models" below).
`model` may also be a `simtte_model` produced by
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md),
which converts a user's own mrgsolve PK/PD model into one usable here –
see
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)'s
own documentation, or
[`sim_tte_ode_models`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
to list the built-in names programmatically. A bare compiled `mrgmod` is
not accepted directly; convert it first.

## Boundary guard

the administrative-censoring contract: The shipped library models latch
`TEVT` (the in-solver event time) only while `SOLVERTIME <= END`, so an
internal solver evaluation past the requested follow-up horizon can
never be recorded as an event time. `sim_tte_ode()` additionally
re-applies the boundary rule on the R side, independent of the in-model
guard: a subject is censored (`sim_status = 0`, `sim_time = end`) if no
event was latched at all, **or** if the latched time is `>= end`. This
is deliberately belt-and-braces: the in-model guard is an internal
robustness measure; the R-side rule above is the documented contract.

## Event-time refinement

The raw in-solver `TEVT` is quantized to the ODE solver's own internal
step grid, which is generally finer than, but not independent of, the
requested output grid and solver tolerance. For every subject classified
as an event, `sim_tte_ode()` therefore refines the final `sim_time`
using the already-validated cumulative-hazard interpolation
`.interpolate_log_survival()` (the same method
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md)
offers as `event_time_method = "log_survival"`), applied **between the
solver's own last pre-crossing evaluation and the crossing evaluation
itself** – an interval one internal solver step wide, captured directly
from the model, not between the two *reported* trajectory rows
bracketing `U`. `TEVT` itself therefore only decides **whether** a
subject had an event before `end`; the reported `sim_time` for an event
subject comes from this bracket interpolation, not from the raw `TEVT`
value.

**This makes refinement accuracy independent of `delta`/ `add`**:
because the bracket is a solver step, not a reporting interval, its
width (and hence the constant-hazard-within-bracket approximation error)
does not grow as the requested output grid gets coarser – confirmed
directly across `delta` of 4, 1, and 0.25 for Weibull
`shape = 0.5, 1, 2, 5` and Gompertz `gamma = +/-`, where an earlier,
grid-based refinement approach had instead gotten *worse* than raw
`TEVT` at coarse `delta` for shapes far from 1. `delta`/`add` therefore
now control only the resolution of the returned `$trajectory` (when
`keep_trajectory = TRUE`) and of the reported-grid *fallback* described
next – not the accuracy of a refined event time.

The reported-grid method is kept only as a fallback, used for a subject
when the in-solver bracket is degenerate (its pre-crossing evaluation
coincides with, or is not earlier than, the crossing itself;
`p11(0) = 1 > U` for every valid `U < 1` rules this out at `t = 0`
specifically, but it is not ruled out later in integration – see the
next paragraph) or unavailable (a user-supplied trajectory missing the
`T_PRE`/`P_PRE`/ `P_POST` capture columns). **One
[`message()`](https://rdrr.io/r/base/message.html) is emitted per
`sim_tte_ode()` call whenever any subject falls back**, naming how many
did; a finer `delta`/`add` improves those subjects' event times
specifically (the fallback is exactly the earlier grid-based method, so
its accuracy still depends on the grid). This message never changes the
returned result and is suppressible with
[`suppressMessages()`](https://rdrr.io/r/base/message.html), like every
guardrail message this function emits (see also "Weibull shape support"
below).

**Two fallback-triggering mechanisms have been characterized so far,
with different remedies**, and the message differentiates them
(classified purely from the already-captured bracket columns, no new
model-side logic): a **Weibull-type** signature (`P_POST` outside
`[0, 1]`, a genuine solver overshoot across one large internal step,
remedied by a finer `delta`/`add`), and an **M-spline-type** signature
(`T_PRE == TEVT` exactly, with both `P_PRE` and `P_POST` legitimate
probabilities – a same-timestamp, multiple-corrector-iteration
degeneracy near a locally steep hazard, remedied by tighter
`rtol`/`atol`, **not** `delta`). A subject matching neither signature
(including every subject when the bracket columns are unavailable at
all) is reported as unclassified, with generic advice. Only the non-zero
buckets are named in the message.

**A very steep hazard at a coarse `delta` can trigger this fallback for
an entire cohort at once**: mrgsolve/lsoda does not take an internal
solver step past the next requested output time, so a coarse `delta` can
force one large internal step across a region where the hazard changes
by orders of magnitude, and a trial/corrector evaluation partway through
that step can be momentarily unphysical – detected and safely handled by
the fallback described above, not a silent wrong answer. Measured
directly for Weibull `shape = 5` and `shape = 10`: fallback engagement
drops from 100% of subjects to 0% between `delta = 1` and `delta = 0.5`
at `shape = 10`. The remedy is exactly what the fallback message already
states: request a finer `delta`/`add`. This is a named, disclosed,
non-blocking regime, not a separate guardrail.

## Weibull shape support

The closed-form
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)'s
Weibull model has no restriction on `shape` (an exact analytical
expression is evaluated at every reported time). This function's
`"weibull"` library model genuinely integrates the Weibull hazard as an
ODE state, which diverges as \\t \to 0^+\\ for `shape < 1` – a real
mathematical property of the hazard, not a floating-point artifact
(integrating this exact hazard as an ODE fails outright for small shapes
without a floor). `weibull_ode.cpp` floors the time value used *inside
the hazard expression only* (never `SOLVERTIME` itself) at `1e-100`,
which was verified directly to eliminate solver
failure/`NaN`/non-monotone `p11` entirely across `shape` from 0.01 to
10, with accuracy at or below mrgsolve's own default solver tolerance
for `shape >= 0.05`, degrading gracefully (larger but still finite,
still monotone bias) below that. **This model always runs (no error) for
any `shape > 0`**; for `shape < 0.05` specifically, prefer
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)`(type = "weibull")`
if exactness near `t = 0` matters more than the in-solver mechanism's
other properties. `sim_tte_ode()` additionally emits one
[`message()`](https://rdrr.io/r/base/message.html) per call when the
effective `shape` (after `param` overrides) is below `0.05`, pointing
back to this section; the simulation itself is unaffected (suppressible
with [`suppressMessages()`](https://rdrr.io/r/base/message.html), like
the refinement-fallback message described above).

## M-spline knot-count variants

`model = "mspline"` evaluates a degree-2 (quadratic) normalized M-spline
basis (Ramsay 1988 scaling, matching
`splines2::mSpline(..., degree = 2, intercept = TRUE)` exactly –
verified to floating-point precision) **inside the compiled model, in
C++**, not as a precomputed covariate the way
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)`(type = "ms")`
does: `HAZ = eta * sum_m coefs[m] * M_m(t)`, a genuinely continuous,
smooth hazard, integrated as a real `$ODE` state like every other
`sim_tte_ode()` model (chosen specifically so the in-solver bracket
refinement above, which assumes a smooth hazard, applies to the M-spline
model too – a piecewise-constant covariate does not give that).

Because the basis is compiled in, only a small, fixed set of
interior-knot counts is shipped (**3, 5, or 7** interior knots, a
moderate, fixed-set compromise) – `length(knots)` selects which shipped
variant is used; any other count errors, naming the supported counts and
pointing to
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)`(type = "ms")`,
which accepts an arbitrary knot count (at the cost of being a
piecewise-constant hazard on its own reporting grid, not a continuous
one). `knots` (interior), `boundary_knots` (`c(lower, upper)`, default
`c(0, end)`), and `coefs` (length `length(knots) + 3`, **non-negative**
– the basis itself is non-negative, so a negative coefficient would make
the hazard negative) together define the same hazard
`sim_tte(type = "ms")` would produce from a
[`splines2::mSpline()`](https://wwenjie.org/splines2/reference/mSpline.html)
basis built with the same knots/degree and `coefs` as its own `coefs`
argument – knots and coefficients are directly portable between the two
functions.

No time floor is needed (unlike `"weibull"`'s divergence at `t = 0`):
the M-spline basis is bounded and finite everywhere on
`[boundary_knots[1], boundary_knots[2]]`, including exactly at the
boundary – verified directly, not assumed. The hazard **is** identically
zero outside that interval, though (a consequence of the clamped-knot
spline construction): `end` (and every `idata$END`) must not exceed
`boundary_knots[2]`, enforced as an error, not a silent extrapolation –
the same restriction `sim_tte(type = "ms")`'s own `end_time > max(time)`
check already applies, for the same reason.

## Simulating from a Cox model

A Cox model specifies \\h(t \mid X) = h_0(t) \exp(X\beta)\\ but leaves
the baseline \\h_0(t)\\ unspecified, so simulating from one requires
choosing a baseline. `model = "mspline"` supplies a flexible one: its
hazard is `exp(mu + lp) * sum_m coefs[m] * M_m(t)`, i.e. \\h_0(t) =
e^{\mu} \sum_m c_m M_m(t)\\ with covariates entering multiplicatively
through `lp = X %*% beta` (see "Covariates and the linear predictor"
below). The simulated data therefore satisfy proportional hazards
exactly, and
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
fitted to them estimates `beta`. See
[`vignette("pkpd-time-to-event", package = "simtte")`](https://csetraynor.github.io/simtte/articles/pkpd-time-to-event.md)
for a worked example recovering `beta` this way.

## PK/PD-linked hazard models

Six library models link the hazard to a PK/PD quantity computed by a
genuine mrgsolve PK/PD backbone, each reusing an unedited
[`mrgsolve::modlib()`](https://mrgsolve.org/docs/reference/modlib.html)
model verbatim and adding only the survival scaffold and the hazard link
itself (the per-model mechanical-edit list is documented at the top of
each `.cpp` file, and is exactly what
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
automates for a model of your own):

- `"pk_hazard"` (backbone `"pk2cmt"`): concentration- driven,
  log-linear: `HAZ = H0 * exp(beta_cp * CP)`.

- `"irm1_hazard"`–`"irm4_hazard"` (backbones `"irm1"`–`"irm4"`, the four
  indirect-response types): response-driven, relative to the subject's
  own pre-dose baseline: `HAZ = H0 * exp(beta_r * (RESP / RESP0 - 1))`,
  where `RESP0 = KIN / KOUT`. The sign of `beta_r` gives the direction
  (positive: a response above baseline raises the hazard; negative: it
  lowers it).

- `"tmdd_hazard"` (backbone `"tmdd"`): linked to `RC`, the drug-target
  complex, not to free drug concentration – chosen because `RC` is the
  quantity unique to target-mediated disposition (a `CP`-based link
  would not differ in kind from `"pk_hazard"`);
  `HAZ = H0 * exp(beta_rc * RC)`.

Every one of the six reduces exactly to the exponential model with
hazard `H0` when its `beta_*` parameter is 0. Dosing and other event
data are supplied via this function's existing `data` argument (mrgsolve
event data, e.g. `evid`/`amt`/ `cmt`/`time` rows) exactly as for any
mrgsolve model; `covariates`/`beta` (see "Covariates and the linear
predictor" above) compose with these models unchanged – `lp` is folded
into `eta` the same way the closed-form models do, wherever a given
backbone itself multiplies a hazard-scale quantity by `eta` (verified
directly, not assumed, for the models where `lp` is wired into the
hazard link this way). TMDD is the intentionally stiff member of this
set (fast binding kinetics next to slow target turnover) and is this
library's designated steepness stress case. An outright ODE solver
failure (not the reported-grid fallback described above, an
mrgsolve/`lsoda` error that stops the call) is possible for
`"tmdd_hazard"` at a sufficiently extreme parameter/dose combination –
far past any shipped default – and is not wrapped with a more
informative message: it means the parameters are too extreme for the
requested tolerance; try tighter `rtol`/`atol` or a smaller dose.

## Covariates and the linear predictor

`covariates`/`beta` generalize
[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)'s
single-covariate `lp_data` mechanism to an arbitrary named set: \\lp(t)
= \sum_k \beta_k X_k(t)\\, computed in R and supplied to the model as an
ordinary time-varying `lp` covariate (last-observation-carried-forward,
the same convention
[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
already uses) – no library model file needs to change to support this,
since every one of them already recomputes `eta = exp(mu + lp)` inside
`$ODE`. Coverage follows the same rule as
[`sim_tte()`](https://csetraynor.github.io/simtte/reference/sim_tte.md)'s
Weibull `lp_data` path: every subject's `covariates` must include an
observation at `time = 0`; the last known value is carried forward for
any remaining follow-up (no M-spline-style "must reach `end`"
requirement, since these are genuine ODE models, not a
piecewise-constant-on-a-fixed-grid closed form). **Limitation**: when
`covariates` has an `ID` column, it must equal `1:nrow(idata)` exactly –
subject-specific covariates for a non-default (non-`1:n`) subject
numbering are not yet supported. **Combining with `data`**: when a
dosing `data` is also supplied, every dosing record's `lp` is filled
with the covariate value already in force for that subject at that time
(the same last-observation-carried-forward rule as above), so a dosing
regimen and time-varying covariates compose freely, with no mrgsolve
warning. `data` already carrying its own `lp` column is an error (which
value should apply at a dosing record would otherwise be ambiguous).

**Legacy form** (`formula` not supplied): every name in `beta` must be a
raw numeric column of `covariates`, used as-is – no factor coding,
interactions, or transforms.

**Formula form**: supplying `formula` builds a design matrix `X` via
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(formula, covariates)`,
and `lp(t) = X(t) %*% beta`; `beta`'s names must then match
`colnames(X)` exactly (an error lists the expected names otherwise).
`covariates` may be **baseline** (no `time` column – equivalent to a
single `time = 0` row per subject, carried forward for the whole
follow-up) or **time-varying** (as in the legacy form). A baseline frame
with no `ID` column is either a **population-level constant**, recycled
to every subject, if it has exactly one row, or **one row per subject**
if it has more than one row – in the latter case, `n` (see above) is
inferred from its row count when not supplied explicitly. An implied
intercept column is dropped automatically, with a
[`message()`](https://rdrr.io/r/base/message.html) naming why (the
model's own baseline hazard parameter, e.g. `H0`/`mu`, already plays
that role, so an intercept in `lp` would be additively unidentifiable
against it) – write `formula = ~ 0 + ...` to suppress that message and
get one design-matrix column per factor level instead of the usual
`k - 1` treatment-contrast columns for the *first* factor term in
`formula` (R's own
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) convention
for keeping a no-intercept design full rank; every other factor term
still gets `k - 1` columns). A separate
[`message()`](https://rdrr.io/r/base/message.html) names which term this
happened for, since it changes what `beta`'s expected names are.

Because `X` only ever needs to be linear in `beta`, not in the
covariates themselves, terms nonlinear in a *covariate* are in scope:
`formula = ~ 0 + age + I(age^2) + log(dose) + splines::ns(weight, 3) + arm`,
or a covariate column already computed as an Emax-shaped transform of a
*known, fixed* `EC50` (e.g. a column `dose / (EC50 + dose)` for a
literal numeric `EC50`), all work directly. Terms nonlinear in a
*parameter* the model itself estimates or a caller varies across
scenarios – an Emax link whose `EC50` is itself unknown or varied – are
out of scope for `lp`/`formula`: that belongs in the model's `HAZ`
expression instead
([`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)`(hazard = "H0 * (CP / (EC50 + CP))", params = list(EC50 = 5))`),
where it composes with PK/PD states, BSV, and dosing directly. The two
express the same kind of relationship, but only the second lets `EC50`
itself be estimated or varied; a fixed `EC50` can be written either way,
and gives the same event-time distribution either way (verified
directly).

## Between-subject variability

`omega` adds between-subject variability (BSV) to the **structural PK/PD
parameters** of the six models above (not to `H0`/`lp`/`beta_*`, the
hazard-scale parameters this package's own scaffold adds): a log-normal
draw, \\\theta_i = \theta\_{pop} \cdot \exp(\eta_i)\\, computed in R
with [`mvgauss`](https://mrgsolve.org/docs/reference/mvgauss.html) and
supplied as ordinary `idata` columns – exactly the mechanism already
used for `U` (no model file needs a declared `$OMEGA` block for this to
work). The valid targets are every structural parameter of the model's
own mrgsolve backbone (see each `*_hazard.cpp` file's own `[PROB]` block
for the exact list, e.g. `pk_hazard`: `CL`, `V2`, `Q`, `V3`, `KA`,
`KA2`, `VMAX`, `KM`).

`omega` may be **named** (matching `dimnames`, any subset of the model's
targets, any order – e.g. BSV on `CL` only) or **unnamed** (dimension
must equal the full target list, applied positionally in the order
documented for that model). Population values come from the *resolved*
parameters (`param` overrides applied first), so `param = list(CL = 3)`
together with `omega` on `"CL"` centers the draw on 3. Supplying both an
`idata` column and an `omega` target with the same name is an error (no
silent override in either direction).

**Coexistence with a user-supplied model** (a compiled `mrgmod` with its
own `$OMEGA`/`ETA()` wiring): `omega` is dispatched on whether the
resolved model already declares a matching block. A model that does (any
hand-written model using the standard `TVCL`/`ETA()` idiom) goes through
[`omat`](https://mrgsolve.org/docs/reference/omega.html) exactly as
before – unaffected by this mechanism. Only a model with **no** declared
block falls through to the `idata` route, and only if it is one of the
named models above; `omega` on any other model (the closed-form-style
library models, or a user model with neither) keeps the existing
informative error naming the missing `$OMEGA` block.

`sigma` is unaffected by any of this and stays exactly what it was:
forwarded to [`smat`](https://mrgsolve.org/docs/reference/sigma.html),
usable only for a model (necessarily user-supplied, for now) that
already declares a `$SIGMA` block – no built-in library model's hazard
depends on a residual-error-perturbed quantity.

## Right censoring

`$events` always carries `sim_reason`, distinguishing `"event"` from
administrative censoring at `end` (`"administrative"`). Supplying
`censoring` adds a second, independent censoring process on top – a
per-subject censoring time drawn from the requested distribution, with
the observed outcome becoming `min(event time, censoring time, end)` and
`sim_reason = "censored"` for a subject pulled earlier by the draw. See
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
for the full spec (the same function this argument calls internally);
this is applied in R after simulation rather than as a second ODE
compartment (out of scope for this release: censoring that depends on a
subject's own simulated PK/PD trajectory, e.g. dropout driven by
toxicity).

## Interval censoring

Supplying `visits` maps the exact (and, if `censoring` was also
supplied, already right-censored) outcome onto an interval, via
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md):
an event is only knowable at the first visit at or after it happened,
giving `(L, R]`; a censored subject is only known event-free through
their last visit at or before their own censoring time, giving
`(L, Inf)`. `sim_time`/`sim_status`/`sim_reason` are unchanged –
`sim_time_left`/`sim_time_right` are added on top. Always applied
**after** any `censoring` (see
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
"Ordering with right censoring" for why the order matters); a
`list(every = , jitter = )` value builds the schedule via
[`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)
first, inside this call's own seeded block – a post-simulation R
mapping, not an in-solver mechanism.

A realistic (missed-visit, dropout, or outcome-reactive) schedule is
**not** a `visits =` form this function accepts directly – build it
standalone, then pass the result as `visits`:


    sched <- visit_schedule(n = 50, every = 4, end = 24, jitter = 1, seed = 1)
    realized <- thin_visits(sched, p_miss = 0.1, seed = 2)
    sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
      n = 50, end = 24, delta = 2, visits = realized, seed = 1)

[`visit_schedule_informative`](https://csetraynor.github.io/simtte/reference/visit_schedule_informative.md)
(which needs `$events`, so it cannot run before simulation at all)
composes the same way, applied to `sim$events` after the fact and fed
back through
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
directly rather than through this function's `visits` argument a second
time.

## Reproducibility

If `idata` does not supply `U`, exactly one `stats::runif(1)` is drawn
per row of `idata` (one row per subject), in row order, *before*
mrgsolve is called. With a fixed `seed` (or a
[`set.seed()`](https://rdrr.io/r/base/Random.html) call by the caller),
results are therefore reproducible for a fixed `idata` row order.

The `idata`-based `omega` draw ("Between-subject variability" above) is
likewise reproducible for a fixed `seed` – verified directly (not
assumed from mrgsolve's own documentation) that
[`mvgauss`](https://mrgsolve.org/docs/reference/mvgauss.html)'s draw,
called without its own `seed` argument, is insensitive to any
[`runif()`](https://rdrr.io/r/stats/Uniform.html)/[`rnorm()`](https://rdrr.io/r/stats/Normal.html)
draws made in between (such as the `U` draw above), so one `seed`
deterministically reproduces both regardless of draw order.

## Parallel simulation

`model` is a character library name or a `simtte_model` every time;
either way, `sim_tte_ode()` builds (or reuses an already-built) compiled
mrgsolve model and
[`loadso`](https://mrgsolve.org/docs/reference/loadso.html)s it before
simulating – always, whether or not you prepare anything yourself. On a
single process this needs no attention. On a
[`makePSOCKcluster`](https://rdrr.io/r/parallel/makeCluster.html) worker
– a separate R process, with no memory of what the parent process
already built – each worker independently paying its own first-call
compile cost is not wrong, just wasteful (harmless for a fast-compiling
built-in model; can be a real cost, up to a couple of seconds, for a
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
one).

[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
builds a model once and returns an object `sim_tte_ode()` never
rebuilds; call it once per worker at cluster initialization, then pass
its return value as `model` in every later `sim_tte_ode()` call on that
worker:


    cl <- parallel::makePSOCKcluster(2)
    parallel::clusterCall(cl, simtte::simtte_prepare_model,
      model = "irm1_hazard")
    results <- parallel::parLapply(cl, seeds, function(s) {
      sim_tte_ode(model = "irm1_hazard", param = list(H0 = 0.01, beta_r = 1),
        n = 50, end = 30, delta = 2, seed = s)
    })
    parallel::stopCluster(cl)

For a
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
model, build it *once* in the parent and broadcast the already-compiled
object instead of having every worker convert/compile it independently:


    tm <- tte_model(my_model, hazard = "H0 * exp(lp + beta * CP)",
      params = list(H0 = 0.01, beta = 0))
    prepared <- simtte_prepare_model(tm)
    parallel::clusterCall(cl, simtte::simtte_prepare_model, model = prepared)
    results <- parallel::parLapply(cl, seeds, function(s) {
      sim_tte_ode(model = prepared, n = 50, end = 30, delta = 2, seed = s)
    })

This is also the only pattern verified race-free under repeated,
fresh-cache stress testing: a worker that only ever calls `loadso()` on
an already-compiled object cannot race (a pure read, never a write),
unlike a worker independently re-triggering its own build against a
directory another process might be writing to at the same time – see the
investigation report cited below for the measurements.

For `model = "mspline"`,
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
only takes `knots` (only its *length* matters – it selects which shipped
variant is compiled); `coefs`/`boundary_knots` are supplied to
`sim_tte_ode()` as usual, not to
[`simtte_prepare_model()`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
(see
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)'s
own documentation for why).

**Fork clusters**
([`mclapply`](https://rdrr.io/r/parallel/mclapply.html), Unix-only): a
forked worker inherits the parent process's memory directly, compiled
model included, so none of the above is needed – simulate immediately,
with no
[`simtte_prepare_model()`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)/`clusterCall()`
step at all.

[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
controls where a model is built (default: a directory under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), i.e. session-scoped
and CRAN-safe); see its own documentation for persisting a compiled
model across R sessions, and
[`simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md)
to remove what accumulates there (e.g. after an R/simtte/mrgsolve
upgrade, for a persistent cache).
`sim_tte_ode()`/[`simtte_prepare_model()`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
never call [`options()`](https://rdrr.io/r/base/options.html) themselves
– `mrgsolve.project`/`mrgsolve.soloc` are passed explicitly, never set
globally on your behalf.

## See also

Choosing a model:
[`sim_tte_ode_models`](https://csetraynor.github.io/simtte/reference/sim_tte_ode_models.md)
lists the built-in names;
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
converts your own mrgsolve model into something `model` accepts.
Post-processing the returned `$events`:
[`add_censoring`](https://csetraynor.github.io/simtte/reference/add_censoring.md)
for independent right censoring (what `censoring` calls internally;
[`censoring_rate_for`](https://csetraynor.github.io/simtte/reference/censoring_rate_for.md)
helps choose its rate) and
[`add_interval_censoring`](https://csetraynor.github.io/simtte/reference/add_interval_censoring.md)
for assessment-time intervals (what `visits` calls; build the schedule
with
[`visit_schedule`](https://csetraynor.github.io/simtte/reference/visit_schedule.md)).
Parallel or multi-process runs:
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
builds a model once per process,
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
says where it is stored, and
[`simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md)
removes it (see "Parallel simulation" above). The older grid-based
engine:
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md)
and
[`sim_tte_df`](https://csetraynor.github.io/simtte/reference/sim_tte_df.md),
which this function does not replace.

## Examples

``` r
# \donttest{
sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
  n = 50, end = 30, delta = 5, seed = 1)
#> Loading model from cache.
head(sim$events)
#>   ID   sim_time sim_status sim_reason
#> 1  1 13.2610805          1      event
#> 2  2  9.8852799          1      event
#> 3  3  5.5712546          1      event
#> 4  4  0.9628203          1      event
#> 5  5 16.0106348          1      event
#> 6  6  1.0715140          1      event
# }
```
