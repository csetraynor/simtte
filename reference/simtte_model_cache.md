# Where simtte builds and caches compiled models

Resolves and returns the directory
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
and
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)
pass as mrgsolve's own `soloc` argument (where a model's compiled shared
object is built and stored) – simtte itself never sets
[`options()`](https://rdrr.io/r/base/options.html) on a user's behalf,
for `mrgsolve`'s own options or this package's `simtte.cache_dir`; this
function only ever *reads* `getOption("simtte.cache_dir")` and appends a
version-qualified subdirectory before returning it.

## Usage

``` r
simtte_model_cache()
```

## Value

Character scalar: an existing, writable directory path.

## Details

**Default** (`simtte.cache_dir` unset): a directory under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), i.e. session-scoped
and automatically cleaned up – CRAN-safe (nothing is ever written
outside a session temp directory unless you opt in below) and the right
choice for any single script or session, including a local parallel
cluster (every worker of a *local* cluster shares the parent's
filesystem, so a model built once under the parent's own cache directory
is directly reachable by every worker's
[`loadso`](https://mrgsolve.org/docs/reference/loadso.html) call – see
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
"Parallel simulation").

**Persistent cache**: set
`options(simtte.cache_dir = tools::R_user_dir("simtte", "cache"))` (the
standard, CRAN-compliant location for a package's own cache data –
`simtte_model_cache` never chooses this for you) to keep a compiled
model's shared object across R restarts, avoiding a repeat compile in a
later session. This matters most for a slow-to-compile
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)-converted
model reused across many scripts, or a genuinely distributed cluster
whose worker nodes do *not* share the parent's filesystem (where every
node needs to be pointed, itself, at the same shared, network- visible
directory).

**Version-key caveat**: whichever directory `simtte.cache_dir` resolves
to, a subdirectory keyed by the current R version, this simtte version,
and the installed mrgsolve version is always appended before it is
returned or created. This is deliberate, not incidental: mrgsolve's own
build-cache key is content/mtime based, not version-aware (verified: its
own `?mread` documentation states plainly that `soloc` is "a directory
location where the model shared object is built and stored", with no
mention of invalidating it across an R or package upgrade) – a
*persistent* cache surviving an R version change or a simtte/mrgsolve
upgrade risks reusing a compiled `.so` built against a different,
ABI-incompatible R/package version. The version key makes that class of
staleness unreachable by construction rather than relying on mrgsolve's
own file-level staleness detection to catch it (which can be defeated by
two builds landing in the same filesystem-timestamp second – see
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)'s
own `name` documentation). The **default**,
[`tempdir()`](https://rdrr.io/r/base/tempfile.html)-based location needs
no such protection on its own (a session cannot upgrade R or a package
mid-run), but the key is applied unconditionally anyway, so there is
exactly one code path to reason about, not two.

## See also

[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md),
which builds into this directory (as mrgsolve's `soloc`).
[`simtte_model_cache_clear`](https://csetraynor.github.io/simtte/reference/simtte_model_cache_clear.md),
to remove what accumulates here.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
and
[`tte_model`](https://csetraynor.github.io/simtte/reference/tte_model.md)
also build into it; see "Parallel simulation" in
[`?sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md).

## Examples

``` r
simtte_model_cache()
#> [1] "/tmp/RtmpRbowND/simtte-model-cache/R4.6.2-simtte2.0.0.9000-mrgsolve2.0.1"
```
