# Clear simtte's compiled-model cache

Removes compiled model artifacts that
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md)/
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md)
built under
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md).
Only ever removes directories *this package itself* created
(version-qualified subdirectories named exactly like
[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)'s
own return value) – never the configured base directory itself, and
never anything under mrgsolve's own default
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) locations that this
package did not create. Safe to call even if `simtte.cache_dir` happens
to be pointed at a directory with other, unrelated content in it (e.g. a
user's home directory): only the recognizably-named version
subdirectories are ever candidates for removal, and each is
double-checked to actually resolve under the configured base before
being deleted.

## Usage

``` r
simtte_model_cache_clear(all = FALSE)
```

## Arguments

- all:

  Logical, default `FALSE`: remove only the *current* version-qualified
  subdirectory – the one
  [`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md)
  would return right now. If `TRUE`, remove every version-qualified
  subdirectory found under the configured base – useful after an R,
  simtte, or mrgsolve upgrade, when a persistent `simtte.cache_dir` has
  accumulated subdirectories for versions no longer in use.

## Value

Character vector of the directory paths actually removed, invisibly
(`character(0)` if there was nothing to remove).

## See also

[`simtte_model_cache`](https://csetraynor.github.io/simtte/reference/simtte_model_cache.md),
which returns the directory this clears.
[`simtte_prepare_model`](https://csetraynor.github.io/simtte/reference/simtte_prepare_model.md),
which rebuilds a model the next time you call it.
[`sim_tte_ode`](https://csetraynor.github.io/simtte/reference/sim_tte_ode.md),
which does the same implicitly.

## Examples

``` r
simtte_model_cache()
#> [1] "/tmp/RtmpRbowND/simtte-model-cache/R4.6.2-simtte2.0.0.9000-mrgsolve2.0.1"
simtte_model_cache_clear()
```
