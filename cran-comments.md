## Submission

This is an update of simtte from 1.0.2 (on CRAN) to 2.0.0. The maintainer
is unchanged.

2.0.0 adds in-solver time-to-event simulation (`sim_tte_ode()`,
`tte_model()`), right and interval censoring, visit schedules, and
helpers for running simulations in parallel workers. The full list is in
`NEWS.md`.

The major version number reflects changes to existing `sim_tte()`
behaviour. They fix bugs in 1.0.2 and change results or checks:

* `time` now sets the output grid. 1.0.2 ignored its spacing.
* The Weibull survival probability is computed in closed form. 1.0.2
  was off by up to 0.16 for `shape < 1` at small times.
* M-spline hazard values are carried forward over the following
  interval, not backward.
* `...` rejects `tgrid`, `obsonly`, `nocb`, `carry_out` and `data`, and
  M-spline inputs and `end_time` are validated more strictly.

Each is listed under "Breaking changes" in `NEWS.md` with the remedy.
No package on CRAN depends on simtte (see below).

## Test environments

* local macOS (aarch64, Apple Silicon), R 4.6.0

## R CMD check results

`R CMD check --as-cran` on the built source, 2.0.0:

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* simtte compiles `mrgsolve` model code at run time. Compiled models are
  written under `tempdir()` by default (`simtte_model_cache()`). A
  persistent location is used only if the user sets
  `options(simtte.cache_dir = )`.

* Spell checkers may flag "pharmacokinetic" and "pharmacodynamic" in the
  DESCRIPTION. Both are standard domain terms and are spelled correctly.

## Downstream dependencies

None. No package on CRAN lists simtte in Depends, Imports, Suggests or
LinkingTo (checked against the CRAN package database on 2026-11-06).
