## Summary

This is a patch release (1.0.2) with bug fixes and a package title/description
update. It has not yet been submitted to CRAN.

### Bug fixes

* Fixed time-column position bug in `sim_tte_df()`: time is now resolved by
  column name, not position.
* Preserved original subject IDs in output (no longer renumbered).
* Corrected `.get_tte()` to use the conventional `<= U` inverse-transform
  boundary rule.
* Added optional `time_var` argument to `sim_tte_df()`.

### Title/Description update

Changed the package Title from "Simulate Time-to-Event Data Using Weibull and
Spline Models" to "Simulate Bespoke Time-to-Event Models Using ODEs", and
rewrote the Description to frame the package around bespoke ODE-based
time-to-event models (built-in Weibull/M-spline hazards plus fully custom
`mrgsolve` ODE systems) rather than Weibull and spline models specifically.
No functional changes accompany this rewording.

## Test environments
* local Ubuntu 20.04.6 LTS, R 4.4.2
* local macOS (aarch64, Apple Silicon), R 4.6.0
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.

  This NOTE is due to an outdated HTML Tidy on the local check environment and
  does not reflect an issue with the package. It does not appear on CRAN's check
  machines or on win-builder.

## Downstream dependencies
None (no reverse dependencies on CRAN).
