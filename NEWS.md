# simtte 1.0.2

* Fixed time-column position bug in `sim_tte_df()`: time is now resolved by
  column name, not position.
* Preserved original subject IDs in output (no longer renumbered).
* Corrected `.get_tte()` to use the conventional `<= U` inverse-transform
  boundary rule.
* Added optional `time_var` argument to `sim_tte_df()`.
* Changed package Title to "Simulate Bespoke Time-to-Event Models Using
  ODEs" and rewrote the Description to frame **simtte** around bespoke
  ODE-based time-to-event models (built-in Weibull/M-spline hazards plus
  fully custom `mrgsolve` ODE systems), rather than Weibull and spline
  models specifically.

# simtte 1.0.1

* Fixed issues from CRAN review:
  * Added `\value` tag to `pipe.Rd` documenting the return value of the pipe operator.
  * Reduced example run times and unwrapped short examples from `\donttest{}` so they run during automated checks; examples now execute in under 5 seconds.

# simtte 1.0.0

* Initial CRAN release.
* Simulate time-to-event datasets using Weibull and M-spline baseline hazard models.
* Inverse transform sampling from cumulative hazard functions via mrgsolve ODE solver.
* Exported functions: `sim_tte()`, `sim_tte_df()`, `explore_pi_tq_surv()`.
* Includes example M-spline dataset (`ms_data`).
* Two vignettes: "Getting Started" and "Advanced Simulation Scenarios".
