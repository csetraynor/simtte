## Resubmission

This is a resubmission addressing CRAN reviewer comments.

### Missing \value tags in .Rd files

All exported functions already had `\value` sections. The `pipe.Rd` file
for the re-exported `%>%` operator was missing `\value`; this has been
added.

### Examples wrapped in \donttest{}

All three exported functions (`sim_tte`, `sim_tte_df`,
`explore_pi_tq_surv`) now have fast runnable examples outside
`\donttest{}`:

- `sim_tte()`: small 5-subject Weibull example (< 5 sec).
- `sim_tte_df()`: pure-R mock data example with no mrgsolve dependency
  (< 1 sec); `\donttest{}` wrapper removed entirely.
- `explore_pi_tq_surv()`: small 5-point pi grid, short time horizon
  (< 5 sec).

Longer illustrative examples are retained inside `\donttest{}`.

## Test environments
* local Ubuntu 20.04.6 LTS, R 4.4.2
* local macOS (aarch64, Apple Silicon), R 4.6.0
* win-builder (R-devel)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Downstream dependencies
None (new package).


