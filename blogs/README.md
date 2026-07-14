# blogs/

Long-form, non-package content (blog posts, articles) about `simtte`. This
folder is excluded from the built package via `.Rbuildignore` (`^blogs$`)
and has no effect on CRAN submission, `R CMD build`, or `R CMD check`.

## Posts

1. [Simulating flexible survival hazards with M-splines in `simtte`](msplines-simulation.md) —
   building a custom M-spline baseline hazard and simulating event times from it.
2. [Simulating random (functional) censoring on top of an M-spline hazard](functional-censoring-simulation.md) —
   reuses the M-spline event times from post 1 and layers independent Weibull
   censoring on top, contrasting it with the administrative censoring used
   previously.
3. [Driving a hazard from a PK/PD ODE system in simtte](pkpd-hazard-ode-simulation.md) —
   couples a pharmacokinetic/pharmacodynamic `mrgsolve` model directly to the
   survival probability as an ODE state (mirroring the package's own built-in
   Weibull/M-spline models), then feeds it into `sim_tte_df()`; reuses the
   Weibull censoring approach from post 2.

## Writing a new post

Posts are written as [Quarto](https://quarto.org) `.qmd` files and rendered
to GitHub-flavored markdown so they display natively on GitHub without
needing Quarto installed to *view* them.

1. Create `blogs/<post-name>.qmd` with a YAML header, e.g.:

   ```yaml
   ---
   title: "Post title"
   subtitle: "One-line description"
   author: "Your Name"
   date: "YYYY-MM-DD"
   format: gfm
   ---
   ```

2. Write prose and runnable R code chunks (` ```{r} `). Avoid `eval: false`
   unless a chunk is illustrative only — code in a post should actually run
   against the current version of `simtte`.

3. Render it:

   ```sh
   quarto render blogs/<post-name>.qmd --to gfm
   ```

   This produces `blogs/<post-name>.md` plus a
   `blogs/<post-name>_files/` folder holding any generated plot images.

4. The `.qmd` is the source of truth — commit both the `.qmd` and the
   rendered `.md` (and its `_files/` folder, so images display on GitHub),
   but never hand-edit the `.md` directly; re-render instead.

## What not to commit

Quarto's rendering cache/state (`.quarto/`, `_freeze/`, `*.knit.md`) is
already listed in `.gitignore` — it's local build state, not needed to view
a rendered post, and shouldn't be committed.
