## Phase G validation: time-varying lp(t) survival functions against an
## independent analytical reference.
##
## simtte's time-varying-lp(t) implementations compute, inside compiled
## mrgsolve models:
##   - weibull_tv (inst/models/weibull_tv.cpp): a closed-form segment
##     summation, H(t) = sum_j exp(mu + lp_j) * (t_j^shape - t_{j-1}^shape)
##   - ms (inst/models/ms.cpp) with lp_data: h(t) = basehaz(t) * exp(mu +
##     lp(t)), integrated by mrgsolve's ODE solver.
##
## Both hazards are piecewise-constant-in-lp(t) by construction. This
## script derives the cumulative hazard at each requested time by
## stats::integrate() -- genuinely independent numerical quadrature over
## the underlying instantaneous hazard function, not the closed-form
## summation formula used inside the compiled models -- and compares
## S(t) = exp(-H(t)) against simtte's actual output. It does not call
## weibull_tv.cpp's or ms.cpp's formulas directly anywhere.
##
## Reproducibility: deterministic given parameters; no seed required.

suppressMessages(library(simtte))

## ---- Independent reference: piecewise lp(t) via LOCF -------------------

## Returns a function lp(t) that is the last-observation-carried-forward
## step function defined by (knot_times, knot_values); this mirrors the
## *documented contract* (nocb = FALSE), not any package source code.
make_locf_fun <- function(knot_times, knot_values) {
    stopifnot(!is.unsorted(knot_times))
    function(t) {
        idx <- findInterval(t, knot_times)
        idx[idx == 0L] <- 1L
        knot_values[idx]
    }
}

## ---- (a) Weibull + time-varying lp(t): independent quadrature ----------
##
## True instantaneous hazard: h(t) = shape * t^(shape - 1) * exp(mu +
## lp(t)). Cumulative hazard H(t) = integral_0^t h(u) du, computed here
## by stats::integrate() (adaptive quadrature, R's own base numerical
## routine) -- an independent numerical method from the closed-form
## segment summation used inside weibull_tv.cpp.
weibull_tv_H_quadrature <- function(t, mu, shape, lp_fun) {
    if (t <= 0) {
        return(0)
    }
    h <- function(u) shape * u^(shape - 1) * exp(mu + lp_fun(u))
    stats::integrate(h, lower = 0, upper = t, rel.tol = 1e-10)$value
}

run_weibull_tv_grid <- function() {
    cases <- list(
        list(mu = -1, shape = 1.5, knot_t = c(0, 1, 2, 3),
            knot_lp = c(0, 1, -0.5, -0.5), times = c(0.5, 1, 1.5, 2, 2.5, 3)),
        list(mu = -2, shape = 0.7, knot_t = c(0, 0.5, 1.5, 4),
            knot_lp = c(-1, 0.5, 2, -2), times = seq(0.1, 4, by = 0.3)),
        list(mu = 0.3, shape = 3, knot_t = c(0, 0.2, 0.9, 2, 5),
            knot_lp = c(0, 0, 1, -1, -1), times = c(0.05, 0.2, 0.9, 2, 3, 5)),
        list(mu = -0.5, shape = 1, knot_t = c(0, 10),
            knot_lp = c(0.25, 0.25), times = seq(0, 10, by = 1)) # constant lp
    )
    results <- lapply(cases, function(cs) {
        lp_fun <- make_locf_fun(cs$knot_t, cs$knot_lp)
        lp_data <- data.frame(ID = 1, time = cs$knot_t, lp = cs$knot_lp)
        out <- simtte:::.sim_surv_df(log_hr = 0, mu = cs$mu,
            shape = cs$shape, type = "weibull", times = cs$times,
            end_time = max(cs$times), lp_data = lp_data)
        H_ref <- vapply(cs$times, weibull_tv_H_quadrature, numeric(1),
            mu = cs$mu, shape = cs$shape, lp_fun = lp_fun)
        S_ref <- exp(-H_ref)
        got <- setNames(out$p11, as.character(out$time))
        data.frame(mu = cs$mu, shape = cs$shape,
            time = cs$times, p11 = unname(got[as.character(cs$times)]),
            p11_ref_quadrature = S_ref,
            abs_err = abs(unname(got[as.character(cs$times)]) - S_ref))
    })
    dplyr::bind_rows(results)
}

## Consistency check: when lp(t) is constant, weibull_tv must reduce to
## simtte's existing closed-form baseline model (inst/models/weibull.cpp,
## unmodified in this phase) -- a second, package-internal but
## structurally independent code path (a different compiled model).
run_weibull_tv_vs_baseline <- function() {
    grid <- seq(0.1, 8, by = 0.4)
    params <- data.frame(mu = c(-1.5, -0.5, 0.2), shape = c(0.6, 1, 2.2),
        lp = c(-0.8, 0, 1.3))
    out_list <- lapply(seq_len(nrow(params)), function(i) {
        mu <- params$mu[i]; shape <- params$shape[i]; lp <- params$lp[i]
        lp_data <- data.frame(ID = 1, time = c(0, grid), lp = lp)
        out_tv <- simtte:::.sim_surv_df(log_hr = 0, mu = mu, shape = shape,
            type = "weibull", times = grid, end_time = max(grid),
            lp_data = lp_data)
        out_baseline <- simtte:::.sim_surv_df(log_hr = lp, mu = mu,
            shape = shape, type = "weibull", times = grid)
        got <- setNames(out_tv$p11, as.character(out_tv$time))
        ref <- setNames(out_baseline$p11, as.character(out_baseline$time))
        data.frame(mu = mu, shape = shape, lp = lp, time = grid,
            p11_tv = unname(got[as.character(grid)]),
            p11_baseline = unname(ref[as.character(grid)]),
            abs_err = abs(unname(got[as.character(grid)]) -
                unname(ref[as.character(grid)])))
    })
    dplyr::bind_rows(out_list)
}

## ---- (b) M-spline + time-varying lp(t): independent quadrature ---------
##
## True instantaneous hazard: h(t) = basehaz(t) * exp(mu + lp(t)), both
## step functions under LOCF. H(t) computed by stats::integrate() over
## that step-function hazard -- independent of the ODE integration
## performed inside ms.cpp.
ms_tv_H_quadrature <- function(t, mu, basehaz_fun, lp_fun) {
    if (t <= 0) {
        return(0)
    }
    h <- function(u) basehaz_fun(u) * exp(mu + lp_fun(u))
    stats::integrate(h, lower = 0, upper = t, rel.tol = 1e-10,
        subdivisions = 500L)$value
}

run_ms_tv_grid <- function() {
    cases <- list(
        list(mu = 0, bh_t = c(0, 1, 2), bh_v = c(1, 2, 4),
            lp_t = c(0, 1, 2), lp_v = c(0, 1, 1), times = c(0, 1, 2),
            end_time = 2),
        list(mu = -0.3, bh_t = c(0, 1, 2, 3), bh_v = c(0.5, 1.5, 0.5, 0.5),
            lp_t = c(0, 0.5, 1.5, 2.5), lp_v = c(-1, 0.5, 1, -0.5),
            times = c(0, 1, 2, 3), end_time = 3),
        list(mu = 0.2, bh_t = c(0, 2, 4), bh_v = c(2, 0.5, 0.5),
            lp_t = c(0, 4), lp_v = c(0, 0), # constant lp
            times = c(0, 2, 4), end_time = 4)
    )
    results <- lapply(cases, function(cs) {
        bh_fun <- make_locf_fun(cs$bh_t, cs$bh_v)
        lp_fun <- make_locf_fun(cs$lp_t, cs$lp_v)
        basehaz <- matrix(bh_fun(cs$times), ncol = 1)
        lp_data <- data.frame(ID = 1, time = cs$lp_t, lp = cs$lp_v)
        out <- simtte:::.sim_surv_df(log_hr = 0, mu = cs$mu, shape = NULL,
            type = "ms", times = cs$times, basehaz = basehaz,
            end_time = cs$end_time, lp_data = lp_data)
        H_ref <- vapply(cs$times, ms_tv_H_quadrature, numeric(1),
            mu = cs$mu, basehaz_fun = bh_fun, lp_fun = lp_fun)
        S_ref <- exp(-H_ref)
        got <- setNames(out$p11, as.character(out$time))
        data.frame(mu = cs$mu, time = cs$times,
            p11 = unname(got[as.character(cs$times)]),
            p11_ref_quadrature = S_ref,
            abs_err = abs(unname(got[as.character(cs$times)]) - S_ref))
    })
    dplyr::bind_rows(results)
}

## ---- Run + report --------------------------------------------------------

cat("Running weibull_tv vs. independent quadrature reference...\n")
weibull_tv_results <- run_weibull_tv_grid()
cat("\n==== weibull_tv vs. stats::integrate() quadrature ====\n")
print(weibull_tv_results)
cat("Max absolute error:", max(weibull_tv_results$abs_err), "\n")

cat("\nRunning weibull_tv vs. the (unmodified) constant-lp baseline model...\n")
weibull_tv_vs_baseline <- run_weibull_tv_vs_baseline()
cat("\n==== weibull_tv vs. inst/models/weibull.cpp (constant lp) ====\n")
cat("Max absolute error:", max(weibull_tv_vs_baseline$abs_err), "\n")

cat("\nRunning M-spline + lp_data vs. independent quadrature reference...\n")
ms_tv_results <- run_ms_tv_grid()
cat("\n==== ms + lp_data vs. stats::integrate() quadrature ====\n")
print(ms_tv_results)
cat("Max absolute error:", max(ms_tv_results$abs_err), "\n")

cat("\n==== Session info ====\n")
cat("R version:", R.version.string, "\n")
cat("mrgsolve version:", as.character(utils::packageVersion("mrgsolve")),
    "\n")
cat("simtte version:", as.character(utils::packageVersion("simtte")), "\n")
