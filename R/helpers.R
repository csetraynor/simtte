#' Get time at a given row index
#' @param simdat Simulation data frame.
#' @param etime Row index.
#' @return Numeric time value.
#' @noRd
.get_time <- function(simdat, etime) {
    simdat[["time"]][etime]
}

#' Get event time index via inverse transform sampling
#' @param U Uniform random variate.
#' @param pcurr Survival probabilities.
#' @return Integer row index or -99 if no event.
#' @noRd
.get_tte <- function(U, pcurr) {
    match(TRUE, unlist(pcurr) <= U, nomatch = -99L)
}

#' Get maximum time from simulation data
#' @param simdat Simulation data frame.
#' @return Numeric time value.
#' @noRd
.get_max_time <- function(simdat) {
    simdat[["time"]][nrow(simdat)]
}

#' Get path to installed model files
#' @return Character path.
#' @noRd
.cfile_dir <- function() {
    system.file("models", package = "simtte")
}

#' Read and cache mrgsolve model
#' @param model Character model name ("weibull", "weibull_tv", or "ms").
#'   \code{"weibull_tv"} is the time-varying-\code{lp} companion to
#'   \code{"weibull"} (see \code{inst/models/weibull_tv.cpp}); the
#'   constant-\code{lp} \code{"weibull"} model itself is unmodified.
#' @return Compiled mrgsolve model object.
#' @noRd
.read_model_static_cache <- function(model) {
    pkg_model_file <- .cfile_dir()
    if (model %in% c("weibull", "weibull_tv", "ms")) {
        mod_surv <- mrgsolve::mread_cache(model = model,
            project = pkg_model_file)
    } else {
        stop("Model '", model, "' must be 'ms', 'weibull', or ",
            "'weibull_tv'.")
    }
    return(mod_surv)
}

#' Validate a numeric simulation time grid
#'
#' Shared validation for the \code{time} argument of \code{\link{sim_tte}}
#' and related internal functions.
#'
#' @param times Numeric vector of candidate time points.
#' @param arg Character. Name to use for \code{times} in error messages.
#' @return \code{times}, invisibly, if valid.
#' @noRd
.validate_time_grid <- function(times, arg = "time") {
    if (!is.numeric(times)) {
        stop("'", arg, "' must be numeric.", call. = FALSE)
    }
    if (length(times) < 1L) {
        stop("'", arg, "' must contain at least one time point.",
            call. = FALSE)
    }
    if (anyNA(times)) {
        stop("'", arg, "' must not contain NA or NaN values.",
            call. = FALSE)
    }
    if (any(!is.finite(times))) {
        stop("'", arg, "' must not contain Inf or -Inf values.",
            call. = FALSE)
    }
    if (any(times < 0)) {
        stop("'", arg, "' must not contain negative values.", call. = FALSE)
    }
    invisible(times)
}

#' Validate an administrative censoring / follow-up horizon
#'
#' \code{end_time = 0} is accepted: it is a coherent, if degenerate,
#' request for a zero-duration follow-up (output grid is just time 0;
#' since every supported survival model has \eqn{S(0) = 1}, every
#' subject is censored at 0). Negative values are rejected as they
#' would violate the non-negative time contract documented in
#' \code{\link{sim_tte}}.
#'
#' @param end_time Numeric scalar.
#' @return \code{end_time}, invisibly, if valid.
#' @noRd
.validate_end_time <- function(end_time) {
    if (!is.numeric(end_time) || length(end_time) != 1L) {
        stop("'end_time' must be a numeric scalar.", call. = FALSE)
    }
    if (!is.finite(end_time)) {
        stop("'end_time' must be finite (not NA, NaN, Inf, or -Inf).",
            call. = FALSE)
    }
    if (end_time < 0) {
        stop("'end_time' must be non-negative.", call. = FALSE)
    }
    invisible(end_time)
}

#' Resolve the deterministic mrgsolve output time grid
#'
#' Determines the exact set of times at which the survival trajectory is
#' reported by \code{mrgsim()}, and hence at which event/censoring times
#' can be resolved by \code{\link{sim_tte_df}}. This is the single
#' source of truth for the public contract documented in
#' \code{\link{sim_tte}}: the output grid is the sorted, de-duplicated
#' union of the user-supplied \code{time} points that do not exceed
#' \code{end_time}, plus \code{end_time} itself (so that administrative
#' censoring always lands exactly on \code{end_time}).
#'
#' For M-spline models the baseline hazard (\code{basehaz}) is only
#' defined on the supplied \code{time} grid, so \code{end_time} cannot
#' extend the grid beyond \code{max(time)}: doing so would require
#' extrapolating the baseline hazard, which this package does not do.
#' Symmetrically, \code{end_time} cannot fall below \code{min(time)}
#' either: a censoring horizon earlier than the first supplied hazard
#' observation is not scientifically well-defined (there is no hazard
#' information to censor against), so it is also rejected;
#' \code{end_time == min(time)} is the earliest accepted value, and
#' yields a one-point output grid at that time. For Weibull models the
#' hazard has a closed form for all \code{t >= 0}, so \code{end_time}
#' may extend (or truncate) the grid freely, all the way down to 0.
#'
#' @param times Numeric vector. Already-validated by the caller is not
#'   required; this function validates it.
#' @param end_time Numeric scalar. Administrative censoring horizon.
#' @param type Character. \code{"weibull"} or \code{"ms"}.
#' @return A sorted, de-duplicated numeric vector: the output time grid.
#' @noRd
.resolve_output_grid <- function(times, end_time, type) {
    .validate_time_grid(times)
    .validate_end_time(end_time)

    grid <- sort(unique(times))

    if (type == "ms" && end_time > max(grid)) {
        stop("'end_time' (", end_time, ") exceeds max(time) (",
            max(grid), "). The M-spline baseline hazard is only ",
            "defined on the supplied 'time' grid and cannot be ",
            "extrapolated beyond it; supply a 'time' vector that ",
            "spans the desired follow-up horizon instead.",
            call. = FALSE)
    }
    if (type == "ms" && end_time < min(grid)) {
        stop("'end_time' (", end_time, ") is earlier than min(time) (",
            min(grid), "). The M-spline baseline hazard trajectory ",
            "must cover the requested follow-up interval; 'end_time' ",
            "cannot precede the first supplied 'time' point.",
            call. = FALSE)
    }

    grid <- grid[grid <= end_time]
    if (length(grid) == 0L || grid[length(grid)] != end_time) {
        grid <- c(grid, end_time)
    }
    grid
}

#' Numerical tolerance for survival-trajectory validation
#'
#' Used to distinguish genuine violations (a non-increasing survival
#' function that actually increases, or probabilities outside [0, 1])
#' from harmless floating-point / ODE-solver noise. mrgsolve's default
#' solver tolerances are \code{rtol = atol = 1e-8}, so deviations at
#' that scale are expected numerical noise rather than a scientifically
#' invalid trajectory.
#'
#' Policy (Phase pre-B hardening): values within this tolerance of
#' \code{[0, 1]} are accepted and then \strong{clamped} exactly to
#' \code{[0, 1]} by \code{\link{.normalize_survival}} before any
#' monotonicity check or event-time selection uses them (values outside
#' the tolerance remain a hard error, unchanged). Clamping, not
#' rejection, is used because \code{sim_tte_df()}'s primary supported
#' use case is arbitrary user-supplied \pkg{mrgsolve} ODE output, where
#' solver noise at this scale is expected and legitimate; strict
#' rejection would make the package brittle for exactly the custom-model
#' workflow it exists to support. Clamping is safe here specifically
#' because \code{stats::runif()} draws \eqn{U} from \eqn{[0, 1)}: for any
#' raw value \eqn{v} with \eqn{|v| \le} \code{.SURV_TOL} of the boundary,
#' the crossing test \eqn{v \le U} and the clamped test
#' \eqn{\mathrm{clamp}(v) \le U} are provably identical for every
#' achievable \eqn{U} (a value already just below 0 satisfies both
#' unconditionally; a value already just above 1 satisfies neither), so
#' clamping cannot change which trajectories are classified as events
#' vs. censoring. Clamping to exactly \code{[0, 1]} is also what makes
#' the eventual \code{-log(S)} transform needed for Phase B interpolation
#' safe: no accepted, normalized value can make \code{-log(S)} negative
#' or undefined, other than the boundary \eqn{S = 0} case (infinite
#' cumulative hazard), which must be handled explicitly by that future
#' code.
#' @noRd
.SURV_TOL <- 1e-8

#' Clamp validated survival probabilities to the closed interval [0, 1]
#'
#' Deterministic normalization applied once, immediately after the
#' \code{[0, 1] +/- .SURV_TOL} range check passes (see \code{.SURV_TOL}
#' for why this is safe and cannot alter event classification). Every
#' downstream consumer (monotonicity validation, event-time selection,
#' and eventually Phase B's \code{-log(S)} interpolation) operates on
#' these normalized values, never on the raw column.
#'
#' @param surv Numeric vector, already validated to lie within
#'   \code{[-.SURV_TOL, 1 + .SURV_TOL]}.
#' @return \code{surv} clamped elementwise to \code{[0, 1]}.
#' @noRd
.normalize_survival <- function(surv) {
    pmin(pmax(surv, 0), 1)
}

#' Validate a resolved M-spline baseline hazard matrix
#'
#' Validates \code{basehaz} (either computed as \code{basis \%*\% coefs}
#' in \code{\link{sim_tte}}, or supplied directly to
#' \code{\link{explore_pi_tq_surv}}) before it reaches
#' \code{mrgsolve::mrgsim()}. This is the single funnel both entry
#' points pass through, so validating here catches invalid hazards
#' regardless of how \code{basehaz} was produced.
#'
#' @param basehaz Numeric matrix, one row per time point, one column
#'   per candidate baseline hazard curve.
#' @param times Numeric vector, the time points corresponding to the
#'   rows of \code{basehaz} (in the order supplied, not yet sorted).
#' @return \code{basehaz}, invisibly, if valid.
#' @noRd
.validate_basehaz <- function(basehaz, times) {
    if (!is.matrix(basehaz) || !is.numeric(basehaz)) {
        stop("'basehaz' (the M-spline baseline hazard, 'basis %*% ",
            "coefs') must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(basehaz) != length(times)) {
        stop("'basehaz' must have one row per element of 'time' ",
            "(nrow(basehaz) = ", nrow(basehaz), ", length(time) = ",
            length(times), ").", call. = FALSE)
    }
    if (any(!is.finite(basehaz))) {
        stop("'basehaz' must not contain NA, NaN, Inf, or -Inf values.",
            call. = FALSE)
    }
    # Floating-point tolerance for the basis %*% coefs matrix product,
    # not a modeling tolerance: a true-zero hazard can compute to a
    # tiny negative residual (~1e-15) purely from floating-point
    # arithmetic. Anything beyond this scale reflects genuinely
    # negative coefficients, which is a real modeling error (a hazard
    # cannot be negative).
    haz_tol <- 1e-10
    if (any(basehaz < -haz_tol)) {
        stop("The M-spline baseline hazard ('basis %*% coefs' / ",
            "'basehaz') must be non-negative at every supplied time; ",
            "the minimum computed value was ", min(basehaz), ". Check ",
            "the sign of 'coefs'.", call. = FALSE)
    }

    dup_times <- unique(times[duplicated(times)])
    if (length(dup_times)) {
        for (tt in dup_times) {
            idx <- which(times == tt)
            for (col in seq_len(ncol(basehaz))) {
                vals <- basehaz[idx, col]
                if (length(unique(vals)) > 1L) {
                    stop("'time' contains duplicated value ", tt,
                        " with conflicting baseline hazard values (",
                        paste(unique(vals), collapse = ", "), "). The ",
                        "hazard at a given time must be unambiguous; ",
                        "remove or reconcile the conflicting rows. ",
                        "(Duplicated times with identical hazard values ",
                        "are accepted.)", call. = FALSE)
                }
            }
        }
    }
    invisible(basehaz)
}

#' mrgsim() arguments controlled internally by simtte
#'
#' These arguments are set explicitly by \code{.sim_surv_df()} to
#' guarantee the output-grid, trajectory, and hazard-carry contracts
#' documented in \code{\link{sim_tte}} / \code{\link{sim_tte_df}}.
#' Supplying any of them via \code{...} would silently override
#' package-controlled behavior; this was verified empirically for each
#' one:
#' \itemize{
#'   \item \code{tgrid}: replaces the package-resolved deterministic
#'     output grid (the whole point of the \code{time}/\code{end_time}
#'     contract).
#'   \item \code{obsonly}: \code{FALSE} reintroduces the internal
#'     bookkeeping/covariate rows into the output, corrupting the
#'     1:1 correspondence between the output grid and \code{time}.
#'   \item \code{nocb}: controls the M-spline baseline-hazard carry
#'     convention (last-observation-carried-forward vs.
#'     next-observation-carried-backward); simtte fixes this to
#'     \code{FALSE} (see the M-spline hazard-carry contract).
#'   \item \code{carry_out} / \code{carry.out}: replacing this list can
#'     drop columns (\code{lp}, \code{basehaz_id}) that
#'     \code{\link{explore_pi_tq_surv}} depends on.
#'   \item \code{data}: \code{.sim_surv_df()} already bundles the
#'     subject/covariate data into the model object passed as the first
#'     argument; a separate \code{data} silently replaces that data set
#'     entirely.
#' }
#' @noRd
.RESERVED_MRGSIM_ARGS <- c("tgrid", "obsonly", "nocb", "carry_out",
    "carry.out", "data")

#' Reject \code{...} arguments that would override package-controlled
#' \code{mrgsim()} behavior.
#'
#' @param dots A named list, typically \code{list(...)} from the
#'   caller.
#' @return \code{TRUE}, invisibly, if no conflicts are found.
#' @noRd
.check_reserved_dots <- function(dots) {
    conflicts <- intersect(names(dots), .RESERVED_MRGSIM_ARGS)
    if (length(conflicts)) {
        stop("The following argument(s) are controlled internally by ",
            "simtte and cannot be overridden via '...': ",
            paste("'", conflicts, "'", sep = "", collapse = ", "),
            ". simtte sets these explicitly to guarantee its documented ",
            "output-grid and trajectory contract (see ?sim_tte); remove ",
            "them from your call.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate a single subject's survival trajectory
#'
#' Enforces the ordering and monotonicity requirements of the custom
#' trajectory contract documented in \code{\link{sim_tte_df}}: within a
#' subject, times must be unique and strictly increasing (no silent
#' re-sorting), and survival must be non-increasing (no silent repair).
#' Column-level checks (type, finiteness, range) are performed once, up
#' front, in \code{sim_tte_df()} itself, before the data is split by
#' subject.
#'
#' @param id Subject identifier, used only for the error message.
#' @param times Numeric vector of times for this subject, in the order
#'   supplied.
#' @param surv Numeric vector of survival probabilities, same order.
#' @return \code{TRUE}, invisibly, if the trajectory is valid;
#'   otherwise throws an informative error.
#' @noRd
.validate_survival_trajectory <- function(id, times, surv) {
    if (length(times) < 1L) {
        stop("Subject '", id, "' has no trajectory points.", call. = FALSE)
    }

    if (anyDuplicated(times)) {
        dup <- times[anyDuplicated(times)]
        stop("Subject '", id, "': 'time' contains a duplicated value (",
            dup, "). Each subject must have exactly one row per ",
            "distinct time so that the survival probability at that ",
            "time is unambiguous; aggregate or remove duplicate times ",
            "before calling sim_tte_df().", call. = FALSE)
    }
    if (is.unsorted(times)) {
        stop("Subject '", id, "': 'time' must be sorted in ascending ",
            "order. sim_tte_df() treats the reported trajectory as the ",
            "output grid used for event-time resolution and does not ",
            "sort it automatically; sort each subject's rows by time ",
            "before calling sim_tte_df() (e.g. dplyr::arrange(dat, ID, ",
            "time)).", call. = FALSE)
    }
    if (length(surv) >= 2L) {
        increase <- which(diff(surv) > .SURV_TOL)
        if (length(increase)) {
            i <- increase[1]
            stop("Subject '", id, "': survival probability increases ",
                "from ", surv[i], " at time ", times[i], " to ",
                surv[i + 1L], " at time ", times[i + 1L], ". A ",
                "survival function must be non-increasing; ",
                "sim_tte_df() does not repair non-monotone ",
                "trajectories automatically.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

#' Canonicalize a time-varying `lp(t)` input to one row per (subject, time)
#'
#' Validates the user-supplied \code{lp_data} argument of
#' \code{\link{sim_tte}} and expands it to the canonical internal
#' representation: \code{data.frame(ID, time, lp)} covering every
#' subject \code{1:n_subjects} (the same internal numbering
#' \code{sim_tte()} already uses for \code{xdata}, i.e.
#' \code{seq_along(pi)}).
#'
#' Two input shapes are accepted:
#' \itemize{
#'   \item no \code{ID} column: a single population-level trajectory,
#'     recycled identically to every subject;
#'   \item an \code{ID} column present: subject-specific trajectories.
#'     \code{sort(unique(lp_data$ID))} must equal \code{1:n_subjects}
#'     exactly (every subject must have a trajectory; no partial
#'     coverage, no silent recycling of a subset).
#' }
#'
#' No sorting or repair is performed; ordering/duplicate/coverage
#' validation happens in \code{\link{.validate_lp_data_trajectories}}
#' and \code{\link{.check_lp_data_coverage}}, called separately once
#' this function has produced the canonical shape.
#'
#' @param lp_data Data frame as documented for \code{\link{sim_tte}}'s
#'   \code{lp_data} argument.
#' @param n_subjects Integer. Number of subjects (\code{length(pi)}).
#' @return \code{data.frame(ID, time, lp)}, one row per supplied
#'   \code{(subject, time)} pair, \code{ID} always in \code{1:n_subjects}.
#' @noRd
.canonicalize_lp_data <- function(lp_data, n_subjects) {
    lp_data <- as.data.frame(lp_data)
    if (anyDuplicated(names(lp_data))) {
        stop("'lp_data' has duplicated column names.", call. = FALSE)
    }
    for (v in c("time", "lp")) {
        if (!v %in% names(lp_data)) {
            stop("Column '", v, "' not found in 'lp_data'.", call. = FALSE)
        }
    }
    if (!is.numeric(lp_data$time)) {
        stop("Column 'time' in 'lp_data' must be numeric.", call. = FALSE)
    }
    if (!is.numeric(lp_data$lp)) {
        stop("Column 'lp' in 'lp_data' must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(lp_data$time))) {
        stop("Column 'time' in 'lp_data' must not contain NA, NaN, Inf, ",
            "or -Inf values.", call. = FALSE)
    }
    if (any(lp_data$time < 0)) {
        stop("Column 'time' in 'lp_data' must not contain negative ",
            "values.", call. = FALSE)
    }
    # lp is a linear predictor, not a hazard: negative and zero values
    # are valid and expected, only non-finite values are rejected.
    if (any(!is.finite(lp_data$lp))) {
        stop("Column 'lp' in 'lp_data' must not contain NA, NaN, Inf, ",
            "or -Inf values.", call. = FALSE)
    }

    if ("ID" %in% names(lp_data)) {
        if (anyNA(lp_data$ID)) {
            stop("Column 'ID' in 'lp_data' must not contain missing ",
                "values.", call. = FALSE)
        }
        supplied_ids <- sort(unique(lp_data$ID))
        expected_ids <- seq_len(n_subjects)
        if (!identical(as.numeric(supplied_ids), as.numeric(expected_ids))) {
            stop("'lp_data$ID' must contain exactly one trajectory per ",
                "subject, with values 1:", n_subjects, " (matching the ",
                "subject ordering of 'pi'); got IDs: ",
                paste(supplied_ids, collapse = ", "), ".", call. = FALSE)
        }
        out <- data.frame(ID = lp_data$ID, time = lp_data$time,
            lp = lp_data$lp)
    } else {
        out <- do.call(rbind, lapply(seq_len(n_subjects), function(i) {
            data.frame(ID = i, time = lp_data$time, lp = lp_data$lp)
        }))
    }
    out
}

#' Validate each subject's `lp(t)` trajectory (ordering, duplicates)
#'
#' Mirrors \code{\link{.validate_survival_trajectory}}'s ordering
#' philosophy for the new `lp(t)` covariate: no silent sorting, no
#' silent deduplication. Duplicate times are accepted only when the
#' \code{lp} value at that time is identical for every duplicate (same
#' policy as \code{\link{.validate_basehaz}}).
#'
#' @param lp_canonical Output of \code{\link{.canonicalize_lp_data}}.
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.validate_lp_data_trajectories <- function(lp_canonical) {
    by_id <- split(lp_canonical, lp_canonical$ID)
    for (id in names(by_id)) {
        times_i <- by_id[[id]]$time
        lp_i <- by_id[[id]]$lp
        if (is.unsorted(times_i)) {
            stop("'lp_data': 'time' must be sorted in ascending order ",
                "for subject '", id, "'. lp_data is treated as the ",
                "covariate-update grid and is not sorted automatically; ",
                "sort each subject's rows by time before calling ",
                "sim_tte().", call. = FALSE)
        }
        dup <- anyDuplicated(times_i)
        if (dup) {
            dup_time <- times_i[dup]
            vals <- lp_i[times_i == dup_time]
            if (length(unique(vals)) > 1L) {
                stop("'lp_data': subject '", id, "' has duplicated time ",
                    dup_time, " with conflicting lp values (",
                    paste(unique(vals), collapse = ", "), "). The lp ",
                    "value at a given time must be unambiguous.",
                    call. = FALSE)
            }
        }
    }
    invisible(TRUE)
}

#' Check `lp_data` covers the requested follow-up horizon
#'
#' Every subject's trajectory must include an observation at time 0
#' (there is no implicit `lp(0)`, mirroring \code{sim_tte_df()}'s
#' refusal to invent an implicit \eqn{S(0)=1} state -- see
#' \code{?sim_tte_df}). For M-spline models, which cannot extrapolate a
#' time-varying input beyond its own supplied grid (the same reasoning
#' already applied to \code{basehaz} in
#' \code{\link{.resolve_output_grid}}), the trajectory must also reach
#' at least \code{end_time}. For Weibull models this upper bound is not
#' required: the last known \code{lp} value is carried forward (LOCF)
#' for any remaining follow-up, exactly as the constant-lp baseline
#' model already allows \code{end_time} to exceed \code{max(time)}
#' freely.
#'
#' @param lp_canonical Output of \code{\link{.canonicalize_lp_data}}.
#' @param end_time Numeric scalar.
#' @param type Character. \code{"weibull"} or \code{"ms"}.
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.check_lp_data_coverage <- function(lp_canonical, end_time, type) {
    by_id <- split(lp_canonical, lp_canonical$ID)
    for (id in names(by_id)) {
        times_i <- by_id[[id]]$time
        if (!(0 %in% times_i)) {
            stop("'lp_data' must include an observation at time = 0 for ",
                "subject '", id, "'; sim_tte() does not invent an ",
                "implicit lp(0).", call. = FALSE)
        }
        if (type == "ms" && max(times_i) < end_time) {
            stop("'lp_data' must cover the requested follow-up horizon ",
                "for subject '", id, "': max(time) (", max(times_i),
                ") is less than end_time (", end_time, "). The M-spline ",
                "model cannot extrapolate a time-varying 'lp' beyond ",
                "its own supplied grid, the same restriction already ",
                "applied to 'basehaz'; supply 'lp_data' covering the ",
                "full follow-up period.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

#' Last-observation-carried-forward lookup
#'
#' For each element of \code{query_times}, returns the value from the
#' largest \code{known_times} entry that is \code{<= query_time} --
#' i.e. the R-side equivalent of mrgsolve's \code{nocb = FALSE}
#' covariate carry rule, used to pre-fill a second, independently-timed
#' covariate onto a merged grid before both are handed to
#' \code{mrgsolve::data_set()} (verified necessary: mrgsolve requires
#' every input row to carry a value for every covariate column; see
#' PHASE_G_REPORT.md). Query times before \code{min(known_times)} use
#' the first known value (backward-fill), matching \code{mrgsim()}'s own
#' \code{filbak} default and avoiding \code{NA} propagation into the
#' compiled model.
#'
#' @param known_times,known_values Numeric vectors, same length, sorted
#'   ascending by \code{known_times} (not re-sorted here -- callers
#'   already have validated, sorted per-subject trajectories).
#' @param query_times Numeric vector of times to evaluate at.
#' @return Numeric vector, same length as \code{query_times}.
#' @noRd
.locf_at <- function(known_times, known_values, query_times) {
    idx <- findInterval(query_times, known_times)
    idx[idx == 0L] <- 1L
    known_values[idx]
}

# ---------------------------------------------------------------------
# sim_tte_ode() Phase 1 helpers (reports/03_implementation_plan.md
# Phase 1; see reports/02_technical_design.md section 2 for the
# in-solver event-detection mechanism these support).
# ---------------------------------------------------------------------

#' Path to installed sim_tte_ode() library models
#' @return Character path.
#' @noRd
.ode_library_dir <- function() {
    system.file("models", "library", package = "simtte")
}

#' Map a public sim_tte_ode() model name to its library file basename
#'
#' Kept as a tiny, explicit lookup (not a naming convention every future
#' file must follow) so the public, user-facing name
#' (\code{\link{sim_tte_library_models}}) can stay short while each
#' shipped file keeps a descriptive, \code{_ode}-suffixed name
#' distinguishing it from the closed-form engine models
#' (\code{weibull.cpp} etc.) that live one directory up. Phase 5's
#' directory-scan dispatch (\code{03_implementation_plan.md}) is
#' expected to replace this map with plain file discovery once more
#' than a handful of names exist.
#' @noRd
.ODE_LIBRARY_FILES <- c(exponential = "exponential_ode")

#' Load (and cache) a bundled sim_tte_ode() library model
#' @param model Character. One of \code{\link{sim_tte_library_models}}.
#' @return Compiled mrgsolve model object.
#' @noRd
.read_ode_library_model <- function(model) {
    file <- unname(.ODE_LIBRARY_FILES[model])
    mrgsolve::mread_cache(model = file, project = .ode_library_dir())
}

#' Validate the sim_tte_ode() model contract
#'
#' A model usable by \code{\link{sim_tte_ode}} must expose a \code{p11}
#' survival compartment and \code{U}/\code{END} parameters (populated
#' either by a \code{$PARAM} default, \code{param =}, or per-subject
#' \code{idata}; all three ultimately require the name to already be a
#' declared model parameter, which is what is checked here). This is an
#' explicit-contract check, not auto-injection (design report section
#' 5.2): a model failing this check is rejected with an informative
#' error, never silently patched. Reused unchanged from Phase 1 onward,
#' including once a user-supplied \code{mrgmod} is accepted (Phase 5).
#'
#' @param mod A compiled mrgsolve model object.
#' @return \code{TRUE}, invisibly, if the contract is satisfied.
#' @noRd
.validate_ode_model_contract <- function(mod) {
    cmts <- names(mrgsolve::init(mod))
    if (!"p11" %in% cmts) {
        stop("This mrgsolve model has no 'p11' compartment. ",
            "sim_tte_ode() requires a survival compartment named 'p11' ",
            "with dxdt_p11 = -p11 * HAZ; see ?sim_tte_ode.", call. = FALSE)
    }
    pars <- names(mrgsolve::param(mod))
    missing_pars <- setdiff(c("U", "END"), pars)
    if (length(missing_pars)) {
        stop("This mrgsolve model is missing required parameter(s): ",
            paste(missing_pars, collapse = ", "), ". sim_tte_ode() ",
            "requires 'U' (the per-subject uniform draw) and 'END' ",
            "(the administrative censoring horizon) to be declared ",
            "$PARAM entries; see ?sim_tte_ode.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Build/validate the idata table for sim_tte_ode()
#'
#' Canonicalizes \code{idata} to a data frame with \code{ID}, \code{U},
#' and \code{END} columns. Missing \code{U} is drawn as one independent
#' \code{stats::runif()} value per row (i.e. per subject, since
#' \code{idata} is one row per \code{ID} by construction) -- "keyed by
#' ID" in the sense that each row's draw is permanently associated with
#' that row's \code{ID} from the moment it is drawn, never reassigned by
#' a later reorder (see the design report section 2.6 and \code{?sim_tte_ode}
#' "Reproducibility"). Missing \code{END} is filled with \code{end}.
#'
#' @param idata \code{NULL} or a data frame with an \code{ID} column.
#' @param n Integer. Number of subjects, used only when \code{idata} is
#'   \code{NULL}.
#' @param end Numeric scalar. Fallback value for a missing \code{END}
#'   column.
#' @return A data frame with columns \code{ID}, \code{U}, \code{END}
#'   (plus any other columns the caller supplied).
#' @noRd
.build_ode_idata <- function(idata, n, end) {
    if (is.null(idata)) {
        if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
            n < 1 || n != round(n)) {
            stop("'n' must be a positive integer scalar when 'idata' is ",
                "not supplied.", call. = FALSE)
        }
        idata <- data.frame(ID = seq_len(n))
    } else {
        idata <- as.data.frame(idata)
        if (!"ID" %in% names(idata)) {
            stop("'idata' must contain an 'ID' column.", call. = FALSE)
        }
        if (anyDuplicated(idata$ID)) {
            stop("'idata' must not contain duplicated 'ID' values.",
                call. = FALSE)
        }
    }
    if (!"U" %in% names(idata)) {
        idata$U <- stats::runif(nrow(idata))
    } else if (any(!is.finite(idata$U)) || any(idata$U < 0) ||
        any(idata$U >= 1)) {
        stop("'idata$U' must contain finite values in [0, 1).",
            call. = FALSE)
    }
    if (!"END" %in% names(idata)) {
        idata$END <- end
    }
    idata
}

#' Resolve sim_tte_ode() event/censoring times from a simulated trajectory
#'
#' Applies the censoring rule and refinement step documented in
#' \code{?sim_tte_ode} ("Boundary guard" / "Event-time refinement") to
#' the \code{mrgsim()} output of a library model following the
#' \code{p11}/\code{TEVT}/\code{event_found}/\code{U}/\code{END}
#' contract (\code{.validate_ode_model_contract()}). Reuses
#' \code{\link{.normalize_survival}}, \code{\link{.get_tte}}, and
#' \code{\link{.interpolate_log_survival}} unchanged rather than
#' reimplementing crossing detection/interpolation -- this is the same
#' machinery \code{\link{sim_tte_df}} already uses for
#' \code{event_time_method = "log_survival"}.
#'
#' @param traj Data frame: \code{mrgsim()} output with columns \code{ID},
#'   \code{time}, \code{p11}, \code{TEVT}, \code{event_found}, \code{U},
#'   \code{END}, one or more rows per subject in ascending \code{time}
#'   order (mrgsolve's own reporting order; re-sorted defensively below).
#' @return Data frame with columns \code{ID}, \code{sim_time},
#'   \code{sim_status}.
#' @noRd
.resolve_ode_events <- function(traj) {
    ids <- unique(traj$ID)
    rows <- lapply(ids, function(id) {
        sub <- traj[traj$ID == id, , drop = FALSE]
        sub <- sub[order(sub$time), ]
        last <- sub[nrow(sub), ]
        u_i <- last$U
        end_i <- last$END

        # Censoring rule (?sim_tte_ode "Boundary guard"): no crossing
        # latched, or the latched time is >= end -- checked here on the
        # R side regardless of the in-model SOLVERTIME <= END guard.
        if (!isTRUE(as.logical(last$event_found)) || last$TEVT >= end_i) {
            return(data.frame(ID = id, sim_time = end_i, sim_status = 0L))
        }

        p11 <- .normalize_survival(sub$p11)
        etime_idx <- .get_tte(u_i, p11)
        if (etime_idx == -99L) {
            # Defensive fallback, not expected to be reached: `end` is
            # always a reported row (.resolve_output_grid()) and p11 is
            # monotonically non-increasing, so a latched pre-`end` event
            # implies some reported row has p11 <= U by the time `end`
            # is reached. Falls back to the raw (already < end_i, by the
            # check above) in-solver estimate.
            sim_time <- last$TEVT
        } else if (etime_idx == 1L) {
            # Crossing already present at the first reported observation
            # -- no earlier point to interpolate from (identical rule to
            # .simulate_survival_id()).
            sim_time <- sub$time[1]
        } else {
            i <- etime_idx - 1L
            sim_time <- .interpolate_log_survival(t_i = sub$time[i],
                t_ip1 = sub$time[etime_idx], s_i = p11[i],
                s_ip1 = p11[etime_idx], u = u_i)
        }
        data.frame(ID = id, sim_time = sim_time, sim_status = 1L)
    })
    dplyr::bind_rows(rows)
}
