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
#'
#' @param subdir Character or \code{NULL}. An optional subdirectory of
#'   \code{inst/models/} (e.g. \code{"library"} for
#'   \code{inst/models/library/}, the \code{\link{sim_tte_ode}} model
#'   library). \code{NULL} (default) returns \code{inst/models/} itself,
#'   the closed-enum internal engine directory used by
#'   \code{\link{.read_model_static_cache}} -- unchanged from every call
#'   site that existed before this parameter was added.
#' @return Character path.
#' @noRd
.cfile_dir <- function(subdir = NULL) {
    if (is.null(subdir)) {
        system.file("models", package = "simtte")
    } else {
        system.file("models", subdir, package = "simtte")
    }
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
            "'weibull_tv'.", call. = FALSE)
    }
    mod_surv
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
#' @param value_cols Character vector of covariate-value column names to
#'   validate/carry through. Default \code{"lp"} reproduces every prior
#'   (Phase G) call site and its behavior exactly: a single covariate
#'   named \code{lp}. \code{\link{sim_tte_ode}}'s multi-covariate
#'   support (Phase 2, \code{03_implementation_plan.md}) is the only
#'   caller that supplies more than one name, generalizing this from
#'   "time-varying lp(t)" to "time-varying named covariates" without a
#'   second, duplicated implementation.
#' @return \code{data.frame(ID, time, <value_cols>)}, one row per
#'   supplied \code{(subject, time)} pair, \code{ID} always in
#'   \code{1:n_subjects}.
#' @noRd
.canonicalize_lp_data <- function(lp_data, n_subjects, value_cols = "lp") {
    lp_data <- as.data.frame(lp_data)
    if (anyDuplicated(names(lp_data))) {
        stop("'lp_data' has duplicated column names.", call. = FALSE)
    }
    for (v in c("time", value_cols)) {
        if (!v %in% names(lp_data)) {
            stop("Column '", v, "' not found in 'lp_data'.", call. = FALSE)
        }
    }
    if (!is.numeric(lp_data$time)) {
        stop("Column 'time' in 'lp_data' must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(lp_data$time))) {
        stop("Column 'time' in 'lp_data' must not contain NA, NaN, Inf, ",
            "or -Inf values.", call. = FALSE)
    }
    if (any(lp_data$time < 0)) {
        stop("Column 'time' in 'lp_data' must not contain negative ",
            "values.", call. = FALSE)
    }
    for (v in value_cols) {
        if (!is.numeric(lp_data[[v]])) {
            stop("Column '", v, "' in 'lp_data' must be numeric.",
                call. = FALSE)
        }
        # A linear predictor/covariate, not a hazard: negative and zero
        # values are valid and expected, only non-finite values rejected.
        if (any(!is.finite(lp_data[[v]]))) {
            stop("Column '", v, "' in 'lp_data' must not contain NA, ",
                "NaN, Inf, or -Inf values.", call. = FALSE)
        }
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
        out <- data.frame(ID = lp_data$ID, time = lp_data$time)
    } else {
        out <- do.call(rbind, lapply(seq_len(n_subjects), function(i) {
            data.frame(ID = i, time = lp_data$time)
        }))
    }
    n_rep <- nrow(out) / nrow(lp_data)
    for (v in value_cols) {
        out[[v]] <- rep(lp_data[[v]], n_rep)
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
#' @param value_cols Character vector of covariate-value column names to
#'   check for time-conflicts. Default \code{"lp"} (see
#'   \code{\link{.canonicalize_lp_data}}'s \code{value_cols} for why
#'   this default reproduces every prior call site unchanged).
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.validate_lp_data_trajectories <- function(lp_canonical, value_cols = "lp") {
    by_id <- split(lp_canonical, lp_canonical$ID)
    for (id in names(by_id)) {
        times_i <- by_id[[id]]$time
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
            for (v in value_cols) {
                vals <- by_id[[id]][[v]][times_i == dup_time]
                if (length(unique(vals)) > 1L) {
                    stop("'lp_data': subject '", id, "' has duplicated ",
                        "time ", dup_time, " with conflicting ", v,
                        " values (", paste(unique(vals), collapse = ", "),
                        "). The ", v, " value at a given time must be ",
                        "unambiguous.", call. = FALSE)
                }
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
#' \code{?sim_tte_df}). For models that cannot extrapolate a
#' time-varying input beyond its own supplied grid (the same reasoning
#' already applied to \code{basehaz} in
#' \code{\link{.resolve_output_grid}}), the trajectory must also reach
#' at least \code{end_time}. Otherwise this upper bound is not
#' required: the last known \code{lp}/covariate value is carried
#' forward (LOCF) for any remaining follow-up, exactly as the
#' constant-lp baseline model already allows \code{end_time} to exceed
#' \code{max(time)} freely.
#'
#' \code{strict_coverage} names the actual behaviour directly (does this
#' model require the trajectory to reach \code{end_time}, yes/no)
#' instead of routing it through a \code{sim_tte()} model-type string,
#' which is what this parameter replaced
#' (\code{reports/06_phase2_report.md} section 8 /
#' \code{reports/04_author_decisions.md} "After the test runbook /
#' Phase 2.5"): the M-spline baseline hazard is the only case requiring
#' strict coverage, but that reason is about *why* M-spline needs it,
#' not a fact this function should re-derive from a type name every
#' future caller has to know the meaning of.
#'
#' @param lp_canonical Output of \code{\link{.canonicalize_lp_data}}.
#' @param end_time Numeric scalar.
#' @param strict_coverage Logical. \code{TRUE} if the trajectory must
#'   reach \code{end_time} (M-spline models); \code{FALSE} if the last
#'   known value may be carried forward instead (Weibull models,
#'   \code{sim_tte_ode()}'s genuine-ODE library models).
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.check_lp_data_coverage <- function(lp_canonical, end_time,
    strict_coverage) {
    by_id <- split(lp_canonical, lp_canonical$ID)
    for (id in names(by_id)) {
        times_i <- by_id[[id]]$time
        if (!(0 %in% times_i)) {
            stop("'lp_data' must include an observation at time = 0 for ",
                "subject '", id, "'; sim_tte() does not invent an ",
                "implicit lp(0).", call. = FALSE)
        }
        if (strict_coverage && max(times_i) < end_time) {
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

#' Map a public sim_tte_ode() model name to its library file basename
#'
#' Kept as a tiny, explicit lookup (not a naming convention every future
#' file must follow) so the public, user-facing model name can stay
#' short while each shipped file keeps a descriptive, \code{_ode}-suffixed
#' name distinguishing it from the closed-form engine models
#' (\code{weibull.cpp} etc.) that live one directory up. Phase 5's
#' directory-scan dispatch (\code{03_implementation_plan.md}) is
#' expected to replace this map with plain file discovery once more
#' than a handful of names exist.
#' @noRd
.ODE_LIBRARY_FILES <- c(exponential = "exponential_ode",
    weibull = "weibull_ode", gompertz = "gompertz_ode",
    pk_hazard = "pk_hazard", irm1_hazard = "irm1_hazard",
    irm2_hazard = "irm2_hazard", irm3_hazard = "irm3_hazard",
    irm4_hazard = "irm4_hazard", tmdd_hazard = "tmdd_hazard")

#' Load (and cache) a bundled sim_tte_ode() library model file by basename
#'
#' The shared primitive behind \code{\link{.read_ode_library_model}}
#' (the \code{.ODE_LIBRARY_FILES}-keyed models) and the M-spline
#' knot-count dispatch (\code{\link{.mspline_file_for}}), which resolves
#' its file basename differently (from \code{length(knots)}, not a 1:1
#' public-name map) but loads it the same way.
#'
#' @param file Character. File basename (no extension) under
#'   \code{inst/models/library/}.
#' @return Compiled mrgsolve model object.
#' @noRd
.read_ode_library_model_file <- function(file) {
    mrgsolve::mread_cache(model = file, project = .cfile_dir("library"))
}

#' Load (and cache) a bundled sim_tte_ode() library model
#' @param model Character. One of \code{names(.ODE_LIBRARY_FILES)}.
#' @return Compiled mrgsolve model object.
#' @noRd
.read_ode_library_model <- function(model) {
    .read_ode_library_model_file(unname(.ODE_LIBRARY_FILES[model]))
}

# ---------------------------------------------------------------------
# sim_tte_ode() Phase 3 M-spline helpers (reports/09_phase3_report.md;
# convention decided there section 1, matching
# splines2::mSpline(..., degree = 2, intercept = TRUE) exactly).
# "mspline" is a single public model name dispatching to one of a fixed
# set of shipped interior-knot-count variants
# (reports/02_technical_design.md section 4 option 1), selected from
# length(knots) -- unlike the 1:1 .ODE_LIBRARY_FILES map the other
# three models use.
# ---------------------------------------------------------------------

#' Fixed degree of every shipped mspline_ode_k*.cpp variant (quadratic
#' M-splines, matching the convention this session pinned down -- see
#' reports/09_phase3_report.md section 1). Not user-configurable: a
#' different degree would need a different set of shipped model files.
#' @noRd
.MSPLINE_DEGREE <- 2L

#' Interior-knot counts shipped as library model variants.
#' @noRd
.MSPLINE_INTERIOR_COUNTS <- c(3L, 5L, 7L)

#' Resolve a length(knots) to its shipped mspline_ode_k*.cpp file
#' @param n_interior Integer. \code{length(knots)}.
#' @return Character file basename.
#' @noRd
.mspline_file_for <- function(n_interior) {
    if (!n_interior %in% .MSPLINE_INTERIOR_COUNTS) {
        stop("sim_tte_ode(model = \"mspline\") ships a fixed set of ",
            "interior-knot-count variants: ",
            paste(.MSPLINE_INTERIOR_COUNTS, collapse = ", "),
            "; got length(knots) = ", n_interior, ". Use one of the ",
            "supported counts, or sim_tte(type = \"ms\") for an ",
            "arbitrary knot count.", call. = FALSE)
    }
    paste0("mspline_ode_k", n_interior)
}

#' Validate sim_tte_ode(model = "mspline")'s knots/boundary_knots/coefs
#'
#' Mirrors the M-spline domain/coefficient constraints a fitted hazard
#' must satisfy: strictly increasing interior knots inside the open
#' boundary-knot interval (the same requirement
#' \code{splines2::mSpline()} itself enforces), and non-negative
#' coefficients (the hazard \code{eta * sum(c_m * M_m(t))} would
#' otherwise go negative, since \code{M_m(t) >= 0} always).
#' \code{end}/\code{idata$END} must not exceed \code{boundary_knots[2]}:
#' the basis (and therefore the hazard) is identically zero beyond it
#' (see \code{mspline_ode_k3.cpp}'s "DOMAIN" note), so silently allowing
#' \code{end} past it would silently stop hazard accrual rather than
#' erroring the way \code{sim_tte(type = "ms")}'s own
#' \code{end_time > max(time)} check already does for the same reason.
#'
#' @param knots Numeric vector of interior knots.
#' @param boundary_knots Numeric length-2 vector.
#' @param coefs Numeric vector, length \code{length(knots) + .MSPLINE_DEGREE + 1}.
#' @param end Numeric scalar, the requested follow-up horizon.
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.validate_mspline_args <- function(knots, boundary_knots, coefs, end) {
    if (!is.numeric(knots) || !length(knots)) {
        stop("sim_tte_ode(model = \"mspline\") requires 'knots' (a ",
            "numeric vector of interior knots).", call. = FALSE)
    }
    # Checked first, ahead of the strict-interior-placement check below:
    # an unsupported knot count is the more useful error to see first,
    # regardless of where the (wrong-count) knots happen to sit relative
    # to boundary_knots.
    .mspline_file_for(length(knots))
    if (!is.numeric(boundary_knots) || length(boundary_knots) != 2L ||
        any(!is.finite(boundary_knots)) || boundary_knots[1] >= boundary_knots[2]) {
        stop("'boundary_knots' must be a finite numeric vector ",
            "c(lower, upper) with lower < upper.", call. = FALSE)
    }
    if (any(!is.finite(knots)) || is.unsorted(knots, strictly = TRUE)) {
        stop("'knots' must be finite and strictly increasing (not ",
            "sorted automatically).", call. = FALSE)
    }
    if (knots[1] <= boundary_knots[1] || knots[length(knots)] >= boundary_knots[2]) {
        stop("'knots' must lie strictly inside (boundary_knots[1], ",
            "boundary_knots[2]) = (", boundary_knots[1], ", ",
            boundary_knots[2], ").", call. = FALSE)
    }
    n_expected_coefs <- length(knots) + .MSPLINE_DEGREE + 1L
    if (!is.numeric(coefs) || length(coefs) != n_expected_coefs) {
        stop("'coefs' must be numeric, length(knots) + ", .MSPLINE_DEGREE,
            " + 1 = ", n_expected_coefs, " (got ", length(coefs), ").",
            call. = FALSE)
    }
    if (any(!is.finite(coefs))) {
        stop("'coefs' must not contain NA, NaN, Inf, or -Inf values.",
            call. = FALSE)
    }
    if (any(coefs < 0)) {
        stop("'coefs' must be non-negative: the M-spline basis is ",
            "itself non-negative, so a negative coefficient would make ",
            "the hazard negative.", call. = FALSE)
    }
    if (end > boundary_knots[2]) {
        stop("'end' (", end, ") exceeds boundary_knots[2] (",
            boundary_knots[2], "). The M-spline basis is identically ",
            "zero beyond the upper boundary knot, so the hazard would ",
            "silently drop to zero there instead of extrapolating; ",
            "supply a wider 'boundary_knots' or a smaller 'end' (the ",
            "same restriction sim_tte(type = \"ms\")'s own ",
            "'end_time > max(time)' check already applies, for the ",
            "same reason).", call. = FALSE)
    }
    invisible(TRUE)
}

#' Build the sim_tte_ode(model = "mspline") $PARAM override list
#'
#' Maps \code{knots}/\code{boundary_knots}/\code{coefs} to the
#' \code{bk_lo}/\code{bk_hi}/\code{k1..kK}/\code{c1..cM} names every
#' \code{mspline_ode_k*.cpp} variant declares.
#'
#' @param knots,boundary_knots,coefs As validated by
#'   \code{\link{.validate_mspline_args}}.
#' @return Named list, mergeable into \code{sim_tte_ode()}'s own
#'   \code{param} argument.
#' @noRd
.build_mspline_param <- function(knots, boundary_knots, coefs) {
    p <- list(bk_lo = boundary_knots[1], bk_hi = boundary_knots[2])
    for (i in seq_along(knots)) p[[paste0("k", i)]] <- knots[i]
    for (i in seq_along(coefs)) p[[paste0("c", i)]] <- coefs[i]
    p
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

#' Build sim_tte_ode() time-varying-covariate rows from named covariates
#'
#' Generalizes \code{sim_tte()}'s single-covariate \code{lp_data}
#' mechanism (Phase F/G, \code{PHASE_F_DESIGN_REPORT.md} section 5) to
#' an arbitrary named set of log-linear covariates plus a coefficient
#' vector: \eqn{lp(t) = \sum_k \beta_k X_k(t)}, combined in R into a
#' single \code{lp} column and handed to \pkg{mrgsolve} as ordinary
#' covariate-update rows -- exactly the mechanism every
#' \code{sim_tte_ode()} library model already supports via its existing
#' \code{lp} parameter (Phase 1; no model-file change is required for
#' covariate support, see \code{reports/06_phase2_report.md}).
#'
#' Reuses \code{\link{.canonicalize_lp_data}},
#' \code{\link{.validate_lp_data_trajectories}}, and
#' \code{\link{.check_lp_data_coverage}} with \code{value_cols =
#' names(beta)} rather than duplicating them for the multi-covariate
#' case.
#'
#' @param covariates Data frame with a \code{time} column, an optional
#'   \code{ID} column, and one numeric column per name in \code{beta}
#'   (same population-level-vs-subject-specific contract as
#'   \code{sim_tte()}'s \code{lp_data}: no \code{ID} column means one
#'   shared trajectory recycled to every subject).
#' @param beta Named numeric vector of log-linear coefficients; every
#'   name must be a column of \code{covariates}.
#' @param n_subjects Integer. Number of subjects (must equal
#'   \code{covariates}' implied subject count and, when \code{covariates}
#'   has an \code{ID} column, that column must equal
#'   \code{1:n_subjects} -- see "Limitations" in \code{?sim_tte_ode}).
#' @param end Numeric scalar. Forwarded to
#'   \code{\link{.check_lp_data_coverage}} as \code{end_time}, always
#'   with \code{strict_coverage = FALSE}, since every
#'   \code{sim_tte_ode()} model is a genuine ODE that carries the last
#'   known covariate value forward, exactly like \code{sim_tte()}'s
#'   Weibull \code{lp_data} path.
#' @param cmt Integer. The compartment index covariate-update rows must
#'   target (an mrgsolve requirement for any model with at least one
#'   real compartment -- \code{cmt = 0} is only valid for a model with
#'   none, e.g. the closed-form \code{weibull.cpp}; see the Phase 1
#'   finding recorded in \code{reports/05_phase1_report.md} section 1).
#'   A parameter of this shared helper, not a second hardcoded copy, so
#'   a future library model with a different compartment layout only
#'   needs to pass its own value.
#' @return Data frame of covariate-update rows:
#'   \code{data.frame(ID, time, lp, evid = 1, amt = 0, cmt)}.
#' @noRd
.build_ode_covariate_rows <- function(covariates, beta, n_subjects, end,
    cmt = 1L) {
    if (is.null(names(beta)) || any(!nzchar(names(beta)))) {
        stop("'beta' must be a named numeric vector; every name must ",
            "match a column of 'covariates'.", call. = FALSE)
    }
    missing_cols <- setdiff(names(beta), names(covariates))
    if (length(missing_cols)) {
        stop("'covariates' is missing column(s) named in 'beta': ",
            paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
    # The "covariates$ID must equal exactly 1:n_subjects" restriction
    # (?sim_tte_ode "Covariates and the linear predictor", "Limitation")
    # is already enforced by .canonicalize_lp_data() itself (identical
    # check, used unchanged for sim_tte()'s own lp_data) -- not
    # duplicated here.
    canonical <- .canonicalize_lp_data(covariates, n_subjects,
        value_cols = names(beta))
    .validate_lp_data_trajectories(canonical, value_cols = names(beta))
    .check_lp_data_coverage(canonical, end, strict_coverage = FALSE)

    canonical$lp <- as.matrix(canonical[names(beta)]) %*% beta[names(beta)]
    data.frame(ID = canonical$ID, time = canonical$time, lp = canonical$lp,
        evid = 1L, amt = 0, cmt = cmt)
}

#' Merge a caller-supplied dosing \code{data} with covariate-update rows
#'
#' \code{sim_tte_ode()} builds \code{cov_rows} (\code{\link{.build_ode_covariate_rows}}/
#' \code{\link{.build_ode_covariate_rows_formula}}) with their own
#' \code{lp} column, then combines them with the caller's own \code{data}
#' (ordinary \pkg{mrgsolve} dosing/event rows, no \code{lp} of their own)
#' for a single \code{mrgsim()} call. Plain \code{dplyr::bind_rows(data,
#' cov_rows)} makes every dosing row's \code{lp} \code{NA} (a column only
#' \code{cov_rows} has), which is not just cosmetic: \pkg{mrgsolve} does
#' not carry a \code{NA} \code{$PARAM} value forward the way it carries a
#' real one (\code{nocb = FALSE} covariate semantics apply only to
#' \emph{present} values) -- once \code{lp} goes \code{NA} at a data-set
#' record, the \code{p11} ODE state itself becomes (and permanently
#' stays) \code{NaN} for the rest of that subject's trajectory, even
#' after a later covariate-update row gives \code{lp} a real value again
#' (verified directly, \code{reports/experiments/29_lp_na_mechanism.R}):
#' a real correctness bug (silently turning a genuine later event into
#' administrative censoring), not only the \code{mrgsolve::valid_data_set()}
#' warning it happens to also raise. A merged data set with dosing rows
#' out of time order relative to \code{cov_rows} (guaranteed whenever
#' \code{data} is appended before \code{cov_rows} without re-sorting, as
#' the pre-fix code did) can additionally make \pkg{mrgsolve} error
#' outright (\code{"the data set is not sorted by time"}) or silently
#' simulate a subject's rows as two disjoint blocks -- both fixed by the
#' same sort applied here.
#'
#' Every non-\code{cov_rows} row's \code{lp} is filled by
#' \code{\link{.locf_at}} against \emph{that row's own subject}'s
#' \code{cov_rows} trajectory (last covariate value in force at that
#' row's \code{time} -- \code{cov_rows} is guaranteed to cover \code{time
#' = 0} for every subject, \code{\link{.check_lp_data_coverage}}), so the
#' merged frame contains no \code{NA} in \code{lp} by construction, and
#' is sorted by \code{ID}/\code{time} so every subject's rows form one
#' non-decreasing-time block.
#'
#' @param data \code{NULL}, or the caller's own dosing/event data frame
#'   (\code{sim_tte_ode()}'s own \code{data} argument); must have an
#'   \code{ID} column when non-\code{NULL} (an unlabeled row's subject,
#'   and hence which \code{lp} trajectory it should carry, would
#'   otherwise be ambiguous -- an explicit error names this instead of
#'   guessing).
#' @param cov_rows Data frame from \code{\link{.build_ode_covariate_rows}}/
#'   \code{\link{.build_ode_covariate_rows_formula}}: one row per
#'   covariate-update record, columns \code{ID}, \code{time}, \code{lp},
#'   \code{evid}, \code{amt}, \code{cmt}.
#' @return \code{data} and \code{cov_rows} combined into one data frame,
#'   sorted by \code{ID} then \code{time}, with no \code{NA} in \code{lp}.
#' @noRd
.merge_ode_covariate_rows <- function(data, cov_rows) {
    if (is.null(data) || !nrow(data)) {
        return(cov_rows[order(cov_rows$ID, cov_rows$time), ])
    }
    data <- as.data.frame(data)
    if ("lp" %in% names(data)) {
        stop("'data' already has an 'lp' column, which conflicts with ",
            "the 'lp' column sim_tte_ode() builds from 'covariates'/",
            "'beta' (or 'formula'): supplying both leaves no way to ",
            "tell which one should apply at a dosing record. Fold ",
            "'data$lp' into 'covariates'/'beta' (or 'formula') instead ",
            "of supplying it directly in 'data', or rename 'data's own ",
            "column if it means something else.", call. = FALSE)
    }
    if (!"ID" %in% names(data)) {
        stop("'data' must have an 'ID' column when combined with ",
            "'covariates'/'beta' (or 'formula'): sim_tte_ode() cannot ",
            "otherwise tell which subject's 'lp' trajectory a dosing ",
            "row in 'data' belongs to.", call. = FALSE)
    }
    cov_by_id <- split(cov_rows, cov_rows$ID)
    unknown <- setdiff(unique(data$ID), cov_rows$ID)
    if (length(unknown)) {
        stop("'data' has ID(s) not present in 'covariates': ",
            paste(unknown, collapse = ", "), ". Every subject dosed ",
            "via 'data' must also have a covariate trajectory.",
            call. = FALSE)
    }
    data$lp <- vapply(seq_len(nrow(data)), function(i) {
        traj <- cov_by_id[[as.character(data$ID[i])]]
        .locf_at(traj$time, traj$lp, data$time[i])
    }, numeric(1))
    merged <- dplyr::bind_rows(data, cov_rows)
    merged[order(merged$ID, merged$time), ]
}

#' Apply an omega/sigma matrix to a sim_tte_ode() model, with an
#' informative error instead of mrgsolve's own cryptic one
#'
#' \code{mrgsolve::omat(mod, matrix)}/\code{smat()} only \emph{update}
#' an already-declared OMEGA/SIGMA block; a model that declares none at
#' all (every \code{sim_tte_ode()} library model as of Phase 4, none of
#' which wires \code{ETA(n)} into any parameter -- verified directly,
#' including against an unmodified \code{mrgsolve::mread()} of
#' \code{pk2cmt} itself, which fails identically) errors with
#' \code{"improper signature: omat"}, which names neither the argument
#' nor the reason. This wraps that call so the user sees a clear
#' explanation instead (\code{reports/10_phase4_report.md} open risks).
#'
#' @param mod A compiled mrgsolve model object.
#' @param value The \code{omega}/\code{sigma} matrix supplied by the
#'   caller.
#' @param arg Character. \code{"omega"} or \code{"sigma"}, for the error
#'   message.
#' @param fn \code{mrgsolve::omat} or \code{mrgsolve::smat}.
#' @return The updated model object.
#' @noRd
.apply_ode_matlist <- function(mod, value, arg, fn) {
    tryCatch(fn(mod, value), error = function(e) {
        stop("Could not apply '", arg, "' to this sim_tte_ode() model: ",
            conditionMessage(e), ". This usually means the model does ",
            "not declare a matching $", toupper(arg), " block: mrgsolve's ",
            "omat()/smat() only update an already-declared block, they ",
            "do not create one from nothing, and none of the built-in ",
            "sim_tte_ode() library models declares one as of this ",
            "release (see ?sim_tte_ode \"Between-subject variability\"). ",
            "'", arg, "' is only usable with a model you supply yourself ",
            "that already declares a $", toupper(arg), " block of the ",
            "same dimension.", call. = FALSE)
    })
}

#' Between-subject variability (BSV) targets for the built-in PK/PD
#' hazard library
#'
#' \code{reports/11_bsv_review.md} option (b): per-subject parameter
#' values via \code{idata} (the same mechanism \code{U}/\code{END}
#' already use), not a declared \code{$OMEGA} block -- verified
#' directly (that session's \code{reports/experiments/}) to need zero
#' model-file changes, unlike the classic \code{TVCL}/\code{ETA()}
#' rename idiom.
#'
#' One entry per PK/PD library model (\code{reports/10_phase4_report.md}),
#' the ordered character vector of every *structural* PK/PD parameter
#' declared in that model's own \code{$PARAM} block(s) -- i.e. the
#' backbone's own parameters, taken from the compiled model itself
#' while writing this list (\code{reports/12_bsv_implementation_report.md}),
#' not from memory. The hazard-side parameters this package's own
#' scaffold adds (\code{H0}, \code{lp}, \code{beta_cp}/\code{beta_r}/
#' \code{beta_rc}, \code{U}, \code{END}) are deliberately excluded --
#' BSV on the hazard's own log-linear terms is not what this mechanism
#' is for (a user wanting that already has \code{lp}/\code{covariates}).
#' \code{exponential}/\code{weibull}/\code{gompertz}/\code{mspline}
#' have no entry: \code{omega}/\code{sigma} on those keeps the existing
#' \code{\link{.apply_ode_matlist}} informative error.
#' @noRd
.ODE_BSV_TARGETS <- list(
    pk_hazard = c("CL", "V2", "Q", "V3", "KA", "KA2", "VMAX", "KM"),
    irm1_hazard = c("CL", "V2", "Q", "V3", "KA", "KA2", "KIN", "KOUT",
        "IC50", "IMAX", "n", "VMAX", "KM"),
    irm2_hazard = c("CL", "V2", "Q", "V3", "KA", "KA2", "KIN", "KOUT",
        "IC50", "IMAX", "n", "VMAX", "KM"),
    irm3_hazard = c("CL", "V2", "Q", "V3", "KA", "KA2", "KIN", "KOUT",
        "EC50", "EMAX", "n", "VMAX", "KM"),
    irm4_hazard = c("CL", "V2", "KA", "KA2", "Q", "V3", "KIN", "KOUT",
        "EC50", "EMAX", "VMAX", "KM", "n"),
    tmdd_hazard = c("KPT", "KTP", "V2", "KA", "KA2", "KEL", "R0", "KDEG",
        "KINT", "KON", "KOFF")
)

#' Draw per-subject BSV parameter values and append them to \code{idata}
#'
#' Implements \code{reports/11_bsv_review.md} option (b): draws ETAs
#' from \code{omega} with \code{\link[mrgsolve]{mvgauss}} (already
#' exported by mrgsolve, no new dependency), transforms them
#' log-normally (\code{param_i = param_pop * exp(eta)}, the convention
#' already implicit in every one of these models' own PK parameters),
#' and appends the result as ordinary \code{idata} columns -- exactly
#' the mechanism \code{\link{.build_ode_idata}} already uses for
#' \code{U}. No model file is read or changed by this function.
#'
#' \strong{Seeding}: draws \code{n = nrow(idata)} rows in one call,
#' without passing \code{seed} to \code{\link[mrgsolve]{mvgauss}}
#' itself, relying on \code{sim_tte_ode()}'s own top-level
#' \code{set.seed(seed)} the way \code{U} already does. Verified
#' directly (\code{reports/experiments/}, this session) that
#' \code{mvgauss()}'s draw, for a fixed \code{seed}, is reproducible
#' and -- surprisingly -- \strong{insensitive to any \code{runif()}/
#' \code{rnorm()} draws made in between} (e.g. the \code{U} draw
#' already made by \code{\link{.build_ode_idata}}): calling it after
#' \code{U} or calling it as the very first draw of the session gives
#' the identical result for the same \code{seed}. This was confirmed
#' empirically, not assumed from mrgsolve's documentation. It does mean
#' one \code{seed} deterministically reproduces \strong{both} \code{U}
#' and the BSV draw regardless of draw order -- the property this
#' function's reproducibility contract relies on -- but two separate
#' \code{mvgauss()} calls in the same session (not done here; this
#' function only ever calls it once, on the full requested matrix) do
#' differ from each other, so this is not "always the same output"
#' more generally, only "insensitive to intervening ordinary R draws".
#'
#' @param omega Numeric square matrix. Either \strong{named}
#'   (\code{dimnames} identical row/column parameter names, any subset
#'   of \code{targets}, any order) or \strong{unnamed} (dimension must
#'   equal \code{length(targets)}; applied positionally, in
#'   \code{targets} order).
#' @param targets Character vector of valid BSV target names for this
#'   model -- \code{.ODE_BSV_TARGETS[[model]]} for a built-in library
#'   model, or a \code{tte_model()}-converted model's own
#'   \code{bsv_targets} (\code{reports/13_converter_design.md} section
#'   5; generalized from a \code{model} character-name lookup to this
#'   plain vector so both sources share one code path).
#' @param model_label Character. Used only in error messages (a
#'   built-in library model's name, or a converted model's own
#'   \code{name}).
#' @param n Integer. Expected number of subjects; checked against
#'   \code{nrow(idata)} defensively (the draw itself uses
#'   \code{nrow(idata)}).
#' @param idata Data frame, already carrying \code{ID}/\code{U}/
#'   \code{END} (\code{\link{.build_ode_idata}}'s output). Must not
#'   already contain a column named after any drawn target.
#' @param param Named list or named numeric vector of the
#'   \strong{already-resolved} population parameter values (i.e.
#'   \code{as.list(mrgsolve::param(mod))} after the caller's own
#'   \code{param =} overrides have been merged in) -- so
#'   \code{param = list(CL = 3)} together with \code{omega} on
#'   \code{"CL"} centers the draw on 3, not the model file's own
#'   default.
#' @return \code{idata} with one new column per drawn target.
#' @noRd
.build_ode_bsv_idata <- function(omega, targets, model_label, n, idata,
    param) {
    if (is.null(targets)) {
        stop("sim_tte_ode(): model = \"", model_label, "\" has no ",
            "between-subject-variability targets. See ?sim_tte_ode ",
            "\"Between-subject variability\".", call. = FALSE)
    }
    if (!is.matrix(omega) || !is.numeric(omega) || nrow(omega) != ncol(omega)) {
        stop("'omega' must be a square numeric matrix. Valid ",
            "between-subject-variability targets for model = \"",
            model_label, "\": ", paste(targets, collapse = ", "), ".",
            call. = FALSE)
    }
    dn <- dimnames(omega)
    named <- !is.null(dn) && !is.null(dn[[1]]) && !is.null(dn[[2]])
    used <- if (named) {
        if (!identical(dn[[1]], dn[[2]])) {
            stop("'omega' row and column names must be identical (same ",
                "parameter names, same order).", call. = FALSE)
        }
        if (anyDuplicated(dn[[1]])) {
            stop("'omega' has duplicated parameter names.", call. = FALSE)
        }
        unknown <- setdiff(dn[[1]], targets)
        if (length(unknown)) {
            stop("'omega' names parameter(s) not in model = \"",
                model_label, "\"'s between-subject-variability targets: ",
                paste(unknown, collapse = ", "), ". Valid targets: ",
                paste(targets, collapse = ", "), ".", call. = FALSE)
        }
        dn[[1]]
    } else if (!is.null(dn)) {
        stop("'omega' must have both row and column names, or neither ",
            "(positional, matching all ", length(targets), " targets ",
            "in order: ", paste(targets, collapse = ", "), ").",
            call. = FALSE)
    } else {
        if (nrow(omega) != length(targets)) {
            stop("An unnamed 'omega' must have dimension ", length(targets),
                " x ", length(targets), " (one row/column per target, in ",
                "order: ", paste(targets, collapse = ", "), ") for model = \"",
                model_label, "\"; got ", nrow(omega), " x ", ncol(omega),
                ". Name 'omega's dimnames to supply a subset instead.",
                call. = FALSE)
        }
        targets
    }
    clash <- intersect(used, names(idata))
    if (length(clash)) {
        stop("'idata' already has column(s) matching between-subject-",
            "variability target(s): ", paste(clash, collapse = ", "),
            ". Supply either 'idata' columns or 'omega' for a given ",
            "parameter, not both.", call. = FALSE)
    }
    param <- as.list(param)
    missing_pop <- setdiff(used, names(param))
    if (length(missing_pop)) {
        stop("No resolved population value for: ",
            paste(missing_pop, collapse = ", "), ".", call. = FALSE)
    }
    if (!isTRUE(n == nrow(idata))) {
        stop("'n' (", n, ") does not match nrow(idata) (", nrow(idata),
            ").", call. = FALSE)
    }

    omega_use <- if (named) omega[used, used, drop = FALSE] else omega
    eta <- mrgsolve::mvgauss(omega_use, n = nrow(idata))
    for (j in seq_along(used)) {
        idata[[used[j]]] <- as.numeric(param[[used[j]]]) * exp(eta[, j])
    }
    idata
}

#' Classify why a subject's event time used the reported-grid fallback
#'
#' Two fallback-triggering mechanisms have been characterized so far,
#' purely from the bracket columns \code{\link{.resolve_ode_events}}
#' already captures -- no new model-side logic
#' (\code{reports/04_author_decisions.md} "After Phase 3" decision 1):
#' \itemize{
#'   \item \strong{weibull_type}: \code{p_post} outside \code{[0, 1]} --
#'     a genuine solver overshoot across one large internal step
#'     (\code{reports/08_phase2_5_report.md} section 4, Weibull
#'     \code{shape >= 5}). Remedy: a finer \code{delta}/\code{add}.
#'   \item \strong{mspline_type}: \code{t_pre == tevt} exactly, with
#'     both \code{p_pre} and \code{p_post} legitimate probabilities --
#'     a same-timestamp, multiple-corrector-iteration degeneracy near a
#'     locally steep hazard (\code{reports/09_phase3_report.md} open
#'     risk 1). Remedy: tighter \code{rtol}/\code{atol}, not
#'     \code{delta} (confirmed flat across \code{delta} in that
#'     report).
#'   \item \strong{unclassified}: neither signature matches (including
#'     when the bracket columns are unavailable at all, i.e.
#'     \code{has_bracket = FALSE} in the caller) -- generic advice only.
#' }
#'
#' @param t_pre,p_pre,tevt,p_post Numeric scalars, the bracket as
#'   captured from the model, or \code{NA_real_} if unavailable.
#' @return Character scalar: one of \code{"weibull_type"},
#'   \code{"mspline_type"}, \code{"unclassified"}.
#' @noRd
.classify_ode_fallback <- function(t_pre, p_pre, tevt, p_post) {
    if (is.finite(p_post) && (p_post < 0 || p_post > 1)) {
        return("weibull_type")
    }
    if (is.finite(t_pre) && is.finite(tevt) && t_pre == tevt &&
        is.finite(p_pre) && p_pre >= 0 && p_pre <= 1 &&
        is.finite(p_post) && p_post >= 0 && p_post <= 1) {
        return("mspline_type")
    }
    "unclassified"
}

#' Refine an event time from the in-solver pre/post-crossing bracket
#'
#' Grid-free counterpart to the old reported-grid refinement
#' (\code{reports/06_phase2_report.md} section 4 /
#' \code{reports/04_author_decisions.md} "After the test runbook /
#' Phase 2.5"): interpolates between the solver's own last pre-crossing
#' evaluation (\code{t_pre}, \code{p_pre}) and the crossing evaluation
#' itself (\code{tevt}, \code{p_post}) -- an interval one internal
#' solver step wide, captured directly from the model
#' (\code{T_PRE}/\code{P_PRE}/\code{P_POST}, see the library \file{.cpp}
#' models), rather than between the two \emph{reported} rows bracketing
#' \code{u}. Because the bracket does not depend on the requested output
#' grid, neither does the refined time.
#'
#' @param t_pre,p_pre,tevt,p_post Numeric scalars, captured from the
#'   model as described above.
#' @param u Numeric scalar. The subject's own uniform draw.
#' @return Numeric scalar, or \code{NA_real_} if the bracket is
#'   degenerate (\code{t_pre >= tevt}, or a non-finite input) -- the
#'   caller falls back to the reported-grid method in that case.
#' @noRd
.refine_ode_event_time_insolver <- function(t_pre, p_pre, tevt, p_post, u) {
    if (!is.finite(t_pre) || !is.finite(p_pre) || !is.finite(p_post) ||
        !is.finite(tevt) || t_pre >= tevt) {
        return(NA_real_)
    }
    .interpolate_log_survival(t_i = t_pre, t_ip1 = tevt, s_i = p_pre,
        s_ip1 = p_post, u = u)
}

#' Resolve sim_tte_ode() event/censoring times from a simulated trajectory
#'
#' Applies the censoring rule and refinement step documented in
#' \code{?sim_tte_ode} ("Boundary guard" / "Event-time refinement") to
#' the \code{mrgsim()} output of a library model following the
#' \code{p11}/\code{TEVT}/\code{event_found}/\code{T_PRE}/\code{P_PRE}/
#' \code{P_POST}/\code{U}/\code{END} contract
#' (\code{.validate_ode_model_contract()}).
#'
#' Refinement is grid-free by default
#' (\code{\link{.refine_ode_event_time_insolver}}): it uses the solver's
#' own pre/post-crossing bracket, not the reported output grid, so its
#' accuracy no longer depends on \code{delta}
#' (\code{reports/08_phase2_5_report.md}). The previous reported-grid
#' method (reusing \code{\link{.get_tte}}/
#' \code{\link{.interpolate_log_survival}} between the two \emph{reported}
#' rows bracketing \code{u} -- the same machinery
#' \code{\link{sim_tte_df}} uses for \code{event_time_method =
#' "log_survival"}) is kept as a fallback for the rare case the in-solver
#' bracket is degenerate (a hand-supplied \code{traj} missing the
#' \code{T_PRE}/\code{P_PRE}/\code{P_POST} columns, or, in principle, a
#' model that never records a genuine pre-crossing evaluation). One
#' \code{message()} is emitted per call naming how many subjects fell
#' back, so a coarse-\code{delta}, poor-refinement outcome is never
#' silent (\code{?sim_tte_ode} "Event-time refinement"; suppressible via
#' \code{suppressMessages()}, never changes the returned result -- see
#' the general guardrail policy in
#' \code{reports/04_author_decisions.md}).
#'
#' @param traj Data frame: \code{mrgsim()} output with columns \code{ID},
#'   \code{time}, \code{p11}, \code{TEVT}, \code{event_found}, \code{U},
#'   \code{END}, one or more rows per subject in ascending \code{time}
#'   order (mrgsolve's own reporting order; re-sorted defensively below).
#'   \code{T_PRE}/\code{P_PRE}/\code{P_POST} are optional (their absence
#'   forces the reported-grid fallback for every event subject) so this
#'   function keeps working against the synthetic/hand-built trajectories
#'   used in the test suite's own defensive-fallback tests.
#' @return Data frame with columns \code{ID}, \code{sim_time},
#'   \code{sim_status}, \code{sim_reason} (\code{"event"} or
#'   \code{"administrative"} -- \code{sim_tte_ode()}'s own \code{censoring
#'   =} argument, when supplied, further updates this column via
#'   \code{\link{add_censoring}}; see \code{reports/16_censoring_design.md}).
#' @noRd
.resolve_ode_events <- function(traj) {
    has_bracket <- all(c("T_PRE", "P_PRE", "P_POST") %in% names(traj))
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
            return(data.frame(ID = id, sim_time = end_i, sim_status = 0L,
                sim_reason = "administrative", used_fallback = FALSE,
                fallback_signature = NA_character_))
        }

        sim_time <- if (has_bracket) {
            .refine_ode_event_time_insolver(t_pre = last$T_PRE,
                p_pre = last$P_PRE, tevt = last$TEVT, p_post = last$P_POST,
                u = u_i)
        } else {
            NA_real_
        }
        used_fallback <- is.na(sim_time)
        fallback_signature <- if (!used_fallback) {
            NA_character_
        } else if (has_bracket) {
            .classify_ode_fallback(t_pre = last$T_PRE, p_pre = last$P_PRE,
                tevt = last$TEVT, p_post = last$P_POST)
        } else {
            "unclassified"
        }

        if (used_fallback) {
            # Fallback: reported-grid refinement (pre-Phase-2.5 method).
            p11 <- .normalize_survival(sub$p11)
            etime_idx <- .get_tte(u_i, p11)
            if (etime_idx == -99L) {
                # Defensive, not expected to be reached: `end` is always
                # a reported row (.resolve_output_grid()) and p11 is
                # monotonically non-increasing, so a latched pre-`end`
                # event implies some reported row has p11 <= U by the
                # time `end` is reached. Falls back to the raw (already
                # < end_i, by the check above) in-solver estimate.
                sim_time <- last$TEVT
            } else if (etime_idx == 1L) {
                # Crossing already present at the first reported
                # observation -- no earlier point to interpolate from.
                sim_time <- sub$time[1]
            } else {
                i <- etime_idx - 1L
                sim_time <- .interpolate_log_survival(t_i = sub$time[i],
                    t_ip1 = sub$time[etime_idx], s_i = p11[i],
                    s_ip1 = p11[etime_idx], u = u_i)
            }
        }
        data.frame(ID = id, sim_time = sim_time, sim_status = 1L,
            sim_reason = "event", used_fallback = used_fallback,
            fallback_signature = fallback_signature)
    })
    out <- dplyr::bind_rows(rows)
    n_fallback <- sum(out$used_fallback)
    if (n_fallback > 0L) {
        # Differentiated fallback message (reports/04_author_decisions.md
        # "After Phase 3" decision 1): classify purely from the already-
        # captured bracket columns (.classify_ode_fallback()), report
        # only the non-zero buckets, each with its own matching remedy.
        n_weibull <- sum(out$fallback_signature == "weibull_type",
            na.rm = TRUE)
        n_mspline <- sum(out$fallback_signature == "mspline_type",
            na.rm = TRUE)
        n_other <- sum(out$fallback_signature == "unclassified",
            na.rm = TRUE)
        parts <- character(0)
        if (n_weibull > 0L) {
            parts <- c(parts, paste0(n_weibull, " subject(s) had an ",
                "unphysical P_POST outside [0, 1] (a genuine solver ",
                "overshoot across one large internal step) -- a finer ",
                "'delta'/'add' resolves this"))
        }
        if (n_mspline > 0L) {
            parts <- c(parts, paste0(n_mspline, " subject(s) had ",
                "T_PRE == TEVT with both P_PRE and P_POST legitimate ",
                "probabilities (a same-timestamp corrector-iteration ",
                "degeneracy near a locally steep hazard) -- tighter ",
                "'rtol'/'atol' resolves this, not 'delta'"))
        }
        if (n_other > 0L) {
            parts <- c(parts, paste0(n_other, " subject(s) fell back ",
                "for an unclassified reason (or the bracket columns ",
                "were unavailable) -- try a finer 'delta'/'add' and/or ",
                "tighter 'rtol'/'atol'"))
        }
        message("sim_tte_ode(): ", n_fallback, " subject(s) used the ",
            "reported-grid refinement fallback: ",
            paste(parts, collapse = "; "), ". See ?sim_tte_ode ",
            "\"Event-time refinement\".")
    }
    out$used_fallback <- NULL
    out$fallback_signature <- NULL
    out
}

#' Emit a guardrail message for a very small Weibull shape
#'
#' \code{weibull_ode.cpp} always runs (no error, no warning) for any
#' \code{shape > 0} (\code{reports/06_phase2_report.md} risk R1), but
#' accuracy near \eqn{t = 0} degrades below \code{shape = 0.05}. This is
#' informational only -- it never changes what runs, matching every
#' other guardrail in this package (\code{reports/04_author_decisions.md}
#' "After the test runbook / Phase 2.5"); suppress with
#' \code{suppressMessages()} if not wanted.
#'
#' @param mod A compiled mrgsolve model object, with \code{param}
#'   already merged (so \code{param(mod)[["shape"]]} reflects any
#'   user override).
#' @param model Character. The public \code{sim_tte_ode()} model name.
#' @return \code{NULL}, invisibly. Called for the \code{message()} side
#'   effect only.
#' @noRd
.check_weibull_shape_guardrail <- function(mod, model) {
    if (identical(model, "weibull")) {
        shape <- mrgsolve::param(mod)[["shape"]]
        if (is.finite(shape) && shape < 0.05) {
            message("sim_tte_ode(): Weibull 'shape' = ", shape, " is ",
                "below 0.05; the in-solver mechanism's accuracy near ",
                "t = 0 degrades below this threshold (see ?sim_tte_ode ",
                "\"Weibull shape support\"). The simulation proceeds ",
                "unchanged; if exactness near t = 0 matters more, use ",
                "sim_tte(type = \"weibull\") instead.")
        }
    }
    invisible(NULL)
}
