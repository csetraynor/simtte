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
#' @param model Character model name ("weibull" or "ms").
#' @return Compiled mrgsolve model object.
#' @noRd
.read_model_static_cache <- function(model) {
    pkg_model_file <- .cfile_dir()
    if (model == "weibull" || model == "ms") {
        mod_surv <- mrgsolve::mread_cache(model = model,
            project = pkg_model_file)
    } else {
        stop("Model '", model, "' must be 'ms' or 'weibull'.")
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
