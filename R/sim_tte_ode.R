#' Simulate Time-to-Event Data via In-Solver Event Detection
#'
#' Phase 1 implementation of the in-solver event-detection mechanism
#' described in \code{reports/02_technical_design.md} section 2: the
#' event/censoring time is located \emph{during} ODE integration (at the
#' solver's own internal evaluations), by drawing a per-subject
#' \eqn{U \sim \mathrm{Uniform}(0, 1)} and integrating the survival
#' probability \code{p11} as a state, latching the first solver time at
#' which \code{p11 <= U}. This is a different mechanism from
#' \code{\link{sim_tte}}/\code{\link{sim_tte_df}}, which resolve event
#' times \emph{after} the fact from an already-materialized trajectory;
#' see \code{reports/01_codebase_audit.md} section 4 and
#' \code{reports/02_technical_design.md} section 1 for why the two are
#' separate code paths, not two modes of one function.
#'
#' Only \code{model = "exponential"} (constant hazard) is implemented in
#' this phase. Weibull, Gompertz, M-spline, and PK/PD-linked models, and
#' an exported model-discovery function analogous to
#' \code{\link{simtte_example_models}}, are later phases
#' (\code{reports/03_implementation_plan.md} Phase 5).
#'
#' @param model Character string naming a bundled library model.
#'   Currently only \code{"exponential"} is implemented. A user-supplied
#'   compiled \code{mrgmod} is not yet accepted (Phase 5).
#' @param param Named list or named numeric vector of population-level
#'   \code{$PARAM} overrides, forwarded to \code{\link[mrgsolve]{param}}.
#'   \code{U} and \code{END} are ordinarily supplied per subject via
#'   \code{idata} (see below), not here.
#' @param omega,sigma Optional matrices, forwarded to
#'   \code{\link[mrgsolve]{omat}}/\code{\link[mrgsolve]{smat}}.
#' @param n Integer. Number of subjects. Ignored (and inferred from
#'   \code{nrow(idata)}) when \code{idata} is supplied.
#' @param end Numeric scalar. Administrative censoring horizon. Written
#'   into the model's \code{END} parameter for every subject that does
#'   not already supply its own \code{END} column in \code{idata}, and
#'   used to guard in-solver event detection against latching a time
#'   past the requested follow-up (see "Boundary guard" below).
#' @param delta Numeric scalar. Reported-grid spacing, used to build the
#'   requested output grid \code{seq(0, end, by = delta)} when \code{add}
#'   is not supplied. Ignored if \code{add} is supplied.
#' @param add Numeric vector of explicit reported times. If \code{NULL}
#'   (default), \code{seq(0, end, by = delta)} is used. Either way,
#'   \code{end} itself is always included in the resolved grid (via
#'   \code{.resolve_output_grid()}), which the censoring rule and
#'   refinement step below both depend on.
#' @param idata \code{NULL} (default) or a data frame with an \code{ID}
#'   column and, optionally, \code{U} and/or \code{END} columns. Missing
#'   \code{U} is drawn as \code{stats::runif(nrow(idata))} (one
#'   independent draw per row/subject; see "Reproducibility" below).
#'   Missing \code{END} is filled with \code{end}. Any other column is
#'   passed through to \pkg{mrgsolve} as an ordinary per-subject
#'   covariate/parameter.
#' @param data Optional dosing/time-varying-covariate event data,
#'   forwarded to \code{\link[mrgsolve]{mrgsim}}'s own \code{data}
#'   argument.
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is
#'   called before \code{U} is drawn.
#' @param keep_trajectory Logical. If \code{TRUE}, the full \code{mrgsim()}
#'   output is retained in the returned object's \code{$trajectory}.
#'   Default \code{FALSE} (the coarse-output/in-solver-detection design
#'   is partly motivated by \emph{not} requiring a dense materialized
#'   trajectory for large cohorts; see
#'   \code{reports/02_technical_design.md} section 5.3).
#' @param ... Forwarded to \code{\link[mrgsolve]{mrgsim}}. \code{tgrid},
#'   \code{obsonly}, \code{nocb}, and \code{carry_out}/\code{carry.out}
#'   are controlled internally to guarantee the contracts documented
#'   here and cannot be overridden; supplying them raises an error
#'   naming the argument (the same \code{.check_reserved_dots()}
#'   mechanism \code{\link{sim_tte}} already uses). This is where
#'   solver-tolerance arguments relevant to in-solver event-time
#'   accuracy (\code{rtol}, \code{atol}, \code{hmax}, ...) are set; see
#'   \code{reports/02_technical_design.md} section 2.4/5.4.
#'
#' @return An object of class \code{"simtte_ode_sim"}, a list with:
#' \describe{
#'   \item{events}{Data frame with columns \code{ID}, \code{sim_time},
#'     \code{sim_status} (1 = event, 0 = censored).}
#'   \item{trajectory}{The full \code{mrgsim()} output data frame, or
#'     \code{NULL} unless \code{keep_trajectory = TRUE}.}
#'   \item{model}{The compiled \code{mrgmod} object used.}
#'   \item{seed}{The \code{seed} argument, as supplied (possibly
#'     \code{NULL}).}
#'   \item{call}{The matched call, for reproducibility bookkeeping.}
#' }
#'
#' @section Boundary guard: the administrative-censoring contract:
#' The shipped library models latch \code{TEVT} (the in-solver event
#' time) only while \code{SOLVERTIME <= END}, so an internal solver
#' evaluation past the requested follow-up horizon can never be recorded
#' as an event time (this guards against the exact defect documented in
#' \code{reports/PHASE_E_THRESHOLD_TRACKING_REPORT.md} section 10 and
#' reproduced/fixed in \code{reports/02_technical_design.md} section
#' 2.5). \code{sim_tte_ode()} additionally re-applies the boundary rule
#' on the R side, independent of the in-model guard: a subject is
#' censored (\code{sim_status = 0}, \code{sim_time = end}) if no event
#' was latched at all, \strong{or} if the latched time is \code{>= end}.
#' This is deliberately belt-and-braces: the in-model guard is an
#' internal robustness measure; the R-side rule above is the documented
#' contract.
#'
#' @section Event-time refinement:
#' The raw in-solver \code{TEVT} is quantized to the ODE solver's own
#' internal step grid, which is generally finer than, but not
#' independent of, the requested output grid and solver tolerance (see
#' \code{reports/02_technical_design.md} section 2.4). For every subject
#' classified as an event, \code{sim_tte_ode()} therefore refines the
#' final \code{sim_time} using the already-validated cumulative-hazard
#' interpolation \code{.interpolate_log_survival()} (the same
#' method \code{\link{sim_tte_df}} offers as
#' \code{event_time_method = "log_survival"}), applied between the two
#' \emph{reported} trajectory rows bracketing that subject's own
#' \code{U} threshold. \code{TEVT} itself therefore only decides
#' \strong{whether} a subject had an event before \code{end}; the
#' reported \code{sim_time} for an event subject comes from this
#' interpolation step, not from the raw \code{TEVT} value (the latter is
#' used only as a defensive fallback in the extremely unlikely case that
#' no reported row's \code{p11} is at or below \code{U}, which should not
#' occur given the grid-construction/monotonicity guarantees above).
#'
#' @section Reproducibility:
#' If \code{idata} does not supply \code{U}, exactly one
#' \code{stats::runif(1)} is drawn per row of \code{idata} (one row per
#' subject), in row order, \emph{before} \pkg{mrgsolve} is called. With a
#' fixed \code{seed} (or a \code{set.seed()} call by the caller), results
#' are therefore reproducible for a fixed \code{idata} row order.
#'
#' @seealso \code{\link{sim_tte}}, \code{\link{sim_tte_df}}, for the
#'   grid-based mechanism this function does not replace.
#' @export
#' @examples
#' \donttest{
#' sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
#'   n = 50, end = 30, delta = 5, seed = 1)
#' head(sim$events)
#' }
sim_tte_ode <- function(model, param = list(), omega = NULL, sigma = NULL,
    n = 1L, end = 100, delta = 1, add = NULL, idata = NULL, data = NULL,
    seed = NULL, keep_trajectory = FALSE, ...) {

    .check_reserved_dots(list(...))
    model <- match.arg(model, choices = names(.ODE_LIBRARY_FILES))
    .validate_end_time(end)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    mod <- .read_ode_library_model(model)
    .validate_ode_model_contract(mod)

    if (!is.null(omega)) {
        mod <- mrgsolve::omat(mod, omega)
    }
    if (!is.null(sigma)) {
        mod <- mrgsolve::smat(mod, sigma)
    }
    if (length(param)) {
        mod <- mrgsolve::param(mod, as.list(param))
    }

    idata <- .build_ode_idata(idata, n = n, end = end)
    grid <- .resolve_output_grid(if (is.null(add)) seq(0, end, by = delta)
        else add, end, type = "weibull")

    out <- as.data.frame(mrgsolve::mrgsim(mod, idata = idata, data = data,
        tgrid = grid, obsonly = TRUE, nocb = FALSE,
        carry_out = c("U", "END"), ...))

    events <- .resolve_ode_events(out)

    result <- list(events = events,
        trajectory = if (isTRUE(keep_trajectory)) out else NULL,
        model = mod, seed = seed, call = match.call())
    class(result) <- "simtte_ode_sim"
    result
}

#' Print an in-solver simulation result
#' @param x A \code{"simtte_ode_sim"} object.
#' @param ... Ignored.
#' @export
#' @method print simtte_ode_sim
print.simtte_ode_sim <- function(x, ...) {
    n_event <- sum(x$events$sim_status == 1L)
    n_cens <- sum(x$events$sim_status == 0L)
    cat("<simtte_ode_sim>", nrow(x$events), "subjects:", n_event,
        "event(s),", n_cens, "censored\n")
    if (!is.null(x$trajectory)) {
        cat(" trajectory kept:", nrow(x$trajectory), "rows\n")
    }
    print(utils::head(x$events))
    invisible(x)
}
