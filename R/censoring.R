#' Apply independent right censoring to a simulated events data frame
#'
#' Adds a second, independent censoring process to an already-simulated
#' time-to-event data frame -- the mechanism described in
#' \code{reports/16_censoring_design.md} (option A): a censoring time
#' \eqn{C_i} is drawn per subject from \code{censoring}, and the observed
#' outcome becomes \eqn{\min(T_i, C_i, end)}, with \code{sim_status = 1}
#' only if \eqn{T_i} (the original \code{sim_time}) is the minimum.
#' Works on the output of \code{\link{sim_tte_ode}}, \code{\link{sim_tte}},
#' or \code{\link{sim_tte_df}} alike -- all three already share the same
#' \code{ID}/\code{sim_time}/\code{sim_status} shape.
#'
#' @param events Data frame with (by default) \code{ID}, \code{sim_time},
#'   \code{sim_status} columns -- column names are configurable via
#'   \code{id_var}/\code{time_var}/\code{status_var} for a data frame
#'   using different names. Any other column is left untouched.
#' @param censoring A list describing the censoring-time distribution,
#'   with a \code{dist} element:
#'   \describe{
#'     \item{\code{"exponential"}}{\code{rate} (positive scalar).
#'       \code{C ~ Exp(rate)}.}
#'     \item{\code{"weibull"}}{\code{shape}, \code{scale} (positive
#'       scalars). \code{C ~ Weibull(shape, scale)}
#'       (\code{\link[stats]{rweibull}} parameterization).}
#'     \item{\code{"uniform"}}{\code{min}, \code{max} (\code{0 <= min <
#'       max}). \code{C ~ Uniform(min, max)} -- the common case for a
#'       uniform-accrual administrative censoring time.}
#'     \item{\code{"user"}}{\code{fun}, a function \code{n -> numeric(n)}
#'       returning \code{n} non-negative, finite censoring times. The
#'       escape hatch for any distribution not covered above.}
#'   }
#' @param end Numeric scalar. Administrative censoring horizon -- the
#'   third term of \code{min(T, C, end)}. Always wins: no observed time
#'   in the result exceeds \code{end}. Required (part of the censoring
#'   formula itself, not inferred from \code{events}).
#' @param id_var,time_var,status_var Character scalars naming the
#'   subject-ID, event/censoring-time, and status columns of
#'   \code{events}. Default \code{"ID"}/\code{"sim_time"}/\code{"sim_status"}
#'   (the column names every \pkg{simtte} simulation function already
#'   uses).
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is
#'   called before \eqn{C} is drawn. Leave \code{NULL} when calling this
#'   from inside a larger already-seeded workflow (this is how
#'   \code{\link{sim_tte_ode}}'s own \code{censoring =} argument uses it,
#'   so that one \code{seed} reproduces \code{U}, any between-subject
#'   draw, and the censoring draw together, in that fixed order -- see
#'   \code{?sim_tte_ode} "Reproducibility").
#'
#' @return \code{events} with \code{time_var}/\code{status_var} updated
#'   in place, plus one new column, \code{sim_reason}: \code{"event"},
#'   \code{"censored"} (the independent draw, or a pre-existing censoring
#'   time from \code{events} itself, was earlier than both the event and
#'   \code{end}), or \code{"administrative"} (the observed time equals
#'   \code{end} and the subject is not an event). Every other column of
#'   \code{events} is returned unchanged.
#'
#' @seealso \code{\link{censoring_rate_for}}, to choose \code{rate}/
#'   \code{scale} for a target censoring fraction; the "Right censoring"
#'   section of \code{vignette("pkpd-time-to-event", package = "simtte")}
#'   for a worked example, including the scale-helper workflow.
#' @export
#' @examples
#' \donttest{
#' sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
#'   n = 100, end = 30, delta = 2, seed = 1)
#' out <- add_censoring(sim$events, censoring = list(dist = "exponential",
#'   rate = 0.02), end = 30, seed = 2)
#' table(out$sim_reason)
#'
#' # Also works directly on the legacy sim_tte() output.
#' lp <- matrix(rnorm(50, 0, 0.5), nrow = 50)
#' legacy <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
#'   time = seq(0, 30, by = 1), type = "weibull", end_time = 30)
#' out2 <- add_censoring(legacy, censoring = list(dist = "uniform",
#'   min = 5, max = 30), end = 30, seed = 3)
#' table(out2$sim_reason)
#' }
add_censoring <- function(events, censoring, end, id_var = "ID",
    time_var = "sim_time", status_var = "sim_status", seed = NULL) {
    events <- as.data.frame(events)
    for (v in c(id_var, time_var, status_var)) {
        if (!v %in% names(events)) {
            stop("Column '", v, "' not found in 'events'.", call. = FALSE)
        }
    }
    if (!.is_finite_scalar(end) || end < 0) {
        stop("'end' must be a non-negative finite numeric scalar.",
            call. = FALSE)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n <- nrow(events)
    c_times <- .draw_censoring_times(censoring, n)

    orig_time <- events[[time_var]]
    orig_status <- events[[status_var]]
    if (!is.numeric(orig_time) || any(!is.finite(orig_time))) {
        stop("'events[[time_var]]' must be finite numeric.", call. = FALSE)
    }

    new_time <- pmin(orig_time, c_times, end)
    # sim_status = 1 only if the *original* event time is the minimum
    # of the three (ties -- T_i == C_i or T_i == end exactly -- favor
    # the event, reports/16_censoring_design.md option A).
    event_survives <- (orig_status == 1L) & (orig_time <= c_times) &
        (orig_time <= end)
    new_status <- as.integer(event_survives)
    reason <- ifelse(new_status == 1L, "event",
        ifelse(new_time == end, "administrative", "censored"))

    events[[time_var]] <- new_time
    events[[status_var]] <- new_status
    events$sim_reason <- reason
    events
}

#' Solve for a censoring-distribution parameter giving a target censoring fraction
#'
#' Users think in terms of "about 20% independently censored before
#' \code{end}", not in the rate/scale units \code{\link{add_censoring}}
#' needs, and the right units depend on the model's own time scale. This
#' solves for the exponential \code{rate} (or, for a fixed Weibull
#' \code{shape}, the Weibull \code{scale}) that would give a target
#' independent-censoring fraction, treating \code{times} -- the event
#' times from an already-run, \strong{uncensored} simulation -- as an
#' empirical sample of the latent event time \eqn{T}: it solves
#' \code{mean(1 - exp(-rate * times)) == target} (exponential) or
#' \code{mean(1 - exp(-(times / scale)^shape)) == target} (Weibull) via
#' \code{\link[stats]{uniroot}}. This is the Monte Carlo estimate of
#' \eqn{P(C < T)} under independence, using the supplied sample in place
#' of \eqn{T}'s true (generally not closed-form once between-subject
#' variability or covariates are in the model) distribution --
#' \strong{approximate, and dependent on the run it is given}: a
#' different (even same-seed-different-\code{n}) uncensored run will
#' give a slightly different answer. Re-check the realized fraction on a
#' large cohort after drawing with the solved parameter; see
#' \code{vignette("pkpd-time-to-event", package = "simtte")} "Right
#' censoring" for a worked before/after check.
#'
#' Only genuine event times are informative here: an
#' administratively-censored subject's true \eqn{T_i} is unknown beyond
#' \code{end}, so it carries no information about where independent
#' censoring would actually bite. Pass \code{times = events$sim_time[events$sim_status
#' == 1]} from an uncensored run, not the full \code{events$sim_time}.
#'
#' @param times Numeric vector of event times (\code{sim_status == 1}
#'   only) from an uncensored simulation. Must be non-negative and
#'   finite.
#' @param target Numeric scalar in \eqn{(0, 1)}. The desired
#'   \eqn{P(C < T)}.
#' @param dist \code{"exponential"} (default) or \code{"weibull"}.
#' @param shape Required (positive numeric scalar), and used only, for
#'   \code{dist = "weibull"}: the Weibull shape is fixed by the caller
#'   and this function solves for \code{scale} alone.
#'
#' @return Numeric scalar: \code{rate} for \code{dist = "exponential"},
#'   \code{scale} for \code{dist = "weibull"}.
#' @seealso \code{\link{add_censoring}}.
#' @export
#' @examples
#' \donttest{
#' sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
#'   n = 500, end = 30, delta = 2, seed = 1)
#' event_times <- sim$events$sim_time[sim$events$sim_status == 1]
#' rate <- censoring_rate_for(event_times, target = 0.2)
#' out <- add_censoring(sim$events, censoring = list(dist = "exponential",
#'   rate = rate), end = 30, seed = 2)
#' mean(out$sim_reason == "censored") # close to 0.2, on the events subset
#' }
censoring_rate_for <- function(times, target,
    dist = c("exponential", "weibull"), shape = NULL) {
    dist <- match.arg(dist)
    if (!is.numeric(times) || length(times) < 1L ||
        any(!is.finite(times)) || any(times < 0)) {
        stop("'times' must be a non-empty, non-negative, finite numeric ",
            "vector (the event times from an uncensored run).",
            call. = FALSE)
    }
    if (!.is_finite_scalar(target) || target <= 0 || target >= 1) {
        stop("'target' must be a numeric scalar strictly between 0 and 1.",
            call. = FALSE)
    }
    if (dist == "exponential") {
        # (Two differently-shaped local closures named the same, 'f',
        # would trip R CMD check's "multiple local function definitions"
        # static-analysis NOTE -- named distinctly instead.)
        rate_gap <- function(rate) mean(1 - exp(-rate * times)) - target
        return(stats::uniroot(rate_gap, interval = c(1e-10, 1e6))$root)
    }
    if (!.is_positive_scalar(shape)) {
        stop("'shape' must be a positive numeric scalar for dist = ",
            "\"weibull\" (censoring_rate_for() solves for 'scale' at a ",
            "fixed 'shape').", call. = FALSE)
    }
    scale_gap <- function(scale) mean(1 - exp(-(times / scale)^shape)) - target
    stats::uniroot(scale_gap, interval = c(1e-6, 1e6 * max(times, 1)))$root
}

#' Is x a single finite numeric value?
#'
#' Shared predicate for the scalar-parameter validation repeated across
#' \code{\link{add_censoring}}/\code{\link{censoring_rate_for}}/
#' \code{\link{.draw_censoring_times}} (\code{end}, \code{target},
#' \code{rate}, \code{shape}, \code{scale}, ...). \code{NULL} is safely
#' rejected (\code{is.numeric(NULL)} is \code{FALSE}, not an error), so
#' callers need no separate \code{is.null()} guard.
#'
#' @param x Value to check.
#' @return Logical scalar.
#' @noRd
.is_finite_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
}

#' Is x a single finite, strictly positive numeric value?
#' @param x Value to check.
#' @return Logical scalar.
#' @noRd
.is_positive_scalar <- function(x) {
    .is_finite_scalar(x) && x > 0
}

#' Draw n independent censoring times from a censoring spec
#'
#' Internal dispatcher shared by \code{\link{add_censoring}}. Kept
#' separate from \code{add_censoring()} itself so the spec-validation/
#' draw logic has one home regardless of caller (only one caller today).
#'
#' @param censoring The user-supplied spec (see \code{?add_censoring}).
#' @param n Integer. Number of draws.
#' @return Numeric vector of length \code{n}, all finite and
#'   non-negative.
#' @noRd
.draw_censoring_times <- function(censoring, n) {
    if (!is.list(censoring) || is.null(censoring$dist)) {
        stop("'censoring' must be a list with a 'dist' element (one of ",
            "\"exponential\", \"weibull\", \"uniform\", \"user\"); see ",
            "?add_censoring.", call. = FALSE)
    }
    dist <- censoring$dist
    times <- switch(dist,
        exponential = {
            rate <- censoring$rate
            if (!.is_positive_scalar(rate)) {
                stop("censoring$rate must be a positive numeric scalar ",
                    "for dist = \"exponential\".", call. = FALSE)
            }
            stats::rexp(n, rate = rate)
        },
        weibull = {
            shp <- censoring$shape
            scl <- censoring$scale
            if (!.is_positive_scalar(shp) || !.is_positive_scalar(scl)) {
                stop("censoring$shape and censoring$scale must be positive ",
                    "numeric scalars for dist = \"weibull\".", call. = FALSE)
            }
            stats::rweibull(n, shape = shp, scale = scl)
        },
        uniform = {
            mn <- censoring$min
            mx <- censoring$max
            if (!.is_finite_scalar(mn) || !.is_finite_scalar(mx) ||
                mn < 0 || mx <= mn) {
                stop("censoring$min and censoring$max must be numeric ",
                    "scalars with 0 <= min < max for dist = \"uniform\".",
                    call. = FALSE)
            }
            stats::runif(n, min = mn, max = mx)
        },
        user = {
            fun <- censoring$fun
            if (is.null(fun) || !is.function(fun)) {
                stop("censoring$fun must be a function(n) for dist = ",
                    "\"user\".", call. = FALSE)
            }
            out <- fun(n)
            if (!is.numeric(out) || length(out) != n) {
                stop("censoring$fun(n) must return a numeric vector of ",
                    "length n (", n, "); got length ", length(out), ".",
                    call. = FALSE)
            }
            out
        },
        stop("censoring$dist must be one of \"exponential\", \"weibull\", ",
            "\"uniform\", \"user\" (got \"", dist, "\").", call. = FALSE)
    )
    if (any(!is.finite(times)) || any(times < 0)) {
        stop("The censoring-time draw contained a non-finite or negative ",
            "value; check 'censoring' (a dist = \"user\" function is the ",
            "most likely cause).", call. = FALSE)
    }
    times
}
