#' Map simulated event/censoring times onto an assessment-time interval
#'
#' Adds interval-censoring bounds to an already-simulated (and, if
#' applicable, already right-censored -- see "Ordering" below) events
#' data frame, the mechanism described in
#' \code{reports/18_interval_censoring_design.md} (option A): given a
#' visit schedule, an event known to have occurred exactly at
#' \code{sim_time} is only \emph{detectable} at the first visit at or
#' after it, and a censored subject is only known event-free through
#' their last visit at or before their own censoring time. Works on the
#' output of \code{\link{sim_tte_ode}}, \code{\link{sim_tte}}, or
#' \code{\link{sim_tte_df}} alike, same as \code{\link{add_censoring}}.
#'
#' @param events Data frame with (by default) \code{ID}, \code{sim_time},
#'   \code{sim_status} columns -- column names are configurable via
#'   \code{id_var}/\code{time_var}/\code{status_var}. Every other column
#'   is left untouched, and so are these three: this function only adds
#'   \code{sim_time_left}/\code{sim_time_right} (see "Return" below).
#' @param visits A visit schedule, either:
#'   \describe{
#'     \item{a numeric vector}{applied to every subject (e.g.
#'       \code{c(0, 4, 8, 12, 24)}). Must start at \code{0} and be
#'       strictly increasing (sorted, no duplicate visit times).}
#'     \item{a data frame}{with \code{ID}/\code{time} columns, for a
#'       per-subject schedule. Every ID in \code{events} must have a
#'       schedule (starting at \code{0}, strictly increasing), and every
#'       ID in \code{visits} must be one that appears in \code{events}.}
#'   }
#'   \code{\link{visit_schedule}} builds the data-frame form from a
#'   fixed spacing plus optional per-subject jitter, rather than this
#'   function accepting a jitter spec itself -- see its own
#'   documentation.
#' @param id_var,time_var,status_var Character scalars naming the
#'   subject-ID, event/censoring-time, and status columns of
#'   \code{events}. Default \code{"ID"}/\code{"sim_time"}/\code{"sim_status"}.
#'
#' @return \code{events} with two new columns:
#' \describe{
#'   \item{sim_time_left}{Left (open) bound of the interval, \code{L}.}
#'   \item{sim_time_right}{Right (closed) bound, \code{R}, or \code{Inf}
#'     for an interval with no known upper bound (a censored subject, or
#'     an event that occurred after the subject's last visit and so was
#'     never confirmed by the schedule -- \code{sim_status} still reads
#'     \code{1} in that case: the exact columns are the simulation's
#'     ground truth, unaffected by this mapping; see "Interval
#'     conventions" in the design report for why this combination is
#'     intentional, not an inconsistency).}
#' }
#' \code{sim_time}/\code{sim_status}/\code{sim_reason} (if present) are
#' returned unchanged -- this function only adds columns, so every
#' existing consumer of those three keeps working. For
#' \code{survival::Surv(time = sim_time_left, time2 = sim_time_right,
#' type = "interval2")}, which expects \code{NA} rather than \code{Inf}
#' for an open upper bound, convert with \code{ifelse(is.infinite(
#' sim_time_right), NA, sim_time_right)} -- see the vignette's
#' "Interval censoring" section for a full worked example.
#'
#' @section Ordering with right censoring:
#' Call this \strong{after} \code{\link{add_censoring}} (or after the
#' exact simulation, if no independent right censoring was requested),
#' not before: the interval mapping reads \code{sim_time}/\code{sim_status}
#' as given, so a dropout applied afterward would leave the
#' already-added interval columns stale and inconsistent with the final
#' exact columns. \code{\link{sim_tte_ode}}'s own \code{censoring}/
#' \code{visits} convenience arguments apply them in this order
#' internally.
#'
#' @seealso \code{\link{visit_schedule}}, to build a per-subject
#'   schedule with fixed spacing and optional jitter;
#'   \code{\link{add_censoring}}, which should run first.
#' @export
#' @examples
#' \donttest{
#' sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
#'   n = 50, end = 24, delta = 2, seed = 1)
#' out <- add_interval_censoring(sim$events, visits = c(0, 4, 8, 12, 16, 20, 24))
#' head(out)
#' }
add_interval_censoring <- function(events, visits, id_var = "ID",
    time_var = "sim_time", status_var = "sim_status") {
    events <- .validate_events_columns(events, id_var, time_var, status_var)
    ids <- events[[id_var]]
    if (anyDuplicated(ids)) {
        stop("'events[[id_var]]' must not contain duplicated IDs (one ",
            "row per subject).", call. = FALSE)
    }
    times <- events[[time_var]]
    statuses <- events[[status_var]]
    if (!is.numeric(times) || any(!is.finite(times))) {
        stop("'events[[time_var]]' must be finite numeric.", call. = FALSE)
    }

    sched <- .normalize_visit_schedule(visits, ids)
    bounds <- mapply(.interval_bounds_for_subject, ids, times, statuses,
        MoreArgs = list(sched = sched))

    events$sim_time_left <- bounds[1, ]
    events$sim_time_right <- bounds[2, ]
    events
}

#' Generate a per-subject visit schedule with fixed spacing and jitter
#'
#' Builds the per-subject data-frame form of \code{visits} that
#' \code{\link{add_interval_censoring}} accepts: assessments every
#' \code{every} time units from \code{0} to \code{end}, optionally
#' perturbed per subject by up to \code{jitter} (a common pattern in
#' trial simulation -- scheduled visits rarely land exactly on the
#' nominal day). The baseline visit (\code{time = 0}) is never
#' jittered, so the result always satisfies
#' \code{add_interval_censoring()}'s "must start at 0" requirement;
#' every other visit's jitter is clamped at \code{0} (never negative)
#' and the result is re-sorted per subject, so a large \code{jitter}
#' relative to \code{every} cannot silently reorder a subject's own
#' visits.
#'
#' @param n Integer. Number of subjects; schedules are generated for
#'   \code{ID = 1:n} (remap the \code{ID} column afterward if your
#'   \code{events} uses different subject identifiers).
#' @param every Positive numeric scalar. Nominal spacing between visits.
#' @param end Non-negative numeric scalar. Last nominal visit time (the
#'   schedule is \code{seq(0, end, by = every)} before jitter).
#' @param jitter Non-negative numeric scalar (default \code{0}, no
#'   jitter). Each non-baseline visit is perturbed by an independent
#'   \code{stats::runif(1, -jitter, jitter)} draw, clamped at \code{0}.
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is
#'   called before any jitter is drawn. Leave \code{NULL} when calling
#'   this from inside a larger already-seeded workflow (this is how
#'   \code{\link{sim_tte_ode}}'s own \code{visits =} argument uses it
#'   when given a jitter spec, so one \code{seed} reproduces \code{U},
#'   any between-subject draw, censoring, and the visit jitter together).
#'
#' @return A data frame with \code{ID}/\code{time} columns, one row per
#'   subject per visit -- directly usable as
#'   \code{\link{add_interval_censoring}}'s \code{visits} argument.
#' @seealso \code{\link{add_interval_censoring}}.
#' @export
#' @examples
#' visit_schedule(n = 5, every = 4, end = 20, jitter = 1, seed = 1)
visit_schedule <- function(n, every, end, jitter = 0, seed = NULL) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
        n != round(n)) {
        stop("'n' must be a positive integer scalar.", call. = FALSE)
    }
    if (!.is_positive_scalar(every)) {
        stop("'every' must be a positive numeric scalar.", call. = FALSE)
    }
    if (!.is_finite_scalar(end) || end < 0) {
        stop("'end' must be a non-negative finite numeric scalar.",
            call. = FALSE)
    }
    if (!.is_finite_scalar(jitter) || jitter < 0) {
        stop("'jitter' must be a non-negative numeric scalar.",
            call. = FALSE)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    base_visits <- seq(0, end, by = every)
    k <- length(base_visits)
    rows <- lapply(seq_len(n), function(i) {
        t <- base_visits
        if (jitter > 0 && k > 1L) {
            t[-1] <- pmax(0, t[-1] + stats::runif(k - 1L, -jitter, jitter))
            t <- sort(unique(t))
        }
        data.frame(ID = i, time = t)
    })
    dplyr::bind_rows(rows)
}

#' Validate a single subject's visit vector
#'
#' Shared by both branches of \code{\link{.normalize_visit_schedule}}.
#'
#' @param v Numeric vector.
#' @return \code{invisible(TRUE)}; called for the validation side
#'   effect (raises an error on a malformed schedule).
#' @noRd
.validate_visit_vector <- function(v) {
    if (!is.numeric(v) || length(v) < 1L || any(!is.finite(v)) ||
        any(v < 0)) {
        stop("A visit schedule must be a non-empty, non-negative, ",
            "finite numeric vector.", call. = FALSE)
    }
    if (is.unsorted(v, strictly = TRUE)) {
        stop("A visit schedule must be strictly increasing (sorted, no ",
            "duplicate visit times).", call. = FALSE)
    }
    if (v[1] != 0) {
        stop("A visit schedule must start at 0 (the baseline ",
            "assessment).", call. = FALSE)
    }
    invisible(TRUE)
}

#' Normalize a user-supplied visits spec to one schedule per subject
#'
#' @param visits See \code{\link{add_interval_censoring}}'s \code{visits}.
#' @param ids The \code{events[[id_var]]} column (defines which
#'   subjects need a schedule, and, for the vector form, how many
#'   copies to make).
#' @return Named list, keyed by \code{as.character(ids)}, of validated
#'   numeric visit vectors, one per subject in \code{ids}.
#' @noRd
.normalize_visit_schedule <- function(visits, ids) {
    ids_chr <- as.character(ids)
    if (is.data.frame(visits)) {
        if (!all(c("ID", "time") %in% names(visits))) {
            stop("'visits' must have 'ID' and 'time' columns when ",
                "supplied as a per-subject data frame.", call. = FALSE)
        }
        visits_ids <- as.character(visits$ID)
        unknown <- setdiff(unique(visits_ids), ids_chr)
        if (length(unknown)) {
            stop("'visits' contains ID(s) not present in 'events': ",
                paste(unknown, collapse = ", "), ".", call. = FALSE)
        }
        missing_sched <- setdiff(ids_chr, unique(visits_ids))
        if (length(missing_sched)) {
            stop("'visits' is missing a schedule for ID(s): ",
                paste(missing_sched, collapse = ", "), ".", call. = FALSE)
        }
        sched <- lapply(split(visits$time, visits_ids), function(v) {
            .validate_visit_vector(v)
            v
        })
        return(sched[ids_chr])
    }
    if (is.numeric(visits)) {
        .validate_visit_vector(visits)
        return(stats::setNames(rep(list(visits), length(ids_chr)), ids_chr))
    }
    stop("'visits' must be a numeric vector (a common schedule) or a ",
        "data frame with 'ID'/'time' columns (a per-subject schedule); ",
        "see ?add_interval_censoring.", call. = FALSE)
}

#' Interval bounds for one subject, given their exact outcome and schedule
#'
#' The two rules from \code{reports/18_interval_censoring_design.md}
#' "Interval conventions": an event is only detectable at the first
#' visit at or after it happened; a censored subject is only known
#' event-free through their last visit at or before their own censoring
#' time. Both default \code{L} to \code{0} when no earlier visit
#' qualifies (an event before the first post-baseline visit, or a
#' subject censored before their first post-baseline visit).
#'
#' @param id Subject ID (used only to look up \code{sched}).
#' @param t Exact \code{sim_time}.
#' @param s Exact \code{sim_status} (\code{1} = event, \code{0} =
#'   censored).
#' @param sched Named list from \code{\link{.normalize_visit_schedule}}.
#' @return Numeric vector of length 2: \code{c(L, R)}.
#' @noRd
.interval_bounds_for_subject <- function(id, t, s, sched) {
    v <- sched[[as.character(id)]]
    if (isTRUE(s == 1L)) {
        ge <- v[v >= t]
        if (length(ge)) {
            r <- min(ge)
            lt <- v[v < r]
            l <- if (length(lt)) max(lt) else 0
            return(c(l, r))
        }
        return(c(max(v), Inf))
    }
    le <- v[v <= t]
    l <- if (length(le)) max(le) else 0
    c(l, Inf)
}
