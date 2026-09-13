#' Map simulated event/censoring times onto an assessment-time interval
#'
#' Adds interval-censoring bounds to an already-simulated (and, if
#' applicable, already right-censored -- see "Ordering" below) events
#' data frame: given a
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
#' An event with \code{sim_time_left == 0} (an event before the first
#' post-baseline visit) \strong{is} a left-censored observation -- the
#' event is known to have occurred before \code{sim_time_right}, exact
#' time unknown -- with no separate mechanism needed; see the vignette's
#' "Left censoring" subsection for a worked example.
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
#' \code{add_interval_censoring()}'s "must start at 0" requirement.
#'
#' \code{jitter} is validated against \code{every} \emph{before} any
#' draw is made: the maximum possible deviation (\code{jitter} for
#' \code{"uniform"}, \code{jitter_trunc * jitter} for \code{"normal"})
#' must be less than \code{every / 2}, which guarantees -- by
#' construction, not by a post-hoc check -- that jittered visits can
#' never cross order or collide (adjacent nominal visits are
#' \code{every} apart; two deviations each smaller than half of that
#' cannot close the gap) or go negative (the first non-baseline visit
#' is at least \code{every / 2} above \code{0}). A \code{jitter}
#' violating this bound is an error, naming the bound, rather than a
#' silently reordered or clamped schedule.
#'
#' @param n Integer. Number of subjects; schedules are generated for
#'   \code{ID = 1:n} (remap the \code{ID} column afterward if your
#'   \code{events} uses different subject identifiers).
#' @param every Positive numeric scalar. Nominal spacing between visits.
#' @param end Non-negative numeric scalar. Last nominal visit time (the
#'   schedule is \code{seq(0, end, by = every)} before jitter).
#' @param jitter Non-negative numeric scalar (default \code{0}, no
#'   jitter). Each non-baseline visit is perturbed by an independent
#'   draw from \code{jitter_dist}: \code{jitter} is the draw's
#'   half-width for \code{"uniform"} (\code{stats::runif(1, -jitter,
#'   jitter)}) or its standard deviation for \code{"normal"}, which is
#'   truncated at \code{+/- jitter_trunc * jitter} (see
#'   \code{jitter_trunc}) and drawn via inverse-CDF
#'   (\code{stats::qnorm()} on a \code{stats::runif()} restricted to the
#'   truncated range) -- one uniform draw per visit either way, so both
#'   distributions consume the seeded stream the same shape of way.
#' @param jitter_dist \code{"uniform"} (default; unchanged behavior and
#'   seeds from before this argument existed) or \code{"normal"}
#'   (truncated, see \code{jitter_trunc}; its seeds changed when the
#'   truncation/bound was added -- see \code{NEWS.md}).
#' @param jitter_trunc Positive numeric scalar, default \code{2}. Only
#'   used for \code{jitter_dist = "normal"}: the truncation half-width,
#'   in standard deviations, of the truncated normal jitter draw.
#'   Ignored for \code{"uniform"}.
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
#' @seealso \code{\link{add_interval_censoring}}, \code{\link{thin_visits}}
#'   and \code{\link{visit_schedule_informative}}, to layer missed
#'   visits/dropout on top of the schedule this function returns.
#' @export
#' @examples
#' visit_schedule(n = 5, every = 4, end = 20, jitter = 1, seed = 1)
#' visit_schedule(n = 5, every = 4, end = 20, jitter = 0.4,
#'   jitter_dist = "normal", seed = 1)
visit_schedule <- function(n, every, end, jitter = 0,
    jitter_dist = c("uniform", "normal"), jitter_trunc = 2, seed = NULL) {
    jitter_dist <- match.arg(jitter_dist)
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
    if (!.is_positive_scalar(jitter_trunc)) {
        stop("'jitter_trunc' must be a positive numeric scalar.",
            call. = FALSE)
    }
    if (jitter > 0) {
        max_dev <- if (jitter_dist == "uniform") jitter else jitter_trunc * jitter
        if (max_dev >= every / 2) {
            stop("'jitter' (", jitter, ") is too large relative to ",
                "'every' (", every, ") for jitter_dist = \"", jitter_dist,
                "\": the maximum possible deviation (", max_dev, ") must ",
                "be less than every / 2 (", every / 2, ") so that visits ",
                "cannot cross order. Use a smaller 'jitter'",
                if (jitter_dist == "normal") " or 'jitter_trunc'" else "",
                ".", call. = FALSE)
        }
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    base_visits <- seq(0, end, by = every)
    k <- length(base_visits)
    rows <- lapply(seq_len(n), function(i) {
        t <- base_visits
        if (jitter > 0 && k > 1L) {
            draws <- if (jitter_dist == "uniform") {
                stats::runif(k - 1L, -jitter, jitter)
            } else {
                # Truncated normal via inverse-CDF: one uniform draw per
                # visit, restricted to the CDF range of [-jitter_trunc,
                # jitter_trunc] standard deviations, then mapped through
                # qnorm() -- deterministic given that single draw, which
                # keeps seed stability the same shape as the uniform
                # branch above (reports/23_phase6c_report.md).
                lo <- stats::pnorm(-jitter_trunc)
                hi <- stats::pnorm(jitter_trunc)
                stats::qnorm(stats::runif(k - 1L, lo, hi), sd = jitter)
            }
            jittered <- t[-1] + draws
            if (is.unsorted(jittered, strictly = TRUE) || any(jittered <= 0)) {
                # Unreachable for either built-in distribution once the
                # jitter/every bound above holds (see this function's
                # details) -- an internal invariant, not a user-facing
                # condition. If this ever fires, the bound check has a
                # bug; it is not meant to be worked around by the caller.
                stop("visit_schedule(): internal invariant violated -- ",
                    "jittered visits were not strictly increasing and ",
                    "positive despite passing the jitter/every bound ",
                    "check. Please report this as a bug.", call. = FALSE)
            }
            t <- c(0, jittered)
        }
        data.frame(ID = i, time = t)
    })
    dplyr::bind_rows(rows)
}

#' Thin a visit schedule with missed visits and dropout (C1)
#'
#' Layers schedule-only randomness on top of a fixed visit schedule --
#' no simulated outcome is read, so this is the "safe", always-
#' non-informative member of the visit-process generators (contrast
#' with \code{\link{visit_schedule_informative}}, which reads the
#' simulated outcome by design). Each non-baseline visit is independently missed
#' with probability \code{p_miss}; independently of that, each subject
#' may also start dropping out: with probability \code{p_dropout}, one
#' of their non-baseline visits is chosen uniformly at random as the
#' dropout onset, and it and every later visit is dropped. The baseline
#' visit (\code{time = 0}) is never dropped by either mechanism.
#'
#' @param visits A per-subject schedule (a data frame with \code{ID}/
#'   \code{time} columns, e.g. from \code{\link{visit_schedule}}), or a
#'   common numeric vector -- in which case \code{ids} is required, and
#'   the vector is expanded to one copy per ID first (via the same
#'   internal validation \code{\link{add_interval_censoring}} itself
#'   uses).
#' @param p_miss,p_dropout Numeric scalars in \code{[0, 1]}, default
#'   \code{0} (no thinning; \code{visits} is returned as given). The
#'   dropout-onset visit is drawn from the subject's \emph{original}
#'   non-baseline visits, independently of any \code{p_miss} draws for
#'   that subject.
#' @param ids Required only when \code{visits} is a numeric vector;
#'   the subject IDs to expand it across (e.g. \code{events$ID}).
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is
#'   called before any draw. Leave \code{NULL} inside a larger
#'   already-seeded workflow.
#'
#' @return A data frame with \code{ID}/\code{time} columns -- a
#'   (possibly shorter, per subject) schedule directly usable as
#'   \code{\link{add_interval_censoring}}'s \code{visits} argument. Must
#'   be applied before \code{\link{add_interval_censoring}}, and
#'   composes with \code{\link{add_censoring}} in either relative order
#'   (they touch disjoint objects: this function only ever reads/writes
#'   a visit schedule, never \code{events}).
#' @seealso \code{\link{visit_schedule}}, \code{\link{visit_schedule_informative}}
#'   (C2, the outcome-reactive counterpart), \code{\link{add_interval_censoring}}.
#' @export
#' @examples
#' sched <- visit_schedule(n = 20, every = 4, end = 24, seed = 1)
#' thinned <- thin_visits(sched, p_miss = 0.1, p_dropout = 0.05, seed = 2)
#' nrow(thinned) <= nrow(sched)
thin_visits <- function(visits, p_miss = 0, p_dropout = 0, ids = NULL,
    seed = NULL) {
    if (!.is_probability_scalar(p_miss)) {
        stop("'p_miss' must be a numeric scalar in [0, 1].", call. = FALSE)
    }
    if (!.is_probability_scalar(p_dropout)) {
        stop("'p_dropout' must be a numeric scalar in [0, 1].", call. = FALSE)
    }
    if (is.numeric(visits)) {
        if (is.null(ids)) {
            stop("'ids' is required when 'visits' is a numeric vector ",
                "(the common schedule is expanded per subject); see ",
                "?thin_visits.", call. = FALSE)
        }
        sched_list <- .normalize_visit_schedule(visits, ids)
        visits <- dplyr::bind_rows(lapply(seq_along(ids), function(i) {
            data.frame(ID = ids[i], time = sched_list[[as.character(ids[i])]])
        }))
    } else if (!is.data.frame(visits) ||
        !all(c("ID", "time") %in% names(visits))) {
        stop("'visits' must be a numeric vector (a common schedule, ",
            "with 'ids') or a data frame with 'ID'/'time' columns (a ",
            "per-subject schedule).", call. = FALSE)
    }
    if (p_miss == 0 && p_dropout == 0) {
        return(visits)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    subj_ids <- unique(visits$ID)
    rows <- lapply(subj_ids, function(id) {
        v <- visits$time[visits$ID == id]
        .validate_visit_vector(v)
        non_baseline <- seq_along(v)[-1]
        keep <- rep(TRUE, length(v))
        if (length(non_baseline)) {
            if (p_miss > 0) {
                miss <- stats::runif(length(non_baseline)) < p_miss
                keep[non_baseline[miss]] <- FALSE
            }
            if (p_dropout > 0 && stats::runif(1) < p_dropout) {
                onset_pos <- floor(stats::runif(1) * length(non_baseline)) + 1L
                onset <- non_baseline[onset_pos]
                keep[non_baseline[non_baseline >= onset]] <- FALSE
            }
        }
        data.frame(ID = id, time = v[keep])
    })
    dplyr::bind_rows(rows)
}

#' Generate an outcome-informative visit schedule (C2)
#'
#' Layers visit-process randomness that \strong{reacts to the simulated
#' event time} on top of a fixed visit schedule: a scheduled
#' visit shortly before a subject's event may be missed at an elevated
#' rate (\code{miss_near_event}, e.g. "too unwell to attend"), and/or an
#' unscheduled extra visit may be added shortly after it
#' (\code{extra_visit_after_event}, e.g. "brought in for follow-up after
#' a reported symptom"). Both are optional and independent of each
#' other; supplying neither is equivalent to \code{\link{thin_visits}}
#' with \code{p_miss = p_miss_base}.
#'
#' @section Informative assessment:
#' \strong{This function's output is not a neutral assessment schedule.}
#' By construction, whether/when a subject is assessed now depends on
#' whether/when they had the simulated event -- exactly the situation
#' most interval-censored analyses (Turnbull NPMLE, standard parametric
#' interval-censored likelihoods) assume does \emph{not} hold when they
#' treat assessment times as ignorable. Using this function's schedule
#' is appropriate for \emph{studying} that bias (e.g. comparing a fit
#' against \code{\link{thin_visits}}'s non-informative schedule at the
#' same nominal miss rate, as the vignette's "Interval censoring"
#' section does), not as a drop-in, general-purpose visit generator. A
#' \code{message()} naming this is emitted on every call.
#'
#' @param visits A per-subject schedule (data frame with \code{ID}/
#'   \code{time}), or a common numeric vector (expanded per subject
#'   using \code{events[[id_var]]} -- no separate \code{ids} argument
#'   needed here, unlike \code{\link{thin_visits}}, since \code{events}
#'   already supplies the ID list).
#' @param events Data frame with (by default) \code{ID}, \code{sim_time},
#'   \code{sim_status} columns, \strong{after} any \code{\link{add_censoring}}
#'   call (see "Pipeline position" below) -- column names configurable
#'   via \code{id_var}/\code{time_var}/\code{status_var}.
#' @param end Numeric scalar. An unscheduled extra visit
#'   (\code{extra_visit_after_event}) landing after \code{end} is
#'   dropped, not added -- the same administrative horizon
#'   \code{\link{add_censoring}} uses.
#' @param miss_near_event \code{NULL} (default, no effect), or
#'   \code{list(window = <positive scalar>, p = <scalar in [0, 1]>)}.
#'   For an event subject (\code{sim_status == 1}), a scheduled visit
#'   \code{v} with \code{t - window <= v < t} (\code{t} the event time)
#'   is missed with probability \code{1 - (1 - p_miss_base) * (1 - p)}
#'   -- \code{p_miss_base}'s background miss chance and \code{p}'s
#'   near-event miss chance are independent, competing risks, so they
#'   add rather than one overriding the other (this collapses exactly
#'   to \code{p} when \code{p_miss_base = 0}, the common case). Has no
#'   effect on a censored subject, or on a visit outside the window --
#'   both use \code{p_miss_base} alone.
#' @param extra_visit_after_event \code{NULL} (default, no effect), or a
#'   distribution spec in \code{\link{add_censoring}}'s \code{censoring}
#'   shape (\code{list(dist = "exponential", rate = ...)},
#'   \code{"weibull"}, \code{"uniform"}, \code{"lognormal"},
#'   \code{"gamma"}, or \code{"user"}). For an event subject only, one
#'   delay is drawn and an unscheduled visit is added at
#'   \code{sim_time + delay}, unless that would exceed \code{end}.
#' @param p_miss_base Numeric scalar in \code{[0, 1]}, default \code{0}.
#'   Background miss probability applied to every non-baseline visit not
#'   otherwise covered by \code{miss_near_event} (i.e. every visit for a
#'   censored subject, and every out-of-window visit for an event
#'   subject) -- the same mechanism as \code{\link{thin_visits}}'s
#'   \code{p_miss}, folded in here so a caller does not need to chain
#'   both functions just to add a uniform background miss rate on top of
#'   the near-event effect.
#' @param id_var,time_var,status_var Character scalars naming the
#'   subject-ID, event/censoring-time, and status columns of
#'   \code{events}. Default \code{"ID"}/\code{"sim_time"}/\code{"sim_status"}.
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is
#'   called before any draw.
#'
#' @return A data frame with \code{ID}/\code{time} columns -- directly
#'   usable as \code{\link{add_interval_censoring}}'s \code{visits}
#'   argument.
#'
#' @section Pipeline position:
#' \code{events} should already reflect any independent right
#' censoring: exact simulation -> \code{\link{add_censoring}} ->
#' \code{visit_schedule_informative()} -> \code{\link{add_interval_censoring}}.
#' Generating the schedule from pre-censoring event times would make the
#' near-event/extra-visit mechanics react to a \code{sim_time} that the
#' final \code{events} no longer has -- the same ordering reasoning as
#' \code{\link{add_interval_censoring}}'s own "Ordering with right
#' censoring" section.
#'
#' @seealso \code{\link{thin_visits}} (C1, the non-informative
#'   counterpart), \code{\link{add_interval_censoring}}.
#' @export
#' @examples
#' \donttest{
#' sim <- sim_tte_ode(model = "exponential", param = list(H0 = 0.05),
#'   n = 50, end = 24, delta = 2, seed = 1)
#' sched <- visit_schedule(n = 50, every = 4, end = 24, seed = 1)
#' realized <- visit_schedule_informative(sched, sim$events, end = 24,
#'   miss_near_event = list(window = 3, p = 0.5), seed = 2)
#' out <- add_interval_censoring(sim$events, visits = realized)
#' head(out)
#' }
visit_schedule_informative <- function(visits, events, end,
    miss_near_event = NULL, extra_visit_after_event = NULL,
    p_miss_base = 0, id_var = "ID", time_var = "sim_time",
    status_var = "sim_status", seed = NULL) {
    events <- .validate_events_columns(events, id_var, time_var, status_var)
    if (!.is_finite_scalar(end) || end < 0) {
        stop("'end' must be a non-negative finite numeric scalar.",
            call. = FALSE)
    }
    if (!.is_probability_scalar(p_miss_base)) {
        stop("'p_miss_base' must be a numeric scalar in [0, 1].",
            call. = FALSE)
    }
    if (!is.null(miss_near_event)) {
        if (!is.list(miss_near_event) || is.null(miss_near_event$window) ||
            is.null(miss_near_event$p)) {
            stop("'miss_near_event' must be a list with 'window' and 'p' ",
                "elements; see ?visit_schedule_informative.", call. = FALSE)
        }
        if (!.is_positive_scalar(miss_near_event$window)) {
            stop("'miss_near_event$window' must be a positive numeric ",
                "scalar.", call. = FALSE)
        }
        if (!.is_probability_scalar(miss_near_event$p)) {
            stop("'miss_near_event$p' must be a numeric scalar in [0, 1].",
                call. = FALSE)
        }
    }

    ids <- events[[id_var]]
    sched <- .normalize_visit_schedule(visits, ids)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    times <- events[[time_var]]
    statuses <- events[[status_var]]
    rows <- mapply(function(id, t, s) {
        v <- sched[[as.character(id)]]
        keep <- rep(TRUE, length(v))
        non_baseline <- seq_along(v)[-1]
        if (length(non_baseline)) {
            p_eff <- rep(p_miss_base, length(non_baseline))
            if (isTRUE(s == 1L) && !is.null(miss_near_event)) {
                w <- miss_near_event$window
                in_window <- v[non_baseline] < t & v[non_baseline] >= (t - w)
                # Additive/competing-risks, not override: a visit is
                # missed if it would have been missed by either
                # mechanism (independent events) -- collapses exactly
                # to miss_near_event$p when p_miss_base = 0.
                p_eff[in_window] <- 1 - (1 - p_miss_base) *
                    (1 - miss_near_event$p)
            }
            if (any(p_eff > 0)) {
                miss <- stats::runif(length(non_baseline)) < p_eff
                keep[non_baseline[miss]] <- FALSE
            }
        }
        v_kept <- v[keep]
        if (isTRUE(s == 1L) && !is.null(extra_visit_after_event)) {
            delay <- .draw_censoring_times(extra_visit_after_event, 1)
            extra_time <- t + delay
            if (extra_time <= end) {
                v_kept <- sort(unique(c(v_kept, extra_time)))
            }
        }
        data.frame(ID = id, time = v_kept)
    }, ids, times, statuses, SIMPLIFY = FALSE)

    message("visit_schedule_informative(): the returned schedule depends ",
        "on the simulated event times; interval-censored analyses ",
        "assuming non-informative assessment times will be biased by ",
        "design (see ?visit_schedule_informative \"Informative ",
        "assessment\").")

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
