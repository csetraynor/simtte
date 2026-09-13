#' Inverse Transform Sampling for a Single Subject
#'
#' Internal function that performs inverse transform sampling on
#' survival probabilities to generate an event time for a single
#' subject.
#'
#' @param simdat Data frame of simulated survival probabilities
#'   for one subject.
#' @param covs Data frame of covariates for one subject.
#' @param id Subject identifier.
#' @param id_var Character string. Name of the ID variable.
#'
#' @return A data frame with event time and status, merged with
#'   covariate data.
#' @noRd
#' @importFrom stats runif
.simulate_survival <- function(simdat, id, covs, id_var) {
    if (missing(id_var)) {
        id_var <- "ID"
    }
    newdat <- .simulate_survival_id(simdat, id, id_var)
    dplyr::left_join(newdat, covs, by = id_var)
}

#' Workhorse for Inverse Transform Sampling
#'
#' Draws a uniform random variate and finds the first reported time at
#' which the survival probability is at or below that value. Validates
#' the trajectory (ordering and monotonicity) before sampling; see
#' \code{\link{sim_tte_df}} for the full custom-trajectory contract.
#'
#' For \code{event_time_method = "grid"}, event-time resolution is
#' grid-based, not interpolated: the returned time is always one of the
#' times present in \code{simdat}, namely the first one at which
#' \code{p11 <= U}. For \code{event_time_method = "log_survival"}, a
#' crossing that occurs strictly after the first reported observation is
#' refined by linear interpolation in cumulative hazard between the two
#' reported points surrounding it (see \code{\link{.interpolate_log_survival}}
#' and the "Event-time interpolation" section of \code{\link{sim_tte_df}});
#' censoring, and a crossing already present at the very first
#' observation (no earlier point to interpolate from), are identical
#' between the two methods. If the trajectory never crosses \code{U}, the
#' subject is censored at its last reported time
#' (\code{simdat[["time"]][nrow(simdat)]}) — i.e. at that subject's own
#' final observation, not necessarily a common study cutoff shared
#' across subjects (see \code{\link{sim_tte_df}}).
#'
#' Exactly one \code{stats::runif(1)} is drawn per subject regardless of
#' \code{event_time_method}.
#'
#' @param simdat Data frame with columns \code{ID}, \code{time}, and
#'   \code{p11} (canonicalized survival probabilities) for one subject,
#'   already in ascending time order.
#' @param id Subject identifier.
#' @param id_var Character string. Name of the ID variable.
#' @param event_time_method Character string. \code{"grid"} (default) or
#'   \code{"log_survival"}. Not validated here (callers are expected to
#'   have already resolved it via \code{match.arg()}).
#'
#' @return A tibble with columns \code{time}, \code{status}, and
#'   the ID variable.
#' @noRd
#' @importFrom stats runif
.simulate_survival_id <- function(simdat, id, id_var,
    event_time_method = "grid") {
    .validate_survival_trajectory(id, simdat[["time"]], simdat[["p11"]])
    u <- stats::runif(1)
    p <- simdat[["p11"]]
    etime <- .get_tte(u, p)
    if (etime != -99) {
        if (event_time_method == "log_survival" && etime > 1L) {
            # Interpolate only within the interval .get_tte() already
            # identified as the crossing interval -- never an
            # independent search. This is what guarantees "grid" and
            # "log_survival" classify every subject identically, and
            # that a flat survival segment elsewhere in the trajectory
            # can never produce a zero-denominator: by construction,
            # simdat[["p11"]][etime - 1] > u >= simdat[["p11"]][etime].
            i <- etime - 1L
            eventtime <- .interpolate_log_survival(
                t_i = simdat[["time"]][i], t_ip1 = simdat[["time"]][etime],
                s_i = p[i], s_ip1 = p[etime], u = u)
        } else {
            eventtime <- .get_time(simdat, etime)
        }
        outdata <- dplyr::tibble(time = eventtime, status = 1,
            ID = id)
    } else {
        eventtime <- .get_max_time(simdat)
        outdata <- dplyr::tibble(time = eventtime, status = 0,
            ID = id)
    }
    outdata
}

#' Interpolate an event time linearly in cumulative hazard
#'
#' Given a crossing interval already identified by \code{\link{.get_tte}}
#' (\code{s_i > u >= s_ip1}, with \code{t_i < t_ip1} adjacent reported
#' times), computes an event time by linear interpolation of cumulative
#' hazard \code{H = -log(S)} between the two endpoints. This is
#' equivalent to assuming the hazard is constant over \code{[t_i, t_ip1]}
#' (see the "Event-time interpolation" section of
#' \code{\link{sim_tte_df}}).
#'
#' \code{s_ip1 == 0} is special-cased: \code{H(t_ip1) = -log(0) = Inf},
#' and a finite target cumulative hazard is only reached in the limit as
#' \code{t -> t_ip1}, so the event time is \code{t_ip1} exactly. Naively
#' evaluating the general formula here would silently compute
#' \code{finite / Inf == 0} (fraction 0, i.e. \code{t_event = t_i}),
#' which is both mathematically wrong and inconsistent with the grid
#' method (which reports \code{t_ip1} whenever survival first reaches
#' 0).
#'
#' The resulting interpolation fraction is clamped to \code{[0, 1]}
#' before use: \code{log()}/\code{exp()} round-trips are not guaranteed
#' exact to the last bit, so an algebraically-boundary fraction (0 or 1,
#' i.e. \code{u} exactly equal to \code{s_i} or \code{s_ip1}) could
#' otherwise evaluate a hair outside \code{[0, 1]}. The clamp guarantees
#' the returned time is always within \code{[t_i, t_ip1]}, and in
#' particular never after the subject's final reported time.
#'
#' @param t_i,t_ip1 Numeric scalars. Times bounding the crossing
#'   interval.
#' @param s_i,s_ip1 Numeric scalars. Survival probabilities at
#'   \code{t_i}/\code{t_ip1}; \code{s_i > u >= s_ip1} by construction of
#'   the caller.
#' @param u Numeric scalar. The uniform draw that triggered this
#'   crossing.
#' @return Numeric scalar. The interpolated event time, always in
#'   \code{[t_i, t_ip1]}.
#' @noRd
.interpolate_log_survival <- function(t_i, t_ip1, s_i, s_ip1, u) {
    if (s_ip1 == 0) {
        return(t_ip1)
    }
    H_i <- -log(s_i)
    H_ip1 <- -log(s_ip1)
    H_u <- -log(u)
    fraction <- (H_u - H_i) / (H_ip1 - H_i)
    fraction <- min(max(fraction, 0), 1)
    t_i + fraction * (t_ip1 - t_i)
}
