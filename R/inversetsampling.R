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
#' Event-time resolution is grid-based, not interpolated: the returned
#' time is always one of the times present in \code{simdat}, namely the
#' first one at which \code{p11 <= U}. If the trajectory never crosses
#' \code{U}, the subject is censored at its last reported time
#' (\code{simdat[["time"]][nrow(simdat)]}) — i.e. at that subject's own
#' final observation, not necessarily a common study cutoff shared
#' across subjects (see \code{\link{sim_tte_df}}).
#'
#' @param simdat Data frame with columns \code{ID}, \code{time}, and
#'   \code{p11} (canonicalized survival probabilities) for one subject,
#'   already in ascending time order.
#' @param id Subject identifier.
#' @param id_var Character string. Name of the ID variable.
#'
#' @return A tibble with columns \code{time}, \code{status}, and
#'   the ID variable.
#' @noRd
#' @importFrom stats runif
.simulate_survival_id <- function(simdat, id, id_var) {
    .validate_survival_trajectory(id, simdat[["time"]], simdat[["p11"]])
    u <- stats::runif(1)
    p <- simdat[["p11"]]
    etime <- .get_tte(u, p)
    if (etime != -99) {
        eventtime <- .get_time(simdat, etime)
        outdata <- dplyr::tibble(time = eventtime, status = 1,
            ID = id)
    } else {
        eventtime <- .get_max_time(simdat)
        outdata <- dplyr::tibble(time = eventtime, status = 0,
            ID = id)
    }
    return(outdata)
}
