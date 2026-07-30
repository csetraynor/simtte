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
#' Draws a uniform random variate and finds the time at which the
#' survival probability crosses below that value.
#'
#' @param simdat Data frame with columns \code{ID}, \code{time}, and
#'   \code{p11} (canonicalized survival probabilities) for one subject.
#' @param id Subject identifier.
#' @param id_var Character string. Name of the ID variable.
#'
#' @return A tibble with columns \code{time}, \code{status}, and
#'   the ID variable.
#' @noRd
#' @importFrom stats runif
.simulate_survival_id <- function(simdat, id, id_var) {
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
