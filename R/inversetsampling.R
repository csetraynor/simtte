#' Inverse transform sampling
#' Internal function to coduct inverse transform sampling
#' specific to use with the simulated survival datasets in \code{sim_surv_df}.
#' @param simdat simulated dataset.
#' @param covs covariate dataset.
#' @param id ID of the patient must be named id if no id_var is passed
#' @param id_var optional id var if id is not named id
simulate_survival <- function(simdat, id, covs, id_var) {
    if (missing(id_var)) {
        id_var <- "ID"
    }
    newdat <- simulate_survival_id(simdat, id, id_var)
    dplyr::left_join(newdat, covs, by = id_var)
}

#' Simulate survival internal
#' Workhorse to conduct inverse transform sampling.
#' @param simdat simulated dataset.
#' @param id ID of the patient must be named id if no id_var is passed
#' @param id_var optional id var if id is not named id
#' @importFrom stats runif
simulate_survival_id <- function(simdat, id, id_var) {
    u <- runif(1)  ## sample u for inverse transform sampling
    p <- simdat$p11
    etime = get_tte(u, p)
    if (etime != -99) {
        eventtime = get_time(simdat, etime)
        outdata <-  dplyr::tibble(time = eventtime, status = 1) %>% dplyr::mutate(ID = id)
    } else {
        eventtime = get_max_time(simdat)
        outdata <-  dplyr::tibble(time = eventtime, status = 0) %>% dplyr::mutate(ID = id)
    }
    return(outdata)
}
