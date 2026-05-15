#' Get time at a given row index
#' @param simdat Simulation data frame.
#' @param etime Row index.
#' @return Numeric time value.
#' @noRd
.get_time <- function(simdat, etime) {
    unlist(simdat[etime, 2])
}

#' Get event time index via inverse transform sampling
#' @param U Uniform random variate.
#' @param pcurr Survival probabilities.
#' @return Integer row index or -99 if no event.
#' @noRd
.get_tte <- function(U, pcurr) {
    match(1, cumsum(pcurr %>% unlist() < U), nomatch = -99L)
}

#' Get maximum time from simulation data
#' @param simdat Simulation data frame.
#' @return Numeric time value.
#' @noRd
.get_max_time <- function(simdat) {
    unlist(simdat[nrow(simdat), 2])
}

#' Get path to installed model files
#' @return Character path.
#' @noRd
.cfile_dir <- function() {
    file.path(path.package("simtte"), "models")
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
