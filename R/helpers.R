## helpers functions
get_time <- function(simdat, etime) {
    unlist(simdat[etime, 2])
}

get_tte <- function(U, pcurr) {
    match(1, cumsum(pcurr %>% unlist() < U), nomatch = -99L)
}

get_max_time <- function(simdat) {
    unlist(simdat[nrow(simdat), 2])
}

# special function to read the static cache build.
read_model_static_cache <- function(model) {
    if (model == "weibull") {
      mod_surv <- suppressMessages(mrgsolve::mread_cache("solve_mrgmodels/wei_prop_haz.cpp"))
    } else if (model == "ms") {
      mod_surv <- suppressMessages(mrgsolve::mread_cache("solve_mrgmodels/ms_prop_haz.cpp"))
    } else {
        stop("Model ", model, " must be ms, weibull.")
    }
    return(mod_surv)
}
