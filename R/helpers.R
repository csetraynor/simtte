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

cfile_dir <- function() {
  file.path(path.package("simtte"), "models")
}

# special function to read the static cache build.
read_model_static_cache <- function(model) {
    pkg_model_file <- cfile_dir()
    if (model == "weibull" | model == "ms") {
      mod_surv <- mrgsolve::mread_cache(model = model, project = pkg_model_file)
    } else {
        stop("Model ", model, " must be 'ms' or 'weibull' ")
    }
    return(mod_surv)
}
