#' Simulate time-to-event
#' Main function of \code{simtte}. Simulates a time-to-event dataset,
#' also known as survival dataset. It can simulate parametric models Weibull and exponential
#' and splines models.
#' @param pi prognostic index in the relative hazard for each individual.
#' Given a matrix of measurements X and a vector of hazard coefficients b,
#' the \code{pi} is given by the matrix multiplication.
#' assuming a linear hazard model. \code{pi} is also known as linear predictor.
#' @param log_pi Is the prognostic index in the log-scale defulat \code{TRUE}.
#' @param mu The intercept of the model for Weibull and ex
#' @param coefs For M-splines and B-slines model the coefficient of each spline.
#' @param basis the basis matrix for splines models.
#' @param time For M-splines the times corresponding to the basis matrix, ignored for the Weibull.
#' @param end_time Administrative censoring time, or \code{max(times)}.
#' @param type Either type Weibull or M-splines.
#' @param ... Parameters passed to mrgsolve mrgsim function.
#' @export
#' @import dplyr
#' @examples
#' \donttest{
#' data("ms_data")
#' mu <- ms_data$mu
#' basis <- ms_data$basis
#' coefs <- ms_data$coefs
#' time <- ms_data$time
#' lp <- matrix(runif(nrow(basis)),  nrow = nrow(basis))
#' wei_sim <- sim_tte(pi = lp, mu = -1, coefs = 1.1, time = time, type = "weibull", obs.only = F, obs.aug = T, delta = 0.5, end_time = 100)
#' ms_sim <- sim_tte(pi = lp, mu = mu, basis = basis, coefs = coefs, time = time, type = "ms")
#' }
sim_tte <- function(pi, log_pi = TRUE, mu = -3, coefs = 0, basis = NULL,
    time = 100, end_time, type = "weibull", ...) {
    ID <- NULL
    if (type == "ms") {
        if (ncol(basis) != length(coefs)) {
            stop("Basis columns ans coefficients have to have same length.")
        }
        basehaz = basis %*% coefs
    }
    if (type == "weibull") {
        if (length(coefs) > 1) {
            stop("Only 1 coefficient (shape) allowed for the Weibull model.")
        }
        shape <- coefs
    }
    if (!log_pi) {
        if (any(log_pi) <= 0) {
            stop("All prognostic index must be positive to transfor to log scale")
        }
        pi <- log(pi)
    }
    data_sim = sim_surv_df(log_hr = pi, mu = mu, basehaz = basehaz,
        type = type, shape = shape, times = time, ...)
    data_splited_id_sim <- data_sim %>% dplyr::group_by(ID) %>%
        dplyr::group_split()
    pi <- as.numeric(pi)
    xdata <- data.frame(ID = seq_along(pi), lp = pi)
    covs <- xdata %>% dplyr::group_by(ID) %>% dplyr::group_split()
    pred_surv <- mapply(FUN = simulate_survival, simdat = data_splited_id_sim,
        id = seq_along(pi), covs = covs, MoreArgs = list(id_var = "ID"),
        SIMPLIFY = FALSE)
    pred_surv <- pred_surv %>% dplyr::bind_rows()
    colnames(pred_surv)[1:2] <- c("sim_time", "sim_status")
    dat <- dplyr::inner_join(pred_surv, xdata, by = c("ID", "lp"))
    return(dat)
}

#' Internal workhorse function to simulate the time-to-event dataframe
#' This is an internal function and is not to be called directly.
#' @param log_hr the prognostic index.
#' @param mu intercept parameter
#' @param shape shape parameter
#' @param type type of model Weibull, splines.
#' @param times times to simulate
#' @param basehaz basehazard calculated.
#' @param end_time the censoring time.
#' @param ... passed to mrgsolve mrgsim.
#' @import mrgsolve
#' @importFrom dplyr %>%
sim_surv_df <- function(log_hr, mu, shape, type, times, basehaz,
    end_time, ...) {
    if (missing(end_time)) {
        end_time = max(times)
    }
    mod_surv <- read_model_static_cache(type)
    if (type == "weibull") {
        ev1 <- as.data.frame(expand.grid(shape, log_hr))
        colnames(ev1) <- c("shape", "lp")
        data_surv <- ev1 %>% dplyr::mutate(ID = seq_len(nrow(ev1)), cmt = 1,
            amt = 0, evid = 1, time = 0, mu = mu, basehaz_id = as.numeric(as.factor(shape)))
    } else {
        ## type m splines
        basehaz_id = rep(seq_len(ncol(basehaz)), each = nrow(basehaz))
        ev1 <- expand.grid(c(basehaz), log_hr) %>% as.data.frame()
        ev1$basehaz_id <- rep(basehaz_id, length(log_hr))
        colnames(ev1) <- c("basehaz", "lp", "basehaz_id")
        ev1$time <- c(rep(times, ncol(basehaz) * length(log_hr)))
        ev1$ID <- c(rep(seq_len(ncol(basehaz) * length(log_hr)),
            each = nrow(basehaz)))
        data_surv <- ev1 %>% dplyr::mutate(cmt = 1, amt = 0, evid = 1,
            mu = mu)
    }
    mod_surv %>% data_set(data_surv) %>% mrgsim(carry.out = c("basehaz_id", "lp", "mu", "shape",
            "basehaz"), end = end_time, ...) %>% as.data.frame()
}


#' Study relation fo hazard ratio and t50
#' The function allows calculating the difference in survival between
#' different log hazard ratios. Useful for understanding distributions
#' in a quantile non-classical way.
#' @param q quantile default to 0.5
#' @param log_hr the prognostic index.
#' @param mu intercept parameter
#' @param shape shape parameter
#' @param type type of model Weibull, splines.
#' @param times times to simulate
#' @param basehaz basehazard calculated.
#' @param end_time the censoring time.
#' @param ... passed to mrgsolve mrgsim.
#' @export
study_hr_t50_surv <- function(q, log_hr, mu, shape, type, times, basehaz,
                              end_time, ...) {
    lp <- p11 <- basehaz_id <- survrel50 <- time <- p11_baseline <-  NULL
    sim_surv <- sim_surv_df(log_hr, mu, shape, type, times, basehaz,
                            end_time, ...)
    t50_baseline <- sim_surv %>% filter(lp == 0) %>% mutate(survrel50 = abs(p11 -
        q)) %>% group_by(basehaz_id) %>% filter(survrel50 ==
        min(survrel50)) %>% ungroup() %>% select(time, basehaz_id) %>%
        distinct()
    ratio_treat <- inner_join(sim_surv, t50_baseline, by = c("time",
        "basehaz_id")) %>% distinct()
    ratio_baseline <- ratio_treat %>% filter(lp == 0) %>% select(basehaz_id,
        p11_baseline = p11)
    left_join(ratio_treat, ratio_baseline, by = "basehaz_id") %>%
        mutate(survdiff_t50base = p11 - p11_baseline) %>% distinct()
}
