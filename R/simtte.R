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
#' data("ms_data")
#' mu <- ms_data$mu
#' basis <- ms_data$basis
#' coefs <- ms_data$coefs
#' time <- ms_data$time
#' lp <- matrix(runif(nrow(basis)),  nrow = nrow(basis))
#' wei_sim <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
#'  time = time, type = "weibull", obs.only = FALSE,
#'  obs.aug = TRUE, delta = 0.5, end_time = 100)
#' ms_sim <- sim_tte(pi = lp, mu = mu, basis = basis,
#'  coefs = coefs, time = time, type = "ms")
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
    dat <- sim_tte_df(data_sim, id_var = "ID", xdata = xdata)
    return(dat)
}

#' Simulate time-to-event dataframe
#' This function can be used to use inverse transform sampling from a
#' cosutmised time-to-event model implemented in **mrgsolve**.
#' One only needs to pass the dataframe and the state variable
#' the probability to remain alive, by default p11.
#' @param dat The mrgsim output dataframe.
#' @param surv_var the state probability to remain with no event , default to p11.
#' @param id_var ID vaiable default to ID.
#' @param xdata additional covariate data to pass to the final dataframe.
#'  Must have the same \code{id_var} as \code{dat}.
#'  @export
sim_tte_df <- function(dat,
                       surv_var = "p11",
                       id_var = "ID",
                       xdata = NULL){
  ID <- NULL
  dat <- as.data.frame(dat)
  if(surv_var != "p11"){
    dat$p11 <- dat[[surv_var]]
  }
  if(id_var != "ID") {
    dat$ID <- dat[[id_var]]
  }
  if(is.null(xdata)){
    xdata <- data.frame(ID = seq_along(unique(dat$ID)))
  } else {
    xdata$ID <- xdata[[id_var]]
  }
  covs <- xdata %>% dplyr::group_by(ID) %>% dplyr::group_split()
  data_splited_id_sim <- dat %>% dplyr::group_by(ID) %>%
    dplyr::group_split()
  pred_surv <- mapply(FUN = simulate_survival, simdat = data_splited_id_sim,
                      id = seq_along(unique(dat$ID)), covs = covs, MoreArgs = list(id_var = "ID"),
                      SIMPLIFY = FALSE)
  pred_surv <- pred_surv %>% dplyr::bind_rows()
  colnames(pred_surv)[1:2] <- c("sim_time", "sim_status")
  dat <- dplyr::inner_join(pred_surv, xdata, by = c(colnames(xdata)))
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


#' Simulate relation of prognostic index and quantile times
#' The function allows calculating the difference in survival between
#' different log hazard ratios. Useful for understanding distributions
#' in a quantile non-classical way.
#' @param q quantile default to 0.5
#' @param pi the prognostic index.
#' @param mu intercept parameter
#' @param shape shape parameter
#' @param type type of model Weibull, splines.
#' @param times times to simulate
#' @param basehaz basehazard calculated.
#' @param end_time the censoring time.
#' @param ... passed to mrgsolve mrgsim.
#' @export
#' @examples
#' data_sim = explore_pi_tq_surv(pi = seq(-3, 3, by = 0.1),
#'  mu = -1,
#' shape = 1.1,
#'  end_time = 200,
#'   type = "weibull")
#' plot(survdiff_tq ~ lp, data = data_sim)
explore_pi_tq_surv <- function(q = 0.5,
                             pi,
                             mu,
                             shape,
                             type,
                             times,
                             basehaz,
                             end_time, ...) {
    lp <- p11 <- basehaz_id <- survrelq <- time <- p11_baseline <-  NULL
    pi <- c(pi)
    sim_surv <- sim_surv_df(pi,
                            mu,
                            shape,
                            type,
                            times,
                            basehaz,
                            end_time,
                            ...)
    tq_baseline <- sim_surv %>% dplyr::filter(lp == 0) %>%
      dplyr::mutate(survrelq = abs(p11 - q)) %>%
      dplyr::group_by(basehaz_id) %>%
      dplyr::filter(survrelq == min(survrelq)) %>%
      dplyr::ungroup() %>%
      dplyr::select(time, basehaz_id) %>%
      dplyr::distinct()
    ratio_treat <- dplyr::inner_join(sim_surv, tq_baseline, by = c("time",
        "basehaz_id")) %>% dplyr::distinct()
    ratio_baseline <- ratio_treat %>% dplyr::filter(lp == 0) %>%
      dplyr::select(basehaz_id, p11_baseline = p11)
    dplyr::left_join(ratio_treat, ratio_baseline, by = "basehaz_id") %>%
      dplyr::mutate(survdiff_tq = p11 - p11_baseline) %>%
      dplyr::distinct()
}
