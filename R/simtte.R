#' Simulate Time-to-Event Data
#'
#' Main function of \pkg{simtte}. Simulates a time-to-event dataset
#' (survival dataset) using either a Weibull parametric model or an
#' M-spline flexible baseline hazard model. Event times are generated
#' via inverse transform sampling from the cumulative hazard function
#' computed by \pkg{mrgsolve}.
#'
#' @param pi Numeric matrix or vector. Prognostic index (linear predictor)
#'   for each individual. Given a covariate matrix X and coefficient
#'   vector b, \code{pi = X \%*\% b}.
#' @param log_pi Logical; is the prognostic index already on the log scale?
#'   Default \code{TRUE}.
#' @param mu Numeric scalar. Intercept parameter of the model. Default
#'   \code{-3}.
#' @param coefs Numeric vector. For M-spline models, the coefficients of
#'   each spline basis function. For Weibull models, the shape parameter
#'   (scalar).
#' @param basis Numeric matrix. Basis matrix for M-spline models.
#'   Ignored for Weibull models. Default \code{NULL}.
#' @param time Numeric vector. For M-spline models, the time points
#'   corresponding to rows of the basis matrix. For Weibull models,
#'   used only to determine \code{end_time} if not supplied. Default
#'   \code{100}.
#' @param end_time Numeric scalar. Administrative censoring time.
#'   Defaults to \code{max(time)}.
#' @param type Character string. Model type: \code{"weibull"} or
#'   \code{"ms"} (M-splines). Default \code{"weibull"}.
#' @param ... Additional arguments passed to
#'   \code{\link[mrgsolve]{mrgsim}}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{sim_time}{Simulated event or censoring time.}
#'   \item{sim_status}{Event indicator (1 = event, 0 = censored).}
#'   \item{ID}{Subject identifier.}
#'   \item{lp}{Log prognostic index (linear predictor).}
#' }
#'
#' @references
#' Bender R, Augustin T, Blettner M (2005). Generating survival times to
#' simulate Cox proportional hazards models. \emph{Statistics in Medicine},
#' 24(11), 1713--1723. \doi{10.1002/sim.2059}
#'
#' Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
#' and proportional-odds models for censored survival data, with application
#' to prognostic modelling and estimation of treatment effects.
#' \emph{Statistics in Medicine}, 21(15), 2175--2197.
#' \doi{10.1002/sim.1203}
#'
#' @export
#' @examples
#' # Fast Weibull example with a small dataset
#' set.seed(1)
#' lp <- matrix(rnorm(5, 0, 0.5), nrow = 5)
#' result <- sim_tte(
#'   pi = lp, mu = -1, coefs = 1.1,
#'   time = seq(0.1, 10, by = 0.5),
#'   type = "weibull", end_time = 10
#' )
#' head(result)
#' \donttest{
#' # Larger examples using bundled ms_data
#' data("ms_data")
#' mu <- ms_data$mu
#' basis <- ms_data$basis
#' coefs <- ms_data$coefs
#' time <- ms_data$time
#' lp <- matrix(runif(nrow(basis)), nrow = nrow(basis))
#' wei_sim <- sim_tte(pi = lp, mu = -1, coefs = 1.1,
#'   time = time, type = "weibull", end_time = 100)
#' ms_sim <- sim_tte(pi = lp, mu = mu, basis = basis,
#'   coefs = coefs, time = time, type = "ms")
#' }
sim_tte <- function(pi, log_pi = TRUE, mu = -3, coefs = 0, basis = NULL,
    time = 100, end_time, type = "weibull", ...) {
    ID <- NULL
    basehaz <- NULL
    shape <- NULL

    if (!is.numeric(pi)) {
        stop("'pi' must be numeric.")
    }
    if (!is.numeric(mu) || length(mu) != 1L) {
        stop("'mu' must be a numeric scalar.")
    }
    if (!is.numeric(coefs)) {
        stop("'coefs' must be numeric.")
    }
    type <- match.arg(type, choices = c("weibull", "ms"))

    if (type == "ms") {
        if (is.null(basis)) {
            stop("'basis' must be provided for M-spline models.")
        }
        if (ncol(basis) != length(coefs)) {
            stop("Basis columns and coefficients must have the same length.")
        }
        basehaz <- basis %*% coefs
    }
    if (type == "weibull") {
        if (length(coefs) != 1L) {
            stop("Only 1 coefficient (shape) allowed for the Weibull model.")
        }
        if (coefs <= 0) {
            stop("Weibull shape parameter must be positive.")
        }
        shape <- coefs
    }
    if (!log_pi) {
        if (any(pi <= 0)) {
            stop("All prognostic index values must be positive for log transform.")
        }
        pi <- log(pi)
    }
    data_sim <- .sim_surv_df(log_hr = pi, mu = mu, basehaz = basehaz,
        type = type, shape = shape, times = time, ...)
    pi <- as.numeric(pi)
    xdata <- data.frame(ID = seq_along(pi), lp = pi)
    dat <- sim_tte_df(data_sim, id_var = "ID", xdata = xdata)
    return(dat)
}

#' Simulate Time-to-Event Data from a Custom mrgsolve Output
#'
#' Applies inverse transform sampling to a data frame produced by an
#' \pkg{mrgsolve} simulation that contains a survival probability column.
#' This function is useful when you have a custom time-to-event model
#' implemented in \pkg{mrgsolve} and want to generate event times.
#'
#' @param dat Data frame. Output from \code{\link[mrgsolve]{mrgsim}}
#'   containing at minimum columns for time, subject ID, and survival
#'   probability. Columns are resolved by name (not position).
#' @param surv_var Character string. Name of the column containing the
#'   survival probability (probability of remaining event-free). Default
#'   \code{"p11"}.
#' @param id_var Character string. Name of the subject ID column.
#'   Default \code{"ID"}.
#' @param time_var Character string. Name of the time column.
#'   Default \code{"time"}.
#' @param xdata Optional data frame of additional covariates to merge
#'   into the output. Must contain a column matching \code{id_var}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{sim_time}{Simulated event or censoring time.}
#'   \item{sim_status}{Event indicator (1 = event, 0 = censored).}
#'   \item{ID}{Subject identifier.}
#' }
#' Plus any additional columns from \code{xdata}.
#'
#' @export
#' @examples
#' # Create a mock survival probability data frame (no mrgsolve required)
#' mock_dat <- data.frame(
#'   ID = rep(1:3, each = 50),
#'   time = rep(seq(0.1, 10, length.out = 50), 3),
#'   p11 = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 3)
#' )
#' result <- sim_tte_df(mock_dat)
#' head(result)
#'
#' # Custom column names
#' mock_dat2 <- data.frame(
#'   subj = rep(101:103, each = 50),
#'   tgrid = rep(seq(0.1, 10, length.out = 50), 3),
#'   surv = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 3)
#' )
#' result2 <- sim_tte_df(mock_dat2, surv_var = "surv",
#'   id_var = "subj", time_var = "tgrid")
#' head(result2)
sim_tte_df <- function(dat,
                       surv_var = "p11",
                       id_var = "ID",
                       time_var = "time",
                       xdata = NULL) {
    dat <- as.data.frame(dat)

    # Validate no duplicated column names
    if (anyDuplicated(names(dat))) {
        stop("'dat' has duplicated column names.", call. = FALSE)
    }
    # Validate required columns exist
    for (v in c(id_var, time_var, surv_var)) {
        if (!v %in% names(dat)) {
            stop("Column '", v, "' not found in 'dat'.", call. = FALSE)
        }
    }
    # Validate types
    if (!is.numeric(dat[[time_var]])) {
        stop("Column '", time_var, "' must be numeric.", call. = FALSE)
    }
    if (!is.numeric(dat[[surv_var]])) {
        stop("Column '", surv_var, "' must be numeric.", call. = FALSE)
    }
    # Validate IDs not missing
    if (anyNA(dat[[id_var]])) {
        stop("Column '", id_var, "' must not contain missing values.",
            call. = FALSE)
    }

    # Canonicalize to internal names: ID, time, p11
    simdat <- data.frame(
        ID = dat[[id_var]],
        time = dat[[time_var]],
        p11 = dat[[surv_var]]
    )

    ids <- unique(simdat$ID)
    data_split <- split(simdat, simdat$ID)

    pred_surv <- lapply(seq_along(ids), function(i) {
        .simulate_survival_id(data_split[[as.character(ids[i])]], ids[i],
            "ID")
    })
    pred_surv <- dplyr::bind_rows(pred_surv)
    colnames(pred_surv)[1:2] <- c("sim_time", "sim_status")

    if (is.null(xdata)) {
        return(pred_surv)
    }

    # Validate xdata
    xdata <- as.data.frame(xdata)
    if (anyDuplicated(names(xdata))) {
        stop("'xdata' has duplicated column names.", call. = FALSE)
    }
    if (!id_var %in% names(xdata)) {
        stop("Column '", id_var, "' not found in 'xdata'.", call. = FALSE)
    }
    xdata$ID <- xdata[[id_var]]
    if (anyDuplicated(xdata$ID)) {
        stop("'xdata' has duplicated IDs.", call. = FALSE)
    }
    dplyr::inner_join(pred_surv, xdata, by = "ID")
}

#' Simulate Survival for Internal Data Frame Construction
#'
#' Internal workhorse function that sets up the mrgsolve input data
#' and runs the ODE simulation to compute survival probabilities
#' over time.
#'
#' @param log_hr Numeric vector. Log hazard ratios (prognostic indices).
#' @param mu Numeric scalar. Intercept parameter.
#' @param shape Numeric scalar. Shape parameter (Weibull only).
#' @param type Character string. Model type.
#' @param times Numeric vector. Time points for simulation.
#' @param basehaz Numeric matrix. Baseline hazard values (M-spline only).
#' @param end_time Numeric scalar. Censoring time.
#' @param ... Additional arguments passed to \code{\link[mrgsolve]{mrgsim}}.
#'
#' @return A data frame of mrgsolve simulation output.
#' @noRd
.sim_surv_df <- function(log_hr, mu, shape, type, times, basehaz,
    end_time, ...) {
    if (missing(end_time)) {
        end_time <- max(times)
    }
    mod_surv <- .read_model_static_cache(type)
    if (type == "weibull") {
        ev1 <- as.data.frame(expand.grid(shape, log_hr))
        colnames(ev1) <- c("shape", "lp")
        data_surv <- ev1 %>% dplyr::mutate(
            ID = seq_len(nrow(ev1)), cmt = 1,
            amt = 0, evid = 1, time = 0, mu = mu,
            basehaz_id = as.numeric(as.factor(shape)))
    } else {
        basehaz_id <- rep(seq_len(ncol(basehaz)), each = nrow(basehaz))
        ev1 <- expand.grid(c(basehaz), log_hr) %>% as.data.frame()
        ev1$basehaz_id <- rep(basehaz_id, length(log_hr))
        colnames(ev1) <- c("basehaz", "lp", "basehaz_id")
        ev1$time <- c(rep(times, ncol(basehaz) * length(log_hr)))
        ev1$ID <- c(rep(seq_len(ncol(basehaz) * length(log_hr)),
            each = nrow(basehaz)))
        data_surv <- ev1 %>% dplyr::mutate(cmt = 1, amt = 0, evid = 1,
            mu = mu)
    }
    out <- mrgsolve::mrgsim(
        mrgsolve::data_set(mod_surv, data_surv),
        carry_out = c("basehaz_id", "lp", "mu", "shape", "basehaz"),
        end = end_time, ...)
    as.data.frame(out)
}


#' Explore Prognostic Index and Quantile Survival Times
#'
#' Calculates the difference in survival probability at a given quantile
#' time across a range of prognostic index values. Useful for understanding
#' how different log hazard ratios affect survival outcomes.
#'
#' @param q Numeric scalar. Survival quantile of interest (default 0.5,
#'   i.e. median survival).
#' @param pi Numeric vector. Prognostic index values to evaluate.
#' @param mu Numeric scalar. Intercept parameter.
#' @param shape Numeric scalar or vector. Shape parameter(s) for the
#'   Weibull model.
#' @param type Character string. Model type: \code{"weibull"} or
#'   \code{"ms"}.
#' @param times Numeric vector. Time points (M-spline models).
#' @param basehaz Numeric matrix. Baseline hazard (M-spline models).
#' @param end_time Numeric scalar. Administrative censoring time.
#' @param ... Additional arguments passed to \code{\link[mrgsolve]{mrgsim}}.
#'
#' @return A data frame including columns \code{lp}, \code{p11},
#'   \code{survdiff_tq}, and model parameters.
#'
#' @export
#' @examples
#' # Small fast example
#' data_sim <- explore_pi_tq_surv(
#'   pi = seq(-1, 1, by = 0.5),
#'   mu = -1,
#'   shape = 1.1,
#'   end_time = 10,
#'   type = "weibull"
#' )
#' head(data_sim)
#' \donttest{
#' # Larger range example
#' data_sim2 <- explore_pi_tq_surv(
#'   pi = seq(-3, 3, by = 0.1),
#'   mu = -1,
#'   shape = 1.1,
#'   end_time = 200,
#'   type = "weibull"
#' )
#' plot(survdiff_tq ~ lp, data = data_sim2)
#' }
explore_pi_tq_surv <- function(q = 0.5,
                               pi,
                               mu,
                               shape,
                               type,
                               times,
                               basehaz,
                               end_time, ...) {
    lp <- p11 <- basehaz_id <- survrelq <- time <- p11_baseline <- NULL
    pi <- c(pi)
    sim_surv <- .sim_surv_df(pi, mu, shape, type, times, basehaz,
        end_time, ...)
    tq_baseline <- sim_surv %>%
        dplyr::filter(lp == 0) %>%
        dplyr::mutate(survrelq = abs(p11 - q)) %>%
        dplyr::group_by(basehaz_id) %>%
        dplyr::filter(survrelq == min(survrelq)) %>%
        dplyr::ungroup() %>%
        dplyr::select(time, basehaz_id) %>%
        dplyr::distinct()
    ratio_treat <- dplyr::inner_join(sim_surv, tq_baseline,
        by = c("time", "basehaz_id")) %>%
        dplyr::distinct()
    ratio_baseline <- ratio_treat %>%
        dplyr::filter(lp == 0) %>%
        dplyr::select(basehaz_id, p11_baseline = p11)
    dplyr::left_join(ratio_treat, ratio_baseline, by = "basehaz_id") %>%
        dplyr::mutate(survdiff_tq = p11 - p11_baseline) %>%
        dplyr::distinct()
}
