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
#'   vector b, \code{pi = X \%*\% b}. Must be non-empty and finite (no
#'   \code{NA}/\code{NaN}/\code{Inf}/\code{-Inf}).
#' @param log_pi Logical; is the prognostic index already on the log scale?
#'   Default \code{TRUE}.
#' @param mu Numeric scalar. Intercept parameter of the model. Must be
#'   finite. Default \code{-3}.
#' @param coefs Numeric vector. For M-spline models, the coefficients of
#'   each spline basis function. For Weibull models, the shape parameter
#'   (scalar), which must be strictly positive. Must be finite in all
#'   cases.
#' @param basis Numeric matrix. Basis matrix for M-spline models, with
#'   one row per element of \code{time} and one column per coefficient
#'   in \code{coefs}. Must be finite. Ignored for Weibull models.
#'   Default \code{NULL}.
#' @param time Numeric vector. The output time grid: the exact set of
#'   times at which the survival trajectory is simulated and reported,
#'   and therefore the set of times at which an event can be detected
#'   (see "Time grid and event-time resolution" below). For M-spline
#'   models, \code{time} must additionally be sorted in non-decreasing
#'   order and have one element per row of \code{basis} (the baseline
#'   hazard is only defined at these points). Must be non-negative,
#'   finite, and contain no \code{NA}/\code{NaN}. Default
#'   \code{seq(0, 100, by = 1)}.
#' @param end_time Numeric scalar. Administrative censoring horizon.
#'   Defaults to \code{max(time)}. Any requested \code{time} points
#'   beyond \code{end_time} are dropped, and \code{end_time} itself is
#'   always added to the output grid so that censoring lands exactly on
#'   it. Must be non-negative; \code{end_time = 0} is accepted (a
#'   degenerate zero-duration follow-up: the output grid is just time 0,
#'   and since \eqn{S(0) = 1} for every supported model, every subject
#'   is censored at 0). For Weibull models \code{end_time} may otherwise
#'   be any non-negative value, including below \code{min(time)} or
#'   above \code{max(time)} (the closed-form hazard is defined for all
#'   \code{t >= 0}). For M-spline models \code{end_time} must satisfy
#'   \code{min(time) <= end_time <= max(time)}: the baseline hazard is
#'   only defined on the supplied \code{time} grid and this package does
#'   not extrapolate it in either direction.
#' @param type Character string. Model type: \code{"weibull"} or
#'   \code{"ms"} (M-splines). Default \code{"weibull"}.
#' @param ... Additional arguments passed to
#'   \code{\link[mrgsolve]{mrgsim}}. \code{tgrid}, \code{obsonly},
#'   \code{nocb}, \code{carry_out}/\code{carry.out}, and \code{data} are
#'   controlled internally to guarantee the contracts documented here
#'   and cannot be overridden; supplying them raises an error naming the
#'   argument.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{sim_time}{Simulated event or censoring time.}
#'   \item{sim_status}{Event indicator (1 = event, 0 = censored).}
#'   \item{ID}{Subject identifier.}
#'   \item{lp}{Log prognostic index (linear predictor).}
#' }
#'
#' @section Time grid and event-time resolution:
#' The output grid used for simulation is \code{sort(unique(time))},
#' restricted to values \code{<= end_time}, with \code{end_time}
#' appended if it is not already present. \pkg{mrgsolve} is instructed
#' to report the survival trajectory \emph{only} at these times (via
#' its \code{tgrid} mechanism), so the grid is deterministic and is not
#' silently replaced by \pkg{mrgsolve}'s own default output schedule.
#'
#' Event times are then resolved by \code{\link{sim_tte_df}} as the
#' first grid point at which the survival probability falls to or below
#' a uniform random draw. No interpolation between grid points is
#' performed: the resolution of the simulated event time is exactly the
#' spacing of \code{time}. A subject whose survival never crosses the
#' threshold on the grid is censored at \code{end_time}. See
#' \code{\link{sim_tte_df}} for the full contract (this is the function
#' \code{sim_tte()} delegates to internally).
#'
#' @section M-spline hazard carry convention:
#' For \code{type = "ms"}, \code{basis \%*\% coefs} gives a baseline
#' hazard value at each element of \code{time}; this is supplied to
#' \pkg{mrgsolve} as a time-varying input, not evaluated as a continuous
#' spline function at arbitrary times. The package fixes the piecewise
#' interpretation explicitly (\pkg{mrgsolve}'s \code{nocb = FALSE}, i.e.
#' last-observation-carried-forward): the hazard value reported at
#' \code{time[i]} applies from \code{time[i]} (inclusive) until
#' \code{time[i + 1]} (exclusive) -- a right-continuous step function in
#' the hazard. Concretely, for knots/hazards \code{time = c(0, 1, 2)},
#' \code{basehaz = c(1, 2, 4)}: \eqn{S(0.5) = \exp(-0.5)},
#' \eqn{S(1) = \exp(-1)}, \eqn{S(1.5) = \exp(-2)},
#' \eqn{S(2) = \exp(-3)} (verified directly against \pkg{mrgsolve}
#' output; this is \strong{not} the \pkg{mrgsolve} default, which is
#' next-observation-carried-backward and would instead give
#' \eqn{S(0.5) = \exp(-1)}).
#'
#' A consequence of this convention: the hazard value associated with
#' the \emph{last} element of \code{time} has no interval after it to
#' apply to, so it never affects the survival trajectory (it would only
#' matter if \code{end_time} extended past \code{max(time)}, which is
#' not permitted; see \code{end_time} above). \code{basis}/\code{coefs}
#' therefore define \code{length(time) - 1} active piecewise-constant
#' hazard segments, not \code{length(time)}.
#'
#' This is a grid-based piecewise-constant hazard, not continuous
#' M-spline evaluation: the package receives precomputed basis values at
#' specific times and never re-evaluates the spline basis at arbitrary
#' times.
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
    time = seq(0, 100, by = 1), end_time, type = "weibull", ...) {
    ID <- NULL
    basehaz <- NULL
    shape <- NULL

    if (!is.numeric(pi)) {
        stop("'pi' must be numeric.")
    }
    if (length(pi) < 1L) {
        stop("'pi' must contain at least one individual (empty cohorts ",
            "are not allowed).")
    }
    if (any(!is.finite(pi))) {
        stop("'pi' must not contain NA, NaN, Inf, or -Inf values.")
    }
    if (!is.numeric(mu) || length(mu) != 1L) {
        stop("'mu' must be a numeric scalar.")
    }
    if (!is.finite(mu)) {
        stop("'mu' must be finite (not NA, NaN, Inf, or -Inf).")
    }
    if (!is.numeric(coefs)) {
        stop("'coefs' must be numeric.")
    }
    if (any(!is.finite(coefs))) {
        stop("'coefs' must not contain NA, NaN, Inf, or -Inf values.")
    }
    .validate_time_grid(time, arg = "time")
    if (!missing(end_time)) {
        .validate_end_time(end_time)
    } else {
        end_time <- max(time)
    }
    type <- match.arg(type, choices = c("weibull", "ms"))

    if (type == "ms") {
        if (is.null(basis)) {
            stop("'basis' must be provided for M-spline models.")
        }
        if (!is.matrix(basis) || !is.numeric(basis)) {
            stop("'basis' must be a numeric matrix.")
        }
        if (ncol(basis) != length(coefs)) {
            stop("Basis columns and coefficients must have the same length.")
        }
        if (nrow(basis) != length(time)) {
            stop("'basis' must have one row per element of 'time' ",
                "(nrow(basis) = ", nrow(basis), ", length(time) = ",
                length(time), ").")
        }
        if (any(!is.finite(basis))) {
            stop("'basis' must not contain NA, NaN, Inf, or -Inf values.")
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
        type = type, shape = shape, times = time, end_time = end_time,
        ...)
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
#' @section Trajectory contract:
#' \code{dat} must contain, for each subject, one ordered survival
#' trajectory satisfying:
#' \itemize{
#'   \item at least one row;
#'   \item \code{time} values that are non-negative, finite, and
#'     strictly increasing within the subject, with no duplicates.
#'     Rows are \strong{not} sorted automatically: unsorted or
#'     duplicated times raise an error rather than being silently
#'     re-ordered or resolved by row order;
#'   \item survival probabilities in \eqn{[0, 1]} (a small numerical
#'     tolerance of \eqn{10^{-8}} is allowed for floating-point/ODE
#'     solver noise: values outside \eqn{[0, 1]} but within this
#'     tolerance are accepted and then deterministically \strong{clamped}
#'     to exactly \eqn{0} or \eqn{1} before monotonicity checking and
#'     event-time selection; values further outside \eqn{[0, 1]} are a
#'     hard error, never silently clamped or dropped). This clamping
#'     cannot change whether a subject is classified as an event or
#'     censored, because \code{stats::runif()} draws from \eqn{[0, 1)}:
#'     a raw value already below 0 and its clamped value of exactly 0
#'     both satisfy \eqn{value \le U} for every achievable \eqn{U}, and
#'     symmetrically a raw value already above 1 and its clamped value
#'     of exactly 1 both fail it;
#'   \item survival that is non-increasing within the subject (checked
#'     on the clamped values). A trajectory that increases by more than
#'     the same \eqn{10^{-8}} tolerance at any step is rejected rather
#'     than silently repaired, because that pattern cannot arise from a
#'     valid survival function and usually indicates a data or model
#'     error; an increase within the tolerance is treated as numerical
#'     noise and accepted as-is (not further adjusted).
#' }
#' A trajectory is \strong{not} required to start at survival = 1 or at
#' time = 0: \code{sim_tte_df()} never invents an implicit \eqn{S(0) =
#' 1} state. Whatever survival value is reported at the first observed
#' time is used exactly as given (see "Censoring and event-time
#' resolution" below for what this implies when that first value is
#' already at or below the sampled threshold). Subjects may have
#' different time ranges; see below for what this means for censoring.
#'
#' @section Censoring and event-time resolution:
#' For each subject, a uniform random variate \eqn{U \sim
#' \mathrm{Uniform}(0, 1)} is drawn, and the event time is the
#' \strong{first reported time} at which survival is at or below
#' \eqn{U}. This is a grid-based lookup, not interpolation: the
#' returned time is always one of the times present in \code{dat} for
#' that subject, never a value between two reported times. If the
#' first reported observation already satisfies \eqn{S(t_1) \le U},
#' that first time is returned as the event time; no earlier crossing
#' is inferred.
#'
#' If a subject's trajectory never crosses \eqn{U}, that subject is
#' censored (\code{sim_status = 0}) at the \strong{last reported time
#' in their own trajectory}, i.e. \code{max(time)} for that subject.
#' Because this censoring time is subject-specific, \code{sim_tte_df()}
#' does \strong{not} inherently represent a common administrative
#' censoring cutoff unless every subject's trajectory happens to share
#' the same final time — each subject's final observation time is
#' effectively that subject's own censoring horizon. There is
#' currently no separate \code{end_time} argument to
#' \code{sim_tte_df()} itself; \code{\link{sim_tte}} provides one for
#' its built-in Weibull/M-spline models by ensuring every subject's
#' generated trajectory is extended to a common \code{end_time} (see
#' \code{?sim_tte}).
#'
#' @section Reproducibility:
#' \code{sim_tte_df()} draws exactly one \code{stats::runif(1)} per
#' subject, in the order subjects are first encountered in \code{dat}
#' (row order), not in ID-sorted order. With a fixed \code{set.seed()},
#' results are therefore reproducible \strong{for a fixed row order} of
#' \code{dat}, but reordering the rows of \code{dat} so that subjects
#' are first encountered in a different sequence changes which random
#' draw each subject receives, even with the same seed and the same set
#' of per-subject trajectories. This is standard, unmodified R RNG
#' behavior (a single shared random stream consumed sequentially), not a
#' bug; \code{sim_tte_df()} does not currently use per-subject or
#' per-ID RNG streams.
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

    if (nrow(dat) == 0L) {
        stop("'dat' has no rows.", call. = FALSE)
    }
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
    if (is.numeric(dat[[id_var]]) && any(!is.finite(dat[[id_var]]))) {
        stop("Column '", id_var, "' must not contain Inf, -Inf, or NaN ",
            "values.", call. = FALSE)
    }
    # Validate time and survival values: no NA/NaN/Inf/-Inf. is.finite()
    # is FALSE for all four, so a single check covers them.
    if (any(!is.finite(dat[[time_var]]))) {
        stop("Column '", time_var, "' must not contain NA, NaN, Inf, or ",
            "-Inf values.", call. = FALSE)
    }
    if (any(!is.finite(dat[[surv_var]]))) {
        stop("Column '", surv_var, "' must not contain NA, NaN, Inf, or ",
            "-Inf values.", call. = FALSE)
    }
    # Time must be non-negative.
    if (any(dat[[time_var]] < 0)) {
        stop("Column '", time_var, "' must not contain negative values.",
            call. = FALSE)
    }
    # Survival must be a probability, in [0, 1] (allowing a small
    # tolerance for floating-point / ODE solver noise; see .SURV_TOL).
    # Values outside the tolerance band are a hard error. Values inside
    # it are accepted, and then deterministically clamped to exactly
    # [0, 1] below via .normalize_survival() before anything else (event
    # selection, monotonicity checking) sees them; see .SURV_TOL for why
    # this clamping is safe and cannot change which trajectories are
    # classified as events vs. censoring.
    if (any(dat[[surv_var]] < -.SURV_TOL) ||
        any(dat[[surv_var]] > 1 + .SURV_TOL)) {
        stop("Column '", surv_var, "' must contain survival ",
            "probabilities in [0, 1] (a tolerance of ", .SURV_TOL,
            " is allowed for floating-point/solver noise).",
            call. = FALSE)
    }

    # Canonicalize to internal names: ID, time, p11. Survival is
    # normalized (clamped to [0, 1]) here, once, so every downstream
    # consumer -- monotonicity validation and event-time selection --
    # operates on the same normalized values, never the raw column.
    simdat <- data.frame(
        ID = dat[[id_var]],
        time = dat[[time_var]],
        p11 = .normalize_survival(dat[[surv_var]])
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
    # Reject `...` arguments that would silently override the
    # output-grid/trajectory/hazard-carry contract this function
    # guarantees (see ?.RESERVED_MRGSIM_ARGS and Phase pre-B audit).
    # Checked before any other work so the error surfaces immediately.
    .check_reserved_dots(list(...))

    if (missing(end_time)) {
        end_time <- max(times)
    }
    if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) {
        stop("'mu' must be a finite numeric scalar.", call. = FALSE)
    }
    if (!is.numeric(log_hr) || length(log_hr) < 1L ||
        any(!is.finite(log_hr))) {
        stop("'pi' (the linear predictor) must be finite numeric and ",
            "contain at least one value.", call. = FALSE)
    }
    if (type == "weibull") {
        if (!is.numeric(shape) || length(shape) < 1L ||
            any(!is.finite(shape))) {
            stop("'shape' must be finite numeric.", call. = FALSE)
        }
        if (any(shape <= 0)) {
            stop("'shape' (Weibull shape parameter) must be strictly ",
                "positive.", call. = FALSE)
        }
    }
    if (type == "ms" && is.unsorted(times)) {
        stop("'time' must be sorted in non-decreasing order for ",
            "M-spline models: each row of 'basis' corresponds ",
            "positionally to the matching element of 'time'.",
            call. = FALSE)
    }
    if (type == "ms") {
        .validate_basehaz(basehaz, times)
    }
    # The reported output grid is resolved once, deterministically, and
    # used verbatim as mrgsolve's `tgrid`. This is the mechanism that
    # makes the public `time` argument of sim_tte() actually control
    # the grid on which event times are resolved (see the "Time grid
    # and event-time resolution" section of ?sim_tte).
    grid <- .resolve_output_grid(times, end_time, type)

    mod_surv <- .read_model_static_cache(type)
    if (type == "weibull") {
        ev1 <- as.data.frame(expand.grid(shape, log_hr))
        colnames(ev1) <- c("shape", "lp")
        data_surv <- ev1 %>% dplyr::mutate(
            ID = seq_len(nrow(ev1)), cmt = 0,
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
        tgrid = grid, obsonly = TRUE,
        # nocb = FALSE (last-observation-carried-forward): the hazard
        # value reported at time t_i applies from t_i until the next
        # knot. This only affects the M-spline model (the Weibull model
        # has no time-varying mrgsolve input); see the "M-spline hazard
        # carry convention" section of ?sim_tte for the full contract
        # and rationale.
        nocb = FALSE, ...)
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
#' @param pi Numeric vector. Prognostic index values to evaluate. Must
#'   be finite.
#' @param mu Numeric scalar. Intercept parameter. Must be finite.
#' @param shape Numeric scalar or vector. Shape parameter(s) for the
#'   Weibull model. Must be finite and strictly positive.
#' @param type Character string. Model type: \code{"weibull"} or
#'   \code{"ms"}.
#' @param times Numeric vector. Output time grid; see the "Time grid and
#'   event-time resolution" and "M-spline hazard carry convention"
#'   sections of \code{\link{sim_tte}} for the exact contract. For
#'   M-spline models this must also be the time points corresponding to
#'   the rows of \code{basehaz}. If missing, defaults to
#'   \code{seq(0, end_time, by = 1)}.
#' @param basehaz Numeric matrix. Baseline hazard (M-spline models), one
#'   row per element of \code{times}. Must be finite and non-negative;
#'   see \code{\link{sim_tte}}'s "M-spline hazard carry convention" for
#'   how it is interpreted between supplied times.
#' @param end_time Numeric scalar. Administrative censoring time. See
#'   the \code{end_time} contract documented in \code{\link{sim_tte}}
#'   (including the M-spline \code{min(times) <= end_time <= max(times)}
#'   requirement).
#' @param ... Additional arguments passed to
#'   \code{\link[mrgsolve]{mrgsim}}. \code{tgrid}, \code{obsonly},
#'   \code{nocb}, \code{carry_out}/\code{carry.out}, and \code{data}
#'   cannot be overridden; see \code{\link{sim_tte}}.
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
    if (missing(times)) {
        # Preserves the pre-Phase-A default output resolution (a unit
        # grid) now that `times` genuinely controls the output grid
        # instead of being silently ignored (see ?sim_tte).
        times <- seq(0, end_time, by = 1)
    }
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
