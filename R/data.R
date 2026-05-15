#' Simulated Data for M-Splines Model Demonstration
#'
#' A list containing parameters for an M-spline survival model,
#' suitable for use with \code{\link{sim_tte}}.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{mu}{Numeric scalar. Intercept parameter.}
#'   \item{basis}{Numeric matrix. M-spline basis matrix with rows
#'     corresponding to time points and columns to basis functions.}
#'   \item{coefs}{Numeric vector. Coefficients for each basis function.}
#'   \item{time}{Numeric vector. Time points corresponding to the rows
#'     of the basis matrix.}
#' }
#' @source Simulated example data.
"ms_data"
