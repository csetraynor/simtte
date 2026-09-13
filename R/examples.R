#' Path to simtte's shipped example mrgsolve models
#'
#' Kept separate from \code{\link{.cfile_dir}} (\code{inst/models/}),
#' which is the private working directory of the internal
#' \code{weibull}/\code{weibull_tv}/\code{ms} engine models loaded by
#' \code{\link{.read_model_static_cache}}. \code{inst/models/examples/}
#' holds only user-facing, editable example models, discovered purely by
#' directory listing (see \code{\link{simtte_example_models}}); nothing
#' under this directory is ever loaded internally by \code{\link{sim_tte}}
#' or \code{\link{sim_tte_df}}.
#' @return Character path, or \code{""} if the package is not installed
#'   (mirrors \code{\link[base]{system.file}}'s own convention).
#' @noRd
.example_models_dir <- function() {
    system.file("models", "examples", package = "simtte")
}

#' List simtte's bundled example mrgsolve models
#'
#' Returns the names of the mechanistic mrgsolve models shipped as
#' user-facing examples/templates in \code{inst/models/examples/}. These
#' are plain, editable \code{.cpp} model source files, not part of the
#' package's internal Weibull/M-spline simulation engine (see
#' \code{\link{sim_tte}}) and not used by it. They demonstrate the
#' \code{dxdt_p11 = -p11 * HAZ} survival co-integration pattern that a
#' user's own custom mrgsolve model is expected to follow when working
#' with \code{\link{sim_tte_df}} directly; see "Shipped example models"
#' below.
#'
#' Discovery is a directory listing of \code{*.cpp} files in
#' \code{inst/models/examples/}, not a maintained registry or catalog:
#' adding a new example model in a future package version requires only
#' adding a new \code{.cpp} file to that directory, with no change to
#' this function.
#'
#' @section Shipped example models:
#' \describe{
#'   \item{\code{pkpd_idr_hazard}}{A one-compartment oral-absorption PK
#'     model coupled to an indirect-response (turnover) mediator model,
#'     in which drug exposure inhibits mediator production and the event
#'     hazard is proportional to the current mediator level. Includes
#'     between-subject variability (\code{$OMEGA}) on clearance, volume,
#'     and potency. This is the model described in the package's
#'     companion manuscript ("ODE-coupled PK/PD hazard"), shipped here as
#'     a standalone, loadable file.}
#'   \item{\code{pkpd_linear_hazard}}{A minimal one-compartment IV-bolus
#'     PK model with a hazard that decreases linearly in drug
#'     concentration from a constant baseline. No between-subject
#'     variability or turnover dynamics. Both the PK and hazard are
#'     linear, so this model has an independent closed-form survival
#'     function, making it a short, readable starting point for a custom
#'     model as well as this package's validation anchor for mechanistic
#'     hazard models (see \code{inst/validation/}).}
#' }
#'
#' @return Character vector of example model names, suitable for passing
#'   to \code{\link{simtte_example_model}}.
#' @seealso \code{\link{simtte_example_model}}, to load one;
#'   \code{\link{sim_tte_df}}, the function these models are for.
#' @export
#' @examples
#' simtte_example_models()
simtte_example_models <- function() {
    dir <- .example_models_dir()
    tools::file_path_sans_ext(list.files(dir, pattern = "\\.cpp$"))
}

#' Load a bundled simtte example mrgsolve model
#'
#' Compiles (or loads from cache) one of simtte's shipped example
#' mechanistic mrgsolve models; see \code{\link{simtte_example_models}}
#' for the list of available names and what each model represents. This
#' is a thin wrapper around \code{\link[mrgsolve]{mread_cache}} pointed
#' at the package's example-model directory; it adds nothing beyond
#' discoverability and an informative error for an unrecognized name.
#'
#' The returned model is a standard \pkg{mrgsolve} model object: use
#' \code{\link[mrgsolve]{mrgsim}} to simulate it, exactly as with any
#' other \pkg{mrgsolve} model. The underlying \code{.cpp} source file can
#' be found under \code{system.file("models", "examples", package =
#' "simtte")} and copied out and modified freely with any \pkg{mrgsolve}
#' tooling; loading it through this function does not require continuing
#' to use this function afterward.
#'
#' @section Intended workflow:
#' These example models are for use with \code{\link{sim_tte_df}}, the
#' package's entry point for custom/mechanistic models with an
#' ODE-integrated survival compartment (see its "Trajectory contract"
#' section); they are not used by, and cannot be passed to,
#' \code{\link{sim_tte}}, which only ever simulates its own built-in
#' Weibull and M-spline models.
#' \preformatted{
#' mod  <- simtte_example_model("pkpd_linear_hazard")
#' data <- mrgsolve::ev(amt = 100, cmt = 1, time = 0)
#' out  <- as.data.frame(mrgsolve::mrgsim(mod, data = data, end = -1,
#'           add = seq(0, 60, by = 0.5), obsonly = TRUE))
#' sim_tte_df(out[, c("ID", "time", "p11")])
#' }
#'
#' @param name Character scalar. One of \code{\link{simtte_example_models}()}.
#' @param ... Additional arguments passed to
#'   \code{\link[mrgsolve]{mread_cache}} (for example \code{quiet = TRUE}).
#'
#' @return A compiled \pkg{mrgsolve} model object.
#' @seealso \code{\link{simtte_example_models}}, to list available names;
#'   \code{\link{sim_tte_df}}, the function these models are for.
#' @export
#' @examples
#' \donttest{
#' mod <- simtte_example_model("pkpd_linear_hazard")
#' }
simtte_example_model <- function(name, ...) {
    if (!is.character(name) || length(name) != 1L || is.na(name)) {
        stop("'name' must be a single character string.", call. = FALSE)
    }
    available <- simtte_example_models()
    if (!name %in% available) {
        stop("Unknown example model '", name, "'. Available example ",
            "models: ", paste(available, collapse = ", "), ".",
            call. = FALSE)
    }
    mrgsolve::mread_cache(model = name, project = .example_models_dir(),
        ...)
}
