#' Where simtte builds and caches compiled models
#'
#' Resolves and returns the directory \code{\link{sim_tte_ode}} and
#' \code{\link{simtte_prepare_model}} pass as \pkg{mrgsolve}'s own
#' \code{soloc} argument (where a model's compiled shared object is
#' built and stored) -- simtte itself never sets \code{options()} on a
#' user's behalf, for \code{mrgsolve}'s own options or this package's
#' \code{simtte.cache_dir}; this function only ever \emph{reads}
#' \code{getOption("simtte.cache_dir")} and appends a version-qualified
#' subdirectory before returning it.
#'
#' \strong{Default} (\code{simtte.cache_dir} unset): a directory under
#' \code{\link[base]{tempdir}()}, i.e. session-scoped and automatically
#' cleaned up -- CRAN-safe (nothing is ever written outside a session
#' temp directory unless you opt in below) and the right choice for any
#' single script or session, including a local \pkg{parallel} cluster
#' (every worker of a \emph{local} cluster shares the parent's
#' filesystem, so a model built once under the parent's own cache
#' directory is directly reachable by every worker's
#' \code{\link[mrgsolve]{loadso}} call -- see \code{?sim_tte_ode}
#' "Parallel simulation").
#'
#' \strong{Persistent cache}: set \code{options(simtte.cache_dir =
#' tools::R_user_dir("simtte", "cache"))} (the standard, CRAN-compliant
#' location for a package's own cache data -- \code{\link{simtte_model_cache}}
#' never chooses this for you) to keep a compiled model's shared object
#' across R restarts, avoiding a repeat compile in a later session. This
#' matters most for a slow-to-compile \code{\link{tte_model}}-converted
#' model reused across many scripts, or a genuinely distributed cluster
#' whose worker nodes do \emph{not} share the parent's filesystem (where
#' every node needs to be pointed, itself, at the same shared, network-
#' visible directory).
#'
#' \strong{Version-key caveat}: whichever directory \code{simtte.cache_dir}
#' resolves to, a subdirectory keyed by the current R version, this
#' \pkg{simtte} version, and the installed \pkg{mrgsolve} version is
#' always appended before it is returned or created. This is deliberate,
#' not incidental: \pkg{mrgsolve}'s own build-cache key is content/mtime
#' based, not version-aware (verified: its own \code{?mread}
#' documentation states plainly that \code{soloc} is "a directory
#' location where the model shared object is built and stored", with no
#' mention of invalidating it across an R or package upgrade) -- a
#' \emph{persistent} cache surviving an R version change or a
#' \pkg{simtte}/\pkg{mrgsolve} upgrade risks reusing a compiled
#' \code{.so} built against a different, ABI-incompatible R/package
#' version. The version key makes that class of staleness unreachable
#' by construction rather than relying on \pkg{mrgsolve}'s own
#' file-level staleness detection to catch it (which a prior session
#' already found can be defeated by two builds landing in the same
#' filesystem-timestamp second -- see \code{\link{tte_model}}'s own
#' \code{name} documentation). The \strong{default}, \code{tempdir()}-based
#' location needs no such protection on its own (a session cannot
#' upgrade R or a package mid-run), but the key is applied unconditionally
#' anyway, so there is exactly one code path to reason about, not two.
#'
#' @return Character scalar: an existing, writable directory path.
#' @seealso \code{\link{simtte_prepare_model}}, which uses this as
#'   \pkg{mrgsolve}'s \code{soloc}. \code{\link{simtte_model_cache_clear}},
#'   to remove what this builds up. \code{?sim_tte_ode} "Parallel
#'   simulation" for the intended cluster-init pattern.
#' @export
#' @examples
#' simtte_model_cache()
simtte_model_cache <- function() {
    base_dir <- getOption("simtte.cache_dir",
        file.path(tempdir(), "simtte-model-cache"))
    dir <- file.path(base_dir, .simtte_cache_version_tag())
    if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE)) {
        stop("simtte_model_cache() could not create its cache directory '",
            dir, "' -- check that '", base_dir, "' (from options(",
            "simtte.cache_dir = ...), or the default location under ",
            "tempdir()) is writable.", call. = FALSE)
    }
    dir
}

#' Current version-qualified cache subdirectory name
#'
#' Shared by \code{\link{simtte_model_cache}} (appends it to build the
#' directory it returns) and \code{\link{simtte_model_cache_clear}}
#' (matches it, and its own regex form, to decide what "current"/"every
#' version" means) -- kept as one function so the two never drift apart.
#'
#' @return Character scalar, e.g. \code{"R4.3.1-simtte1.1.0-mrgsolve2.0.1"}.
#' @noRd
.simtte_cache_version_tag <- function() {
    # ("Version-key caveat", simtte_model_cache()'s own docs above):
    # current R version, installed simtte version, installed mrgsolve
    # version.
    paste0("R", getRversion(), "-simtte", utils::packageVersion("simtte"),
        "-mrgsolve", utils::packageVersion("mrgsolve"))
}

#' Validate and normalize \code{sim_tte_ode()}/\code{simtte_prepare_model()}'s
#' \code{model} argument, \emph{before} unwrapping a prepared object
#'
#' Shared by both public entry points so the "character library name or
#' \code{simtte_model}" contract, and its error message, is defined
#' exactly once.
#'
#' @param model As supplied by the caller (already unwrapped if it was
#'   a \code{simtte_prepared_model} -- see \code{\link{.load_ode_model}}).
#' @return List: \code{model} (the input, or the \code{match.arg()}-resolved
#'   library name), \code{is_converted} (logical).
#' @noRd
.resolve_ode_model_spec <- function(model) {
    is_converted <- inherits(model, "simtte_model")
    if (!is_converted) {
        if (!is.character(model)) {
            stop("'model' must be a character library-model name (see ",
                "?sim_tte_ode) or a simtte_model produced by ",
                "tte_model(); a bare compiled mrgsolve model is not ",
                "accepted directly -- convert it first with tte_model() ",
                "(see ?tte_model). An object returned by ",
                "simtte_prepare_model() is also accepted.", call. = FALSE)
        }
        model <- match.arg(model, choices = c(names(.ODE_LIBRARY_FILES), "mspline"))
    }
    list(model = model, is_converted = is_converted)
}

#' Build (if needed) and \code{loadso()} an ODE model -- the one loading
#' path shared by \code{\link{sim_tte_ode}} and
#' \code{\link{simtte_prepare_model}}
#'
#' \code{\link[mrgsolve]{loadso}} is called unconditionally, even for a
#' model just compiled in this same process (verified cheap and safe to
#' repeat: ~0.2ms for 20 repeated calls, reports/experiments/30_worker_loading_test.R) --
#' this is the fix for simttepower's own open risk that a built-in
#' model reached via \code{sim_tte_ode()} on a worker node has no
#' \code{loadso()} safety net the way a \code{\link{tte_model}}-converted
#' one does: now every path gets one, unconditionally, not just the
#' \code{tte_model()} one.
#'
#' @param model Already validated/normalized by
#'   \code{\link{.resolve_ode_model_spec}}: a library name, \code{"mspline"},
#'   or a \code{simtte_model}.
#' @param is_converted As returned by \code{\link{.resolve_ode_model_spec}}.
#' @param knots Only used when \code{model = "mspline"}, to select the
#'   shipped knot-count variant (\code{\link{.mspline_file_for}}).
#' @return Compiled, \code{loadso()}-ed \code{mrgmod}.
#' @noRd
.load_ode_model <- function(model, is_converted, knots) {
    mod <- if (is_converted) {
        model$mod
    } else if (identical(model, "mspline")) {
        .read_ode_library_model_file(.mspline_file_for(length(knots)))
    } else {
        .read_ode_library_model(model)
    }
    mrgsolve::loadso(mod)
    mod
}

#' Prepare a model for simulation on this process, once
#'
#' Builds (or, if this process has already built the exact same model
#' before -- \code{\link[mrgsolve]{mread_cache}}'s/\code{\link[mrgsolve]{mcode_cache}}'s
#' own in-session cache -- reuses) and \code{\link[mrgsolve]{loadso}}s a
#' model, returning it wrapped in a class \code{\link{sim_tte_ode}}
#' recognizes and never rebuilds. \code{\link{sim_tte_ode}} calls this
#' function internally for every \code{model} it is given (a character
#' library name, a \code{\link{tte_model}}-converted \code{simtte_model},
#' or an object this function already returned) -- there is exactly one
#' model-loading path, whether or not a caller prepares anything ahead
#' of time.
#'
#' \strong{The point of calling this yourself}: a \pkg{parallel} PSOCK
#' worker is a separate R process with no memory of what the parent
#' process has already built or loaded. Calling \code{sim_tte_ode()}
#' directly, unprepared, on each worker is not wrong -- every path
#' through it now calls \code{\link[mrgsolve]{loadso}} -- but it does
#' mean each worker independently pays its own full compile cost the
#' first time (harmless for a fast-compiling built-in model; measured
#' up to ~1.5s for a \code{tte_model()}-converted one, "Parallel
#' simulation" section of \code{?sim_tte_ode}). Calling
#' \code{simtte_prepare_model()} once per worker at cluster
#' initialization avoids that redundant work \emph{and} is the
#' pattern verified race-free (see this file's own top-of-file note,
#' and \code{reports/experiments/30_worker_race_test.R}/
#' \code{30_worker_race_test2.R}): a worker receiving an
#' already-compiled model object (via ordinary \pkg{parallel}
#' serialization -- \code{parallel::clusterCall()}/\code{clusterExport()}
#' both work) and only ever calling \code{loadso()} on it never races,
#' unlike a worker independently re-triggering \code{mcode_cache()}/
#' \code{mread_cache()} against a directory another process might be
#' writing to at the same time.
#'
#' Idempotent: calling this on its own return value (exactly what
#' happens when the same \code{model =} argument is broadcast to every
#' worker via \code{clusterCall()}, each worker then calling this
#' function on it again) only re-\code{loadso()}s (cheap, safe to
#' repeat) and returns the object unchanged -- never a second build.
#'
#' @param model A character library-model name (as
#'   \code{\link{sim_tte_ode}}'s own \code{model} argument), a
#'   \code{simtte_model} from \code{\link{tte_model}}, or an object this
#'   function already returned.
#' @param knots Only used, and required, when \code{model = "mspline"}:
#'   its \emph{length} selects which of the three shipped knot-count
#'   variants (\code{?sim_tte_ode} "M-spline knot-count variants") gets
#'   compiled -- the knot \emph{values} themselves, like
#'   \code{coefs}/\code{boundary_knots}, only ever affect \code{$PARAM}
#'   values, applied later by \code{\link{sim_tte_ode}} itself, never
#'   which compiled variant is loaded (checked directly: this function's
#'   only use of \code{knots} is \code{length(knots)}). \code{coefs}/
#'   \code{boundary_knots} are therefore not parameters of this function
#'   at all -- supplying either raises R's own "unused argument" error,
#'   naming it, rather than silently accepting a value that would not
#'   affect what gets built.
#' @return An object of class \code{"simtte_prepared_model"}. Pass it
#'   directly as \code{\link{sim_tte_ode}}'s own \code{model} argument.
#' @seealso \code{\link{simtte_model_cache}}, where the compiled model
#'   is stored (and \code{\link{simtte_model_cache_clear}}, to remove
#'   it). \code{?sim_tte_ode} "Parallel simulation" for the full
#'   cluster-init pattern.
#' @export
#' @examples
#' \donttest{
#' prepared <- simtte_prepare_model("exponential")
#' sim <- sim_tte_ode(model = prepared, param = list(H0 = 0.1), n = 10,
#'   end = 20, delta = 2, seed = 1)
#' head(sim$events)
#' }
simtte_prepare_model <- function(model, knots = NULL) {
    if (inherits(model, "simtte_prepared_model")) {
        mrgsolve::loadso(model$mod)
        return(model)
    }
    spec <- .resolve_ode_model_spec(model)
    if (identical(spec$model, "mspline") && is.null(knots)) {
        stop("model = \"mspline\" requires 'knots' (its length selects ",
            "which compiled variant to load; 'coefs'/'boundary_knots' ",
            "are not needed here -- they only set $PARAM values, applied ",
            "later by sim_tte_ode() itself); see ?sim_tte_ode \"M-spline ",
            "knot-count variants\".", call. = FALSE)
    }
    mod <- .load_ode_model(spec$model, spec$is_converted, knots)
    structure(list(spec = spec$model, mod = mod,
        mspline_n_knots = if (identical(spec$model, "mspline")) length(knots)
            else NULL), class = "simtte_prepared_model")
}

#' Clear simtte's compiled-model cache
#'
#' Removes compiled model artifacts that \code{\link{simtte_prepare_model}}/
#' \code{\link{sim_tte_ode}} built under \code{\link{simtte_model_cache}}.
#' Only ever removes directories \emph{this package itself} created
#' (version-qualified subdirectories named exactly like
#' \code{\link{simtte_model_cache}}'s own return value) -- never the
#' configured base directory itself, and never anything under
#' \pkg{mrgsolve}'s own default \code{tempdir()} locations that this
#' package did not create. Safe to call even if \code{simtte.cache_dir}
#' happens to be pointed at a directory with other, unrelated content in
#' it (e.g. a user's home directory): only the recognizably-named
#' version subdirectories are ever candidates for removal, and each is
#' double-checked to actually resolve under the configured base before
#' being deleted.
#'
#' @param all Logical, default \code{FALSE}: remove only the
#'   \emph{current} version-qualified subdirectory -- the one
#'   \code{\link{simtte_model_cache}} would return right now. If
#'   \code{TRUE}, remove every version-qualified subdirectory found
#'   under the configured base -- useful after an R, simtte, or
#'   mrgsolve upgrade, when a persistent \code{simtte.cache_dir} has
#'   accumulated subdirectories for versions no longer in use.
#' @return Character vector of the directory paths actually removed,
#'   invisibly (\code{character(0)} if there was nothing to remove).
#' @seealso \code{\link{simtte_model_cache}}, which this clears.
#' @export
#' @examples
#' simtte_model_cache()
#' simtte_model_cache_clear()
simtte_model_cache_clear <- function(all = FALSE) {
    base_dir <- getOption("simtte.cache_dir",
        file.path(tempdir(), "simtte-model-cache"))
    targets <- if (isTRUE(all)) {
        # Only ever a candidate for deletion if its name is shaped like
        # .simtte_cache_version_tag()'s own output -- conservative on
        # purpose, in case the configured base has unrelated content.
        children <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
        version_tag_re <- "^R[0-9][0-9.-]*-simtte[0-9][0-9.-]*-mrgsolve[0-9][0-9.-]*$"
        children[grepl(version_tag_re, basename(children))]
    } else {
        file.path(base_dir, .simtte_cache_version_tag())
    }
    targets <- targets[dir.exists(targets)]
    base_norm <- normalizePath(base_dir, mustWork = FALSE)
    for (target in targets) {
        target_norm <- normalizePath(target, mustWork = FALSE)
        if (!startsWith(target_norm, paste0(base_norm, .Platform$file.sep))) {
            stop("simtte_model_cache_clear() refused to remove '", target,
                "': it does not resolve to a subdirectory of the ",
                "configured cache base ('", base_dir, "') -- this should ",
                "not happen; please report it.", call. = FALSE)
        }
        unlink(target, recursive = TRUE)
    }
    invisible(targets)
}

#' Print a prepared model
#' @param x A \code{"simtte_prepared_model"} object.
#' @param ... Ignored.
#' @export
#' @method print simtte_prepared_model
print.simtte_prepared_model <- function(x, ...) {
    label <- if (inherits(x$spec, "simtte_model")) x$spec$name else x$spec
    cat("<simtte_prepared_model>", label, "-- ready to simulate in this",
        "process (see ?simtte_prepare_model)\n")
    invisible(x)
}
