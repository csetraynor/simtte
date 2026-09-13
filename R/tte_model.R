#' Convert a PK/PD mrgsolve model into a `sim_tte_ode()`-ready model
#'
#' Writes the survival scaffold (\code{p11}, \code{U}, \code{END}, the
#' \code{TEVT}/\code{event_found}/\code{T_PRE}/\code{P_PRE}/\code{P_POST}
#' in-solver event-detection latch -- \code{\link{sim_tte_ode}}'s own
#' mechanism, \code{reports/02_technical_design.md} section 2) onto a
#' user's own \pkg{mrgsolve} PK/PD model, exactly the mechanical edit
#' list every shipped \code{*_hazard.cpp} library model already
#' documents by hand at the top of its \code{[PROB]} block
#' (\code{reports/10_phase4_report.md}) -- generalized here to any
#' model with an \code{$ODE}/\code{$DES} block, not just the six
#' shipped backbones. See \code{reports/13_converter_design.md} for the
#' full design (block-locating/splicing mechanics, name-collision
#' policy, why a hand-written \code{RESP0}-style helper variable is
#' never needed) and \code{reports/14_phase5a_report.md} for what was
#' built and the self-consistency check against the shipped library.
#'
#' Nothing in the input model is edited in place except its
#' \code{$GLOBAL}/\code{$MAIN} blocks (both singleton in \pkg{mrgsolve};
#' every other addition -- a new \code{$PARAM}, \code{$INIT}, an
#' \code{$ODE} block carrying only the hazard link, a \code{$CAPTURE}
#' block -- is appended as a brand-new block at the end of the file,
#' never spliced into the model's own \code{$ODE}/\code{$DES} content.
#'
#' @param model One of: a file path to a model file; a character vector
#'   of model code (either one string with embedded newlines, or
#'   already one line per element -- both accepted); or a compiled
#'   \code{mrgmod} (its own source is read back from \code{mod@code}).
#'   Must declare an \code{$ODE} or \code{$DES} block -- a
#'   \code{$PKMODEL}/\code{$PRED}-only model is rejected, since the
#'   survival compartment has to be co-integrated in the same ODE
#'   system; use the ODE-form equivalent from
#'   \code{\link[mrgsolve]{modlib}} instead (e.g. \code{"pk2cmt"}, not
#'   \code{"pk2"}). A \code{simtte_model} (this function's own return
#'   class) is rejected -- convert the original model again with a
#'   different \code{hazard}/\code{params} instead of converting an
#'   already-converted one.
#' @param hazard Character scalar: one C++ expression for the
#'   instantaneous hazard \code{HAZ}, written in terms of any name
#'   already in scope inside the model's own \code{$ODE}/\code{$DES}
#'   (states, \code{$GLOBAL} \code{#define}s, \code{$MAIN} locals,
#'   \code{$PARAM} values) plus whatever names \code{params} adds --
#'   e.g. \code{"H0 * exp(lp + beta_cp * CP)"}. Include \code{lp} in
#'   the expression if \code{\link{sim_tte_ode}}'s \code{covariates}/
#'   \code{beta} mechanism should affect this model; \code{lp} is
#'   always declared (default 0) whether or not the expression uses it.
#' @param params Named list, the hazard's own parameters and their
#'   population-default values (e.g. \code{list(H0 = 0.01, beta_cp =
#'   0)}). Must not name \code{lp}, \code{U}, or \code{END} -- those
#'   three are always added automatically (defaults \code{0}, \code{0},
#'   \code{1e9}, matching every shipped \code{*_hazard.cpp} model) and
#'   are rejected if present here.
#' @param name Character scalar or \code{NULL} (default): passed to
#'   \code{\link[mrgsolve]{mcode_cache}} as the model name (caching and
#'   \code{print()} only). Left at \code{NULL}, a name is derived from
#'   an MD5 hash of the converted code, so two different conversions
#'   never share a build directory/model name even if called back to
#'   back -- \code{mcode_cache()} is normally content-hash aware under
#'   a repeated name (verified,
#'   \code{reports/experiments/09_mcode_cache_test.R}), but two builds
#'   under the same name landing in the same filesystem-timestamp
#'   second can defeat its underlying \code{make}-based rebuild check
#'   (found directly, \code{reports/experiments/14_mcode_cache_staleness_test.R}
#'   -- exactly the pattern an automated test loop hits, converting many
#'   small variants back to back); hashing the default name sidesteps
#'   this rather than relying on \code{mcode_cache()}'s own detection.
#'   A caller-supplied \code{name} does not get this treatment (it is
#'   used verbatim, for predictable caching across sessions) -- avoid
#'   reusing one \code{name} for genuinely different \code{hazard}/
#'   \code{params} content called in quick succession.
#' @param bsv_targets Character vector or \code{NULL} (default): the
#'   parameter names \code{\link{sim_tte_ode}}'s \code{omega} argument
#'   may target via its per-subject-\code{idata} route
#'   (\code{reports/11_bsv_review.md} option (b)). Ignored (with a
#'   \code{message()}) if \code{model} already declares its own
#'   \code{$OMEGA} block -- that model already has a working BSV route
#'   (\code{\link[mrgsolve]{omat}}), and the two never combine for the
#'   same model (\code{reports/11_bsv_review.md} section 4).
#'
#' @return An object of class \code{"simtte_model"}, a list with
#'   \code{code} (the converted model source, printable directly with
#'   \code{cat(x$code, sep = "\n")}), \code{mod} (the compiled
#'   \code{mrgmod}), \code{hazard}, \code{hazard_params} (the fully
#'   resolved parameter list, including the always-added \code{lp}/
#'   \code{U}/\code{END}), \code{bsv_targets}, \code{name}, and
#'   \code{source} (how \code{model} was supplied: \code{"code"},
#'   \code{"file"}, or \code{"mrgmod"}). Pass the whole object as
#'   \code{\link{sim_tte_ode}}'s \code{model} argument.
#'
#' @seealso \code{\link{sim_tte_ode}}, which accepts this function's
#'   return value directly as its \code{model} argument.
#' @export
#' @examples
#' \donttest{
#' # mrgsolve::modlib() itself needs mrgsolve *attached* (library(mrgsolve)),
#' # not just installed; mrgsolve::mread() with an explicit 'project' reaches
#' # the same internal model library without that requirement.
#' irm1 <- mrgsolve::mread(model = "irm1",
#'   project = system.file("models", package = "mrgsolve"), compile = FALSE)
#' tm <- tte_model(irm1,
#'   hazard = "H0 * exp(lp + beta_r * (RESP / (KIN / KOUT) - 1.0))",
#'   params = list(H0 = 0.01, beta_r = 0),
#'   bsv_targets = c("CL", "V2", "KIN", "KOUT"))
#' tm
#' sim <- sim_tte_ode(model = tm, n = 20, end = 30, delta = 5, seed = 1)
#' head(sim$events)
#' }
tte_model <- function(model, hazard, params = list(), name = NULL,
    bsv_targets = NULL) {

    if (inherits(model, "simtte_model")) {
        stop("'model' is already a simtte_model (the output of an ",
            "earlier tte_model() call). Converting a simtte_model ",
            "again is not supported -- call tte_model() on the ",
            "original, unconverted model with the new 'hazard'/",
            "'params' instead.", call. = FALSE)
    }
    if (!is.character(hazard) || length(hazard) != 1L || !nzchar(hazard)) {
        stop("'hazard' must be a single non-empty character string ",
            "(a C++ expression for HAZ).", call. = FALSE)
    }
    params <- as.list(params)
    reserved_param_names <- c("lp", "U", "END")
    always_added <- intersect(names(params), reserved_param_names)
    if (length(always_added)) {
        stop("'params' must not name ",
            paste(sQuote(always_added, q = FALSE), collapse = ", "),
            " -- 'lp', 'U', and 'END' are always added automatically ",
            "by tte_model() (defaults 0, 0, 1e9).", call. = FALSE)
    }

    norm <- .tte_model_normalize_code(model)
    code <- norm$code

    blocks <- mrgsolve::modelparse(code, drop_blank = FALSE,
        keep_mapping = TRUE)
    block_names <- names(blocks)
    .tte_model_check_ode_precondition(block_names)

    .tte_model_reserved_name_check(code, params)

    all_params <- c(params, list(lp = 0, U = 0, END = 1e9))

    ops <- .tte_model_build_ops(code = code, block_names = block_names,
        starts = attr(blocks, "start"), hazard = hazard,
        all_params = all_params)
    new_code <- .tte_model_splice(code, ops)

    name <- if (is.null(name)) .tte_model_default_name(new_code) else name
    mod <- tryCatch(mrgsolve::mcode_cache(model = name, code = new_code),
        error = function(e) {
            stop("tte_model(): the converted model failed to compile: ",
                conditionMessage(e), "\n\n--- converted model code ---\n",
                paste(new_code, collapse = "\n"), call. = FALSE)
        })
    .validate_ode_model_contract(mod)

    has_omega <- nrow(mrgsolve::omat(mod, make = TRUE)) > 0
    if (has_omega && !is.null(bsv_targets)) {
        message("tte_model(): 'model' already declares a $OMEGA block; ",
            "'bsv_targets' is ignored and sim_tte_ode()'s 'omega' will ",
            "be applied via mrgsolve::omat() instead (see ?sim_tte_ode ",
            "\"Between-subject variability\", \"Coexistence with a ",
            "user-supplied model\").")
    }

    structure(list(code = new_code, mod = mod, hazard = hazard,
        hazard_params = all_params, bsv_targets = bsv_targets,
        name = name, source = norm$source), class = "simtte_model")
}

#' Print a converted PK/PD model
#' @param x A \code{"simtte_model"} object.
#' @param ... Ignored.
#' @export
#' @method print simtte_model
print.simtte_model <- function(x, ...) {
    cat("<simtte_model> \"", x$name, "\" (source: ", x$source, ")\n",
        sep = "")
    cat(" hazard: HAZ =", x$hazard, "\n")
    cat(" parameters:\n")
    for (nm in names(x$hazard_params)) {
        cat("  ", nm, "=", x$hazard_params[[nm]], "\n")
    }
    has_omega <- nrow(mrgsolve::omat(x$mod, make = TRUE)) > 0
    if (has_omega) {
        cat(" between-subject variability: via the model's own declared",
            "$OMEGA block (mrgsolve::omat())\n")
    } else if (!is.null(x$bsv_targets)) {
        cat(" between-subject variability: via idata, targets:",
            paste(x$bsv_targets, collapse = ", "), "\n")
    } else {
        cat(" between-subject variability: none configured\n")
    }
    cat(" $code: ", length(x$code), " lines -- cat(x$code, sep = \"\\n\") ",
        "to view\n", sep = "")
    invisible(x)
}

#' Normalize a tte_model() `model` argument to model code text
#'
#' A compiled \code{mrgmod}'s own source is read back from \code{mod@code}
#' (verified directly to round-trip through \code{\link[mrgsolve]{mcode}},
#' \code{reports/experiments/04_model_source_accessor_test.R}) -- code is
#' the canonical representation this function converts, since conversion
#' is a text rewrite (\code{reports/04_author_decisions.md} "After the
#' BSV implementation" decision 3).
#'
#' @param model As documented for \code{\link{tte_model}}.
#' @return A list: \code{code} (character vector, one line per element)
#'   and \code{source} (\code{"file"}, \code{"mrgmod"}, or \code{"code"}).
#' @noRd
.tte_model_normalize_code <- function(model) {
    if (methods::is(model, "mrgmod")) {
        return(list(code = model@code, source = "mrgmod"))
    }
    if (!is.character(model) || !length(model)) {
        stop("'model' must be a file path, a character vector of model ",
            "code, or a compiled mrgmod; see ?tte_model.", call. = FALSE)
    }
    if (length(model) == 1L && !grepl("\n", model, fixed = TRUE) &&
        file.exists(model)) {
        return(list(code = readLines(model, warn = FALSE), source = "file"))
    }
    code <- if (any(grepl("\n", model, fixed = TRUE))) {
        unlist(strsplit(model, "\n", fixed = TRUE), use.names = FALSE)
    } else {
        model
    }
    list(code = code, source = "code")
}

#' Reject a model with no $ODE/$DES block (tte_model() hard precondition)
#'
#' The survival compartment (\code{p11}) must be co-integrated in the
#' same ODE system as the rest of the model
#' (\code{reports/04_author_decisions.md} "After the BSV implementation"
#' decision 4); a \code{$PKMODEL}/\code{$PRED}-only model has no such
#' system to append to.
#'
#' @param block_names Character vector, \code{names()} of
#'   \code{\link[mrgsolve]{modelparse}}'s return value.
#' @return \code{TRUE}, invisibly, if the precondition is satisfied.
#' @noRd
.tte_model_check_ode_precondition <- function(block_names) {
    if (any(c("ODE", "DES") %in% block_names)) {
        return(invisible(TRUE))
    }
    if (any(c("PKMODEL", "PRED") %in% block_names)) {
        stop("tte_model(): this model has a $PKMODEL/$PRED block but no ",
            "$ODE/$DES block. The survival compartment (p11) must be ",
            "co-integrated in the same ODE system as the rest of the ",
            "model, which a $PKMODEL/$PRED-only model does not have; ",
            "use the ODE-form equivalent from mrgsolve::modlib() ",
            "instead (e.g. \"pk2cmt\" instead of \"pk2\"/\"pk2iv\", ",
            "\"irm1\" already has one; see ?mrgsolve::modlib).",
            call. = FALSE)
    }
    stop("tte_model(): this model has no $ODE/$DES block. ",
        "tte_model() adds the survival compartment (p11) as an ODE ",
        "state alongside the model's own equations, which requires ",
        "one to already exist; see ?tte_model.", call. = FALSE)
}

#' Derive a default mcode_cache() model name from the converted code's
#' own content, so two different conversions never share a build
#' directory even if called back to back within the same wall-clock
#' second (\code{reports/experiments/14_mcode_cache_staleness_test.R}
#' -- see \code{\link{tte_model}}'s own \code{name} argument doc for
#' why this is needed rather than relying on
#' \code{\link[mrgsolve]{mcode_cache}}'s own content-hash detection).
#' \code{tools::md5sum()} (base-priority \pkg{tools}, no new dependency)
#' needs a file, so the code is written to a tempfile first.
#' @param code Character vector, the fully assembled converted model code.
#' @return Character scalar, \code{"simtte_tte_model_<12 hex chars>"}.
#' @noRd
.tte_model_default_name <- function(code) {
    tmp <- tempfile(fileext = ".cpp")
    on.exit(unlink(tmp))
    writeLines(code, tmp)
    digest <- unname(tools::md5sum(tmp))
    paste0("simtte_tte_model_", substr(digest, 1L, 12L))
}

#' Names tte_model()'s own scaffold always introduces
#'
#' The nine names every shipped \code{*_hazard.cpp} model's own
#' scaffold uses (\code{p11}, the in-solver latch state/flags) plus
#' \code{lp} -- disclosed addition, not in the literal kickoff-prompt
#' list, but \code{lp} is, like \code{U}/\code{END}, always injected by
#' \code{\link{tte_model}} itself, so a backbone that already declares
#' its own \code{lp} needs the same clash error the other names get
#' (\code{reports/13_converter_design.md} section 3).
#' @noRd
.TTE_MODEL_RESERVED_NAMES <- c("p11", "U", "END", "HAZ", "TEVT",
    "event_found", "T_PRE", "P_PRE", "P_POST", "lp")

#' Check tte_model()'s reserved/params names against the input model's
#' own parameter and compartment names, without compiling C++
#'
#' \code{mrgsolve::mread(..., compile = FALSE)} still fully parses
#' \code{$PARAM}/\code{$CMT}/annotations and returns real
#' \code{param()}/\code{init()} names (verified directly,
#' \code{reports/experiments/08_compile_false_test.R}), so this check
#' costs no C++ compilation. \code{project = tempdir()} is passed
#' explicitly: \code{mread()}'s own default \code{project} is
#' \code{getwd()}, not \code{tempdir()} (unlike \code{mcode()}/
#' \code{mcode_cache()}), and it writes a \code{.cpp} copy of
#' \code{code} into \code{project} even with \code{compile = FALSE} --
#' found via a real \code{R CMD check} NOTE ("non-standard file ...
#' simtte_tte_model_precheck.cpp") before this argument was added.
#'
#' @param code Character vector, the (unedited) input model code.
#' @param params Named list, the caller's \code{tte_model(params = )}.
#' @return \code{TRUE}, invisibly, if no clash is found.
#' @noRd
.tte_model_reserved_name_check <- function(code, params) {
    base <- tryCatch(mrgsolve::mread(model = "simtte_tte_model_precheck",
        code = code, project = tempdir(), compile = FALSE),
        error = function(e) {
            stop("tte_model(): the supplied model failed to parse: ",
                conditionMessage(e), call. = FALSE)
        })
    existing <- c(names(mrgsolve::param(base)), names(mrgsolve::init(base)))

    clash <- intersect(.TTE_MODEL_RESERVED_NAMES, existing)
    if (length(clash)) {
        stop("tte_model(): the supplied model already declares ",
            paste(clash, collapse = ", "), ", which tte_model()'s own ",
            "survival scaffold needs to add itself; rename the ",
            "conflicting parameter(s)/compartment(s) in the source ",
            "model before converting it.", call. = FALSE)
    }
    clash2 <- intersect(names(params), existing)
    if (length(clash2)) {
        stop("tte_model(): 'params' names ",
            paste(clash2, collapse = ", "),
            ", which the supplied model already declares as a ",
            "parameter or compartment; choose different hazard-",
            "parameter name(s).", call. = FALSE)
    }
    invisible(TRUE)
}

#' Build the survival scaffold text fragments (verbatim across every
#' tte_model() conversion; only the $ODE half's HAZ line is
#' model-specific)
#' @param hazard Character scalar, the HAZ expression.
#' @return A named list of character vectors: \code{global}, \code{main},
#'   \code{ode}.
#' @noRd
.tte_model_scaffold_fragments <- function(hazard) {
    list(
        global = c(
            "  // -- BEGIN simtte survival scaffolding (tte_model()) ----",
            "  static int    event_found = 0;",
            "  static double TEVT        = 0.0;",
            "  static double T_PRE       = 0.0;",
            "  static double P_PRE       = 1.0;",
            "  static double P_POST      = 0.0;",
            "  // HAZ is declared here (not as an $ODE-local) so that a",
            "  // brand-new, separate $ODE block can still assign it and",
            "  // have $CAPTURE see it -- mrgsolve's autodec only hoists",
            "  // an $ODE-local into shared storage when it is declared",
            "  // in the SAME $ODE occurrence as the rest of the model's",
            "  // own equations, verified directly not to happen across",
            "  // two separate $ODE (or $ODE + $DES) occurrences (see",
            "  // reports/experiments/10_two_ode_blocks_capture_test.R).",
            "  static double HAZ         = 0.0;",
            "  // -- END simtte survival scaffolding --------------------"),
        main = c(
            "  // -- BEGIN simtte survival scaffolding (tte_model()) ----",
            "  if (NEWIND <= 1) {",
            "    event_found = 0;",
            "    TEVT = 0.0;",
            "    T_PRE = 0.0;",
            "    P_PRE = 1.0;",
            "    P_POST = 0.0;",
            "  }",
            "  // -- END simtte survival scaffolding --------------------"),
        ode = c(
            "  // -- BEGIN simtte survival scaffolding (tte_model()) ----",
            paste0("  HAZ = ", hazard, ";"),
            "  dxdt_p11 = -p11 * HAZ;",
            paste0("  if (!event_found && p11 > U && SOLVERTIME <= END ",
                "&& SOLVERTIME >= T_PRE) {"),
            "    T_PRE = SOLVERTIME;",
            "    P_PRE = p11;",
            "  }",
            "  if (!event_found && p11 <= U && SOLVERTIME <= END) {",
            "    event_found = 1;",
            "    TEVT = SOLVERTIME;",
            "    P_POST = p11;",
            "  }",
            "  // -- END simtte survival scaffolding --------------------")
    )
}

#' Build the (anchor, lines) splice operations for one tte_model() call
#'
#' \code{$GLOBAL}/\code{$MAIN} are mrgsolve singleton blocks
#' (\code{mrgsolve:::block_list_single}, verified directly) so an
#' existing occurrence is edited in place (scaffold lines appended at
#' its own end); an absent one is appended as a brand-new block instead.
#' Every other addition (\code{$ODE}, \code{$PARAM}, \code{$INIT},
#' \code{$CAPTURE}) is always a brand-new block, appended at the file's
#' end, never spliced into the model's own content
#' (\code{reports/13_converter_design.md} section 2).
#'
#' @param code Character vector, the (unedited) input model code.
#' @param block_names,starts \code{names()}/\code{attr(, "start")} of
#'   \code{mrgsolve::modelparse(code, drop_blank = FALSE, keep_mapping
#'   = TRUE)}.
#' @param hazard Character scalar.
#' @param all_params Named list, the fully resolved parameter set
#'   (caller's \code{params} plus \code{lp}/\code{U}/\code{END}).
#' @return A list of \code{list(anchor = integer, lines = character)},
#'   in the order they should appear once spliced (used as the
#'   within-anchor concatenation order by \code{\link{.tte_model_splice}}).
#' @noRd
.tte_model_build_ops <- function(code, block_names, starts, hazard,
    all_params) {
    n <- length(code)
    ends <- c(starts[-1] - 1L, n)
    frag <- .tte_model_scaffold_fragments(hazard)

    # GLOBAL/MAIN are singleton: an existing occurrence is edited in
    # place (anchored at its own end line); an absent one is appended
    # as a brand-new block (anchored at the file's end, with its own
    # header). Everything else this function adds is always brand-new,
    # so those four blocks share the same file-end anchor and are
    # combined into a single splice op below -- .tte_model_splice()
    # would merge them into one anyway (same anchor), so building them
    # as one avoids four near-identical list() calls for no difference
    # in the result.
    .singleton_op <- function(name, header, lines) {
        idx <- which(block_names == name)
        if (!length(idx)) {
            list(anchor = n, lines = c("", header, lines))
        } else {
            list(anchor = max(ends[idx]), lines = lines)
        }
    }

    param_line <- paste("$PARAM", paste(sprintf("%s = %s",
        names(all_params),
        vapply(all_params, format, character(1L), scientific = FALSE,
            trim = TRUE)), collapse = ", "))

    list(
        .singleton_op("GLOBAL", "$GLOBAL", frag$global),
        .singleton_op("MAIN", "$MAIN", frag$main),
        list(anchor = n, lines = c("", "$ODE", frag$ode, "", param_line,
            "", "$INIT", "  p11 = 1", "", "$CAPTURE",
            "  TEVT event_found T_PRE P_PRE P_POST HAZ"))
    )
}

#' Apply a list of (anchor, lines) splice operations to a model code vector
#'
#' Operations sharing the same anchor are concatenated (in the order
#' supplied) into a single insertion so they land in that same relative
#' order; distinct anchors are then applied in descending order so an
#' insertion at a higher line number never invalidates a not-yet-applied
#' anchor at a lower one -- every anchor is computed once, up front,
#' against the original (unmodified) \code{code}
#' (\code{reports/13_converter_design.md} section 2).
#'
#' @param code Character vector, the original model code.
#' @param ops List of \code{list(anchor, lines)}, as built by
#'   \code{\link{.tte_model_build_ops}}.
#' @return Character vector, the spliced model code.
#' @noRd
.tte_model_splice <- function(code, ops) {
    anchors <- vapply(ops, `[[`, integer(1L), "anchor")
    grouped <- split(ops, anchors)
    for (a in sort(unique(anchors), decreasing = TRUE)) {
        lines <- unlist(lapply(grouped[[as.character(a)]], `[[`, "lines"),
            use.names = FALSE)
        code <- append(code, lines, after = a)
    }
    code
}
