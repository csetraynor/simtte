# Formula-based linear-predictor interface for sim_tte_ode()'s
# 'covariates'/'beta' mechanism ("Before Phase 7: covariate interface",
# reports/04_author_decisions.md; reports/25_covariate_interface_review.md).
#
# This file adds exactly one thing: a way to turn (formula, data) into
# a data.frame whose columns are already named to match a coefficient
# vector, so it can be handed to the existing, unmodified
# .build_ode_covariate_rows() (R/helpers.R) -- see the review's section
# 6 for why no change to that pipeline was needed. Everything about
# canonicalization, per-subject validation, coverage, and the
# evid = 1 / amt = 0 / cmt covariate-update row shape stays exactly
# what it already was.

#' Build a linear-predictor design matrix from a formula
#'
#' Ordinary \code{stats::model.matrix()}, with two departures from its
#' defaults, both load-bearing (see
#' \code{reports/25_covariate_interface_review.md} section 2):
#' \code{na.action = stats::na.fail} (the default, \code{na.omit},
#' silently drops rows and would misalign the design matrix against the
#' \code{ID}/\code{time} columns carried alongside it), and an implied
#' intercept column is dropped after being built (not avoided by
#' forcing \code{~ 0 + ...} from the start), matching
#' \code{survival::coxph()}'s own convention -- see "Recommendation (b)"
#' in the review.
#'
#' @param formula A formula, e.g. \code{~ age + arm}.
#' @param data A data frame containing every variable \code{formula}
#'   references.
#' @return A numeric matrix (\code{stats::model.matrix()} output, minus
#'   any \code{"(Intercept)"} column).
#' @noRd
.build_lp_design_matrix <- function(formula, data) {
    if (!inherits(formula, "formula")) {
        stop("'formula' must be a formula, e.g. ~ age + arm; see ",
            "?sim_tte_ode \"Covariates and the linear predictor\".",
            call. = FALSE)
    }
    data <- as.data.frame(data)
    mf <- tryCatch(
        stats::model.frame(formula, data = data, na.action = stats::na.fail),
        error = function(e) {
            stop("Could not build the design matrix from 'formula': ",
                conditionMessage(e), call. = FALSE)
        })
    X <- stats::model.matrix(formula, data = mf)
    if ("(Intercept)" %in% colnames(X)) {
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        message("sim_tte_ode(): dropped the intercept column implied by ",
            "'formula' -- the model's own baseline hazard parameter ",
            "(e.g. 'H0'/'mu') already plays that role, so an intercept ",
            "in 'lp' would be additively unidentifiable against it. ",
            "Write 'formula = ~ 0 + ...' to build the design matrix ",
            "without this message. See ?sim_tte_ode \"Covariates and ",
            "the linear predictor\".")
    } else {
        full_level_factors <- .full_level_factor_terms(formula, mf, X)
        if (length(full_level_factors)) {
            message("sim_tte_ode(): 'formula' has no intercept, so ",
                if (length(full_level_factors) > 1L) "factor(s) " else "factor ",
                paste(sprintf("'%s'", full_level_factors), collapse = ", "),
                " expanded to one design-matrix column per level, not ",
                "the usual reference-level contrast (this is the ",
                "intended effect of writing '~ 0 + ...', not an error). ",
                "See ?sim_tte_ode \"Covariates and the linear predictor\".")
        }
    }
    if (!ncol(X)) {
        stop("'formula' produces an empty design matrix (no covariate ",
            "columns remain after removing the intercept); supply at ",
            "least one covariate term.", call. = FALSE)
    }
    X
}

#' Which main-effect factor terms got a full set of per-level columns
#'
#' With no intercept, R's own \code{stats::model.matrix()} gives the
#' \emph{first} factor term a full-rank set of per-level dummy columns
#' (\code{nlevels(x)} columns, not the usual \code{nlevels(x) - 1}
#' reference-level contrast) and keeps ordinary contrast coding for
#' every other factor term, to keep the whole design full rank --
#' verified directly, not assumed. Named here so
#' \code{\link{.build_lp_design_matrix}} can message about it (decision
#' recorded under "After the covariate interface",
#' \code{reports/04_author_decisions.md} decision 2) instead of leaving
#' a silent difference between "the usual `k - 1` columns" and "all `k`
#' columns" for a user to puzzle out from `colnames(X)` alone.
#'
#' @param formula The formula \code{X} was built from.
#' @param mf The \code{stats::model.frame()} used to build \code{X}.
#' @param X The design matrix (with its \code{"assign"} attribute
#'   intact, i.e. straight from \code{stats::model.matrix()}).
#' @return Character vector of variable names (possibly empty) that got
#'   full-level coding.
#' @noRd
.full_level_factor_terms <- function(formula, mf, X) {
    term_labels <- attr(stats::terms(formula, data = mf), "term.labels")
    assign <- attr(X, "assign")
    out <- character(0)
    for (i in seq_along(term_labels)) {
        v <- term_labels[i]
        # Skips interactions/transforms (e.g. "age:arm", "I(age^2)"):
        # only a bare column name is ever a factor term worth naming
        # here.
        if (!v %in% names(mf)) {
            next
        }
        col <- mf[[v]]
        if (!(is.factor(col) || is.character(col))) {
            next
        }
        if (sum(assign == i) == nlevels(as.factor(col))) {
            out <- c(out, v)
        }
    }
    out
}

#' Validate that `beta`'s names match a design matrix's columns exactly
#'
#' Two-directional, unlike \code{.build_ode_covariate_rows()}'s own
#' (legacy, unchanged) one-directional check -- see
#' \code{reports/25_covariate_interface_review.md} section 2 for why an
#' extra, unmatched \code{beta} name is worth catching explicitly for
#' the formula path (its names come from \code{model.matrix()}, not
#' from the user's own column names, so a mismatch is more likely to be
#' a real mistake than a deliberate partial-coefficient set).
#'
#' @param beta Named numeric vector.
#' @param X Design matrix from \code{\link{.build_lp_design_matrix}}.
#' @return \code{TRUE}, invisibly, if valid.
#' @noRd
.validate_beta_matches_design <- function(beta, X) {
    if (!is.numeric(beta) || is.null(names(beta)) || any(!nzchar(names(beta)))) {
        stop("'beta' must be a named numeric vector when 'formula' is ",
            "supplied; every name must match a column of the design ",
            "matrix built from 'formula'. Expected coefficient name(s): ",
            paste(colnames(X), collapse = ", "), ".", call. = FALSE)
    }
    expected <- colnames(X)
    missing <- setdiff(expected, names(beta))
    extra <- setdiff(names(beta), expected)
    if (length(missing) || length(extra)) {
        stop("'beta' does not match the design matrix built from ",
            "'formula'. Expected coefficient name(s): ",
            paste(expected, collapse = ", "), ".",
            if (length(missing)) paste0(" Missing from 'beta': ",
                paste(missing, collapse = ", "), ".") else "",
            if (length(extra)) paste0(" Not a column of the design ",
                "matrix: ", paste(extra, collapse = ", "), ".") else "",
            call. = FALSE)
    }
    invisible(TRUE)
}

#' Build sim_tte_ode() covariate-update rows from a formula + data frame
#'
#' The formula-based generalization of \code{covariates}/\code{beta}
#' (\code{reports/25_covariate_interface_review.md}): builds a design
#' matrix via \code{\link{.build_lp_design_matrix}}, matches \code{beta}
#' against it via \code{\link{.validate_beta_matches_design}}, and hands
#' the result to \code{\link{.build_ode_covariate_rows}} unchanged --
#' every canonicalization/validation/coverage rule and the covariate-
#' update row shape itself are exactly what the legacy (non-formula)
#' path already uses.
#'
#' A \code{data} frame with no \code{time} column is treated as
#' baseline-only: equivalent to a single \code{time = 0} observation per
#' subject, carried forward for the whole follow-up by the same
#' last-observation-carried-forward rule \code{.check_lp_data_coverage()}
#' already applies (its coverage requirement is exactly "an observation
#' at time = 0", so this needs no change to that function -- see the
#' review section 2). Unlike the time-varying shape's own "no ID
#' column" convention (one shared trajectory, recycled to everyone),
#' a baseline frame with no \code{ID} column is one row \emph{per
#' subject} when it has more than one row, or a population-level
#' constant recycled to every subject when it has exactly one row (the
#' decision recorded under "After the covariate interface",
#' \code{reports/04_author_decisions.md}) -- either way, \code{ID} is
#' filled positionally before \code{\link{.canonicalize_lp_data}} ever
#' sees it. \code{n_subjects} itself, when \code{formula} is used with a
#' multi-row, no-\code{ID} baseline frame and the caller did not set
#' \code{sim_tte_ode()}'s own \code{n} explicitly, is inferred from this
#' frame's row count one level up, by
#' \code{\link{.resolve_n_for_baseline_covariates}} -- not here, since
#' this function only ever sees the already-resolved \code{n_subjects}.
#'
#' @param formula A formula, e.g. \code{~ age + arm}.
#' @param data Covariate data frame: baseline (one row per subject, no
#'   \code{time} column) or time-varying (\code{ID}, \code{time},
#'   covariate columns) -- see \code{?sim_tte_ode} "Covariates and the
#'   linear predictor".
#' @param beta Named numeric vector matching \code{colnames()} of the
#'   design matrix built from \code{formula}.
#' @param n_subjects,end,cmt Forwarded to
#'   \code{\link{.build_ode_covariate_rows}} unchanged.
#' @return Data frame of covariate-update rows, exactly as
#'   \code{\link{.build_ode_covariate_rows}} returns.
#' @noRd
.build_ode_covariate_rows_formula <- function(formula, data, beta,
    n_subjects, end, cmt = 1L) {
    data <- as.data.frame(data)
    if (!"time" %in% names(data)) {
        if (!"ID" %in% names(data)) {
            # decision ("After the covariate interface",
            # reports/04_author_decisions.md): one row with no ID is a
            # population-level constant, recycled to every subject; more
            # than one row is already one row per subject.
            # rep(..., length.out = n_rows) collapses both into the same
            # recycle-then-number step instead of two separate branches.
            n_rows <- if (nrow(data) == 1L) n_subjects else nrow(data)
            data <- data[rep(seq_len(nrow(data)), length.out = n_rows), ,
                drop = FALSE]
            data$ID <- seq_len(n_rows)
        }
        data$time <- 0
    }
    X <- .build_lp_design_matrix(formula, data)
    .validate_beta_matches_design(beta, X)

    # check.names = FALSE: model.matrix() column names for an
    # interaction term (e.g. "age:armB") or a transform (e.g.
    # "I(age^2)") must survive unmangled, so they still match
    # colnames(X) (and therefore names(beta)) exactly once turned into
    # a data.frame.
    design <- as.data.frame(X, check.names = FALSE)
    design$time <- data$time
    if ("ID" %in% names(data)) {
        design$ID <- data$ID
    }

    .build_ode_covariate_rows(design, beta[colnames(X)], n_subjects, end,
        cmt = cmt)
}

#' Resolve `sim_tte_ode()`'s `n` against a baseline covariates frame
#'
#' Decision recorded under "After the covariate interface"
#' (\code{reports/04_author_decisions.md}): a baseline (no \code{time}
#' column) \code{covariates} frame with more than one row and no own
#' \code{ID} column defines the subject count directly -- inferred from
#' \code{nrow(covariates)} when the caller left \code{n} at its default,
#' checked for an exact match when the caller set \code{n} explicitly.
#' A one-row frame is a population-level constant (\code{n} is
#' untouched -- \code{\link{.build_ode_covariate_rows_formula}} recycles
#' the single row to whatever \code{n} already resolved to). Only
#' relevant when \code{idata} is \code{NULL}: once \code{idata} exists,
#' its own row count already \emph{is} \code{n_subjects}, and the
#' existing "must equal 1:n_subjects exactly" check downstream does the
#' equivalent job for that case.
#'
#' @param covariates,formula,idata As supplied to \code{sim_tte_ode()}.
#' @param n \code{sim_tte_ode()}'s own \code{n}, as supplied or
#'   defaulted.
#' @param n_explicit Logical: did the caller supply \code{n}
#'   themselves (\code{!missing(n)} in \code{sim_tte_ode()})?
#' @return The resolved \code{n}.
#' @noRd
.resolve_n_for_baseline_covariates <- function(covariates, formula, idata,
    n, n_explicit) {
    if (!is.null(idata) || is.null(covariates) || is.null(formula)) {
        return(n)
    }
    covariates <- as.data.frame(covariates)
    if ("time" %in% names(covariates) || "ID" %in% names(covariates) ||
        nrow(covariates) <= 1L) {
        return(n)
    }
    n_cov <- nrow(covariates)
    if (n_explicit && n != n_cov) {
        stop("'n' (", n, ") does not match the number of rows in the ",
            "baseline 'covariates' data frame (", n_cov, "); a baseline ",
            "covariate frame (no 'time' column, no 'ID' column) with ",
            "more than one row defines the subject count directly. ",
            "Supply a matching 'n', or omit 'n' and let it be inferred.",
            call. = FALSE)
    }
    n_cov
}
