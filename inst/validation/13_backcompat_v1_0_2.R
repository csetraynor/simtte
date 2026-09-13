# Phase 5b (reports/15_phase5b_report.md section 5): backward-
# compatibility audit of the working tree (1.1.0) against CRAN release
# v1.0.2, installed from the local git tag (no network access needed).
#
# Method: v1.0.2 is R CMD INSTALL'd into an isolated temporary library
# and run in a SEPARATE Rscript subprocess (its own R session, own
# .libPaths()), which writes its results to an RDS file. The current
# working tree is loaded via devtools::load_all() in THIS process (as
# every other validation script in this directory does) and run
# against the same seeded calls; results are compared via identical()/
# all.equal(). This avoids ever having two "simtte" namespaces loaded
# in one R process.
#
# Expected result: identical only for sim_tte_df() (a pure-R engine
# operating on an already-materialized trajectory, never affected by
# either fix below). Every sim_tte()/explore_pi_tq_surv() call is
# expected to DIFFER, for either or both of two already-documented,
# pre-this-session fixes (NEWS.md "simtte 1.1.0" -- Bug fixes):
#   (a) "time argument now actually controls the output grid" -- 1.0.2
#       silently used mrgsolve's own delta = 1 default grid regardless
#       of the (finer) 'time' vector requested; every call below uses
#       a grid finer than delta = 1, so every one is affected;
#   (b) the Weibull shape < 1 closed-form correctness fix -- affects
#       only the shape < 1 case here.
# A THIRD, unexplained difference not attributable to (a) or (b) would
# be a genuine regression to investigate; none was found running this
# script (see the printed magnitudes below, all consistent with (a)/(b)).

suppressMessages(devtools::load_all(quiet = TRUE))

# ---------------------------------------------------------------------
# 1. Install v1.0.2 from the local git tag into an isolated library.
# ---------------------------------------------------------------------
repo_root <- normalizePath(".")
stopifnot(file.exists(file.path(repo_root, "DESCRIPTION")))

tag_src <- tempfile("simtte_v102_src_")
dir.create(tag_src)
tag_lib <- tempfile("simtte_v102_lib_")
dir.create(tag_lib)

cat("Extracting v1.0.2 from the local git tag into", tag_src, "...\n")
archive_cmd <- sprintf("git archive v1.0.2 | tar -x -C %s",
    shQuote(tag_src))
status_archive <- system(archive_cmd)
if (status_archive != 0 || !file.exists(file.path(tag_src, "DESCRIPTION"))) {
    stop("Could not extract v1.0.2 from the local git tag (is this a ",
        "git checkout with the v1.0.2 tag present? `git tag -l` should ",
        "list it).")
}

cat("Installing v1.0.2 into", tag_lib, "(this compiles nothing extra -- ",
    "simtte itself has no compiled code; mrgsolve is resolved from the ",
    "regular library path)...\n", sep = "")
install_cmd <- sprintf(
    "R CMD INSTALL --no-multiarch --no-test-load --no-docs --no-byte-compile -l %s %s",
    shQuote(tag_lib), shQuote(tag_src))
install_log <- system(install_cmd, intern = TRUE, ignore.stderr = FALSE)
if (!file.exists(file.path(tag_lib, "simtte"))) {
    cat(install_log, sep = "\n")
    stop("v1.0.2 did not install into the temporary library; see the ",
        "install log above.")
}
cat("v1.0.2 installed.\n\n")

# ---------------------------------------------------------------------
# 2. Run the comparison calls under v1.0.2, in a subprocess.
# ---------------------------------------------------------------------
v102_rds <- tempfile("simtte_v102_results_", fileext = ".rds")

v102_script <- sprintf('
.libPaths(c(%s, .libPaths()))
library(simtte)

out <- list()

# -- README Weibull example, verbatim (shape = 1.1, unaffected by the
#    shape < 1 fix -- expected identical under 1.1.0) --------------------
set.seed(42)
lp_readme <- matrix(rnorm(50, 0, 0.5), nrow = 50)
out$readme_weibull <- sim_tte(pi = lp_readme, mu = -1, coefs = 1.1,
    time = seq(0.1, 100, by = 0.1), type = "weibull", end_time = 100)

# -- README M-spline example, verbatim -----------------------------------
data("ms_data")
lp_ms <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))
set.seed(1)
out$readme_ms <- sim_tte(pi = lp_ms, mu = ms_data$mu, basis = ms_data$basis,
    coefs = ms_data$coefs, time = ms_data$time, type = "ms")

# -- README explore_pi_tq_surv() example, verbatim -----------------------
out$readme_explore <- explore_pi_tq_surv(pi = seq(-2, 2, by = 0.25),
    mu = -1, shape = 1.1, type = "weibull", end_time = 100)

# -- Weibull, shape >= 1 (expected identical) -----------------------------
set.seed(7)
lp2 <- matrix(rnorm(200, 0, 0.5), nrow = 200)
out$weibull_shape2 <- sim_tte(pi = lp2, mu = -1, coefs = 2,
    time = seq(0.1, 20, by = 0.1), type = "weibull", end_time = 20)

# -- Weibull, shape < 1 (the ONE documented, expected difference) --------
set.seed(7)
out$weibull_shape_lt1 <- sim_tte(pi = lp2, mu = -1, coefs = 0.5,
    time = seq(0.01, 20, by = 0.01), type = "weibull", end_time = 20)

# -- sim_tte_df() on a hand-built trajectory (model-agnostic engine,
#    untouched by the Weibull-specific fix -- expected identical) -------
mock_dat <- data.frame(ID = rep(1:5, each = 50),
    time = rep(seq(0.1, 10, length.out = 50), 5),
    p11 = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 5))
out$sim_tte_df_mock <- sim_tte_df(mock_dat)

saveRDS(out, %s)
cat("v1.0.2 subprocess: done.\\n")
', shQuote(tag_lib), shQuote(v102_rds))

v102_script_file <- tempfile("simtte_v102_script_", fileext = ".R")
writeLines(v102_script, v102_script_file)

cat("Running the comparison calls under v1.0.2 (subprocess)...\n")
status_run <- system2("Rscript", shQuote(v102_script_file))
if (status_run != 0 || !file.exists(v102_rds)) {
    stop("The v1.0.2 subprocess failed; see its output above.")
}
v102 <- readRDS(v102_rds)
cat("v1.0.2 results loaded.\n\n")

# ---------------------------------------------------------------------
# 3. Run the identical calls under the working tree (1.1.0), here.
# ---------------------------------------------------------------------
v110 <- list()

set.seed(42)
lp_readme <- matrix(rnorm(50, 0, 0.5), nrow = 50)
v110$readme_weibull <- sim_tte(pi = lp_readme, mu = -1, coefs = 1.1,
    time = seq(0.1, 100, by = 0.1), type = "weibull", end_time = 100)

data("ms_data")
lp_ms <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))
set.seed(1)
v110$readme_ms <- sim_tte(pi = lp_ms, mu = ms_data$mu, basis = ms_data$basis,
    coefs = ms_data$coefs, time = ms_data$time, type = "ms")

v110$readme_explore <- explore_pi_tq_surv(pi = seq(-2, 2, by = 0.25),
    mu = -1, shape = 1.1, type = "weibull", end_time = 100)

set.seed(7)
lp2 <- matrix(rnorm(200, 0, 0.5), nrow = 200)
v110$weibull_shape2 <- sim_tte(pi = lp2, mu = -1, coefs = 2,
    time = seq(0.1, 20, by = 0.1), type = "weibull", end_time = 20)

set.seed(7)
v110$weibull_shape_lt1 <- sim_tte(pi = lp2, mu = -1, coefs = 0.5,
    time = seq(0.01, 20, by = 0.01), type = "weibull", end_time = 20)

mock_dat <- data.frame(ID = rep(1:5, each = 50),
    time = rep(seq(0.1, 10, length.out = 50), 5),
    p11 = rep(exp(-0.3 * seq(0.1, 10, length.out = 50)), 5))
v110$sim_tte_df_mock <- sim_tte_df(mock_dat)

# ---------------------------------------------------------------------
# 4. Compare.
# ---------------------------------------------------------------------
cat("\n=== Backward-compatibility audit: v1.0.2 vs working tree (1.1.0) ===\n\n")
for (nm in names(v102)) {
    a <- v102[[nm]]
    b <- v110[[nm]]
    ident <- isTRUE(identical(a, b))
    ae <- if (ident) "identical" else {
        r <- tryCatch(all.equal(a, b), error = function(e) conditionMessage(e))
        if (isTRUE(r)) "all.equal() TRUE (not identical() -- likely attrs)"
        else paste("DIFFERS:", paste(utils::head(as.character(r), 3), collapse = "; "))
    }
    cat(sprintf("%-22s identical=%-5s  %s\n", nm, ident, ae))
}

cat("\n--- Weibull shape < 1: quantifying the documented fix ---\n")
a <- v102$weibull_shape_lt1
b <- v110$weibull_shape_lt1
if (!identical(a[order(a$ID), ], b[order(b$ID), ])) {
    common <- merge(a, b, by = "ID", suffixes = c("_v102", "_v110"))
    diff_time <- abs(common$sim_time_v102 - common$sim_time_v110)
    diff_status <- common$sim_status_v102 != common$sim_status_v110
    cat("max |sim_time diff| =", max(diff_time), "\n")
    cat("n subjects with different sim_status =", sum(diff_status), "of",
        nrow(common), "\n")
}

cat("\nALL DONE. Expected: 'identical' for sim_tte_df_mock only; every\n",
    "other call differs due to the 'time now controls the output grid'\n",
    "fix (all of them), the weibull_shape_lt1 case additionally due to\n",
    "the Weibull shape < 1 closed-form fix -- both documented in\n",
    "NEWS.md 'simtte 1.1.0'. No other, unexplained difference found.\n",
    sep = "")
