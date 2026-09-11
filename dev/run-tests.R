#!/usr/bin/env Rscript
# Targeted test runner for simtte. See reports/07_test_runbook.md for
# the full guide; this script's own --help covers the mechanics.
#
#   Rscript dev/run-tests.R                  # fast default: everything, slow tests skipped
#   Rscript dev/run-tests.R ode              # all test-sim-tte-ode-* files
#   Rscript dev/run-tests.R weibull          # one file by name (exponential/weibull/gompertz/covariates)
#   Rscript dev/run-tests.R legacy           # everything except the new API
#   Rscript dev/run-tests.R ode --slow       # sets SIMTTE_SLOW_TESTS=true
#   Rscript dev/run-tests.R --check          # rcmdcheck --as-cran --no-manual, slow tests on
#
# No CLI framework: commandArgs() parsed directly. devtools + testthat
# (already used throughout package development; not declared in
# DESCRIPTION since dev/ is a development-only, .Rbuildignore'd
# directory, same status as inst/validation/).

args <- commandArgs(trailingOnly = TRUE)
slow <- "--slow" %in% args
check <- "--check" %in% args
group <- setdiff(args, c("--slow", "--check"))
group <- if (length(group)) group[1] else "all"

if (slow || check) {
    Sys.setenv(SIMTTE_SLOW_TESTS = "true")
}
Sys.setenv(NOT_CRAN = "true")

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root (Rscript dev/run-tests.R ...), ",
        "not from inside dev/.", call. = FALSE)
}

if (check) {
    cat("Running R CMD check --as-cran --no-manual (SIMTTE_SLOW_TESTS=true)...\n")
    res <- rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"),
        error_on = "never")
    print(res)
    quit(status = if (length(res$errors)) 1L else 0L)
}

suppressMessages(devtools::load_all(quiet = TRUE))

test_dir <- "tests/testthat"
all_files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)
ode_files <- grep("^test-sim-tte-ode-", basename(all_files), value = TRUE)
ode_files <- file.path(test_dir, ode_files)

files <- switch(group,
    all = all_files,
    ode = ode_files,
    legacy = setdiff(all_files, ode_files),
    exponential = ,
    weibull = ,
    gompertz = ,
    covariates = file.path(test_dir, paste0("test-sim-tte-ode-", group, ".R")),
    stop("Unknown group '", group, "'. Use one of: all, ode, legacy, ",
        "exponential, weibull, gompertz, covariates.", call. = FALSE)
)
files <- files[file.exists(files)]
if (!length(files)) {
    stop("No test files matched group '", group, "'.", call. = FALSE)
}

cat("Group:", group, " Slow tests:", slow, " Files:", length(files), "\n\n")

total_fail <- 0L
for (f in files) {
    t0 <- Sys.time()
    res <- testthat::test_file(f, reporter = "silent")
    elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 2)
    df <- as.data.frame(res)
    n_pass <- sum(df$passed)
    n_fail <- sum(df$failed)
    n_skip <- sum(df$skipped)
    total_fail <- total_fail + n_fail
    status <- if (n_fail > 0) "FAIL" else "ok"
    cat(sprintf("[%s] %-42s pass=%-4d fail=%-3d skip=%-3d %6.2fs\n",
        status, basename(f), n_pass, n_fail, n_skip, elapsed))
}

cat("\n", if (total_fail > 0) "FAILED" else "ALL OK", "\n", sep = "")
# quit(status = if (total_fail > 0) 1L else 0L)
