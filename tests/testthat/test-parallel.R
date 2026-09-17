# "simttepower feedback 2: worker model loading"
# (reports/04_author_decisions.md; reports/30_simttepower_feedback_worker_loading.md):
# simtte_model_cache()/simtte_prepare_model() give simttepower (and any
# other parallel consumer) a supported, in-package way to build a model
# once and reuse it across PSOCK workers, instead of reaching for
# mrgsolve::mcode_cache()/loadso() directly with its own retry loop.
# Fast tests only -- see test-slow-parallel.R for the real 2-worker
# PSOCK cluster coverage (SIMTTE_SLOW_TESTS=true, skip_on_cran()).

# ---- 1. simtte_model_cache() ----

test_that("simtte_model_cache() defaults to a directory under tempdir()", {
    dir <- simtte_model_cache()
    expect_true(dir.exists(dir))
    expect_true(startsWith(normalizePath(dir), normalizePath(tempdir())))
})

test_that("simtte_model_cache()'s directory name encodes R/simtte/mrgsolve versions", {
    dir <- simtte_model_cache()
    tag <- basename(dir)
    expect_match(tag, as.character(getRversion()), fixed = TRUE)
    expect_match(tag, as.character(utils::packageVersion("simtte")), fixed = TRUE)
    expect_match(tag, as.character(utils::packageVersion("mrgsolve")), fixed = TRUE)
})

test_that("simtte_model_cache() honors a custom simtte.cache_dir option, still version-keyed", {
    custom <- file.path(tempdir(), "my-custom-simtte-cache")
    withr_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = withr_opt))
    options(simtte.cache_dir = custom)
    dir <- simtte_model_cache()
    expect_true(startsWith(normalizePath(dir), normalizePath(custom)))
})

# ---- 2. simtte_prepare_model(): identical output, no rebuild needed ----

test_that("sim_tte_ode() with a prepared built-in model matches an unprepared call exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    prepared <- simtte_prepare_model("exponential")
    expect_s3_class(prepared, "simtte_prepared_model")
    sim1 <- sim_tte_ode(model = prepared, param = list(H0 = 0.1), n = 10,
        end = 20, delta = 2, seed = 1)
    sim2 <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 10, end = 20, delta = 2, seed = 1)
    expect_identical(sim1$events, sim2$events)
})

test_that("sim_tte_ode() with a prepared tte_model()-converted model matches an unprepared call exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    tm <- .tte_model_fixture()
    prepared <- simtte_prepare_model(tm)
    expect_s3_class(prepared, "simtte_prepared_model")
    sim1 <- sim_tte_ode(model = prepared, n = 5, end = 20, delta = 2, seed = 1)
    sim2 <- sim_tte_ode(model = tm, n = 5, end = 20, delta = 2, seed = 1)
    expect_identical(sim1$events, sim2$events)
})

test_that("sim_tte_ode() with a prepared mspline model matches an unprepared call exactly", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    knots <- c(5, 10, 15)
    coefs <- rep(1, 6)
    prepared <- simtte_prepare_model("mspline", knots = knots)
    sim1 <- sim_tte_ode(model = prepared, knots = knots, coefs = coefs,
        param = list(mu = -1), n = 5, end = 20, delta = 2, seed = 1)
    sim2 <- sim_tte_ode(model = "mspline", knots = knots, coefs = coefs,
        param = list(mu = -1), n = 5, end = 20, delta = 2, seed = 1)
    expect_identical(sim1$events, sim2$events)
})

test_that("simtte_prepare_model() is idempotent: preparing an already-prepared object round-trips", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    prepared <- simtte_prepare_model("exponential")
    prepared2 <- simtte_prepare_model(prepared)
    expect_s3_class(prepared2, "simtte_prepared_model")
    expect_identical(prepared2$spec, prepared$spec)
    sim <- sim_tte_ode(model = prepared2, param = list(H0 = 0.1), n = 5,
        end = 10, delta = 2, seed = 1)
    expect_equal(nrow(sim$events), 5)
})

test_that("an mspline prepared model rejects a mismatched 'knots' count at simulate time", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    prepared <- simtte_prepare_model("mspline", knots = c(5, 10, 15))
    expect_error(
        sim_tte_ode(model = prepared, knots = c(3, 6, 9, 12, 15),
            coefs = rep(1, 8), param = list(mu = -1), n = 5, end = 20,
            delta = 2, seed = 1),
        "must match")
})

test_that("simtte_prepare_model() rejects an invalid 'model' with the same message sim_tte_ode() uses", {
    expect_error(simtte_prepare_model(list(a = 1)),
        "must be a character library-model name")
    expect_error(simtte_prepare_model(42),
        "must be a character library-model name")
})

test_that("simtte_prepare_model(\"mspline\") without 'knots' errors clearly", {
    expect_error(simtte_prepare_model("mspline"), "requires 'knots'")
})

test_that("simtte_prepare_model() has no 'coefs'/'boundary_knots' parameter -- supplying either errors, naming it", {
    expect_error(
        simtte_prepare_model("mspline", knots = c(1, 2, 3), coefs = rep(1, 6)),
        "unused argument")
    expect_error(
        simtte_prepare_model("mspline", knots = c(1, 2, 3),
            boundary_knots = c(0, 10)),
        "unused argument")
})

test_that("print.simtte_prepared_model() names the underlying model", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    prepared <- simtte_prepare_model("weibull")
    expect_output(print(prepared), "exponential|weibull")
})

# ---- 3. No global mrgsolve option is ever touched ----

test_that("simtte_prepare_model()/sim_tte_ode()/simtte_model_cache() never set mrgsolve.project/soloc", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    before_project <- getOption("mrgsolve.project")
    before_soloc <- getOption("mrgsolve.soloc")

    simtte_model_cache()
    prepared <- simtte_prepare_model("exponential")
    sim_tte_ode(model = prepared, param = list(H0 = 0.1), n = 3, end = 10,
        delta = 2, seed = 1)
    sim_tte_ode(model = "exponential", param = list(H0 = 0.1), n = 3,
        end = 10, delta = 2, seed = 1)
    tm <- .tte_model_fixture()
    simtte_prepare_model(tm)

    expect_identical(getOption("mrgsolve.project"), before_project)
    expect_identical(getOption("mrgsolve.soloc"), before_soloc)
})

# ---- 4. simtte_model_cache_clear() ----

test_that("simtte_model_cache_clear() removes only the current version subdirectory by default", {
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    base <- file.path(tempdir(), "simtte-cache-clear-test-1")
    options(simtte.cache_dir = base)

    current <- simtte_model_cache()
    sibling <- file.path(base, "R0.0.0-simtte0.0.0-mrgsolve0.0.0")
    dir.create(sibling, recursive = TRUE)

    removed <- simtte_model_cache_clear()
    expect_identical(removed, current)
    expect_false(dir.exists(current))
    expect_true(dir.exists(sibling))
})

test_that("simtte_model_cache_clear(all = TRUE) removes every version subdirectory, not the base", {
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    base <- file.path(tempdir(), "simtte-cache-clear-test-2")
    options(simtte.cache_dir = base)

    current <- simtte_model_cache()
    sibling <- file.path(base, "R0.0.0-simtte0.0.0-mrgsolve0.0.0")
    dir.create(sibling, recursive = TRUE)
    unrelated <- file.path(base, "not-a-version-dir")
    dir.create(unrelated, recursive = TRUE)

    removed <- simtte_model_cache_clear(all = TRUE)
    expect_setequal(removed, c(current, sibling))
    expect_false(dir.exists(current))
    expect_false(dir.exists(sibling))
    expect_true(dir.exists(base))
    expect_true(dir.exists(unrelated))
})

test_that("simtte_model_cache_clear() refuses to remove a target outside the resolved cache base", {
    skip_on_os("windows")
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    base <- file.path(tempdir(), "simtte-cache-clear-test-3")
    dir.create(base, recursive = TRUE)
    options(simtte.cache_dir = base)

    outside <- file.path(tempdir(), "simtte-cache-clear-test-3-outside")
    dir.create(outside, recursive = TRUE)
    escape_file <- file.path(outside, "sentinel")
    file.create(escape_file)

    link <- file.path(base, "R9-simtte9-mrgsolve9")
    ok <- suppressWarnings(file.symlink(outside, link))
    skip_if_not(isTRUE(ok), "cannot create symlinks in this environment")

    expect_error(simtte_model_cache_clear(all = TRUE),
        "refused to remove")
    expect_true(file.exists(escape_file))
})

test_that("simtte_model_cache_clear() is idempotent on an empty/nonexistent cache", {
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    options(simtte.cache_dir = file.path(tempdir(),
        "simtte-cache-clear-test-empty"))

    expect_identical(simtte_model_cache_clear(), character(0))
    expect_identical(simtte_model_cache_clear(all = TRUE), character(0))
})

test_that("after simtte_model_cache_clear(), the next simtte_prepare_model() rebuilds and simulates identically", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    options(simtte.cache_dir = file.path(tempdir(),
        "simtte-cache-clear-test-rebuild"))

    prepared1 <- simtte_prepare_model("exponential")
    sim1 <- sim_tte_ode(model = prepared1, param = list(H0 = 0.1), n = 5,
        end = 10, delta = 2, seed = 1)

    cache_dir <- simtte_model_cache()
    removed <- simtte_model_cache_clear()
    expect_identical(removed, cache_dir)
    expect_false(dir.exists(cache_dir))

    prepared2 <- simtte_prepare_model("exponential")
    sim2 <- sim_tte_ode(model = prepared2, param = list(H0 = 0.1), n = 5,
        end = 10, delta = 2, seed = 1)
    expect_identical(sim1$events, sim2$events)
})

# ---- 5. simtte_model_cache(): a read-only base gives a clear error ----

test_that("simtte_model_cache() errors clearly (not via a downstream mrgsolve failure) when its base is not writable", {
    skip_on_os(c("windows", "solaris"))
    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    base <- file.path(tempdir(), "simtte-cache-readonly-test")
    dir.create(base, recursive = TRUE)
    Sys.chmod(base, mode = "0500")
    on.exit(Sys.chmod(base, mode = "0700"), add = TRUE)

    probe <- isTRUE(suppressWarnings(dir.create(file.path(base, "probe"))))
    skip_if(probe, "cannot create a non-writable directory in this environment (e.g. running as root)")

    options(simtte.cache_dir = file.path(base, "nested"))
    expect_error(simtte_model_cache(), "writable")
})
