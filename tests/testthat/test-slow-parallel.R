# Real 2-worker PSOCK cluster coverage for "simttepower feedback 2:
# worker model loading" (see test-parallel.R for the fast, no-cluster
# tests). CRAN policy: no test may spawn more than 2 workers.
#
# A PSOCK worker is a separate process and cannot see this session's
# own devtools::load_all()-dev-loaded simtte (find.package("simtte")
# resolves differently there -- verified directly) -- every test file
# in this suite runs under load_all() (dev/run-tests.R's own
# mechanism), so every test below re-loads simtte on each worker from
# the same source path before using it. A no-op under a properly
# installed simtte (R CMD check's own test run: find.package() then
# resolves to the installed library, which has a Meta/package.rds, and
# every worker's ordinary lazy-loading already finds the right one).
.simtte_pkg_path <- find.package("simtte")
.simtte_is_dev_loaded <- !file.exists(file.path(.simtte_pkg_path, "Meta",
    "package.rds"))
.load_simtte_on_worker <- function(pkg_path, dev_loaded) {
    if (dev_loaded) {
        pkgload::load_all(pkg_path, quiet = TRUE)
    }
    invisible(NULL)
}

test_that("2-worker PSOCK cluster: prepare-at-init + parLapply(sim_tte_ode()), built-in model, repeated on a fresh cache, no retry, matches sequential", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    if (.simtte_is_dev_loaded) skip_if_not_installed("pkgload")
    skip_if_not_slow()

    seeds <- 1:6
    seq_events <- lapply(seeds, function(s) {
        sim_tte_ode(model = "exponential", param = list(H0 = 0.1), n = 10,
            end = 20, delta = 2, seed = s)$events
    })

    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    for (rep in 1:3) {
        # A fresh cache dir each repetition -- no leftover build from a
        # previous iteration to accidentally paper over a race.
        options(simtte.cache_dir = file.path(tempdir(),
            paste0("slow-parallel-builtin-", rep)))

        cl <- parallel::makePSOCKcluster(2)
        parallel::clusterCall(cl, .load_simtte_on_worker,
            pkg_path = .simtte_pkg_path, dev_loaded = .simtte_is_dev_loaded)
        parallel::clusterCall(cl, simtte::simtte_prepare_model,
            model = "exponential")
        par_events <- parallel::parLapply(cl, seeds, function(s) {
            simtte::sim_tte_ode(model = "exponential",
                param = list(H0 = 0.1), n = 10, end = 20, delta = 2,
                seed = s)$events
        })
        parallel::stopCluster(cl)

        for (i in seq_along(seeds)) {
            expect_identical(par_events[[i]], seq_events[[i]],
                label = paste0("rep=", rep, " seed=", seeds[i]))
        }
    }
})

test_that("2-worker PSOCK cluster: prepare-at-init + parLapply(sim_tte_ode()), tte_model()-converted, repeated on a fresh cache, no retry, matches sequential", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    if (.simtte_is_dev_loaded) skip_if_not_installed("pkgload")
    skip_if_not_slow()

    seeds <- 1:6
    seq_events <- lapply(seeds, function(s) {
        tm <- .tte_model_fixture()
        sim_tte_ode(model = tm, n = 10, end = 20, delta = 2, seed = s)$events
    })

    old_opt <- getOption("simtte.cache_dir")
    on.exit(options(simtte.cache_dir = old_opt))
    for (rep in 1:3) {
        options(simtte.cache_dir = file.path(tempdir(),
            paste0("slow-parallel-converted-", rep)))

        tm <- .tte_model_fixture()
        cl <- parallel::makePSOCKcluster(2)
        parallel::clusterCall(cl, .load_simtte_on_worker,
            pkg_path = .simtte_pkg_path, dev_loaded = .simtte_is_dev_loaded)
        prepared <- simtte::simtte_prepare_model(tm)
        parallel::clusterCall(cl, simtte::simtte_prepare_model,
            model = prepared)
        par_events <- parallel::parLapply(cl, seeds, function(s, prepared) {
            simtte::sim_tte_ode(model = prepared, n = 10, end = 20,
                delta = 2, seed = s)$events
        }, prepared = prepared)
        parallel::stopCluster(cl)

        for (i in seq_along(seeds)) {
            expect_identical(par_events[[i]], seq_events[[i]],
                label = paste0("rep=", rep, " seed=", seeds[i]))
        }
    }
})

# ---- Genuinely separate processes ("simttepower feedback 3" Task A,
# reports/31_simttepower_feedback_process_loading.md) ----
#
# The PSOCK-cluster tests above cover workers that receive an
# already-compiled model *object* via ordinary parallel serialization
# -- a different scenario from a chunked/array-job task that reads a
# tte_model()-converted model back from disk in a process that shares
# no memory, and no still-running parent, with whatever process built
# it. callr::r() spawns a genuinely separate R process (its own
# tempdir(), no shared state) to exercise exactly that.

.simtte_process_boundary_build <- function(pkg_path, dev_loaded, lib_file) {
    if (dev_loaded) {
        pkgload::load_all(pkg_path, quiet = TRUE)
    }
    irm1 <- mrgsolve::mread(model = "pk2cmt",
        project = system.file("models", package = "mrgsolve"),
        compile = FALSE)
    tm <- simtte::tte_model(irm1, hazard = "H0 * exp(lp + beta_cp * CP)",
        params = list(H0 = 0.01, beta_cp = 0.02))
    saveRDS(tm, lib_file)
    invisible(NULL)
}

.simtte_process_boundary_use <- function(pkg_path, dev_loaded, lib_file, seed) {
    if (dev_loaded) {
        pkgload::load_all(pkg_path, quiet = TRUE)
    }
    before_project <- getOption("mrgsolve.project")
    before_soloc <- getOption("mrgsolve.soloc")
    tm <- readRDS(lib_file)
    prepared <- simtte::simtte_prepare_model(tm)
    sim <- simtte::sim_tte_ode(model = prepared, n = 10, end = 20, delta = 2,
        seed = seed)
    list(events = sim$events,
        project_untouched = identical(getOption("mrgsolve.project"),
            before_project),
        soloc_untouched = identical(getOption("mrgsolve.soloc"),
            before_soloc))
}

test_that("a tte_model()-converted model survives a genuinely separate process (not a live cluster) with no [loadso] failure, and matches a single-process run at the same seed", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_installed("callr")
    if (.simtte_is_dev_loaded) skip_if_not_installed("pkgload")
    skip_if_not_slow()

    lib_file <- tempfile(fileext = ".rds")
    callr::r(.simtte_process_boundary_build,
        args = list(pkg_path = .simtte_pkg_path,
            dev_loaded = .simtte_is_dev_loaded, lib_file = lib_file))

    result <- callr::r(.simtte_process_boundary_use,
        args = list(pkg_path = .simtte_pkg_path,
            dev_loaded = .simtte_is_dev_loaded, lib_file = lib_file, seed = 1))

    expect_true(result$project_untouched)
    expect_true(result$soloc_untouched)

    irm1 <- mrgsolve::mread(model = "pk2cmt",
        project = system.file("models", package = "mrgsolve"),
        compile = FALSE)
    tm_local <- tte_model(irm1, hazard = "H0 * exp(lp + beta_cp * CP)",
        params = list(H0 = 0.01, beta_cp = 0.02))
    sim_local <- sim_tte_ode(model = tm_local, n = 10, end = 20, delta = 2,
        seed = 1)
    expect_identical(result$events, sim_local$events)
})

test_that("a built-in library model works the same way in a genuinely separate process, matching a single-process run", {
    skip_on_cran()
    skip_if_not_installed("mrgsolve")
    skip_if_not_installed("callr")
    if (.simtte_is_dev_loaded) skip_if_not_installed("pkgload")
    skip_if_not_slow()

    use_step <- function(pkg_path, dev_loaded, seed) {
        if (dev_loaded) {
            pkgload::load_all(pkg_path, quiet = TRUE)
        }
        simtte::sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
            n = 10, end = 20, delta = 2, seed = seed)$events
    }
    result <- callr::r(use_step, args = list(pkg_path = .simtte_pkg_path,
        dev_loaded = .simtte_is_dev_loaded, seed = 1))

    local <- sim_tte_ode(model = "exponential", param = list(H0 = 0.1),
        n = 10, end = 20, delta = 2, seed = 1)$events
    expect_identical(result, local)
})
