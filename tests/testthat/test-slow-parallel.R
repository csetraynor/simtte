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
