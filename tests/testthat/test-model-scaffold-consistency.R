# Pre-Phase-7 review (reports/24_pre_phase7_review.md): the survival
# scaffold (GLOBAL statics, MAIN reset, ODE detection block) is
# supposed to be identical, by design, across every shipped
# inst/models/library/*.cpp model -- documented by hand in
# exponential_ode.cpp's own comments, but never previously verified by
# a test. This scripts that verification so any future edit to one
# model's scaffold that isn't mirrored everywhere else fails a test
# instead of silently drifting.

.extract_scaffold <- function(lines) {
    global <- trimws(grep(
        "^\\s*static (int|double) (event_found|TEVT|T_PRE|P_PRE|P_POST)\\s*=",
        lines, value = TRUE))
    main_start <- grep("if \\(NEWIND <= 1\\) \\{", lines)
    main <- if (length(main_start)) {
        rel_end <- which(grepl("^\\s*\\}\\s*$",
            lines[(main_start[1] + 1L):length(lines)]))[1]
        trimws(lines[main_start[1]:(main_start[1] + rel_end)])
    } else {
        character(0)
    }
    dxdt <- trimws(grep("dxdt_p11 = -p11 \\* HAZ;", lines, value = TRUE))
    pre_start <- grep("if \\(!event_found && p11 > U", lines)
    post_start <- grep("if \\(!event_found && p11 <= U", lines)
    ode_extra <- if (length(pre_start) && length(post_start)) {
        trimws(lines[pre_start[1]:(post_start[1] + 3L)])
    } else {
        character(0)
    }
    list(global = global, main = main, dxdt = dxdt, ode_extra = ode_extra)
}

test_that("the survival scaffold is byte-identical (whitespace-trimmed) across every shipped ODE library model", {
    dir <- system.file("models", "library", package = "simtte")
    files <- list.files(dir, pattern = "\\.cpp$", full.names = TRUE)
    expect_true(length(files) >= 10) # sanity: the glob actually found files

    scaffolds <- lapply(files, function(f) .extract_scaffold(readLines(f)))
    names(scaffolds) <- basename(files)
    ref <- scaffolds[[1]]

    for (nm in names(scaffolds)) {
        info <- paste0(nm, " vs. ", names(scaffolds)[1])
        expect_identical(scaffolds[[nm]]$global, ref$global, info = info)
        expect_identical(scaffolds[[nm]]$main, ref$main, info = info)
        expect_identical(scaffolds[[nm]]$dxdt, ref$dxdt, info = info)
        expect_identical(scaffolds[[nm]]$ode_extra, ref$ode_extra, info = info)
        # Every one of these blocks is non-empty for every file -- a
        # silent regex-miss (e.g. a renamed variable) would otherwise
        # pass the identical() checks above vacuously (all-empty).
        expect_true(length(scaffolds[[nm]]$global) > 0, info = info)
        expect_true(length(scaffolds[[nm]]$main) > 0, info = info)
        expect_true(length(scaffolds[[nm]]$dxdt) > 0, info = info)
        expect_true(length(scaffolds[[nm]]$ode_extra) > 0, info = info)
    }
})

test_that("every shipped ODE library model documents the scaffold with BEGIN/END markers", {
    dir <- system.file("models", "library", package = "simtte")
    files <- list.files(dir, pattern = "\\.cpp$", full.names = TRUE)
    for (f in files) {
        lines <- readLines(f)
        expect_true(any(grepl("BEGIN simtte survival scaffolding", lines)),
            info = basename(f))
        expect_true(any(grepl("END simtte survival scaffolding", lines)),
            info = basename(f))
    }
})
