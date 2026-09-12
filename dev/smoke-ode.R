#!/usr/bin/env Rscript
# Interactive counterpart to inst/validation/08_.../09_...: eyeball one
# sim_tte_ode() model quickly, from the console or the shell. Not a
# replacement for the validation scripts (which are the citable,
# comprehensive evidence for the Phase 1/2 reports) -- this is for a
# quick "does this still look right" check while iterating.
#
# Console:
#   model <- "weibull"; shape <- 0.5
#   source("dev/smoke-ode.R")
#
# Shell:
#   Rscript dev/smoke-ode.R model=weibull shape=0.5 n=1000 delta=0.5

# ---- Parameters: top-of-file variables, overridable either by setting
# them before source()-ing (console) or as key=value shell arguments.
# `inherits = FALSE` matters here: `end` (and, in principle, others)
# would otherwise match a base-package function of the same name via
# the search path, so `exists("end")` alone is silently always TRUE. ----
.smoke_default <- function(name, value) {
    if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        assign(name, value, envir = .GlobalEnv)
    }
}
.smoke_default("model", "exponential")
.smoke_default("n", 500)
.smoke_default("delta", 1)
.smoke_default("end", 20)
.smoke_default("seed", 1)
.smoke_default("H0", 0.1)
.smoke_default("mu", -1)
.smoke_default("shape", 1.5)
.smoke_default("gamma", 0.1)
.smoke_default("lp", 0)

if (!interactive()) {
    for (a in commandArgs(trailingOnly = TRUE)) {
        kv <- strsplit(a, "=", fixed = TRUE)[[1]]
        if (length(kv) == 2) {
            val <- suppressWarnings(as.numeric(kv[2]))
            assign(kv[1], if (is.na(val)) kv[2] else val, envir = .GlobalEnv)
        }
    }
}

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root (Rscript dev/smoke-ode.R ...), ",
        "not from inside dev/.", call. = FALSE)
}
suppressMessages(devtools::load_all(quiet = TRUE))

param <- switch(model,
    exponential = list(H0 = H0, lp = lp),
    weibull = list(mu = mu, shape = shape, lp = lp),
    gompertz = list(mu = mu, gamma = gamma, lp = lp),
    stop("Unknown model '", model, "'. Use exponential, weibull, or gompertz."))

analytic_S <- switch(model,
    exponential = { eta <- H0 * exp(lp); function(t) exp(-eta * t) },
    weibull = { eta <- exp(mu + lp); function(t) exp(-eta * t^shape) },
    gompertz = {
        eta <- exp(mu + lp)
        function(t) exp(-(eta / gamma) * (exp(gamma * t) - 1))
    })
analytic_T <- switch(model,
    exponential = { eta <- H0 * exp(lp); function(u) -log(u) / eta },
    weibull = {
        eta <- exp(mu + lp)
        function(u) (-log(u) / eta)^(1 / shape)
    },
    gompertz = {
        eta <- exp(mu + lp)
        function(u) suppressWarnings(log(1 - (log(u) * gamma) / eta) / gamma)
    })

cat("model:", model, " param:", paste(names(param), unlist(param), sep = "=",
    collapse = ", "), "\n")
cat("n:", n, " delta:", delta, " end:", end, " seed:", seed, "\n\n")

sim <- sim_tte_ode(model = model, param = param, n = n, end = end,
    delta = delta, keep_trajectory = TRUE, seed = seed)

cat("---- head(events) ----\n")
print(head(sim$events))

cat("\n---- censoring fraction ----\n")
cat(round(mean(sim$events$sim_status == 0), 4), "\n")

cat("\n---- empirical vs analytical P(event by t) at quantile checkpoints ----\n")
# t_j chosen so P(event by t_j) = prob exactly, by construction of
# analytic_T() (the same approach test-sim-tte-ode-weibull.R uses to
# avoid saturated, uninformative checkpoints -- see its own comment).
probs <- c(0.2, 0.4, 0.6, 0.8)
t_j <- vapply(probs, function(p) analytic_T(1 - p), numeric(1))
p_event_empirical <- vapply(t_j, function(t) {
    mean(sim$events$sim_time <= t & sim$events$sim_status == 1)
}, numeric(1))
print(data.frame(prob_target = probs, t = round(t_j, 3),
    empirical_P_event = round(p_event_empirical, 4)))

cat("\n---- raw TEVT vs refined sim_time, mean abs error vs analytical quantile ----\n")
traj <- sim$trajectory
last <- traj[!duplicated(traj$ID, fromLast = TRUE), ]
u_by_id <- last$U[match(sim$events$ID, last$ID)]
t_true <- analytic_T(u_by_id)
is_event <- sim$events$sim_status == 1L & !is.na(t_true) & t_true < end
raw_tevt <- last$TEVT[match(sim$events$ID, last$ID)]
err_raw <- abs(raw_tevt[is_event] - t_true[is_event])
err_refined <- abs(sim$events$sim_time[is_event] - t_true[is_event])
cat("n compared:", sum(is_event), "\n")
cat("raw      mean abs err:", round(mean(err_raw), 5), "\n")
cat("refined  mean abs err:", round(mean(err_refined), 5), "\n")
