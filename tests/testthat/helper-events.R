# Shared fixture for the censoring/interval-censoring/visit-process test
# files (Phase 6/6b/6c). Builds a minimal events-shaped data frame
# (ID, sim_time, sim_status) without needing a real sim_tte_ode() run --
# used by test-censoring.R, test-interval-censoring.R, and
# test-visit-process.R alike, all three of which operate on this exact
# column shape regardless of which engine produced it (pre-Phase-6c,
# each file defined this identically; consolidated here per the
# pre-Phase-7 review, reports/24_pre_phase7_review.md).
.fake_events <- function(sim_time, sim_status, id = seq_along(sim_time)) {
    data.frame(ID = id, sim_time = sim_time, sim_status = sim_status)
}
