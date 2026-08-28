# Mirrors the generation in data-raw/example-network.R; the drift guard
# below fails if the shipped .rda files stop matching that script (e.g.
# after editing the design there without regenerating), forcing a
# deliberate regeneration instead of silent drift.
regenerate_example_data <- function() {
  truth <- c(
    atorvastatin = -0.30,
    fluvastatin  = -0.17,
    lovastatin   = -0.20,
    pravastatin  = -0.22,
    rosuvastatin = -0.32,
    simvastatin  = -0.26
  )
  design <- data.frame(
    study_id   = c("rct_ator_simv_1", "rct_ator_simv_2", "rct_ator_prav_1",
                   "rct_ator_prav_2", "rct_rosu_ator_1", "rct_rosu_ator_2",
                   "rct_rosu_simv_1", "rct_simv_prav_1", "rct_prav_lova_1",
                   "rct_lova_fluv_1", "rct_simv_fluv_1", "rct_prav_fluv_1"),
    target     = c("atorvastatin", "atorvastatin", "atorvastatin",
                   "atorvastatin", "rosuvastatin", "rosuvastatin",
                   "rosuvastatin", "simvastatin", "pravastatin",
                   "lovastatin", "simvastatin", "pravastatin"),
    comparator = c("simvastatin", "simvastatin", "pravastatin",
                   "pravastatin", "atorvastatin", "atorvastatin",
                   "simvastatin", "pravastatin", "lovastatin",
                   "fluvastatin", "fluvastatin", "fluvastatin"),
    std_error  = c(0.06, 0.09, 0.05, 0.10, 0.07, 0.11,
                   0.08, 0.07, 0.08, 0.10, 0.09, 0.12)
  )
  anchor_design <- data.frame(
    study_id  = c("rct_simv_plac", "rct_prav_plac", "rct_rosu_plac"),
    drug      = c("simvastatin", "pravastatin", "rosuvastatin"),
    reference = "placebo",
    std_error = c(0.06, 0.07, 0.08)
  )
  tau <- 0.02
  set.seed(12)

  delta <- truth[design$target] - truth[design$comparator]
  comparisons <- data.frame(
    study_id   = design$study_id,
    target     = design$target,
    comparator = design$comparator,
    estimate   = round(unname(stats::rnorm(
      nrow(design), delta, sqrt(design$std_error^2 + tau^2)
    )), 3),
    std_error  = design$std_error
  )
  anchors <- data.frame(
    study_id  = anchor_design$study_id,
    drug      = anchor_design$drug,
    reference = anchor_design$reference,
    estimate  = round(unname(stats::rnorm(
      nrow(anchor_design), truth[anchor_design$drug],
      anchor_design$std_error
    )), 3),
    std_error = anchor_design$std_error
  )
  truth_table <- data.frame(
    drug  = c("placebo", names(truth)),
    theta = c(0, unname(truth))
  )
  list(comparisons = comparisons, anchors = anchors, truth = truth_table)
}

test_that("the packaged example network is valid and single-component", {
  de <- direct_effect_network(example_comparisons,
                              anchors = example_anchors,
                              effect_measure = "HR")

  expect_identical(length(de$treatments), 6L)
  expect_identical(nrow(de$comparisons), 12L)
  expect_identical(nrow(de$anchors), 3L)
  expect_identical(max(de$components), 1L)
  expect_true(check_connectivity(de)$absolute_identifiable)
})

test_that("the shipped truth is well-formed and consistent with the data", {
  expect_identical(names(example_truth), c("drug", "theta"))
  expect_identical(nrow(example_truth), 7L)
  expect_identical(example_truth$drug[1], "placebo")
  expect_identical(example_truth$theta[1], 0)

  de <- direct_effect_network(example_comparisons,
                              anchors = example_anchors,
                              effect_measure = "HR")
  expect_setequal(example_truth$drug, c("placebo", de$treatments))
})

test_that("the shipped data matches its generator (drift guard)", {
  regenerated <- regenerate_example_data()
  expect_identical(example_comparisons, regenerated$comparisons)
  expect_identical(example_anchors, regenerated$anchors)
  expect_identical(example_truth, regenerated$truth)
})

test_that("the example workflow recovers the shipped truth", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(example_comparisons,
                              anchors = example_anchors,
                              effect_measure = "HR")
  absolute <- anchor_surface(fit_surface(de, engine = "netmeta"))
  expect_identical(unique(absolute$effects$reference), "placebo")
  expect_identical(nrow(absolute$effects), 6L)

  # Deterministic given the shipped data, so the bounds can sit close
  # to the observed values (RMSE 0.062, coverage 6/6, rank cor 1.0)
  # while still failing on any real regression in fitting or anchoring.
  truth <- example_truth$theta[match(absolute$effects$drug,
                                     example_truth$drug)]
  errors <- absolute$effects$estimate - truth
  covered <- absolute$effects$lower <= truth &
    truth <= absolute$effects$upper

  expect_lt(max(abs(errors)), 0.15)
  expect_lt(sqrt(mean(errors^2)), 0.08)
  expect_gte(mean(covered), 0.8)
  expect_gt(cor(absolute$effects$estimate, truth, method = "spearman"),
            0.9)
})
