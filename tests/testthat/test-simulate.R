test_that("simulation is reproducible under a seed", {
  a <- simulate_direct_effect_network(n_drugs = 6, n_comparisons = 12,
                                      n_anchors = 2, heterogeneity = 0,
                                      seed = 42)
  b <- simulate_direct_effect_network(n_drugs = 6, n_comparisons = 12,
                                      n_anchors = 2, heterogeneity = 0,
                                      seed = 42)
  c <- simulate_direct_effect_network(n_drugs = 6, n_comparisons = 12,
                                      n_anchors = 2, heterogeneity = 0,
                                      seed = 43)

  expect_identical(a$comparisons, b$comparisons)
  expect_identical(a$anchors, b$anchors)
  expect_identical(a$truth, b$truth)
  expect_false(identical(a$comparisons$estimate, c$comparisons$estimate))
})

test_that("simulation does not disturb the global random-number state", {
  set.seed(99)
  before <- .Random.seed
  simulate_direct_effect_network(n_drugs = 4, n_comparisons = 6, seed = 7)
  expect_identical(.Random.seed, before)
})

test_that("generated comparison graphs are connected by construction", {
  for (seed in 1:5) {
    simulation <- simulate_direct_effect_network(
      n_drugs = 12, n_comparisons = 15, n_anchors = 0,
      heterogeneity = 0, seed = seed
    )
    expect_identical(max(simulation$network$components), 1L)
  }
})

test_that("simulation output is network-ready and retains truth", {
  simulation <- simulate_direct_effect_network(
    n_drugs = 8, n_comparisons = 20, n_anchors = 3,
    heterogeneity = 0.05, seed = 11
  )

  expect_s3_class(simulation$network, "directeffect_network")
  expect_identical(nrow(simulation$comparisons), 20L)
  expect_identical(nrow(simulation$anchors), 3L)
  expect_identical(unique(simulation$anchors$reference), "placebo")
  expect_identical(anyDuplicated(simulation$anchors$drug), 0L)

  expect_identical(simulation$truth$drug[1], "placebo")
  expect_identical(simulation$truth$theta[1], 0)
  expect_setequal(simulation$truth$drug,
                  c("placebo", simulation$network$treatments))
  expect_identical(simulation$tau, 0.05)
})

test_that("impossible simulation requests fail loudly", {
  expect_error(
    simulate_direct_effect_network(n_drugs = 10, n_comparisons = 5),
    "connected"
  )
  expect_error(
    simulate_direct_effect_network(n_drugs = 4, n_comparisons = 8,
                                   n_anchors = 5),
    "n_anchors"
  )
  expect_error(
    simulate_direct_effect_network(n_drugs = 1, n_comparisons = 5),
    "n_drugs"
  )
  expect_error(
    simulate_direct_effect_network(heterogeneity = -0.1),
    "heterogeneity"
  )
})

test_that("a low-noise simulation is recovered with near-nominal accuracy", {
  skip_if_not_installed("netmeta")

  simulation <- simulate_direct_effect_network(
    n_drugs = 40, n_comparisons = 250, n_anchors = 0,
    heterogeneity = 0, seed = 2026,
    se_range = c(0.02, 0.06)
  )
  fit <- fit_surface(simulation$network, engine = "netmeta")
  recovery <- validate_recovery(fit, simulation)

  expect_lt(abs(recovery$bias), 0.02)
  expect_lt(recovery$rmse, 0.05)
  expect_gte(recovery$coverage, 0.85)
  expect_gt(recovery$rank_correlation, 0.95)
  expect_identical(recovery$n_drugs, 39L)
})

test_that("recovery validation works on anchored fits too", {
  skip_if_not_installed("netmeta")

  simulation <- simulate_direct_effect_network(
    n_drugs = 12, n_comparisons = 40, n_anchors = 3,
    heterogeneity = 0, seed = 7, se_range = c(0.02, 0.06)
  )
  absolute <- anchor_surface(
    fit_surface(simulation$network, engine = "netmeta")
  )
  recovery <- validate_recovery(absolute, simulation)

  # Anchored fits compare directly against the absolute truth: no drug
  # is pinned, so all twelve contribute, and nothing degenerates to NaN.
  expect_identical(recovery$n_drugs, 12L)
  expect_identical(recovery$reference, "placebo")
  expect_false(is.nan(recovery$bias))
  expect_lt(abs(recovery$bias), 0.05)
  expect_lt(recovery$rmse, 0.1)
  expect_gt(recovery$rank_correlation, 0.9)
})

test_that("validate_recovery validates its inputs", {
  simulation <- simulate_direct_effect_network(
    n_drugs = 4, n_comparisons = 8, n_anchors = 0, seed = 3
  )
  expect_error(validate_recovery(list(), simulation), "directeffect_fit")

  skip_if_not_installed("netmeta")
  fit <- fit_surface(simulation$network, engine = "netmeta")
  expect_error(validate_recovery(fit, list()), "truth")

  other <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  other_fit <- fit_surface(other, engine = "netmeta")
  expect_error(validate_recovery(other_fit, simulation),
               "simulation truth does not")
})
