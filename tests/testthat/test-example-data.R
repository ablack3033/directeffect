test_that("the packaged example network is valid and single-component", {
  de <- direct_effect_network(example_comparisons,
                              anchors = example_anchors,
                              effect_measure = "HR")

  expect_identical(length(de$treatments), 8L)
  expect_identical(nrow(de$comparisons), 16L)
  expect_identical(nrow(de$anchors), 2L)
  expect_identical(max(de$components), 1L)
  expect_true(check_connectivity(de)$absolute_identifiable)
})

test_that("the example network runs the core workflow", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(example_comparisons,
                              anchors = example_anchors,
                              effect_measure = "HR")
  absolute <- anchor_surface(fit_surface(de, engine = "netmeta"))
  expect_identical(unique(absolute$effects$reference), "placebo")
  expect_identical(nrow(absolute$effects), 8L)
})
