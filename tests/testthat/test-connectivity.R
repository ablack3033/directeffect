# Two disconnected components: {A, B, C} with an anchor on C, and {X, Y}
# with no anchor.
two_component_network <- function() {
  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3", "S4"),
    target     = c("A", "A", "B", "X"),
    comparator = c("B", "C", "C", "Y"),
    estimate   = c(0.02, 0.30, 0.25, 0.10),
    std_error  = c(0.07, 0.09, 0.08, 0.05)
  )
  anchors <- data.frame(
    study_id  = "RCT1",
    drug      = "C",
    reference = "placebo",
    estimate  = 0.18,
    std_error = 0.04
  )
  direct_effect_network(comparisons, anchors = anchors, effect_measure = "HR")
}

test_that("check_connectivity reports each component correctly", {
  report <- check_connectivity(two_component_network())

  expect_s3_class(report, "directeffect_connectivity")
  expect_identical(nrow(report), 2L)

  drugs <- attr(report, "drugs")
  abc <- which(vapply(drugs, function(d) "A" %in% d, logical(1)))
  xy <- which(vapply(drugs, function(d) "X" %in% d, logical(1)))

  expect_identical(drugs[[abc]], c("A", "B", "C"))
  expect_identical(drugs[[xy]], c("X", "Y"))
  expect_identical(report$n_drugs[abc], 3L)
  expect_identical(report$n_comparisons[abc], 3L)
  expect_identical(report$n_anchors[abc], 1L)
  expect_true(report$absolute_identifiable[abc])
  expect_identical(report$n_drugs[xy], 2L)
  expect_identical(report$n_comparisons[xy], 1L)
  expect_identical(report$n_anchors[xy], 0L)
  expect_false(report$absolute_identifiable[xy])
})

test_that("a fully connected anchored network is absolutely identifiable", {
  comparisons <- data.frame(
    study_id   = c("S1", "S2"),
    target     = c("A", "B"),
    comparator = c("B", "C"),
    estimate   = c(0.1, 0.2),
    std_error  = c(0.05, 0.05)
  )
  anchors <- data.frame(
    study_id  = "RCT1",
    drug      = "A",
    reference = "placebo",
    estimate  = 0.5,
    std_error = 0.1
  )
  de <- direct_effect_network(comparisons, anchors = anchors,
                              effect_measure = "HR")
  report <- check_connectivity(de)
  expect_identical(nrow(report), 1L)
  expect_true(report$absolute_identifiable)
})

test_that("an unanchored component prints as relative-only", {
  report <- check_connectivity(two_component_network())
  expect_output(print(report), "relative effects identifiable only")
  expect_output(print(report), "absolute effects identifiable")
  expect_output(print(report), "Component 1:")
  expect_output(print(report), "Component 2:")
  expect_output(print(report), "3 drugs")
  expect_output(print(report), "1 absolute anchor\\b")
})

test_that("a network with no anchors at all reports zero anchors", {
  comparisons <- data.frame(
    study_id   = c("S1", "S2"),
    target     = c("A", "B"),
    comparator = c("B", "C"),
    estimate   = c(0.1, 0.2),
    std_error  = c(0.05, 0.05)
  )
  report <- check_connectivity(
    direct_effect_network(comparisons, effect_measure = "HR")
  )
  expect_identical(report$n_anchors, 0L)
  expect_false(report$absolute_identifiable)
})

test_that("check_connectivity rejects non-network input", {
  expect_error(check_connectivity(data.frame()), "directeffect_network")
})
