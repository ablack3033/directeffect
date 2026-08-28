# The spec's 3-drug example tables live in helper-directeffect.R.

test_that("constructor builds a validated network from the spec example", {
  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")

  expect_s3_class(de, "directeffect_network")
  expect_identical(de$treatments, c("A", "B", "C"))
  expect_identical(nrow(de$comparisons), 3L)
  expect_identical(nrow(de$anchors), 1L)
  expect_identical(de$effect_measure, "HR")
  expect_identical(max(de$components), 1L)
  expect_named(de$components)
  expect_true(igraph::is_igraph(de$graph))
  expect_identical(sort(igraph::V(de$graph)$name), de$treatments)
  expect_null(de$fit)
})

test_that("anchors are optional and extra metadata columns are preserved", {
  comparisons <- spec_comparisons()
  comparisons$database <- "claims_db"
  comparisons$n_target <- c(1000L, 1200L, 900L)

  de <- direct_effect_network(comparisons, effect_measure = "OR")
  expect_null(de$anchors)
  expect_true(all(c("database", "n_target") %in% names(de$comparisons)))
})

test_that("construction fails on missing required columns", {
  comparisons <- spec_comparisons()
  comparisons$std_error <- NULL
  expect_error(direct_effect_network(comparisons), "std_error")

  anchors <- spec_anchors()
  anchors$reference <- NULL
  expect_error(direct_effect_network(spec_comparisons(), anchors = anchors),
               "reference")
})

test_that("construction fails on non-positive standard errors", {
  comparisons <- spec_comparisons()
  comparisons$std_error[2] <- 0
  expect_error(direct_effect_network(comparisons), "positive")

  anchors <- spec_anchors()
  anchors$std_error <- -0.1
  expect_error(direct_effect_network(spec_comparisons(), anchors = anchors),
               "positive")
})

test_that("construction fails on self-comparisons", {
  comparisons <- spec_comparisons()
  comparisons$comparator[1] <- "A"
  expect_error(direct_effect_network(comparisons), "self-comparison")
})

test_that("construction fails on missing drug names", {
  comparisons <- spec_comparisons()
  comparisons$target[3] <- NA_character_
  expect_error(direct_effect_network(comparisons), "missing drug names")

  anchors <- spec_anchors()
  anchors$drug <- NA_character_
  expect_error(direct_effect_network(spec_comparisons(), anchors = anchors),
               "missing drug names")
})

test_that("construction fails on anchors for drugs absent from comparisons", {
  anchors <- spec_anchors()
  anchors$drug <- "Z"
  expect_error(direct_effect_network(spec_comparisons(), anchors = anchors),
               "appear in no comparison")
})

test_that("construction fails on zero-row or non-data-frame inputs", {
  expect_error(direct_effect_network(spec_comparisons()[0, ]),
               "at least one comparison")
  expect_error(direct_effect_network(list(study_id = "S1")),
               "must be a data frame")
  expect_error(
    direct_effect_network(spec_comparisons(), anchors = "not a table"),
    "must be a data frame"
  )
})

test_that("zero-row anchors are treated as no anchors", {
  de <- direct_effect_network(spec_comparisons(),
                              anchors = spec_anchors()[0, ],
                              effect_measure = "HR")
  expect_null(de$anchors)
})

test_that("construction fails when an anchor's drug equals its reference", {
  anchors <- spec_anchors()
  anchors$reference <- "C"
  expect_error(direct_effect_network(spec_comparisons(), anchors = anchors),
               "`drug` equals `reference`")
})

test_that("construction fails on non-numeric or missing estimates", {
  comparisons <- spec_comparisons()
  comparisons$estimate <- as.character(comparisons$estimate)
  expect_error(direct_effect_network(comparisons),
               "estimate.*numeric")

  comparisons <- spec_comparisons()
  comparisons$estimate[2] <- NA_real_
  expect_error(direct_effect_network(comparisons),
               "no missing values")
})

test_that("construction fails on an invalid effect measure", {
  expect_error(direct_effect_network(spec_comparisons(), effect_measure = NA),
               "effect_measure")
  expect_error(
    direct_effect_network(spec_comparisons(), effect_measure = c("HR", "OR")),
    "effect_measure"
  )
})

test_that("differing estimand metadata triggers a warning", {
  comparisons <- spec_comparisons()
  comparisons$population <- c("adults", "adults", "children")
  expect_warning(direct_effect_network(comparisons), "population")

  consistent <- spec_comparisons()
  consistent$population <- "adults"
  expect_no_warning(direct_effect_network(consistent))
})

test_that("print method summarises the network", {
  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")
  expect_output(print(de), "directeffect_network")
  expect_output(print(de), "Effect measure: HR")
  expect_output(print(de), "Drugs:\\s+3")
  expect_output(print(de), "Comparisons:\\s+3")
  expect_output(print(de), "Anchors:\\s+1")
  expect_output(print(de), "Components:\\s+1")
})
