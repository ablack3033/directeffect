skip_if_not_installed("ggraph")
skip_if_not_installed("ggplot2")

# ggplot_build() forces layout and aesthetic evaluation, so a plot that
# would fail at render time fails here.
expect_renders <- function(plot) {
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_true(length(built$data) > 0)
  invisible(built)
}

test_that("plot_network works on a bare network, no fit required", {
  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")
  expect_renders(plot_network(de))
})

test_that("optional encodings can be toggled", {
  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")

  expect_renders(plot_network(de, weight_edges = FALSE, size_nodes = FALSE,
                              shape_anchors = FALSE))
  expect_renders(plot_network(de, label_edges = TRUE))
  expect_renders(plot_network(
    de,
    drug_classes = data.frame(drug = c("A", "B"), class = "statin")
  ))
})

test_that("plot_network renders multi-component networks", {
  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "B", "X"),
    comparator = c("B", "C", "Y"),
    estimate   = c(0.1, 0.2, -0.1),
    std_error  = c(0.05, 0.05, 0.05)
  )
  de <- direct_effect_network(comparisons, effect_measure = "HR")
  expect_identical(max(de$components), 2L)
  expect_renders(plot_network(de))
})

test_that("plot_network renders networks with repeat comparisons", {
  comparisons <- rbind(
    spec_comparisons(),
    data.frame(study_id = "S4", target = "A", comparator = "B",
               estimate = 0.05, std_error = 0.08)
  )
  de <- direct_effect_network(comparisons, effect_measure = "HR")
  expect_renders(plot_network(de))
})

test_that("plot_network validates its input", {
  expect_error(plot_network(data.frame()), "directeffect_network")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  expect_error(plot_network(de, drug_classes = data.frame(x = 1)),
               "drug.*class")
})
