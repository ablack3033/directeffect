skip_if_not_installed("ggplot2")
skip_if_not_installed("netmeta")

surface_fit <- function() {
  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_surface(de, engine = "netmeta")
}

anchored_fit <- function() {
  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")
  anchor_surface(fit_surface(de, engine = "netmeta"))
}

expect_renders <- function(plot) {
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_true(length(built$data) > 0)
  invisible(built)
}

test_that("surface plot renders on both scales with labeled axes", {
  fit <- surface_fit()

  log_plot <- plot_effect_surface(fit, scale = "log")
  expect_renders(log_plot)
  expect_match(log_plot$labels$x, "log HR")

  natural_plot <- plot_effect_surface(fit, scale = "natural")
  expect_renders(natural_plot)
  expect_match(natural_plot$labels$x, "HR")
})

test_that("the natural scale back-transforms so log 0 becomes HR = 1", {
  fit <- surface_fit()
  natural_plot <- plot_effect_surface(fit, scale = "natural")

  expect_equal(natural_plot$data$estimate, exp(fit$effects$estimate))
  expect_equal(natural_plot$data$lower, exp(fit$effects$lower))
  # The reference drug sits at log 0, i.e. HR 1.
  ref <- natural_plot$data[natural_plot$data$drug == "A", ]
  expect_equal(ref$estimate, 1)
})

test_that("unanchored fits are clearly labeled as relative", {
  plot <- plot_effect_surface(surface_fit())
  expect_match(plot$labels$subtitle, "arbitrary reference")
  expect_match(plot$labels$subtitle, "not absolute")
})

test_that("anchored fits show the placebo sea-level line", {
  plot <- plot_effect_surface(anchored_fit())
  expect_match(plot$labels$subtitle, "placebo")

  built <- expect_renders(plot)
  vline <- built$data[[1]]
  expect_equal(vline$xintercept, 0)
  expect_identical(unique(vline$linetype), "solid")
})

test_that("the sea-level figure renders for anchored fits only", {
  plot <- plot_sea_level(anchored_fit())
  expect_renders(plot)
  expect_match(plot$labels$y, "log HR")

  natural <- plot_sea_level(anchored_fit(), scale = "natural")
  expect_renders(natural)
  expect_equal(natural$data$estimate,
               exp(anchored_fit()$effects$estimate))

  expect_error(plot_sea_level(surface_fit()), "no sea level")
})

test_that("surface plots are engine-agnostic", {
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), anchors = spec_anchors(),
                              effect_measure = "HR")
  stan_surface <- fit_surface(de, engine = "stan", seed = 1)
  expect_renders(plot_effect_surface(stan_surface))

  stan_anchored <- anchor_surface(stan_surface, seed = 1)
  expect_renders(plot_effect_surface(stan_anchored))
  expect_renders(plot_sea_level(stan_anchored))
})

test_that("surface plots validate their input", {
  expect_error(plot_effect_surface(data.frame()), "directeffect_fit")
  expect_error(plot_sea_level(data.frame()), "directeffect_fit")
})
