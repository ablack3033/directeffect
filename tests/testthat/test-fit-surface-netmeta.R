# The hand oracle wls_surface() and the spec example tables live in
# helper-directeffect.R.

test_that("netmeta surface matches the hand computation on the spec example", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")

  expect_s3_class(fit, "directeffect_fit")
  expect_identical(
    names(fit$effects),
    c("drug", "estimate", "std_error", "lower", "upper", "scale",
      "reference", "engine")
  )
  expect_identical(fit$effects$drug, c("A", "B", "C"))
  expect_identical(unique(fit$effects$scale), "log")
  expect_identical(unique(fit$effects$reference), "A")
  expect_identical(unique(fit$effects$engine), "netmeta")

  oracle <- wls_surface(spec_comparisons(), c("A", "B", "C"), "A")
  free <- fit$effects[fit$effects$drug != "A", ]
  expect_equal(free$estimate, unname(oracle$estimate[free$drug]),
               tolerance = 1e-6)
  expect_equal(free$std_error, unname(oracle$std_error[free$drug]),
               tolerance = 1e-6)

  ref_row <- fit$effects[fit$effects$drug == "A", ]
  expect_equal(ref_row$estimate, 0)
  expect_equal(ref_row$std_error, 0)

  expect_equal(free$lower, free$estimate - qnorm(0.975) * free$std_error,
               tolerance = 1e-6)
  expect_equal(free$upper, free$estimate + qnorm(0.975) * free$std_error,
               tolerance = 1e-6)
})

test_that("the fit carries the full component contract", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")

  expect_identical(
    names(fit),
    c("effects", "comparisons", "anchors", "heterogeneity", "diagnostics",
      "engine", "engine_fit", "network")
  )
  expect_identical(fit$engine, "netmeta")
  expect_identical(fit$comparisons, de$comparisons)
  expect_null(fit$anchors)
  expect_true(all(c("Q", "df", "p_value", "tau", "I2") %in%
                    names(fit$heterogeneity)))
  expect_s3_class(fit$network, "directeffect_network")
  # engine_fit is the sole escape hatch for engine internals
  expect_s3_class(fit$engine_fit, "netmeta")
})

test_that("anchors are deliberately ignored by surface fitting", {
  skip_if_not_installed("netmeta")

  bare <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  anchored <- direct_effect_network(spec_comparisons(),
                                    anchors = spec_anchors(),
                                    effect_measure = "HR")

  fit_bare <- fit_surface(bare, engine = "netmeta")
  fit_anchored <- fit_surface(anchored, engine = "netmeta")

  expect_equal(fit_anchored$effects, fit_bare$effects)
  expect_identical(nrow(fit_anchored$anchors), 1L)
})

test_that("estimated differences are invariant to the reference choice", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_a <- fit_surface(de, engine = "netmeta", reference = "A")
  fit_c <- fit_surface(de, engine = "netmeta", reference = "C")

  diff_a <- with(fit_a$effects,
                 estimate[drug == "A"] - estimate[drug == "B"])
  diff_c <- with(fit_c$effects,
                 estimate[drug == "A"] - estimate[drug == "B"])
  expect_equal(diff_a, diff_c, tolerance = 1e-8)
  expect_identical(unique(fit_c$effects$reference), "C")
  expect_equal(fit_c$effects$estimate[fit_c$effects$drug == "C"], 0)
})

test_that("the heterogeneity block matches the hand-computed Q", {
  skip_if_not_installed("netmeta")

  # Exactly consistent comparisons: the surface reproduces every edge,
  # so the weighted residual sum of squares Q is 0 on df = 1
  # (3 comparisons - 2 free effects).
  consistent <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.0, 0.4, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
  de <- direct_effect_network(consistent, effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")
  expect_equal(fit$heterogeneity$Q, 0, tolerance = 1e-10)
  expect_identical(as.integer(fit$heterogeneity$df), 1L)
  expect_equal(fit$heterogeneity$tau, 0, tolerance = 1e-6)

  # An inconsistent edge produces exactly the Q of the WLS residuals.
  inconsistent <- consistent
  inconsistent$estimate[3] <- 0.9
  de2 <- direct_effect_network(inconsistent, effect_measure = "HR")
  fit2 <- fit_surface(de2, engine = "netmeta")

  oracle <- wls_surface(inconsistent, c("A", "B", "C"), "A")
  theta <- c(A = 0, oracle$estimate)
  predicted <- theta[inconsistent$target] - theta[inconsistent$comparator]
  q_hand <- sum(((inconsistent$estimate - predicted) /
                   inconsistent$std_error)^2)
  expect_equal(fit2$heterogeneity$Q, q_hand, tolerance = 1e-8)
  expect_gt(fit2$heterogeneity$I2, 0)
})

test_that("fitting a multi-component network fails informatively", {
  skip_if_not_installed("netmeta")

  comparisons <- data.frame(
    study_id   = c("S1", "S2"),
    target     = c("A", "X"),
    comparator = c("B", "Y"),
    estimate   = c(0.1, 0.2),
    std_error  = c(0.05, 0.05)
  )
  de <- direct_effect_network(comparisons, effect_measure = "HR")
  expect_error(fit_surface(de, engine = "netmeta"),
               "check_connectivity")
  expect_error(fit_surface(de, engine = "netmeta"), "2 connected components")
})

test_that("fit_surface validates its inputs", {
  expect_error(fit_surface(data.frame()), "directeffect_network")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  expect_error(fit_surface(de, engine = "netmeta", reference = "Z"),
               "must be one of the network's treatments")
})

test_that("the fit print method shows engine and reference", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")
  expect_output(print(fit), "directeffect_fit")
  expect_output(print(fit), "Engine:\\s+netmeta")
  expect_output(print(fit), "Reference:\\s+A")
})
