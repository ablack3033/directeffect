# The design's mandatory deterministic example. True effects:
# placebo = 0, A = 0.7, B = 0.7, C = 0.3. The A-B comparison is null,
# yet neither drug has a null direct effect — the central property the
# package exists to protect.
mandatory_network <- function(anchor_se = 0.04) {
  comparisons <- data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.0, 0.4, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
  anchors <- data.frame(
    study_id  = "RCT1",
    drug      = "C",
    reference = "placebo",
    estimate  = 0.3,
    std_error = anchor_se
  )
  direct_effect_network(comparisons, anchors = anchors,
                        effect_measure = "HR")
}

test_that("mandatory example: netmeta recovers A = 0.7, B = 0.7, C = 0.3", {
  skip_if_not_installed("netmeta")

  de <- mandatory_network()
  absolute <- anchor_surface(fit_surface(de, engine = "netmeta"))
  recovered <- stats::setNames(absolute$effects$estimate,
                               absolute$effects$drug)

  expect_equal(recovered[["A"]], 0.7, tolerance = 1e-6)
  expect_equal(recovered[["B"]], 0.7, tolerance = 1e-6)
  expect_equal(recovered[["C"]], 0.3, tolerance = 1e-6)

  # The null A-B comparison must not produce null direct effects.
  expect_gt(recovered[["A"]], 0.5)
  expect_gt(recovered[["B"]], 0.5)
})

test_that("mandatory example: stan recovers A = 0.7, B = 0.7, C = 0.3", {
  skip_if_not_installed("rstan")

  de <- mandatory_network()
  surface <- fit_surface(de, engine = "stan", seed = 20260828)
  absolute <- anchor_surface(surface, seed = 20260828)
  recovered <- stats::setNames(absolute$effects$estimate,
                               absolute$effects$drug)

  expect_equal(recovered[["A"]], 0.7, tolerance = 0.03)
  expect_equal(recovered[["B"]], 0.7, tolerance = 0.03)
  expect_equal(recovered[["C"]], 0.3, tolerance = 0.03)

  expect_gt(recovered[["A"]], 0.5)
  expect_gt(recovered[["B"]], 0.5)
})

test_that("anchored fits keep the contract with reference = placebo", {
  skip_if_not_installed("netmeta")

  absolute <- anchor_surface(fit_surface(mandatory_network(),
                                         engine = "netmeta"))

  expect_s3_class(absolute, "directeffect_fit")
  expect_identical(
    names(absolute),
    c("effects", "comparisons", "anchors", "heterogeneity", "diagnostics",
      "engine", "engine_fit", "network")
  )
  expect_identical(
    names(absolute$effects),
    c("drug", "estimate", "std_error", "lower", "upper", "scale",
      "reference", "engine")
  )
  expect_identical(unique(absolute$effects$reference), "placebo")
  expect_identical(unique(absolute$effects$engine), "netmeta")
  expect_identical(nrow(absolute$anchors), 1L)
  expect_output(print(absolute), "sea level")
})

test_that("anchored stan fits keep the Bayesian columns and sample every drug", {
  skip_if_not_installed("rstan")

  surface <- fit_surface(mandatory_network(), engine = "stan",
                         seed = 20260828)
  absolute <- anchor_surface(surface, seed = 20260828)

  expect_identical(
    names(absolute$effects),
    c("drug", "estimate", "std_error", "lower", "upper", "scale",
      "reference", "engine",
      "median", "mean", "sd", "q025", "q975",
      "rhat", "ess_bulk", "ess_tail")
  )
  expect_identical(unique(absolute$effects$reference), "placebo")
  # No pinned row: with anchors there is no arbitrary constraint, so
  # every drug is sampled and has real convergence diagnostics.
  expect_true(all(absolute$effects$std_error > 0))
  expect_true(all(!is.na(absolute$effects$rhat)))
})

test_that("a surface without anchors is never silently positioned", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")
  expect_error(anchor_surface(fit), "never\\s+silently picks a sea level")
})

test_that("anchors passed explicitly override the network's anchors", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "netmeta")

  anchors <- data.frame(
    study_id  = "RCT9",
    drug      = "C",
    reference = "placebo",
    estimate  = 0.25,
    std_error = 0.05
  )
  absolute <- anchor_surface(fit, anchors = anchors)
  expect_identical(absolute$anchors$study_id, "RCT9")
  expect_identical(unique(absolute$effects$reference), "placebo")

  bad <- anchors
  bad$drug <- "Z"
  expect_error(anchor_surface(fit, anchors = bad), "appear in no comparison")
})

test_that("multiple anchors combine by precision, checked by hand", {
  skip_if_not_installed("netmeta")

  fit <- fit_surface(mandatory_network(), engine = "netmeta")
  # Two anchors that disagree: C proposes a higher surface than A.
  anchors <- data.frame(
    study_id  = c("RCT1", "RCT2"),
    drug      = c("C", "A"),
    reference = "placebo",
    estimate  = c(0.30, 0.60),
    std_error = c(0.04, 0.10)
  )
  absolute <- anchor_surface(fit, anchors = anchors)

  # Hand computation of the frequentist sea level: each anchor proposes
  # offset a_m - theta_surface, weighted by 1 / (a_se^2 + surface_se^2);
  # every absolute effect is surface + offset with the offset variance
  # added.
  surface <- fit$effects
  position <- match(anchors$drug, surface$drug)
  proposed <- anchors$estimate - surface$estimate[position]
  weight <- 1 / (anchors$std_error^2 + surface$std_error[position]^2)
  offset <- sum(weight * proposed) / sum(weight)
  offset_var <- 1 / sum(weight)

  expect_equal(absolute$effects$estimate,
               surface$estimate + offset, tolerance = 1e-10)
  expect_equal(absolute$effects$std_error,
               sqrt(surface$std_error^2 + offset_var), tolerance = 1e-10)

  # The disagreeing second anchor pulls the surface up, and the more
  # precise anchor dominates: the offset sits between the proposals,
  # nearer the tighter one.
  expect_true(offset > min(proposed) && offset < max(proposed))
  expect_lt(abs(offset - proposed[1]), abs(offset - proposed[2]))
})

test_that("anchor uncertainty propagates into absolute intervals", {
  skip_if_not_installed("netmeta")

  tight <- anchor_surface(fit_surface(mandatory_network(anchor_se = 0.04),
                                      engine = "netmeta"))
  loose <- anchor_surface(fit_surface(mandatory_network(anchor_se = 0.4),
                                      engine = "netmeta"))

  tight_width <- tight$effects$upper - tight$effects$lower
  loose_width <- loose$effects$upper - loose$effects$lower
  expect_true(all(loose_width > tight_width))

  # The anchor does not pin C exactly: C's absolute interval has width.
  c_row <- tight$effects[tight$effects$drug == "C", ]
  expect_gt(c_row$upper - c_row$lower, 0)
})

test_that("stan anchor uncertainty also propagates", {
  skip_if_not_installed("rstan")

  tight <- anchor_surface(
    fit_surface(mandatory_network(anchor_se = 0.04), engine = "stan",
                seed = 1),
    seed = 1
  )
  loose <- anchor_surface(
    fit_surface(mandatory_network(anchor_se = 0.4), engine = "stan",
                seed = 1),
    seed = 1
  )

  tight_width <- tight$effects$upper - tight$effects$lower
  loose_width <- loose$effects$upper - loose$effects$lower
  expect_true(all(loose_width > tight_width))
})

test_that("anchor_surface validates its input", {
  expect_error(anchor_surface(data.frame()), "directeffect_fit")
})
