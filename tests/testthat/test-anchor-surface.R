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
    c("effects", "covariance", "comparisons", "anchors", "heterogeneity",
      "diagnostics", "engine", "engine_fit", "network")
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

test_that("anchoring matches the joint-GLS oracle on the mandatory example", {
  skip_if_not_installed("netmeta")

  de <- mandatory_network()
  absolute <- anchor_surface(fit_surface(de, engine = "netmeta"))

  # Independent one-stage oracle: anchor rows appended to the weighted
  # comparison design, normal equations solved directly. With a single
  # anchor the two-stage pipeline is exactly this joint GLS.
  oracle <- joint_gls(de$comparisons, de$anchors, c("A", "B", "C"))
  expect_equal(absolute$effects$estimate, unname(oracle$estimate),
               tolerance = 1e-8)
  expect_equal(absolute$effects$std_error, unname(oracle$std_error),
               tolerance = 1e-8)
  expect_equal(absolute$covariance, oracle$covariance, tolerance = 1e-8)

  # Closed form: the anchored drug's surface contribution cancels, so
  # its absolute SE is its anchor's SE, exactly.
  expect_equal(absolute$effects$std_error[absolute$effects$drug == "C"],
               0.04, tolerance = 1e-12)
})

test_that("disagreeing anchors: the GLS offset is the joint-GLS reference row", {
  skip_if_not_installed("netmeta")

  fit <- fit_surface(mandatory_network(), engine = "netmeta")
  # Two anchors that disagree: C proposes a higher surface than A does,
  # and A's anchor is less precise.
  anchors <- data.frame(
    study_id  = c("RCT1", "RCT2"),
    drug      = c("C", "A"),
    reference = "placebo",
    estimate  = c(0.30, 0.60),
    std_error = c(0.04, 0.10)
  )
  absolute <- anchor_surface(fit, anchors = anchors)
  oracle <- joint_gls(fit$comparisons, anchors, c("A", "B", "C"))
  surface <- fit$effects

  # The offset is a single location shift: relative structure untouched.
  shift <- absolute$effects$estimate - surface$estimate
  expect_equal(shift, rep(shift[1], 3), tolerance = 1e-10)

  # Profiling the surface out of the joint model is exact for the
  # location, so the offset and its variance equal the joint GLS at the
  # reference drug ("A", whose surface position is pinned at 0).
  expect_equal(absolute$effects$estimate[surface$drug == "A"],
               unname(oracle$estimate["A"]), tolerance = 1e-10)
  expect_equal(absolute$effects$std_error[surface$drug == "A"],
               unname(oracle$std_error["A"]), tolerance = 1e-10)

  # Elsewhere the location-shift estimator is honestly (slightly) less
  # efficient than the joint GLS, which also lets disagreeing anchors
  # update the relative structure: its SEs bound ours from below.
  expect_true(all(absolute$effects$std_error >=
                    unname(oracle$std_error) - 1e-10))

  # The disagreement is reconciled between the proposals, nearer the
  # more precise anchor.
  proposed <- anchors$estimate -
    surface$estimate[match(anchors$drug, surface$drug)]
  offset <- shift[1]
  expect_true(offset > min(proposed) && offset < max(proposed))
  expect_lt(abs(offset - proposed[1]), abs(offset - proposed[2]))
})

# Build a surface fit from the hand WLS oracle — proven equal to the
# netmeta surface to 1e-6 elsewhere in this suite — so Monte Carlo
# loops over the anchoring stage stay fast.
mc_surface_fit <- function(estimates, de) {
  comparisons <- de$comparisons
  comparisons$estimate <- estimates
  oracle <- wls_surface(comparisons, c("A", "B", "C"), "A")
  effects <- data.frame(
    drug = c("A", "B", "C"),
    estimate = c(0, unname(oracle$estimate[c("B", "C")])),
    std_error = c(0, unname(oracle$std_error[c("B", "C")])),
    lower = NA_real_, upper = NA_real_,
    scale = "log", reference = "A", engine = "netmeta"
  )
  directeffect:::new_directeffect_fit(
    effects = effects, covariance = oracle$covariance,
    heterogeneity = list(), engine = "netmeta", engine_fit = NULL,
    network = de
  )
}

# Replicate the two-stage frequentist path under the mandatory truth
# (A = 0.7, B = 0.7, C = 0.3), redrawing comparisons and anchors, and
# return per-drug estimates, coverage indicators, and the reported SEs
# (data-independent here: fixed design, known variances).
mc_anchoring <- function(anchors, reps) {
  truth <- c(A = 0.7, B = 0.7, C = 0.3)
  de <- mandatory_network()
  estimates <- matrix(NA_real_, reps, 3,
                      dimnames = list(NULL, names(truth)))
  covered <- matrix(NA, reps, 3, dimnames = list(NULL, names(truth)))
  reported_se <- NULL
  for (r in seq_len(reps)) {
    y <- stats::rnorm(3, c(0, 0.4, 0.4), 0.05)
    a <- anchors
    a$estimate <- stats::rnorm(nrow(a), truth[a$drug], a$std_error)
    absolute <- anchor_surface(mc_surface_fit(y, de), anchors = a)
    estimates[r, ] <- absolute$effects$estimate
    covered[r, ] <- absolute$effects$lower <= truth &
      truth <= absolute$effects$upper
    reported_se <- absolute$effects$std_error
  }
  list(estimates = estimates, covered = covered,
       reported_se = reported_se)
}

test_that("anchored 95% intervals attain nominal coverage (Monte Carlo)", {
  # The estimator is linear in Gaussian data with known variances, so
  # nominal coverage is exact in distribution and the empirical rate is
  # Binomial(400, 0.95) per drug: 4 binomial-MC standard deviations
  # (4 * sqrt(0.95 * 0.05 / 400) = 0.044) bounds the test at a
  # false-failure rate of ~6e-5 per drug — no ±10-point tolerance.
  set.seed(20260828)
  mc <- mc_anchoring(mandatory_network()$anchors, reps = 400)

  expect_true(all(abs(colMeans(mc$covered) - 0.95) < 0.044))
  # Reported SEs equal the estimator's actual sampling variability,
  # within ~4 MC standard errors of a sample SD at n = 400 (~15%).
  expect_true(all(abs(apply(mc$estimates, 2, sd) / mc$reported_se - 1)
                  < 0.15))
})

test_that("multi-anchor SEs equal the sampling variability (Monte Carlo)", {
  # Two anchors of different precision: the GLS weights come from the
  # full proposal covariance. Non-circular check of the reported SEs
  # against the estimator's empirical sampling SD.
  anchors <- data.frame(
    study_id  = c("RCT1", "RCT2"),
    drug      = c("C", "A"),
    reference = "placebo",
    estimate  = c(0.30, 0.70),
    std_error = c(0.04, 0.10)
  )
  set.seed(20260828)
  mc <- mc_anchoring(anchors, reps = 400)

  expect_true(all(abs(colMeans(mc$covered) - 0.95) < 0.044))
  expect_true(all(abs(apply(mc$estimates, 2, sd) / mc$reported_se - 1)
                  < 0.15))
})

test_that("engines agree on anchored standard errors (mandatory example)", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("rstan")

  de <- mandatory_network()
  nm <- anchor_surface(fit_surface(de, engine = "netmeta"))
  st <- anchor_surface(
    fit_surface(de, engine = "stan", seed = 20260828, iter = 4000),
    seed = 20260828, iter = 4000
  )
  st_effects <- st$effects[match(nm$effects$drug, st$effects$drug), ]

  # Posterior SDs are ~0.04-0.06 with an MCSE of ~1e-3 at these sampler
  # settings; 0.005 leaves a wide margin yet fails loudly on the old
  # mis-propagation (C: 0.070 reported vs 0.040 true).
  expect_lt(max(abs(nm$effects$std_error - st_effects$std_error)), 0.005)
  expect_lt(max(abs(nm$effects$estimate - st_effects$estimate)), 0.02)
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

  skip_if_not_installed("netmeta")
  legacy <- fit_surface(mandatory_network(), engine = "netmeta")
  legacy$covariance <- NULL
  expect_error(anchor_surface(legacy), "covariance")
})
