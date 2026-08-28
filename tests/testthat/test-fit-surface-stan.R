# The Stan engine implements the identical common-effect model, so the
# same WLS hand oracle applies — up to MCMC noise and the negligible
# shrinkage from the deliberately weak normal(0, 5) prior.

spec_stan_fit <- function() {
  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_surface(de, engine = "stan", seed = 20260828)
}

test_that("stan surface matches the hand computation within MCMC tolerance", {
  skip_if_not_installed("rstan")

  fit <- spec_stan_fit()
  oracle <- wls_surface(spec_comparisons(), c("A", "B", "C"), "A")

  free <- fit$effects[fit$effects$drug != "A", ]
  expect_equal(free$estimate, unname(oracle$estimate[free$drug]),
               tolerance = 0.02)
  expect_equal(free$std_error, unname(oracle$std_error[free$drug]),
               tolerance = 0.1)

  ref_row <- fit$effects[fit$effects$drug == "A", ]
  expect_identical(ref_row$estimate, 0)
  expect_identical(ref_row$std_error, 0)
})

test_that("stan fit honours the shared contract with Bayesian extensions", {
  skip_if_not_installed("rstan")

  fit <- spec_stan_fit()

  expect_s3_class(fit, "directeffect_fit")
  expect_identical(
    names(fit),
    c("effects", "covariance", "comparisons", "anchors", "heterogeneity",
      "diagnostics", "engine", "engine_fit", "network")
  )
  expect_identical(
    names(fit$effects),
    c("drug", "estimate", "std_error", "lower", "upper", "scale",
      "reference", "engine",
      "median", "mean", "sd", "q025", "q975",
      "rhat", "ess_bulk", "ess_tail")
  )
  expect_identical(fit$effects$drug, c("A", "B", "C"))
  expect_identical(unique(fit$effects$engine), "stan")
  expect_identical(unique(fit$effects$reference), "A")
  expect_identical(unique(fit$effects$scale), "log")
  expect_s4_class(fit$engine_fit, "stanfit")
})

test_that("convergence diagnostics are reported and healthy runs are quiet", {
  skip_if_not_installed("rstan")

  fit <- spec_stan_fit()
  sampled <- fit$effects[fit$effects$drug != "A", ]
  expect_true(all(sampled$rhat < 1.01))
  expect_true(all(sampled$ess_bulk > 400))
  expect_true(all(sampled$ess_tail > 400))

  ref_row <- fit$effects[fit$effects$drug == "A", ]
  expect_true(is.na(ref_row$rhat))
  expect_true(is.na(ref_row$ess_bulk))
  expect_true(is.na(ref_row$ess_tail))
})

test_that("the stan covariance matches the WLS oracle and its own std_error", {
  skip_if_not_installed("rstan")

  fit <- spec_stan_fit()
  oracle <- wls_surface(spec_comparisons(), c("A", "B", "C"), "A")

  expect_identical(dimnames(fit$covariance),
                   list(c("A", "B", "C"), c("A", "B", "C")))
  # Posterior covariance vs the frequentist oracle: same likelihood, so
  # agreement within Monte Carlo error (entries are ~2e-3; at these
  # sampler settings their MCSE is a few percent of that).
  free <- c("B", "C")
  expect_equal(fit$covariance[free, free], oracle$covariance[free, free],
               tolerance = 0.1)
  # The pinned reference's row and column are exact zeros.
  expect_identical(unname(fit$covariance["A", ]), rep(0, 3))
  expect_identical(unname(fit$covariance[, "A"]), rep(0, 3))
  # Diagonal consistent with the reported std_error on the same draws.
  expect_equal(sqrt(diag(fit$covariance)),
               stats::setNames(fit$effects$std_error, fit$effects$drug),
               tolerance = 1e-6)
})

test_that("the stan engine refuses multi-arm networks, naming the way out", {
  skip_if_not_installed("rstan")

  de <- suppressWarnings(
    direct_effect_network(three_arm_comparisons(), effect_measure = "HR")
  )
  expect_error(fit_surface(de, engine = "stan"), "multi-arm")
  expect_error(fit_surface(de, engine = "stan"), "independent")
  expect_error(fit_surface(de, engine = "stan"), "netmeta")
})

test_that("an unseeded stan fit leaves the global RNG state untouched", {
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  set.seed(99)
  before <- .Random.seed
  suppressWarnings(fit_surface(de, engine = "stan", chains = 1, iter = 300))
  expect_identical(.Random.seed, before)

  # The anchored path (a refit with the anchored model) as well.
  anchored_de <- direct_effect_network(spec_comparisons(),
                                       anchors = spec_anchors(),
                                       effect_measure = "HR")
  surface <- suppressWarnings(
    fit_surface(anchored_de, engine = "stan", chains = 1, iter = 300)
  )
  set.seed(99)
  before <- .Random.seed
  suppressWarnings(anchor_surface(surface, chains = 1, iter = 300))
  expect_identical(.Random.seed, before)
})

test_that("the same explicit seed reproduces the same posterior summaries", {
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_a <- suppressWarnings(
    fit_surface(de, engine = "stan", chains = 2, iter = 500, seed = 123)
  )
  fit_b <- suppressWarnings(
    fit_surface(de, engine = "stan", chains = 2, iter = 500, seed = 123)
  )
  expect_equal(fit_a$effects, fit_b$effects)
  expect_equal(fit_a$covariance, fit_b$covariance)
})

test_that("a run with too few draws warns about convergence", {
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  warnings <- collect_warnings(
    fit_surface(de, engine = "stan", chains = 2, iter = 60,
                warmup = 30, seed = 1)
  )
  expect_true(any(grepl("convergence diagnostics indicate trouble",
                        warnings)))
})

test_that("a non-default reference is respected by the stan engine", {
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit <- fit_surface(de, engine = "stan", reference = "C", seed = 20260828)

  expect_identical(unique(fit$effects$reference), "C")
  expect_identical(fit$effects$drug, c("A", "B", "C"))
  expect_identical(fit$effects$estimate[fit$effects$drug == "C"], 0)

  oracle <- wls_surface(spec_comparisons(), c("A", "B", "C"), "C")
  free <- fit$effects[fit$effects$drug != "C", ]
  expect_equal(free$estimate, unname(oracle$estimate[free$drug]),
               tolerance = 0.02)
})
