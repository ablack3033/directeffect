# The v0.1 goal, continuously enforced: the same direct-effect surface
# reconstructed independently by netmeta and Stan. Simulated under the
# shared model's own assumptions (Gaussian estimates, known standard
# errors, no heterogeneity) with weak priors, the engines must agree at
# every drug within EPSILON on the log scale. EPSILON = 0.02 covers
# Monte Carlo error (posterior-mean MCSE is ~0.002 at these settings)
# and the negligible shrinkage of the normal(0, 5) prior with an order
# of magnitude to spare, while still failing on any real modelling
# discrepancy.
EPSILON <- 0.02

test_that("netmeta and stan reconstruct the same surface (CI equivalence)", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("rstan")

  simulation <- simulate_direct_effect_network(
    n_drugs = 8, n_comparisons = 24, n_anchors = 0,
    heterogeneity = 0, seed = 20260828
  )
  fit_nm <- fit_surface(simulation$network, engine = "netmeta")
  fit_st <- fit_surface(simulation$network, engine = "stan",
                        seed = 20260828, iter = 4000)

  comparison <- compare_engines(fit_nm, fit_st)
  expect_lt(max(abs(comparison$difference)), EPSILON)
})

test_that("comparison table has the ticket schema and is order-agnostic", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("rstan")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_nm <- fit_surface(de, engine = "netmeta")
  fit_st <- fit_surface(de, engine = "stan", seed = 1)

  comparison <- compare_engines(fit_nm, fit_st)
  expect_identical(
    names(comparison),
    c("drug", "netmeta", "stan_mean", "difference",
      "standardized_difference")
  )
  expect_identical(comparison$drug, fit_nm$effects$drug)
  expect_equal(comparison$difference,
               comparison$stan_mean - comparison$netmeta)
  expect_true(is.na(
    comparison$standardized_difference[comparison$drug == "A"]
  ))

  swapped <- compare_engines(fit_st, fit_nm)
  expect_identical(comparison, swapped)
})

test_that("mismatched fits are rejected", {
  skip_if_not_installed("netmeta")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  fit_nm <- fit_surface(de, engine = "netmeta")

  expect_error(compare_engines(fit_nm, fit_nm), "one netmeta fit and one stan")
  expect_error(compare_engines(fit_nm, data.frame()), "directeffect_fit")

  skip_if_not_installed("rstan")
  fit_st_c <- fit_surface(de, engine = "stan", reference = "C", seed = 1)
  expect_error(compare_engines(fit_nm, fit_st_c), "different references")

  other <- direct_effect_network(
    data.frame(study_id = "S1", target = "X", comparator = "Y",
               estimate = 0.1, std_error = 0.05),
    effect_measure = "HR"
  )
  fit_other <- fit_surface(other, engine = "netmeta")
  skip_if_not_installed("rstan")
  fit_st <- fit_surface(de, engine = "stan", seed = 1)
  expect_error(compare_engines(fit_other, fit_st), "different drugs")
})

test_that("engine comparison plot is a ggplot with the identity line", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("rstan")
  skip_if_not_installed("ggplot2")

  de <- direct_effect_network(spec_comparisons(), effect_measure = "HR")
  plot <- plot_engine_comparison(
    fit_surface(de, engine = "netmeta"),
    fit_surface(de, engine = "stan", seed = 1)
  )
  expect_s3_class(plot, "ggplot")
})
