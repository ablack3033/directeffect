# The mandatory example's comparisons are exactly consistent
# (A-B = 0, A-C = 0.4, B-C = 0.4); flipping B-C to 0.9 injects an
# inconsistency no surface can absorb.
consistent_comparisons <- function() {
  data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.0, 0.4, 0.4),
    std_error  = c(0.05, 0.05, 0.05)
  )
}

inconsistent_comparisons <- function() {
  comparisons <- consistent_comparisons()
  comparisons$estimate[3] <- 0.9
  comparisons
}

netmeta_fit <- function(comparisons) {
  de <- direct_effect_network(comparisons, effect_measure = "HR")
  fit_surface(de, engine = "netmeta")
}

test_that("edge residuals have the ticket schema and vanish on consistent data", {
  skip_if_not_installed("netmeta")

  residuals <- edge_residuals(netmeta_fit(consistent_comparisons()))

  expect_identical(
    names(residuals),
    c("target", "comparator", "observed", "predicted", "residual",
      "standardized_residual", "leverage")
  )
  expect_identical(nrow(residuals), 3L)
  expect_equal(residuals$residual, rep(0, 3), tolerance = 1e-8)
  expect_equal(residuals$standardized_residual, rep(0, 3),
               tolerance = 1e-6)
  # Each edge of the equal-precision triangle absorbs 2/3 of its own
  # variance into its prediction.
  expect_equal(residuals$leverage, rep(2 / 3, 3), tolerance = 1e-8)
})

test_that("standardized residuals match the hat-matrix oracle", {
  skip_if_not_installed("netmeta")

  comparisons <- asymmetric_comparisons()
  fit <- netmeta_fit(comparisons)
  residuals <- edge_residuals(fit)

  # Independent oracle from the weighted normal equations: the hat
  # matrix H = X (X'WX)^-1 X' W gives leverage as diag(H) and the
  # residual covariance as (I - H) V (I - H)'.
  oracle <- wls_surface(comparisons, c("A", "B", "C", "D"), "A")
  X <- oracle$design
  W <- oracle$weights
  V <- diag(comparisons$std_error^2)
  H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  residual_cov <- (diag(5) - H) %*% V %*% t(diag(5) - H)

  expect_equal(residuals$leverage, unname(diag(H)), tolerance = 1e-8)
  expect_equal(residuals$standardized_residual,
               residuals$residual / sqrt(diag(residual_cov)),
               tolerance = 1e-8)
})

test_that("a bridge comparison is NA with leverage 1, never 'consistent'", {
  skip_if_not_installed("netmeta")

  bridge_residuals <- function(cd_estimate) {
    comparisons <- data.frame(
      study_id   = c("S1", "S2", "S3", "S4"),
      target     = c("A", "A", "B", "C"),
      comparator = c("B", "C", "C", "D"),
      estimate   = c(0.1, 0.3, 0.25, cd_estimate),
      std_error  = c(0.05, 0.05, 0.05, 0.05)
    )
    edge_residuals(netmeta_fit(comparisons))
  }

  # Whatever the bridge's estimate, it reports leverage 1 and NA —
  # uncorroborated, not perfectly consistent (the old standardization
  # returned exactly 0 here).
  for (cd_estimate in c(0.2, 5.0)) {
    residuals <- bridge_residuals(cd_estimate)
    bridge_row <- residuals[residuals$target == "C" &
                              residuals$comparator == "D", ]
    expect_equal(bridge_row$leverage, 1, tolerance = 1e-8)
    expect_true(is.na(bridge_row$standardized_residual))
    # The cycle edges stay checkable.
    expect_true(all(!is.na(
      residuals$standardized_residual[residuals$comparator != "D"]
    )))
  }
})

test_that("diagnostics announce the independence assumption on multi-arm input", {
  skip_if_not_installed("netmeta")

  de <- suppressWarnings(
    direct_effect_network(three_arm_comparisons(), effect_measure = "HR")
  )
  fit <- fit_surface(de, engine = "netmeta")

  expect_warning(edge_residuals(fit), "independent comparisons")
  expect_warning(cycle_consistency(fit), "independent comparisons")

  # One comparison per study: no such warning.
  quiet_fit <- netmeta_fit(consistent_comparisons())
  expect_no_warning(edge_residuals(quiet_fit))
  expect_no_warning(cycle_consistency(quiet_fit))
})

test_that("cycle sums are ~ 0 on a consistent network", {
  skip_if_not_installed("netmeta")

  cycles <- cycle_consistency(netmeta_fit(consistent_comparisons()))

  expect_identical(
    names(cycles),
    c("cycle", "n_edges", "inconsistency", "std_error", "z")
  )
  expect_identical(nrow(cycles), 1L)
  expect_identical(cycles$n_edges, 3L)
  expect_equal(cycles$inconsistency, 0, tolerance = 1e-10)
})

test_that("an injected inconsistent edge is flagged by both diagnostics", {
  skip_if_not_installed("netmeta")

  fit <- netmeta_fit(inconsistent_comparisons())

  residuals <- edge_residuals(fit)
  expect_gt(max(abs(residuals$standardized_residual)), 1.96)

  cycles <- cycle_consistency(fit)
  expect_gt(max(abs(cycles$z)), 1.96)
  expect_equal(abs(cycles$inconsistency), 0.5, tolerance = 1e-10)
})

test_that("repeat comparisons of a pair are pooled by precision", {
  skip_if_not_installed("netmeta")

  comparisons <- rbind(
    consistent_comparisons(),
    data.frame(study_id = "S4", target = "C", comparator = "B",
               estimate = -0.4, std_error = 0.05)
  )
  cycles <- cycle_consistency(netmeta_fit(comparisons))

  # Still exactly one basis cycle; the reversed duplicate of B-C pools
  # into the same edge and the network stays consistent.
  expect_identical(nrow(cycles), 1L)
  expect_equal(cycles$inconsistency, 0, tolerance = 1e-10)
})

test_that("a network without cycles yields an empty cycle report", {
  skip_if_not_installed("netmeta")

  tree <- data.frame(
    study_id   = c("S1", "S2"),
    target     = c("A", "B"),
    comparator = c("B", "C"),
    estimate   = c(0.1, 0.2),
    std_error  = c(0.05, 0.05)
  )
  cycles <- cycle_consistency(netmeta_fit(tree))
  expect_identical(nrow(cycles), 0L)
  expect_identical(
    names(cycles),
    c("cycle", "n_edges", "inconsistency", "std_error", "z")
  )
})

test_that("diagnostics consume only the fit contract, on either engine", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("rstan")

  de <- direct_effect_network(consistent_comparisons(),
                              effect_measure = "HR")
  fit_nm <- fit_surface(de, engine = "netmeta")
  fit_stan <- fit_surface(de, engine = "stan", seed = 20260828)

  res_nm <- edge_residuals(fit_nm)
  res_stan <- edge_residuals(fit_stan)
  expect_identical(names(res_nm), names(res_stan))
  expect_equal(res_stan$predicted, res_nm$predicted, tolerance = 0.05)

  # Cycle consistency uses only observed data, so the engines agree
  # exactly.
  expect_equal(cycle_consistency(fit_stan), cycle_consistency(fit_nm))
})

test_that("cycle consistency plot is a ggplot", {
  skip_if_not_installed("netmeta")
  skip_if_not_installed("ggplot2")

  plot <- plot_cycle_consistency(netmeta_fit(consistent_comparisons()))
  expect_s3_class(plot, "ggplot")
})

test_that("diagnostics validate their input", {
  expect_error(edge_residuals(data.frame()), "directeffect_fit")
  expect_error(cycle_consistency(data.frame()), "directeffect_fit")

  skip_if_not_installed("netmeta")
  legacy <- netmeta_fit(consistent_comparisons())
  legacy$covariance <- NULL
  expect_error(edge_residuals(legacy), "covariance")
})
