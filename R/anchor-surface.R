#' Position a fitted surface against placebo (sea level)
#'
#' Applies absolute placebo anchors to a fitted relative surface as a
#' distinct second step, turning relative positions into absolute direct
#' effects against `placebo = 0`. Anchors retain their uncertainty: they
#' pull the surface to the right height without pinning any drug
#' exactly, and their standard errors propagate into every absolute
#' interval.
#'
#' The Bayesian path refits the network with the anchored Stan model, in
#' which no arbitrary identification constraint exists — the anchors
#' determine the absolute location. The frequentist path estimates the
#' single location offset that best reconciles the fitted surface with
#' the anchors by generalized least squares over the anchors' proposals,
#' whose covariance combines the anchor variances with the surface
#' covariance of the anchored drugs. Every absolute variance is
#' `Var(theta_d) + Var(offset) + 2 Cov(theta_d, offset)`, computed from
#' the fit's full surface covariance, so the absolute standard errors
#' equal the actual sampling variability of the estimator. With a single
#' anchor the anchored drug's surface contribution cancels exactly: its
#' absolute standard error is its anchor's standard error, and the whole
#' fit equals the joint generalized-least-squares solution of
#' comparisons plus anchor rows — the same answer the anchored Stan
#' model gives. With several disagreeing anchors the offset and its
#' variance still equal that joint solution at the surface's reference
#' drug, but the joint (Bayesian) fit also lets the anchors update the
#' drugs' relative positions, which a location shift by design does
#' not; the frequentist path then reports the honest — slightly larger
#' — sampling variance of its own estimator.
#'
#' A fit whose network component has no anchor cannot be positioned:
#' `anchor_surface()` refuses rather than silently picking a sea level.
#'
#' @param fit A `directeffect_fit` from [fit_surface()].
#' @param anchors Absolute estimates to anchor with (same schema as in
#'   [direct_effect_network()]). Defaults to the anchors already attached
#'   to the fit's network.
#' @param ... Passed to `rstan::sampling()` when re-fitting with the
#'   anchored Stan model (`chains`, `iter`, `seed`, ...). Unused by the
#'   frequentist path. As in [fit_surface()], when no `seed` is supplied
#'   the sampler is seeded from the wall clock and process id and the
#'   caller's global random-number state is left untouched.
#'
#' @return A `directeffect_fit` whose effects are absolute:
#'   `reference` is `"placebo"` and each estimate is the drug's direct
#'   effect versus placebo on the log scale.
#'
#' @seealso [directeffect_formats] for the explicit anchors input
#'   schema and the effects table schema.
#'
#' @examples
#' comparisons <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(0.0, 0.4, 0.4),
#'   std_error  = c(0.05, 0.05, 0.05)
#' )
#' anchors <- data.frame(
#'   study_id  = "RCT1",
#'   drug      = "C",
#'   reference = "placebo",
#'   estimate  = 0.3,
#'   std_error = 0.04
#' )
#' de <- direct_effect_network(comparisons, anchors = anchors,
#'                             effect_measure = "HR")
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   surface <- fit_surface(de, engine = "netmeta")
#'   absolute <- anchor_surface(surface)
#'   absolute$effects
#' }
#' @export
anchor_surface <- function(fit, anchors = NULL, ...) {
  assert_directeffect_fit(fit)

  de <- fit$network
  if (is.null(anchors)) {
    anchors <- de$anchors
  } else {
    anchors <- validate_anchors_table(anchors, de$comparisons)
  }
  if (is.null(anchors) || nrow(anchors) == 0) {
    stop("This network component has no absolute anchor, so its surface ",
         "cannot be positioned against placebo; the package never ",
         "silently picks a sea level. Attach anchors to the network or ",
         "pass them via `anchors`.", call. = FALSE)
  }

  switch(fit$engine,
    netmeta = anchor_surface_offset(fit, anchors),
    stan = anchor_surface_stan(fit, anchors, ...),
    stop("Unknown engine: ", fit$engine, call. = FALSE)
  )
}

# Frequentist sea level: a single location offset, with its full
# covariance propagated. Each anchor proposes offset
# a_m - theta_surface_drug(m); the proposals' covariance is
# diag(a_se^2) + Sigma[anchored, anchored] where Sigma is the surface
# covariance, and the offset is their generalized-least-squares mean.
# Every absolute effect is surface + offset, with variance
# Var(theta_d) + Var(offset) + 2 Cov(theta_d, offset); the covariance
# term is what makes the algebra of the evidence honest — with a single
# anchor the anchored drug's surface contribution cancels exactly and
# its absolute SE is the anchor's SE. Profiling the surface out of the
# joint model shows the offset and its variance equal the joint GLS of
# comparisons plus anchor rows at the reference drug; with one anchor
# the whole fit is the joint GLS.
anchor_surface_offset <- function(fit, anchors) {
  assert_fit_covariance(fit, "anchor_surface()")
  effects <- fit$effects
  covariance <- fit$covariance
  n_drugs <- nrow(effects)
  anchored_drugs <- anchors$drug
  position <- match(anchored_drugs, effects$drug)

  proposed <- anchors$estimate - effects$estimate[position]
  proposal_cov <- diag(anchors$std_error^2, nrow = nrow(anchors)) +
    covariance[anchored_drugs, anchored_drugs, drop = FALSE]
  weight <- solve(proposal_cov, rep(1, nrow(anchors)))
  offset <- sum(weight * proposed) / sum(weight)
  offset_var <- 1 / sum(weight)
  # Cov(theta_d, offset) for every drug d: the offset subtracts the
  # weighted anchored surface positions, hence the minus sign.
  offset_cov <- -drop(t(covariance[anchored_drugs, , drop = FALSE]) %*%
                        (weight / sum(weight)))

  absolute_cov <- covariance +
    outer(offset_cov, rep(1, n_drugs)) +
    outer(rep(1, n_drugs), offset_cov) +
    offset_var
  dimnames(absolute_cov) <- list(effects$drug, effects$drug)

  z <- stats::qnorm(0.975)
  absolute <- data.frame(
    drug = effects$drug,
    estimate = effects$estimate + offset,
    std_error = sqrt(diag(absolute_cov)),
    scale = "log",
    reference = "placebo",
    engine = fit$engine,
    row.names = NULL
  )
  absolute$lower <- absolute$estimate - z * absolute$std_error
  absolute$upper <- absolute$estimate + z * absolute$std_error
  absolute <- absolute[, c("drug", "estimate", "std_error", "lower",
                           "upper", "scale", "reference", "engine")]

  anchored <- new_directeffect_fit(
    effects = absolute,
    covariance = absolute_cov,
    heterogeneity = fit$heterogeneity,
    engine = fit$engine,
    engine_fit = fit$engine_fit,
    network = fit$network
  )
  anchored$anchors <- anchors
  anchored
}
