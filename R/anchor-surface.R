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
#' single precision-weighted location offset that best reconciles the
#' fitted surface with the anchors, propagating both surface and anchor
#' uncertainty.
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
#'   frequentist path.
#'
#' @return A `directeffect_fit` whose effects are absolute:
#'   `reference` is `"placebo"` and each estimate is the drug's direct
#'   effect versus placebo on the log scale.
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

# Frequentist sea level: a single location offset. Each anchor proposes
# offset a_m - theta_surface_drug(m) with variance a_se_m^2 +
# surface_se_drug(m)^2; the estimate is their precision-weighted mean.
# Every absolute effect is then surface + offset, with the offset's
# variance added to the surface variance.
anchor_surface_offset <- function(fit, anchors) {
  effects <- fit$effects
  position <- match(anchors$drug, effects$drug)

  proposed <- anchors$estimate - effects$estimate[position]
  variance <- anchors$std_error^2 + effects$std_error[position]^2
  weight <- 1 / variance
  offset <- sum(weight * proposed) / sum(weight)
  offset_var <- 1 / sum(weight)

  z <- stats::qnorm(0.975)
  absolute <- data.frame(
    drug = effects$drug,
    estimate = effects$estimate + offset,
    std_error = sqrt(effects$std_error^2 + offset_var),
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
    heterogeneity = fit$heterogeneity,
    engine = fit$engine,
    engine_fit = fit$engine_fit,
    network = fit$network
  )
  anchored$anchors <- anchors
  anchored
}
