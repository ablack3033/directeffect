#' Plot the one-dimensional direct-effect surface
#'
#' Renders the fitted surface for a single outcome as it really is —
#' one-dimensional: each drug's position from harmful to beneficial with
#' its uncertainty interval, ordered by estimate. No 2-D layout pretends
#' to represent causal distance. For an anchored fit the placebo = 0
#' sea-level line is drawn and effects read as absolute; for a surface
#' fit the plot is labeled as relative to the arbitrary reference.
#'
#' @param fit A `directeffect_fit` from [fit_surface()] or
#'   [anchor_surface()]; both engines plot identically.
#' @param scale `"log"` draws effects on the log scale (0 = reference /
#'   sea level); `"natural"` back-transforms to the natural scale of a
#'   multiplicative measure on a logarithmic axis (1 = reference / sea
#'   level; log 0 corresponds to e.g. HR = 1).
#'
#' @return A ggplot object.
#'
#' @examples
#' comparisons <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(0.0, 0.4, 0.4),
#'   std_error  = c(0.05, 0.05, 0.05)
#' )
#' de <- direct_effect_network(comparisons, effect_measure = "HR")
#' if (requireNamespace("netmeta", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_effect_surface(fit_surface(de, engine = "netmeta"))
#' }
#' @export
plot_effect_surface <- function(fit, scale = c("log", "natural")) {
  scale <- match.arg(scale)
  effects <- surface_plot_data(fit, scale)
  anchored <- identical(effects$reference[1], "placebo")
  measure <- fit$network$effect_measure
  baseline <- if (scale == "log") 0 else 1

  labels <- surface_plot_labels(fit, scale, anchored, measure)

  plot <- ggplot2::ggplot(
    effects,
    ggplot2::aes(x = .data$estimate,
                 y = stats::reorder(.data$drug, .data$estimate))
  ) +
    ggplot2::geom_vline(
      xintercept = baseline,
      linetype = if (anchored) "solid" else "dashed",
      colour = if (anchored) "steelblue" else "grey55"
    ) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$lower,
                                          xmax = .data$upper)) +
    ggplot2::labs(
      x = labels$axis,
      y = NULL,
      title = "Direct-effect surface",
      subtitle = labels$subtitle
    )
  if (scale == "natural") {
    plot <- plot + ggplot2::scale_x_log10()
  }
  plot
}

#' Plot the surface positioned against sea level
#'
#' The signature figure of the package's conceptual decomposition: the
#' relative structure inferred from comparative evidence, positioned
#' against the placebo = 0 sea-level line by the absolute anchors. Drugs
#' spread along the horizontal axis in order of effect; their heights
#' are the absolute direct effects. Rendering stays one-dimensional in
#' the effect — the horizontal ordering is presentation, not geometry.
#' Only an anchored fit can be drawn: a relative surface has no sea
#' level.
#'
#' @param fit An anchored `directeffect_fit` from [anchor_surface()].
#' @param scale As in [plot_effect_surface()].
#'
#' @return A ggplot object.
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
#' if (requireNamespace("netmeta", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_sea_level(anchor_surface(fit_surface(de, engine = "netmeta")))
#' }
#' @export
plot_sea_level <- function(fit, scale = c("log", "natural")) {
  scale <- match.arg(scale)
  effects <- surface_plot_data(fit, scale)
  if (!identical(effects$reference[1], "placebo")) {
    stop("This fit is a relative surface (reference \"",
         effects$reference[1], "\"), so it has no sea level to plot ",
         "against. Anchor it first with `anchor_surface()`.",
         call. = FALSE)
  }
  measure <- fit$network$effect_measure
  baseline <- if (scale == "log") 0 else 1
  labels <- surface_plot_labels(fit, scale, anchored = TRUE, measure)

  plot <- ggplot2::ggplot(
    effects,
    ggplot2::aes(x = stats::reorder(.data$drug, .data$estimate),
                 y = .data$estimate)
  ) +
    ggplot2::geom_hline(yintercept = baseline, colour = "steelblue") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper)) +
    # After the point layer: an annotation with numeric x placed earlier
    # would type the discrete drug axis as continuous.
    ggplot2::annotate(
      "text",
      x = 1, y = baseline,
      label = if (scale == "log") "sea level: placebo = 0" else
        paste0("sea level: placebo (", measure, " = 1)"),
      vjust = -0.6, hjust = 0, colour = "steelblue", size = 3
    ) +
    ggplot2::labs(
      x = NULL,
      y = labels$axis,
      title = "Direct-effect surface and sea level",
      subtitle = labels$subtitle
    )
  if (scale == "natural") {
    plot <- plot + ggplot2::scale_y_log10()
  }
  plot
}

surface_plot_data <- function(fit, scale) {
  assert_directeffect_fit(fit)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The surface plots require the ggplot2 package.", call. = FALSE)
  }
  effects <- fit$effects[, c("drug", "estimate", "lower", "upper",
                             "reference")]
  if (scale == "natural") {
    effects$estimate <- exp(effects$estimate)
    effects$lower <- exp(effects$lower)
    effects$upper <- exp(effects$upper)
  }
  effects
}

surface_plot_labels <- function(fit, scale, anchored, measure) {
  axis <- if (scale == "log") {
    paste0("Direct effect, log ", measure)
  } else {
    paste0("Direct effect, ", measure, " (log-spaced axis)")
  }
  subtitle <- if (anchored) {
    paste0("Absolute direct effects vs placebo = ",
           if (scale == "log") "0" else paste0("1 (", measure, ")"),
           "; engine: ", fit$engine)
  } else {
    paste0("Relative surface: positions vs arbitrary reference \"",
           fit$effects$reference[1], "\" = ",
           if (scale == "log") "0" else "1",
           " (not absolute effects); engine: ", fit$engine)
  }
  list(axis = axis, subtitle = subtitle)
}
