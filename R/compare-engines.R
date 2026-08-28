#' Compare the two engines' reconstructions of the same surface
#'
#' The package's central v0.1 claim is that netmeta and Stan reconstruct
#' the same direct-effect surface independently. `compare_engines()`
#' makes that a routine check: it lines the two fits up per drug and
#' reports the difference and standardized difference between the
#' engines' estimates. The fits may be passed in either order; they must
#' cover the same drugs and use the same reference so the comparison is
#' apples-to-apples (refit with the same `reference`, or anchor both,
#' otherwise).
#'
#' @param fit_a,fit_b Two `directeffect_fit` objects for the same
#'   network — one from the `"netmeta"` engine and one from `"stan"`,
#'   in either order.
#'
#' @return A data frame with one row per drug: `drug`, `netmeta` (the
#'   frequentist estimate), `stan_mean` (the posterior mean),
#'   `difference` (`stan_mean - netmeta`), and `standardized_difference`
#'   (difference over the combined standard error; `NA` for the
#'   reference drug, whose estimate is exact 0 in both engines).
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
#'     requireNamespace("rstan", quietly = TRUE) && interactive()) {
#'   comparison <- compare_engines(
#'     fit_surface(de, engine = "netmeta"),
#'     fit_surface(de, engine = "stan")
#'   )
#'   comparison
#' }
#' @export
compare_engines <- function(fit_a, fit_b) {
  assert_directeffect_fit(fit_a)
  assert_directeffect_fit(fit_b)

  engines <- c(fit_a$engine, fit_b$engine)
  if (!setequal(engines, c("netmeta", "stan"))) {
    stop("`compare_engines()` needs one netmeta fit and one stan fit; ",
         "got engines ", paste0("\"", engines, "\"", collapse = " and "),
         ".", call. = FALSE)
  }
  nm <- if (fit_a$engine == "netmeta") fit_a else fit_b
  st <- if (fit_a$engine == "stan") fit_a else fit_b

  if (!setequal(nm$effects$drug, st$effects$drug)) {
    stop("The two fits cover different drugs; fit both engines on the ",
         "same network before comparing.", call. = FALSE)
  }
  ref_nm <- nm$effects$reference[1]
  ref_st <- st$effects$reference[1]
  if (!identical(ref_nm, ref_st)) {
    stop("The two fits use different references (\"", ref_nm,
         "\" vs \"", ref_st, "\"), so their effects are not the same ",
         "quantity. Refit with the same `reference`, or anchor both ",
         "fits, before comparing.", call. = FALSE)
  }

  st_effects <- st$effects[match(nm$effects$drug, st$effects$drug), ]
  difference <- st_effects$estimate - nm$effects$estimate
  combined_se <- sqrt(nm$effects$std_error^2 + st_effects$std_error^2)
  standardized <- ifelse(combined_se > 0, difference / combined_se,
                         NA_real_)

  data.frame(
    drug = nm$effects$drug,
    netmeta = nm$effects$estimate,
    stan_mean = st_effects$estimate,
    difference = difference,
    standardized_difference = standardized,
    row.names = NULL
  )
}

#' @rdname compare_engines
#'
#' @return `plot_engine_comparison()` returns a ggplot of the netmeta
#'   estimates against the Stan posterior means with the identity line
#'   y = x; points on the line mean the engines agree.
#' @export
plot_engine_comparison <- function(fit_a, fit_b) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_engine_comparison()` requires the ggplot2 package.",
         call. = FALSE)
  }
  comparison <- compare_engines(fit_a, fit_b)
  ggplot2::ggplot(
    comparison,
    ggplot2::aes(x = .data$netmeta, y = .data$stan_mean)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    ggplot2::geom_point() +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "netmeta estimate (log scale)",
      y = "Stan posterior mean (log scale)",
      title = "Engine comparison",
      subtitle = "Points on the identity line: both engines reconstruct the same surface"
    )
}
