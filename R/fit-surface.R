#' Fit the relative direct-effect surface
#'
#' Estimates the relative positions of all drugs in a connected
#' direct-effect network from the comparative estimates alone. Anchors
#' are deliberately ignored at this stage: the surface answers "are my
#' comparative estimates internally coherent, and what relative structure
#' do they imply?" — not where that structure sits in absolute terms.
#' Identification uses an arbitrary reference drug fixed at 0; use
#' `anchor_surface()` to position the surface absolutely.
#'
#' @param de A `directeffect_network` created by
#'   [direct_effect_network()]. Must be a single connected component; fit
#'   components separately (see [check_connectivity()]) otherwise.
#' @param engine `"netmeta"` (frequentist network meta-analysis) or
#'   `"stan"` (Bayesian, identical likelihood). Both engines fit the same
#'   common-effect model `y_k ~ N(theta_target - theta_comparator,
#'   se_k^2)` and return the same fit contract.
#' @param reference Drug fixed at 0 for identification. Defaults to the
#'   first treatment alphabetically. The choice is arbitrary and does not
#'   affect any estimated difference between drugs.
#'
#' @return An object of class `directeffect_fit` with components
#'   `effects` (tidy per-drug table: `drug`, `estimate`, `std_error`,
#'   `lower`, `upper`, `scale`, `reference`, `engine`), `comparisons`,
#'   `anchors`, `heterogeneity`, `diagnostics`, `engine`, `engine_fit`
#'   (the raw engine object — the only place engine internals appear),
#'   and `network`.
#'
#' @examples
#' comparisons <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(log(1.02), log(1.34), log(1.29)),
#'   std_error  = c(0.07, 0.09, 0.08)
#' )
#' de <- direct_effect_network(comparisons, effect_measure = "HR")
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   fit <- fit_surface(de, engine = "netmeta")
#'   fit$effects
#' }
#' @export
fit_surface <- function(de, engine = c("netmeta", "stan"),
                        reference = NULL) {
  if (!inherits(de, "directeffect_network")) {
    stop("`de` must be a `directeffect_network` created by ",
         "`direct_effect_network()`.", call. = FALSE)
  }
  engine <- match.arg(engine)

  n_components <- max(de$components)
  if (n_components > 1) {
    stop("The network has ", n_components, " connected components; ",
         "relative effects are only defined within a component. ",
         "Run `check_connectivity(de)` to see the components, then fit ",
         "each connected sub-network separately.", call. = FALSE)
  }

  if (is.null(reference)) {
    reference <- de$treatments[1]
  }
  if (!is.character(reference) || length(reference) != 1 ||
      !reference %in% de$treatments) {
    stop("`reference` must be one of the network's treatments: ",
         paste0("\"", de$treatments, "\"", collapse = ", "), ".",
         call. = FALSE)
  }

  switch(engine,
    netmeta = fit_surface_netmeta(de, reference),
    stan = stop("The Stan engine is not implemented yet; ",
                "use `engine = \"netmeta\"`.", call. = FALSE)
  )
}

# Assemble the engine-agnostic fit object. Every engine adapter funnels
# through this constructor so the contract stays in one place.
new_directeffect_fit <- function(effects, heterogeneity, engine,
                                 engine_fit, network) {
  structure(
    list(
      effects = effects,
      comparisons = network$comparisons,
      anchors = network$anchors,
      heterogeneity = heterogeneity,
      diagnostics = list(),
      engine = engine,
      engine_fit = engine_fit,
      network = network
    ),
    class = "directeffect_fit"
  )
}

#' @export
print.directeffect_fit <- function(x, digits = 3, ...) {
  reference <- x$effects$reference[1]
  cat("<directeffect_fit>\n")
  cat("  Engine:         ", x$engine, "\n", sep = "")
  cat("  Effect measure: ", x$network$effect_measure, " (log scale)\n",
      sep = "")
  cat("  Reference:      ", reference,
      " (arbitrary; surface is relative)\n", sep = "")
  cat("\n")
  shown <- x$effects[, c("drug", "estimate", "std_error", "lower", "upper")]
  shown[-1] <- lapply(shown[-1], round, digits = digits)
  print.data.frame(shown, row.names = FALSE)
  invisible(x)
}
