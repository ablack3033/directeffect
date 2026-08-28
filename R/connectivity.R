#' Report per-component identifiability of a direct-effect network
#'
#' Comparative estimates identify only differences within a connected
#' component; a component with no absolute anchor supports relative
#' effects only, and the package never implies such a component can be
#' positioned absolutely. `check_connectivity()` reports, for each
#' connected component, the drug count, comparison count, anchor count,
#' and whether absolute effects are identifiable there.
#'
#' @param de A `directeffect_network` created by
#'   [direct_effect_network()].
#'
#' @return A data frame of class `directeffect_connectivity` with one row
#'   per connected component and columns `component`, `n_drugs`,
#'   `n_comparisons`, `n_anchors`, and `absolute_identifiable`. The
#'   `drugs` attribute holds the drug names per component. Printing
#'   renders the identifiability report.
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
#' check_connectivity(de)
#' @export
check_connectivity <- function(de) {
  if (!inherits(de, "directeffect_network")) {
    stop("`de` must be a `directeffect_network` created by ",
         "`direct_effect_network()`.", call. = FALSE)
  }

  membership <- de$components
  ids <- sort(unique(membership))

  drugs_by_component <- lapply(ids, function(id) {
    sort(names(membership)[membership == id])
  })
  comparison_component <- membership[de$comparisons$target]
  anchor_component <- if (is.null(de$anchors)) {
    integer(0)
  } else {
    membership[de$anchors$drug]
  }

  report <- data.frame(
    component = ids,
    n_drugs = vapply(drugs_by_component, length, integer(1)),
    n_comparisons = vapply(ids, function(id) {
      sum(comparison_component == id)
    }, integer(1)),
    n_anchors = vapply(ids, function(id) {
      sum(anchor_component == id)
    }, integer(1))
  )
  report$absolute_identifiable <- report$n_anchors > 0
  attr(report, "drugs") <- drugs_by_component
  class(report) <- c("directeffect_connectivity", "data.frame")
  report
}

#' @export
print.directeffect_connectivity <- function(x, ...) {
  for (i in seq_len(nrow(x))) {
    cat("Component ", x$component[i], ":\n", sep = "")
    cat("  ", x$n_drugs[i], " drug", if (x$n_drugs[i] != 1) "s" else "",
        "\n", sep = "")
    cat("  ", x$n_comparisons[i], " comparison",
        if (x$n_comparisons[i] != 1) "s" else "", "\n", sep = "")
    cat("  ", x$n_anchors[i], " absolute anchor",
        if (x$n_anchors[i] != 1) "s" else "", "\n", sep = "")
    cat("  ", if (x$absolute_identifiable[i]) {
      "absolute effects identifiable"
    } else {
      "relative effects identifiable only"
    }, "\n", sep = "")
  }
  invisible(x)
}
