#' Edge residuals: does each comparison agree with the fitted surface?
#'
#' For every observed comparison, reports the observed effect, the
#' effect predicted from the fitted surface (`theta_target -
#' theta_comparator`), their difference, and that difference divided by
#' the comparison's standard error. Large standardized residuals mark
#' comparisons that conflict with the surface implied by the rest of the
#' network. Consumes only the fit contract, so every engine is treated
#' identically. Rows align with `fit$comparisons`.
#'
#' @param fit A `directeffect_fit` from [fit_surface()].
#'
#' @return A data frame with columns `target`, `comparator`, `observed`,
#'   `predicted`, `residual`, and `standardized_residual`, one row per
#'   comparison in `fit$comparisons` order.
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
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   fit <- fit_surface(de, engine = "netmeta")
#'   edge_residuals(fit)
#' }
#' @export
edge_residuals <- function(fit) {
  if (!inherits(fit, "directeffect_fit")) {
    stop("`fit` must be a `directeffect_fit` from `fit_surface()`.",
         call. = FALSE)
  }

  comparisons <- fit$comparisons
  theta <- stats::setNames(fit$effects$estimate, fit$effects$drug)
  predicted <- unname(theta[comparisons$target] -
                        theta[comparisons$comparator])
  residual <- comparisons$estimate - predicted

  data.frame(
    target = comparisons$target,
    comparator = comparisons$comparator,
    observed = comparisons$estimate,
    predicted = predicted,
    residual = residual,
    standardized_residual = residual / comparisons$std_error,
    row.names = NULL
  )
}

#' Cycle consistency: do the comparisons agree around closed loops?
#'
#' Sums the observed effects around each cycle in a cycle basis of the
#' evidence network. Under perfect consistency every cycle sum is 0;
#' a cycle whose standardized sum is large contains comparisons that
#' cannot all be right. Repeat comparisons of the same drug pair are
#' pooled by precision first, and the basis is built from a spanning
#' tree — every cycle of the network is a combination of basis cycles,
#' so no cycle enumeration is ever needed.
#'
#' @param fit A `directeffect_fit` from [fit_surface()].
#'
#' @return A data frame with one row per basis cycle: `cycle` (the drug
#'   walk, e.g. `"A - B - C - A"`), `n_edges`, `inconsistency` (the
#'   signed cycle sum on the log scale), `std_error`, and `z`. Zero rows
#'   when the network has no cycles.
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
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   fit <- fit_surface(de, engine = "netmeta")
#'   cycle_consistency(fit)
#' }
#' @export
cycle_consistency <- function(fit) {
  if (!inherits(fit, "directeffect_fit")) {
    stop("`fit` must be a `directeffect_fit` from `fit_surface()`.",
         call. = FALSE)
  }

  pooled <- pool_parallel_edges(fit$comparisons)
  treatments <- fit$network$treatments

  empty <- data.frame(
    cycle = character(0), n_edges = integer(0),
    inconsistency = numeric(0), std_error = numeric(0), z = numeric(0)
  )
  if (nrow(pooled) < length(treatments)) {
    return(empty)
  }

  graph <- igraph::graph_from_data_frame(
    pooled[, c("target", "comparator")],
    directed = FALSE,
    vertices = data.frame(name = treatments)
  )
  tree <- igraph::mst(graph)
  tree_keys <- apply(igraph::as_edgelist(tree), 1, function(edge) {
    edge_key(edge[1], edge[2])
  })
  chords <- pooled[!pooled$key %in% tree_keys, ]
  if (nrow(chords) == 0) {
    return(empty)
  }

  cycles <- lapply(seq_len(nrow(chords)), function(i) {
    path <- igraph::shortest_paths(
      tree,
      from = chords$target[i], to = chords$comparator[i],
      output = "vpath"
    )$vpath[[1]]
    walk <- c(names(path), chords$target[i])
    cycle_row(walk, pooled)
  })
  cycles <- do.call(rbind, cycles)
  rownames(cycles) <- NULL
  cycles
}

#' @rdname cycle_consistency
#'
#' @return `plot_cycle_consistency()` returns a ggplot showing each
#'   basis cycle's standardized inconsistency with reference lines at
#'   z = -1.96, 0, and 1.96.
#' @importFrom rlang .data
#' @export
plot_cycle_consistency <- function(fit) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_cycle_consistency()` requires the ggplot2 package.",
         call. = FALSE)
  }
  cycles <- cycle_consistency(fit)
  ggplot2::ggplot(cycles, ggplot2::aes(x = .data$z, y = .data$cycle)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Standardized cycle inconsistency (z)",
      y = NULL,
      title = "Cycle consistency",
      subtitle = "Under a coherent surface every cycle sum is ~ 0"
    )
}

# One row per unordered drug pair: the precision-weighted pooled
# estimate oriented from the alphabetically first drug to the second.
pool_parallel_edges <- function(comparisons) {
  first <- pmin(comparisons$target, comparisons$comparator)
  second <- pmax(comparisons$target, comparisons$comparator)
  oriented <- ifelse(comparisons$target == first, 1, -1) *
    comparisons$estimate
  weight <- 1 / comparisons$std_error^2
  key <- vapply(seq_along(first), function(i) {
    edge_key(first[i], second[i])
  }, character(1))

  pooled <- lapply(split(seq_along(key), key), function(rows) {
    data.frame(
      key = key[rows[1]],
      target = first[rows[1]],
      comparator = second[rows[1]],
      estimate = sum(weight[rows] * oriented[rows]) / sum(weight[rows]),
      variance = 1 / sum(weight[rows])
    )
  })
  pooled <- do.call(rbind, pooled)
  rownames(pooled) <- NULL
  pooled
}

edge_key <- function(a, b) {
  paste(pmin(a, b), pmax(a, b), sep = "\r")
}

# Signed sum of pooled observed effects along a closed walk, with its
# standard error from the pooled edge variances.
cycle_row <- function(walk, pooled) {
  total <- 0
  variance <- 0
  for (i in seq_len(length(walk) - 1)) {
    from <- walk[i]
    to <- walk[i + 1]
    edge <- pooled[pooled$key == edge_key(from, to), ]
    sign <- if (edge$target == from) 1 else -1
    total <- total + sign * edge$estimate
    variance <- variance + edge$variance
  }
  std_error <- sqrt(variance)
  data.frame(
    cycle = paste(walk, collapse = " - "),
    n_edges = length(walk) - 1L,
    inconsistency = total,
    std_error = std_error,
    z = total / std_error
  )
}
