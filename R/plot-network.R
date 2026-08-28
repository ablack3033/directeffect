#' Plot the evidence network
#'
#' Draws the direct-effect network — drugs as nodes, comparative
#' estimates as edges — from a bare network object, before any model is
#' fitted. Optional visual encodings: edge thickness for precision,
#' node size for the number of comparisons a drug appears in, node
#' shape for the presence of an absolute anchor, node colour for drug
#' class, and edge labels showing the observed effects. Repeat
#' comparisons of the same pair draw as separate fanned edges, and
#' multi-component networks render with each component laid out
#' separately.
#'
#' @param de A `directeffect_network` created by
#'   [direct_effect_network()].
#' @param weight_edges Encode each comparison's precision (`1 /
#'   std_error^2`) as edge thickness.
#' @param size_nodes Encode the number of comparisons a drug appears in
#'   as node size.
#' @param shape_anchors Encode the presence of an absolute anchor as
#'   node shape (triangle = anchored, circle = not).
#' @param label_edges Label each edge with its observed effect
#'   (log scale, 2 decimals). Off by default.
#' @param drug_classes Optional data frame with columns `drug` and
#'   `class` to colour nodes by drug class. Drugs not listed are shown
#'   as unclassified.
#'
#' @return A ggplot (ggraph) object.
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
#' if (requireNamespace("ggraph", quietly = TRUE)) {
#'   plot_network(de)
#' }
#' @export
plot_network <- function(de, weight_edges = TRUE, size_nodes = TRUE,
                         shape_anchors = TRUE, label_edges = FALSE,
                         drug_classes = NULL) {
  if (!inherits(de, "directeffect_network")) {
    stop("`de` must be a `directeffect_network` created by ",
         "`direct_effect_network()`.", call. = FALSE)
  }
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_network()` requires the ggraph and ggplot2 packages.",
         call. = FALSE)
  }

  nodes <- network_node_data(de)
  graph <- de$graph
  position <- match(igraph::V(graph)$name, nodes$drug)
  graph <- igraph::set_vertex_attr(graph, "n_comparisons",
                                   value = nodes$n_comparisons[position])
  graph <- igraph::set_vertex_attr(
    graph, "anchored",
    value = ifelse(nodes$anchored[position], "anchored", "not anchored")
  )
  graph <- igraph::set_edge_attr(
    graph, "precision",
    value = 1 / igraph::E(graph)$std_error^2
  )
  graph <- igraph::set_edge_attr(
    graph, "effect_label",
    value = sprintf("%.2f", igraph::E(graph)$estimate)
  )
  if (!is.null(drug_classes)) {
    if (!is.data.frame(drug_classes) ||
        !all(c("drug", "class") %in% names(drug_classes))) {
      stop("`drug_classes` must be a data frame with columns `drug` ",
           "and `class`.", call. = FALSE)
    }
    class_of <- drug_classes$class[match(igraph::V(graph)$name,
                                         drug_classes$drug)]
    class_of[is.na(class_of)] <- "unclassified"
    graph <- igraph::set_vertex_attr(graph, "class", value = class_of)
  }

  # Aesthetic mappings are assembled as quoted expressions so each
  # encoding can be toggled independently.
  edge_mappings <- list()
  if (weight_edges) edge_mappings$edge_width <- quote(precision)
  if (label_edges) edge_mappings$label <- quote(effect_label)
  edge_aes <- do.call(ggplot2::aes, edge_mappings)
  edge_layer <- if (label_edges) {
    ggraph::geom_edge_fan(edge_aes, alpha = 0.5, colour = "grey40",
                          angle_calc = "along",
                          label_dodge = grid::unit(3, "mm"))
  } else {
    ggraph::geom_edge_fan(edge_aes, alpha = 0.5, colour = "grey40")
  }

  node_mappings <- list()
  if (size_nodes) node_mappings$size <- quote(n_comparisons)
  if (shape_anchors) node_mappings$shape <- quote(anchored)
  if (!is.null(drug_classes)) node_mappings$colour <- quote(class)
  node_aes <- do.call(ggplot2::aes, node_mappings)

  plot <- ggraph::ggraph(graph, layout = "stress") +
    edge_layer +
    ggraph::geom_node_point(node_aes) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                           repel = TRUE, size = 3) +
    ggraph::scale_edge_width(range = c(0.3, 2),
                             name = "Precision (1 / se\u00b2)") +
    ggplot2::scale_size(range = c(2, 7), name = "Comparisons") +
    ggplot2::scale_shape_manual(
      values = c(anchored = 17, `not anchored` = 16),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Direct-effect evidence network",
      subtitle = paste0("Effect measure: ", de$effect_measure,
                        " (log scale)")
    ) +
    ggraph::theme_graph(base_family = "")
  plot
}
