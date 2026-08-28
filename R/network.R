#' Construct a direct-effect network
#'
#' Builds a validated evidence network from a table of comparative
#' (drug-versus-drug) effect estimates, optionally alongside a table of
#' absolute (placebo-anchored) estimates. The returned object defines the
#' evidence network — it contains no fitted model. Estimates for
#' multiplicative effect measures (HR, RR, OR) must be supplied on the log
#' scale.
#'
#' @param comparisons A data frame with one row per estimated
#'   drug-versus-drug effect. Required columns: `study_id`, `target`,
#'   `comparator`, `estimate`, `std_error`. Additional metadata columns
#'   (e.g. `database`, `population`, `design`, `outcome`) are preserved.
#' @param anchors Optional data frame of absolute (placebo-controlled)
#'   estimates. Required columns: `study_id`, `drug`, `reference`
#'   (normally `"placebo"`), `estimate`, `std_error`. Additional columns
#'   are preserved.
#' @param effect_measure A single string naming the effect measure the
#'   estimates are on, e.g. `"HR"`, `"RR"`, or `"OR"`.
#'
#' @return An object of class `directeffect_network` with components
#'   `comparisons`, `anchors`, `treatments`, `graph` (an igraph graph),
#'   `components` (integer component membership named by treatment),
#'   `effect_measure`, and `metadata`.
#'
#' @examples
#' comparisons <- data.frame(
#'   study_id   = c("S1", "S2", "S3"),
#'   target     = c("A", "A", "B"),
#'   comparator = c("B", "C", "C"),
#'   estimate   = c(log(1.02), log(1.34), log(1.29)),
#'   std_error  = c(0.07, 0.09, 0.08)
#' )
#' anchors <- data.frame(
#'   study_id  = "RCT1",
#'   drug      = "C",
#'   reference = "placebo",
#'   estimate  = log(1.20),
#'   std_error = 0.04
#' )
#' de <- direct_effect_network(comparisons, anchors = anchors,
#'                             effect_measure = "HR")
#' de
#' @export
direct_effect_network <- function(comparisons, anchors = NULL,
                                  effect_measure = "HR") {
  comparisons <- validate_comparisons(comparisons)
  anchors <- validate_anchors_table(anchors, comparisons)

  if (!is.character(effect_measure) || length(effect_measure) != 1 ||
      is.na(effect_measure) || !nzchar(effect_measure)) {
    stop("`effect_measure` must be a single non-missing string, ",
         "e.g. \"HR\", \"RR\", or \"OR\".", call. = FALSE)
  }

  warn_on_estimand_differences(comparisons)

  treatments <- sort(unique(c(comparisons$target, comparisons$comparator)))
  graph <- build_network_graph(comparisons, treatments)
  components <- igraph::components(graph)$membership
  storage.mode(components) <- "integer"

  structure(
    list(
      comparisons = comparisons,
      anchors = anchors,
      treatments = treatments,
      graph = graph,
      components = components,
      effect_measure = effect_measure,
      metadata = list(
        estimand_columns = intersect(estimand_columns(), names(comparisons))
      )
    ),
    class = "directeffect_network"
  )
}

#' @export
print.directeffect_network <- function(x, ...) {
  n_components <- max(x$components)
  n_anchors <- if (is.null(x$anchors)) 0L else nrow(x$anchors)
  cat("<directeffect_network>\n")
  cat("  Effect measure: ", x$effect_measure, " (log scale)\n", sep = "")
  cat("  Drugs:          ", length(x$treatments), "\n", sep = "")
  cat("  Comparisons:    ", nrow(x$comparisons), "\n", sep = "")
  cat("  Anchors:        ", n_anchors, "\n", sep = "")
  cat("  Components:     ", n_components, "\n", sep = "")
  cat("Use check_connectivity() for the per-component identifiability report.\n")
  invisible(x)
}

# The columns that describe the estimand a comparison targets. Differing
# values across comparisons signal that statistical connectivity may not
# imply causal transportability.
estimand_columns <- function() {
  c("population", "time_at_risk", "outcome_definition", "database",
    "study_design", "design", "outcome", "estimand")
}

validate_comparisons <- function(comparisons) {
  if (!is.data.frame(comparisons)) {
    stop("`comparisons` must be a data frame.", call. = FALSE)
  }
  required <- c("study_id", "target", "comparator", "estimate", "std_error")
  missing <- setdiff(required, names(comparisons))
  if (length(missing) > 0) {
    stop("`comparisons` is missing required column",
         if (length(missing) > 1) "s" else "", ": ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(comparisons) == 0) {
    stop("`comparisons` must contain at least one comparison.", call. = FALSE)
  }

  comparisons$target <- as.character(comparisons$target)
  comparisons$comparator <- as.character(comparisons$comparator)

  bad_drug <- is.na(comparisons$target) | !nzchar(comparisons$target) |
    is.na(comparisons$comparator) | !nzchar(comparisons$comparator)
  if (any(bad_drug)) {
    stop("`comparisons` has missing drug names in `target`/`comparator` ",
         "on row", if (sum(bad_drug) > 1) "s" else "", " ",
         paste(which(bad_drug), collapse = ", "), ".", call. = FALSE)
  }

  self <- comparisons$target == comparisons$comparator
  if (any(self)) {
    stop("`comparisons` contains self-comparisons (target == comparator) ",
         "on row", if (sum(self) > 1) "s" else "", " ",
         paste(which(self), collapse = ", "), ".", call. = FALSE)
  }

  if (!is.numeric(comparisons$estimate) || anyNA(comparisons$estimate)) {
    stop("`comparisons$estimate` must be numeric with no missing values.",
         call. = FALSE)
  }
  bad_se <- !is.numeric(comparisons$std_error) | is.na(comparisons$std_error) |
    comparisons$std_error <= 0
  if (any(bad_se)) {
    stop("`comparisons$std_error` must be positive on every row; row",
         if (sum(bad_se) > 1) "s" else "", " ",
         paste(which(bad_se), collapse = ", "), " ",
         if (sum(bad_se) > 1) "are" else "is", " not.", call. = FALSE)
  }

  comparisons
}

validate_anchors_table <- function(anchors, comparisons) {
  if (is.null(anchors)) {
    return(NULL)
  }
  if (!is.data.frame(anchors)) {
    stop("`anchors` must be a data frame or NULL.", call. = FALSE)
  }
  required <- c("study_id", "drug", "reference", "estimate", "std_error")
  missing <- setdiff(required, names(anchors))
  if (length(missing) > 0) {
    stop("`anchors` is missing required column",
         if (length(missing) > 1) "s" else "", ": ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(anchors) == 0) {
    return(NULL)
  }

  anchors$drug <- as.character(anchors$drug)
  anchors$reference <- as.character(anchors$reference)

  bad_drug <- is.na(anchors$drug) | !nzchar(anchors$drug) |
    is.na(anchors$reference) | !nzchar(anchors$reference)
  if (any(bad_drug)) {
    stop("`anchors` has missing drug names in `drug`/`reference` on row",
         if (sum(bad_drug) > 1) "s" else "", " ",
         paste(which(bad_drug), collapse = ", "), ".", call. = FALSE)
  }

  self <- anchors$drug == anchors$reference
  if (any(self)) {
    stop("`anchors` contains rows where `drug` equals `reference` on row",
         if (sum(self) > 1) "s" else "", " ",
         paste(which(self), collapse = ", "), ".", call. = FALSE)
  }

  if (!is.numeric(anchors$estimate) || anyNA(anchors$estimate)) {
    stop("`anchors$estimate` must be numeric with no missing values.",
         call. = FALSE)
  }
  bad_se <- !is.numeric(anchors$std_error) | is.na(anchors$std_error) |
    anchors$std_error <= 0
  if (any(bad_se)) {
    stop("`anchors$std_error` must be positive on every row; row",
         if (sum(bad_se) > 1) "s" else "", " ",
         paste(which(bad_se), collapse = ", "), " ",
         if (sum(bad_se) > 1) "are" else "is", " not.", call. = FALSE)
  }

  network_drugs <- unique(c(comparisons$target, comparisons$comparator))
  orphan <- setdiff(anchors$drug, network_drugs)
  if (length(orphan) > 0) {
    stop("`anchors` refers to drug", if (length(orphan) > 1) "s" else "", " ",
         paste0("\"", orphan, "\"", collapse = ", "),
         " that appear in no comparison. Anchors position drugs that are ",
         "part of the comparative network; a drug appearing only in ",
         "anchors is likely a data error.", call. = FALSE)
  }

  anchors
}

warn_on_estimand_differences <- function(comparisons) {
  present <- intersect(estimand_columns(), names(comparisons))
  differing <- present[vapply(present, function(col) {
    values <- unique(comparisons[[col]])
    values <- values[!is.na(values)]
    length(values) > 1
  }, logical(1))]
  if (length(differing) > 0) {
    warning("Comparisons differ in estimand metadata: ",
            paste0("`", differing, "`", collapse = ", "),
            ". Statistical connectivity does not establish causal ",
            "transportability; check that these comparisons target a ",
            "common estimand before interpreting a combined surface.",
            call. = FALSE)
  }
  invisible(differing)
}

build_network_graph <- function(comparisons, treatments) {
  edges <- data.frame(
    from = comparisons$target,
    to = comparisons$comparator,
    study_id = comparisons$study_id,
    estimate = comparisons$estimate,
    std_error = comparisons$std_error
  )
  igraph::graph_from_data_frame(edges, directed = FALSE,
                                vertices = data.frame(name = treatments))
}
