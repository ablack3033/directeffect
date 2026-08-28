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
#' @seealso [directeffect_formats] for the explicit input and output
#'   format reference, including validation rules.
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
  warn_on_multiarm_studies(comparisons)

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
  require_estimate_table(
    comparisons, "comparisons",
    c("study_id", "target", "comparator", "estimate", "std_error")
  )
  if (nrow(comparisons) == 0) {
    stop("`comparisons` must contain at least one comparison.", call. = FALSE)
  }

  comparisons$target <- as.character(comparisons$target)
  comparisons$comparator <- as.character(comparisons$comparator)
  require_drug_names(comparisons, "comparisons",
                     c("target", "comparator"))

  self <- comparisons$target == comparisons$comparator
  if (any(self)) {
    stop("`comparisons` contains self-comparisons (target == comparator) ",
         "on row", if (sum(self) > 1) "s" else "", " ",
         paste(which(self), collapse = ", "), ".", call. = FALSE)
  }

  require_estimate_values(comparisons, "comparisons")
  comparisons
}

validate_anchors_table <- function(anchors, comparisons) {
  if (is.null(anchors)) {
    return(NULL)
  }
  require_estimate_table(
    anchors, "anchors",
    c("study_id", "drug", "reference", "estimate", "std_error")
  )
  if (nrow(anchors) == 0) {
    return(NULL)
  }

  anchors$drug <- as.character(anchors$drug)
  anchors$reference <- as.character(anchors$reference)
  require_drug_names(anchors, "anchors", c("drug", "reference"))

  self <- anchors$drug == anchors$reference
  if (any(self)) {
    stop("`anchors` contains rows where `drug` equals `reference` on row",
         if (sum(self) > 1) "s" else "", " ",
         paste(which(self), collapse = ", "), ".", call. = FALSE)
  }

  require_estimate_values(anchors, "anchors")

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

# Shared table validation: both input tables are "estimate tables" with
# a study id, two drug-name columns, an estimate, and a standard error.

require_estimate_table <- function(table, name, required) {
  if (!is.data.frame(table)) {
    stop("`", name, "` must be a data frame.", call. = FALSE)
  }
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("`", name, "` is missing required column",
         if (length(missing) > 1) "s" else "", ": ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  invisible(table)
}

require_drug_names <- function(table, name, columns) {
  bad <- Reduce(`|`, lapply(columns, function(column) {
    is.na(table[[column]]) | !nzchar(table[[column]])
  }))
  if (any(bad)) {
    stop("`", name, "` has missing drug names in ",
         paste0("`", columns, "`", collapse = "/"), " on row",
         if (sum(bad) > 1) "s" else "", " ",
         paste(which(bad), collapse = ", "), ".", call. = FALSE)
  }
  invisible(table)
}

require_estimate_values <- function(table, name) {
  if (!is.numeric(table$estimate) || anyNA(table$estimate)) {
    stop("`", name, "$estimate` must be numeric with no missing values.",
         call. = FALSE)
  }
  bad_se <- !is.numeric(table$std_error) | is.na(table$std_error) |
    table$std_error <= 0
  if (any(bad_se)) {
    stop("`", name, "$std_error` must be positive on every row; row",
         if (sum(bad_se) > 1) "s" else "", " ",
         paste(which(bad_se), collapse = ", "), " ",
         if (sum(bad_se) > 1) "are" else "is", " not.", call. = FALSE)
  }
  invisible(table)
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

# Studies contributing more than one comparison row: multi-arm trials,
# whose within-study contrasts share arms and are therefore correlated.
multiarm_studies <- function(comparisons) {
  counts <- table(comparisons$study_id)
  names(counts)[counts > 1]
}

# "study \"T1\" contributes" / "studies \"T1\", \"T2\" contribute" —
# the subject clause shared by every multi-arm warning and refusal.
multiarm_clause <- function(studies) {
  paste0("stud", if (length(studies) > 1) "ies " else "y ",
         paste0("\"", studies, "\"", collapse = ", "),
         " contribute", if (length(studies) > 1) "" else "s")
}

# Multi-arm evidence is territory the two engines treat differently:
# netmeta models the within-study covariance correctly, while the Stan
# likelihood (which assumes independent rows) refuses to fit it. Say so
# at construction, before anything is fitted.
warn_on_multiarm_studies <- function(comparisons) {
  studies <- multiarm_studies(comparisons)
  if (length(studies) > 0) {
    clause <- multiarm_clause(studies)
    warning(toupper(substring(clause, 1, 1)), substring(clause, 2),
            " more than one comparison (a multi-arm trial). The ",
            "engines treat such input differently: the \"netmeta\" ",
            "engine models the correlation between contrasts sharing ",
            "an arm correctly, while the \"stan\" engine's ",
            "independence likelihood does not and will refuse to fit ",
            "this network. Diagnostics that assume independent ",
            "comparisons warn separately.", call. = FALSE)
  }
  invisible(studies)
}

# Per-drug attributes derived from the network: how many comparisons
# each drug appears in, and whether it carries an absolute anchor.
# Graph logic stays here with the network object, not in plotting code.
network_node_data <- function(de) {
  counts <- table(c(de$comparisons$target, de$comparisons$comparator))
  anchored_drugs <- if (is.null(de$anchors)) character(0) else de$anchors$drug
  data.frame(
    drug = de$treatments,
    n_comparisons = as.integer(counts[de$treatments]),
    anchored = de$treatments %in% anchored_drugs,
    row.names = NULL
  )
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
