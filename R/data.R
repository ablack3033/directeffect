#' Example comparative-effect estimates
#'
#' A simulated evidence base of eight fictional drugs (abrelor, belvane,
#' cordexa, delmiro, envarin, folzane, gravix, harlent) compared
#' head-to-head: sixteen active-comparator estimates on the log hazard
#' ratio scale. Together with [example_anchors] it powers the package
#' vignettes. Generated deterministically by `data-raw/example-network.R`
#' with mild between-study heterogeneity (tau = 0.02); the generating
#' truth is recorded in that script.
#'
#' @format A data frame with 16 rows and 5 columns:
#' \describe{
#'   \item{study_id}{Study identifier.}
#'   \item{target, comparator}{The two drugs compared.}
#'   \item{estimate}{Estimated log hazard ratio, target vs comparator.}
#'   \item{std_error}{Its standard error.}
#' }
#' @source Simulated; see `data-raw/example-network.R`.
"example_comparisons"

#' Example placebo anchors
#'
#' Two simulated placebo-controlled estimates for drugs in
#' [example_comparisons], on the log hazard ratio scale. Anchors carry
#' the absolute information that positions the relative surface against
#' placebo = 0.
#'
#' @format A data frame with 2 rows and 5 columns:
#' \describe{
#'   \item{study_id}{Trial identifier.}
#'   \item{drug}{The anchored drug.}
#'   \item{reference}{Always `"placebo"`.}
#'   \item{estimate}{Estimated log hazard ratio, drug vs placebo.}
#'   \item{std_error}{Its standard error.}
#' }
#' @source Simulated; see `data-raw/example-network.R`.
"example_anchors"
