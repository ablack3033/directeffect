#' Example head-to-head comparisons: a simulated statin network
#'
#' A realistic but **entirely simulated** evidence base modeled on the
#' statin trial literature: twelve head-to-head estimates among six
#' statins (atorvastatin, fluvastatin, lovastatin, pravastatin,
#' rosuvastatin, simvastatin) for a major-cardiovascular-events outcome,
#' on the log hazard ratio scale. The drug names are real; every number
#' is fake — drawn from the invented truth in [example_truth] — and must
#' not be interpreted clinically. Together with [example_anchors] it
#' powers the package vignettes. Generated deterministically by
#' `data-raw/example-network.R` with mild between-study heterogeneity
#' (tau = 0.02); the test suite verifies the shipped data matches that
#' generator and that the workflow recovers [example_truth].
#'
#' @format A data frame with 12 rows and 5 columns, in the comparisons
#'   input format documented in [directeffect_formats]:
#' \describe{
#'   \item{study_id}{Simulated trial identifier.}
#'   \item{target, comparator}{The two statins compared.}
#'   \item{estimate}{Simulated log hazard ratio, target vs comparator.}
#'   \item{std_error}{Its standard error.}
#' }
#' @source Simulated; see `data-raw/example-network.R`. The generating
#'   truth ships as [example_truth].
#' @seealso [directeffect_formats] for the full format reference.
"example_comparisons"

#' Example placebo anchors: simulated landmark trials
#'
#' Three simulated placebo-controlled estimates for statins in
#' [example_comparisons], on the log hazard ratio scale — the pattern of
#' the real statin literature, where a few landmark placebo trials
#' carry the absolute information while most comparisons are
#' head-to-head. Anchors position the relative surface against
#' placebo = 0. Like the comparisons, the numbers are **simulated**, not
#' real trial results.
#'
#' @format A data frame with 3 rows and 5 columns, in the anchors input
#'   format documented in [directeffect_formats]:
#' \describe{
#'   \item{study_id}{Simulated trial identifier.}
#'   \item{drug}{The anchored statin.}
#'   \item{reference}{Always `"placebo"`.}
#'   \item{estimate}{Simulated log hazard ratio, drug vs placebo.}
#'   \item{std_error}{Its standard error.}
#' }
#' @source Simulated; see `data-raw/example-network.R`. The generating
#'   truth ships as [example_truth].
#' @seealso [directeffect_formats] for the full format reference.
"example_anchors"

#' Generating truth of the example network
#'
#' The true direct effects the example data were simulated from, on the
#' log hazard ratio scale with `placebo = 0`. These values are
#' **inventions** calibrated to the magnitude of the statin literature
#' (hazard ratios roughly 0.7–0.85 versus placebo), not estimates of any
#' real drug's effect. Shipping the truth makes the example honest:
#' recovery of these values by `fit_surface()` + `anchor_surface()` on
#' [example_comparisons] and [example_anchors] is verified in the test
#' suite, not assumed.
#'
#' @format A data frame with 7 rows (placebo plus six statins) and
#'   2 columns:
#' \describe{
#'   \item{drug}{Drug name; the first row is `"placebo"`.}
#'   \item{theta}{True direct effect versus placebo, log hazard ratio
#'     scale (placebo row is exactly 0).}
#' }
#' @source Simulated; see `data-raw/example-network.R`.
"example_truth"
