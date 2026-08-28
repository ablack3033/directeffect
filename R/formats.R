#' Data formats: every input and output table, explicitly
#'
#' The package moves data through a small set of flat tables with fixed
#' schemas. This page is the single explicit reference for all of them:
#' the two input tables accepted by [direct_effect_network()] and the
#' tables every fitting and diagnostic function returns.
#'
#' @section Scale convention:
#' All estimates enter and leave the package on the **log scale** for
#' multiplicative effect measures (HR, RR, OR): pass `log(1.02)`, not
#' `1.02`. The package never log-transforms for you; `effect_measure`
#' only labels output. Plots can back-transform via `scale = "natural"`.
#'
#' @section The comparisons table (input):
#' One row per estimated drug-versus-drug effect (the same pair may
#' appear in any number of rows). Required columns:
#' \describe{
#'   \item{`study_id`}{Character (or coercible). Identifier of the study
#'     or analysis the estimate came from. Repeats are allowed.}
#'   \item{`target`}{Character. The drug the estimate is *for* (the
#'     numerator of a hazard ratio).}
#'   \item{`comparator`}{Character. The drug compared against (the
#'     denominator). Must differ from `target` on every row.}
#'   \item{`estimate`}{Numeric, no missing values. The estimated effect
#'     of `target` versus `comparator` on the log scale.}
#'   \item{`std_error`}{Numeric, strictly positive on every row. The
#'     standard error of `estimate` on the log scale.}
#' }
#' Any additional column is preserved untouched through the pipeline.
#' The estimand-describing columns — `population`, `time_at_risk`,
#' `outcome_definition`, `database`, `study_design`, `design`,
#' `outcome`, `estimand` — additionally trigger a warning at
#' construction when their values differ across comparisons, because
#' statistical connectivity does not establish causal transportability.
#'
#' Construction fails loudly on: a missing required column, zero rows,
#' a self-comparison (`target == comparator`), a missing or empty drug
#' name, a missing estimate, or a non-positive standard error.
#'
#' @section The anchors table (input):
#' Optional. One row per absolute (placebo-controlled) estimate:
#' \describe{
#'   \item{`study_id`}{Character (or coercible). Trial identifier.}
#'   \item{`drug`}{Character. The anchored drug. Must appear in at
#'     least one comparison — a drug appearing only in anchors is
#'     rejected as a likely data error.}
#'   \item{`reference`}{Character, normally `"placebo"`. What the drug
#'     was compared against; defines the zero of the absolute scale.}
#'   \item{`estimate`}{Numeric, no missing values. The estimated effect
#'     of `drug` versus `reference` on the log scale.}
#'   \item{`std_error`}{Numeric, strictly positive. Its standard error.}
#' }
#'
#' @section The effects table (output, `fit$effects`):
#' Every fit — either engine, anchored or not — carries a tidy per-drug
#' effects table with these columns, in this order:
#' \describe{
#'   \item{`drug`}{Character. One row per drug, in network
#'     (alphabetical) order.}
#'   \item{`estimate`}{Numeric. The drug's position on the log scale:
#'     relative to `reference` for a surface fit, absolute versus
#'     placebo for an anchored fit.}
#'   \item{`std_error`}{Numeric. Standard error (netmeta) or posterior
#'     standard deviation (Stan).}
#'   \item{`lower`, `upper`}{Numeric. 95% confidence limits (netmeta)
#'     or central 95% posterior interval (Stan).}
#'   \item{`scale`}{Character, always `"log"`.}
#'   \item{`reference`}{Character. The drug fixed at 0 for a surface
#'     fit; `"placebo"` for an anchored fit.}
#'   \item{`engine`}{Character, `"netmeta"` or `"stan"`.}
#' }
#' In an unanchored fit the reference drug's row is pinned exactly:
#' `estimate`, `std_error`, `lower`, and `upper` are all 0. Stan fits
#' append the posterior summary columns `median`, `mean`, `sd`, `q025`,
#' `q975` and the convergence diagnostics `rhat`, `ess_bulk`,
#' `ess_tail` (`NA` on a pinned reference row, which is a constant, not
#' a sampled quantity).
#'
#' @section The other fit components (output):
#' A `directeffect_fit` also carries `comparisons` and `anchors` (the
#' input tables, unchanged), `heterogeneity` (netmeta: list `Q`, `df`,
#' `p_value`, `tau`, `I2`; Stan common-effect: list `model`, `tau`),
#' `diagnostics` (list, reserved), `engine`, `engine_fit` (the raw
#' netmeta or stanfit object — the only place engine internals appear),
#' and `network` (the `directeffect_network` that was fitted).
#'
#' @section Diagnostic tables (output):
#' \describe{
#'   \item{[check_connectivity()]}{One row per connected component:
#'     `component`, `n_drugs`, `n_comparisons`, `n_anchors`,
#'     `absolute_identifiable`.}
#'   \item{[edge_residuals()]}{One row per comparison, aligned with
#'     `fit$comparisons`: `target`, `comparator`, `observed`,
#'     `predicted`, `residual`, `standardized_residual`.}
#'   \item{[cycle_consistency()]}{One row per basis cycle: `cycle`,
#'     `n_edges`, `inconsistency`, `std_error`, `z`. Zero rows when the
#'     network has no cycles.}
#'   \item{[compare_engines()]}{One row per drug: `drug`, `netmeta`,
#'     `stan_mean`, `difference`, `standardized_difference` (`NA` where
#'     both engines pin the reference exactly).}
#'   \item{[validate_recovery()]}{A list: `bias`, `rmse`, `coverage`,
#'     `rank_correlation`, `n_drugs`, `reference`.}
#' }
#'
#' @name directeffect_formats
#' @aliases data_formats
NULL
