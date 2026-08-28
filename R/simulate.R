#' Simulate a direct-effect network with known truth
#'
#' Generates true direct effects, a connected comparison graph, standard
#' errors, observed comparative estimates, and placebo anchor estimates.
#' The returned object retains the truth so recovery can be measured
#' rather than assumed — see [validate_recovery()].
#'
#' The generative model matches the package's fitting model: observed
#' comparisons are `y ~ N(theta_target - theta_comparator, se^2 + tau^2)`
#' where `tau` is `heterogeneity` (0 gives exact model match), and
#' anchors are `a ~ N(theta_drug, se^2)` against `reference = "placebo"`
#' at 0. The comparison graph is connected by construction: a random
#' spanning tree first, then the remaining comparisons between random
#' drug pairs (repeat comparisons of the same pair are allowed, as in
#' real evidence networks).
#'
#' @param n_drugs Number of drugs.
#' @param n_comparisons Number of comparative estimates; at least
#'   `n_drugs - 1` so the network can be connected.
#' @param n_anchors Number of placebo anchors, each on a distinct drug
#'   (at most `n_drugs`).
#' @param heterogeneity Between-study standard deviation `tau` added to
#'   every comparison's sampling variance. 0 matches the v0.1
#'   common-effect model exactly.
#' @param seed Seed for reproducibility. The global random-number state
#'   is restored afterwards.
#' @param effect_sd Standard deviation of the true direct effects around
#'   0 on the log scale.
#' @param se_range Range the comparison standard errors are drawn from,
#'   uniformly.
#' @param anchor_se_range Range the anchor standard errors are drawn
#'   from, uniformly.
#'
#' @return A list with components `network` (a ready-to-fit
#'   `directeffect_network`), `comparisons`, `anchors`, `truth` (data
#'   frame of `drug`, `theta`, with the placebo row at 0), `tau`, and
#'   `seed`.
#'
#' @examples
#' simulation <- simulate_direct_effect_network(
#'   n_drugs = 5, n_comparisons = 10, n_anchors = 1,
#'   heterogeneity = 0, seed = 1
#' )
#' simulation$truth
#' @export
simulate_direct_effect_network <- function(n_drugs = 20,
                                           n_comparisons = 100,
                                           n_anchors = 3,
                                           heterogeneity = 0.1,
                                           seed = 1,
                                           effect_sd = 0.5,
                                           se_range = c(0.05, 0.15),
                                           anchor_se_range = c(0.03, 0.1)) {
  if (n_drugs < 2) {
    stop("`n_drugs` must be at least 2.", call. = FALSE)
  }
  if (n_comparisons < n_drugs - 1) {
    stop("`n_comparisons` must be at least `n_drugs - 1` (", n_drugs - 1,
         ") so the comparison graph can be connected.", call. = FALSE)
  }
  if (n_anchors < 0 || n_anchors > n_drugs) {
    stop("`n_anchors` must be between 0 and `n_drugs`.", call. = FALSE)
  }
  if (heterogeneity < 0) {
    stop("`heterogeneity` must be non-negative.", call. = FALSE)
  }

  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
  }
  set.seed(seed)

  drugs <- sprintf("drug_%02d", seq_len(n_drugs))
  theta <- stats::rnorm(n_drugs, mean = 0, sd = effect_sd)
  names(theta) <- drugs

  # Spanning tree: each drug after the first compares against a random
  # earlier drug, guaranteeing one connected component; then the rest.
  tree_target <- seq(2, n_drugs)
  tree_comparator <- vapply(tree_target, function(k) {
    sample.int(k - 1, 1)
  }, integer(1))
  n_extra <- n_comparisons - (n_drugs - 1)
  extra <- t(vapply(seq_len(n_extra), function(i) {
    sample.int(n_drugs, 2)
  }, integer(2)))
  target_idx <- c(tree_target, extra[, 1])
  comparator_idx <- c(tree_comparator, extra[, 2])

  se <- stats::runif(n_comparisons, se_range[1], se_range[2])
  delta <- theta[target_idx] - theta[comparator_idx]
  y <- stats::rnorm(n_comparisons, delta,
                    sqrt(se^2 + heterogeneity^2))

  comparisons <- data.frame(
    study_id = sprintf("sim_%03d", seq_len(n_comparisons)),
    target = drugs[target_idx],
    comparator = drugs[comparator_idx],
    estimate = unname(y),
    std_error = se
  )

  anchors <- NULL
  if (n_anchors > 0) {
    anchored <- sample(drugs, n_anchors)
    anchor_se <- stats::runif(n_anchors, anchor_se_range[1],
                              anchor_se_range[2])
    anchors <- data.frame(
      study_id = sprintf("sim_anchor_%02d", seq_len(n_anchors)),
      drug = anchored,
      reference = "placebo",
      estimate = stats::rnorm(n_anchors, theta[anchored], anchor_se),
      std_error = anchor_se
    )
  }

  network <- direct_effect_network(comparisons, anchors = anchors,
                                   effect_measure = "HR")

  truth <- data.frame(
    drug = c("placebo", drugs),
    theta = c(0, unname(theta))
  )

  list(
    network = network,
    comparisons = comparisons,
    anchors = anchors,
    truth = truth,
    tau = heterogeneity,
    seed = seed
  )
}

#' Measure how well a fit recovers simulation truth
#'
#' Compares a fit against the known true effects of the simulation that
#' generated its data. Works from the fit contract only, so it treats
#' every engine identically. For a surface fit the truth is re-centred
#' at the fit's arbitrary reference drug before comparison, and the
#' reference row itself (estimated as exactly 0 by construction) is
#' excluded from bias, RMSE, and coverage. For an anchored fit
#' (`reference = "placebo"`) the truth is already on the placebo = 0
#' scale, so it is compared directly and every drug contributes.
#'
#' @param fit A `directeffect_fit` from [fit_surface()], fitted to
#'   `simulation$network`.
#' @param simulation The result of [simulate_direct_effect_network()].
#'
#' @return A list with `bias` (mean error), `rmse`, `coverage`
#'   (proportion of intervals containing the true value),
#'   `rank_correlation` (Spearman correlation of estimated and true
#'   effect ordering), `n_drugs` (drugs contributing to the summaries),
#'   and `reference`.
#'
#' @examples
#' simulation <- simulate_direct_effect_network(
#'   n_drugs = 5, n_comparisons = 12, n_anchors = 0,
#'   heterogeneity = 0, seed = 1
#' )
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'   fit <- fit_surface(simulation$network, engine = "netmeta")
#'   validate_recovery(fit, simulation)
#' }
#' @export
validate_recovery <- function(fit, simulation) {
  assert_directeffect_fit(fit)
  if (!is.list(simulation) || is.null(simulation$truth)) {
    stop("`simulation` must be the result of ",
         "`simulate_direct_effect_network()` (it has no `truth`).",
         call. = FALSE)
  }

  effects <- fit$effects
  truth <- simulation$truth
  missing <- setdiff(effects$drug, truth$drug)
  if (length(missing) > 0) {
    stop("The fit contains drug",
         if (length(missing) > 1) "s" else "", " ",
         paste0("\"", missing, "\"", collapse = ", "),
         " that the simulation truth does not; was this fit produced ",
         "from `simulation$network`?", call. = FALSE)
  }

  reference <- effects$reference[1]
  theta_true <- truth$theta[match(effects$drug, truth$drug)]
  if (identical(reference, "placebo")) {
    # Anchored fit: truth is already absolute on the placebo = 0 scale,
    # and no drug is pinned, so every drug contributes.
    free <- rep(TRUE, nrow(effects))
  } else {
    theta_true <- theta_true - theta_true[effects$drug == reference]
    free <- effects$drug != reference
  }
  errors <- effects$estimate[free] - theta_true[free]
  covered <- effects$lower[free] <= theta_true[free] &
    theta_true[free] <= effects$upper[free]

  list(
    bias = mean(errors),
    rmse = sqrt(mean(errors^2)),
    coverage = mean(covered),
    rank_correlation = stats::cor(effects$estimate, theta_true,
                                  method = "spearman"),
    n_drugs = sum(free),
    reference = reference
  )
}
