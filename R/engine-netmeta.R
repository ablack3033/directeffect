# Frequentist surface engine: wraps netmeta. All knowledge of netmeta's
# interface and object layout lives in this file; nothing outside it may
# touch a netmeta object except through `engine_fit`.

fit_surface_netmeta <- function(de, reference) {
  if (!requireNamespace("netmeta", quietly = TRUE)) {
    stop("The \"netmeta\" engine requires the netmeta package. ",
         "Install it with install.packages(\"netmeta\").", call. = FALSE)
  }

  comparisons <- de$comparisons
  nm <- netmeta::netmeta(
    TE = comparisons$estimate,
    seTE = comparisons$std_error,
    treat1 = comparisons$target,
    treat2 = comparisons$comparator,
    studlab = comparisons$study_id,
    sm = de$effect_measure,
    common = TRUE,
    random = TRUE,
    reference.group = reference
  )

  drugs <- de$treatments
  effects <- data.frame(
    drug = drugs,
    estimate = nm$TE.common[drugs, reference],
    std_error = nm$seTE.common[drugs, reference],
    lower = nm$lower.common[drugs, reference],
    upper = nm$upper.common[drugs, reference],
    scale = "log",
    reference = reference,
    engine = "netmeta",
    row.names = NULL
  )

  heterogeneity <- list(
    Q = nm$Q,
    df = nm$df.Q,
    p_value = nm$pval.Q,
    tau = nm$tau,
    I2 = nm$I2
  )

  new_directeffect_fit(
    effects = effects,
    covariance = netmeta_surface_covariance(nm, drugs, reference),
    heterogeneity = heterogeneity,
    engine = "netmeta",
    engine_fit = nm,
    network = de
  )
}

# Common-effect covariance of every drug's estimate versus the
# reference, from netmeta's Moore-Penrose pseudoinverse of the weighted
# Laplacian: Cov(theta_i - theta_r, theta_j - theta_r) =
# L+_ij - L+_ir - L+_rj + L+_rr. Unlike Cov.common (indexed by
# comparison pair labels, with orientation-dependent signs), L+ is
# indexed by treatment, so no sign bookkeeping is needed, and it
# reflects netmeta's multi-arm variance adjustments.
netmeta_surface_covariance <- function(nm, drugs, reference) {
  lp <- nm$Lplus.matrix.common[drugs, drugs, drop = FALSE]
  covariance <- lp -
    outer(lp[, reference], rep(1, length(drugs))) -
    outer(rep(1, length(drugs)), lp[reference, ]) +
    lp[reference, reference]
  # The reference is pinned at exactly 0, so its row and column are
  # exact zeros by the contract, not merely zero to machine precision.
  covariance[reference, ] <- 0
  covariance[, reference] <- 0
  dimnames(covariance) <- list(drugs, drugs)
  covariance
}
