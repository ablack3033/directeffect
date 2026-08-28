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
    heterogeneity = heterogeneity,
    engine = "netmeta",
    engine_fit = nm,
    network = de
  )
}
