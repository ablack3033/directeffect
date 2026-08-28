# Bayesian surface engine: wraps rstan. All knowledge of Stan's interface
# and object layout lives in this file; nothing outside it may touch a
# stanfit object except through `engine_fit`. rstan was chosen over
# cmdstanr because it is installable as a self-contained package (no
# CmdStan toolchain download at setup time), which keeps CI and offline
# installation simple.

# Compiled Stan models, cached per session so repeated fits pay the
# compilation cost once.
stan_models <- new.env(parent = emptyenv())

compiled_stan_model <- function(name) {
  if (is.null(stan_models[[name]])) {
    path <- system.file("stan", paste0(name, ".stan"),
                        package = "directeffect")
    stan_models[[name]] <- rstan::stan_model(file = path)
  }
  stan_models[[name]]
}

fit_surface_stan <- function(de, reference, chains = 4, iter = 2000,
                             seed = sample.int(.Machine$integer.max, 1),
                             refresh = 0, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The \"stan\" engine requires the rstan package. ",
         "Install it with install.packages(\"rstan\").", call. = FALSE)
  }

  # The Stan program fixes theta[1] = 0, so order drugs with the chosen
  # reference first; effects are reported back in network order.
  drugs <- c(reference, setdiff(de$treatments, reference))
  index <- stats::setNames(seq_along(drugs), drugs)

  comparisons <- de$comparisons
  stan_data <- list(
    N = nrow(comparisons),
    K = length(drugs),
    target_idx = unname(index[comparisons$target]),
    comparator_idx = unname(index[comparisons$comparator]),
    y = comparisons$estimate,
    se = comparisons$std_error
  )

  sf <- rstan::sampling(
    compiled_stan_model("surface"),
    data = stan_data,
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = refresh,
    ...
  )

  sims <- rstan::extract(sf, pars = "theta", permuted = FALSE,
                         inc_warmup = FALSE)
  posterior <- as.data.frame(rstan::monitor(sims, warmup = 0,
                                            print = FALSE))

  effects <- data.frame(
    drug = drugs,
    estimate = posterior$mean,
    std_error = posterior$sd,
    lower = posterior$`2.5%`,
    upper = posterior$`97.5%`,
    scale = "log",
    reference = reference,
    engine = "stan",
    median = posterior$`50%`,
    mean = posterior$mean,
    sd = posterior$sd,
    q025 = posterior$`2.5%`,
    q975 = posterior$`97.5%`,
    rhat = posterior$Rhat,
    ess_bulk = posterior$Bulk_ESS,
    ess_tail = posterior$Tail_ESS,
    row.names = NULL
  )

  # The reference is a constant, not a sampled quantity: pin it at exactly
  # 0 and mark its convergence diagnostics as not applicable.
  ref_row <- effects$drug == reference
  effects[ref_row, c("estimate", "std_error", "lower", "upper",
                     "median", "mean", "sd", "q025", "q975")] <- 0
  effects[ref_row, c("rhat", "ess_bulk", "ess_tail")] <- NA_real_

  effects <- effects[match(de$treatments, effects$drug), ]
  rownames(effects) <- NULL

  warn_on_convergence(effects)

  new_directeffect_fit(
    effects = effects,
    heterogeneity = list(model = "common", tau = 0),
    engine = "stan",
    engine_fit = sf,
    network = de
  )
}

warn_on_convergence <- function(effects, rhat_max = 1.01, ess_min = 400) {
  sampled <- effects[!is.na(effects$rhat), ]
  bad <- sampled$rhat > rhat_max | sampled$ess_bulk < ess_min |
    sampled$ess_tail < ess_min
  if (any(bad)) {
    warning("MCMC convergence diagnostics indicate trouble for drug",
            if (sum(bad) > 1) "s" else "", " ",
            paste0("\"", sampled$drug[bad], "\"", collapse = ", "),
            " (max rhat ", round(max(sampled$rhat), 3),
            ", min ess ", round(min(pmin(sampled$ess_bulk,
                                         sampled$ess_tail))),
            "). Increase `iter` or inspect `engine_fit` before trusting ",
            "these estimates.", call. = FALSE)
  }
  invisible(effects)
}
