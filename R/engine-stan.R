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
  stan_data <- stan_comparison_data(de$comparisons, index)

  sf <- rstan::sampling(
    compiled_stan_model("surface"),
    data = stan_data,
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = refresh,
    ...
  )

  effects <- stan_effects_table(sf, drugs, reference)

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

# Bayesian sea level: refit with the anchored Stan model, in which the
# anchors (not an arbitrary constraint) identify the absolute location
# and their uncertainty flows into every drug's posterior.
anchor_surface_stan <- function(fit, anchors, chains = 4, iter = 2000,
                                seed = sample.int(.Machine$integer.max, 1),
                                refresh = 0, ...) {
  de <- fit$network
  drugs <- de$treatments
  index <- stats::setNames(seq_along(drugs), drugs)
  stan_data <- c(
    stan_comparison_data(de$comparisons, index),
    list(
      M = nrow(anchors),
      anchor_idx = as.array(unname(index[anchors$drug])),
      a = as.array(anchors$estimate),
      a_se = as.array(anchors$std_error)
    )
  )

  sf <- rstan::sampling(
    compiled_stan_model("anchored_surface"),
    data = stan_data,
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = refresh,
    ...
  )

  effects <- stan_effects_table(sf, drugs, reference = "placebo")
  warn_on_convergence(effects)

  anchored <- new_directeffect_fit(
    effects = effects,
    heterogeneity = list(model = "common", tau = 0),
    engine = "stan",
    engine_fit = sf,
    network = de
  )
  anchored$anchors <- anchors
  anchored
}

# The comparison block of data shared by both Stan programs. as.array
# keeps length-1 inputs dimensioned as Stan arrays.
stan_comparison_data <- function(comparisons, index) {
  list(
    N = nrow(comparisons),
    K = length(index),
    target_idx = as.array(unname(index[comparisons$target])),
    comparator_idx = as.array(unname(index[comparisons$comparator])),
    y = as.array(comparisons$estimate),
    se = as.array(comparisons$std_error)
  )
}

stan_effects_table <- function(sf, drugs, reference) {
  sims <- rstan::extract(sf, pars = "theta", permuted = FALSE,
                         inc_warmup = FALSE)
  posterior <- as.data.frame(rstan::monitor(sims, warmup = 0,
                                            print = FALSE))
  data.frame(
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
