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
                             seed = NULL, refresh = 0, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The \"stan\" engine requires the rstan package. ",
         "Install it with install.packages(\"rstan\").", call. = FALSE)
  }
  refuse_multiarm_stan(de$comparisons)

  # The Stan program fixes theta[1] = 0, so order drugs with the chosen
  # reference first; effects are reported back in network order.
  drugs <- c(reference, setdiff(de$treatments, reference))
  index <- stats::setNames(seq_along(drugs), drugs)
  stan_data <- stan_comparison_data(de$comparisons, index)

  sf <- sampling_preserving_rng(
    compiled_stan_model("surface"),
    data = stan_data,
    chains = chains,
    iter = iter,
    seed = stan_sampling_seed(seed),
    refresh = refresh,
    ...
  )

  effects <- stan_effects_table(sf, drugs, reference)
  covariance <- stan_theta_covariance(sf, drugs)

  # The reference is a constant, not a sampled quantity: pin it at exactly
  # 0 and mark its convergence diagnostics as not applicable.
  ref_row <- effects$drug == reference
  effects[ref_row, c("estimate", "std_error", "lower", "upper",
                     "median", "mean", "sd", "q025", "q975")] <- 0
  effects[ref_row, c("rhat", "ess_bulk", "ess_tail")] <- NA_real_
  covariance[reference, ] <- 0
  covariance[, reference] <- 0

  effects <- effects[match(de$treatments, effects$drug), ]
  rownames(effects) <- NULL
  covariance <- covariance[de$treatments, de$treatments]

  warn_on_convergence(effects)

  new_directeffect_fit(
    effects = effects,
    covariance = covariance,
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
                                seed = NULL, refresh = 0, ...) {
  de <- fit$network
  refuse_multiarm_stan(de$comparisons)
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

  sf <- sampling_preserving_rng(
    compiled_stan_model("anchored_surface"),
    data = stan_data,
    chains = chains,
    iter = iter,
    seed = stan_sampling_seed(seed),
    refresh = refresh,
    ...
  )

  effects <- stan_effects_table(sf, drugs, reference = "placebo")
  warn_on_convergence(effects)

  anchored <- new_directeffect_fit(
    effects = effects,
    covariance = stan_theta_covariance(sf, drugs),
    heterogeneity = list(model = "common", tau = 0),
    engine = "stan",
    engine_fit = sf,
    network = de
  )
  anchored$anchors <- anchors
  anchored
}

# The Stan likelihood treats every comparison row as independent, which
# is wrong for several rows from one study (a multi-arm trial): their
# shared arms correlate the estimates, and ignoring that overstates
# precision. Fail loudly rather than silently approximate.
refuse_multiarm_stan <- function(comparisons) {
  studies <- multiarm_studies(comparisons)
  if (length(studies) > 0) {
    stop("The \"stan\" engine cannot fit this network: ",
         multiarm_clause(studies),
         " more than one comparison (a multi-arm trial), and the Stan ",
         "likelihood treats comparison rows as independent, which ",
         "overstates precision for correlated within-study contrasts. ",
         "Use `engine = \"netmeta\"`, which models the within-study ",
         "covariance of multi-arm trials correctly.", call. = FALSE)
  }
  invisible(comparisons)
}

# When the caller supplies no seed, seed the sampler from the wall
# clock and process id instead of drawing from R's global RNG stream:
# fitting a model must never consume or disturb the caller's
# random-number state (the standard the simulation module already
# meets).
stan_sampling_seed <- function(seed) {
  if (!is.null(seed)) {
    return(seed)
  }
  as.integer((as.numeric(Sys.time()) * 1000 + Sys.getpid()) %%
               (.Machine$integer.max - 1)) + 1L
}

# rstan::sampling() modifies .Random.seed even when given an explicit
# seed; snapshot and restore the caller's state around the call (and
# remove any state the sampler created where none existed) so fitting
# is invisible to the caller's RNG stream.
sampling_preserving_rng <- function(...) {
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())))
  }
  rstan::sampling(...)
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

# Posterior covariance of theta across all retained draws, labelled by
# drug in the Stan program's drug order (the caller reorders).
stan_theta_covariance <- function(sf, drugs) {
  sims <- rstan::extract(sf, pars = "theta", permuted = FALSE,
                         inc_warmup = FALSE)
  draws <- apply(sims, 3, c)
  covariance <- stats::cov(draws)
  dimnames(covariance) <- list(drugs, drugs)
  covariance
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
