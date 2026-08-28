# Stan-side demonstrations for the stats review of directeffect.
# 1. Anchored SEs: the package's own Bayesian engine vs the frequentist
#    offset path on the mandatory example.
# 2. Multi-arm study: netmeta vs Stan on a single 3-arm trial — the
#    "identical likelihood" claim.
# 3. Global RNG state perturbed by fit_surface(engine = "stan") defaults.

# Run from the repository root.
for (f in list.files("R", full.names = TRUE)) source(f)

# Compile from the source tree instead of an installed package.
compiled_stan_model <- function(name) {
  if (is.null(stan_models[[name]])) {
    stan_models[[name]] <- rstan::stan_model(
      file = file.path("inst/stan", paste0(name, ".stan")))
  }
  stan_models[[name]]
}

options(mc.cores = 1)

cat("=== Demo 1: anchored SEs, Stan engine, mandatory example ===\n")
comparisons <- data.frame(
  study_id = c("S1", "S2", "S3"),
  target = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate = c(0.0, 0.4, 0.4),
  std_error = c(0.05, 0.05, 0.05))
anchors <- data.frame(
  study_id = "RCT1", drug = "C", reference = "placebo",
  estimate = 0.3, std_error = 0.04)
de <- direct_effect_network(comparisons, anchors = anchors,
                            effect_measure = "HR")

surface_st <- fit_surface(de, engine = "stan", seed = 20260828,
                          iter = 8000)
absolute_st <- anchor_surface(surface_st, seed = 20260828, iter = 8000)
cat("Stan anchored effects (posterior sd = correct SE):\n")
print(absolute_st$effects[, c("drug", "estimate", "std_error")],
      digits = 4)

surface_nm <- fit_surface(de, engine = "netmeta")
absolute_nm <- anchor_surface(surface_nm)
cat("\nnetmeta-path anchored effects (offset formula):\n")
print(absolute_nm$effects[, c("drug", "estimate", "std_error")],
      digits = 4)

cat("\n=== Demo 2: one 3-arm trial, netmeta vs Stan ===\n")
# Single study S1 with arms A, B, C: consistent contrasts
# (AB + BC = AC), equal contrast SEs 0.1 (arm variance 0.005 each).
multi <- data.frame(
  study_id = "S1",
  target = c("B", "C", "C"),
  comparator = c("A", "A", "B"),
  estimate = c(0.1, 0.3, 0.2),
  std_error = c(0.1, 0.1, 0.1))
de_m <- direct_effect_network(multi, effect_measure = "HR")

fit_m_nm <- fit_surface(de_m, engine = "netmeta")
fit_m_st <- fit_surface(de_m, engine = "stan", seed = 1, iter = 8000)

cat("netmeta (multi-arm aware) effects:\n")
print(fit_m_nm$effects[, c("drug", "estimate", "std_error")], digits = 4)
cat("Stan (rows independent) effects:\n")
print(fit_m_st$effects[, c("drug", "estimate", "std_error")], digits = 4)
cat("compare_engines() difference (what the equivalence test checks):\n")
print(compare_engines(fit_m_nm, fit_m_st), digits = 4)
cat("SE ratio stan/netmeta per drug (never checked anywhere):\n")
print(fit_m_st$effects$std_error / fit_m_nm$effects$std_error, digits = 4)
cat("sqrt(2/3) for reference:", sqrt(2 / 3), "\n")

cat("\n=== Demo 3: global RNG state after unseeded stan fit ===\n")
set.seed(7); expected <- rnorm(1)
set.seed(7)
invisible(fit_surface(de, engine = "stan", chains = 1, iter = 400))
after <- rnorm(1)
cat("rnorm(1) after set.seed(7):                 ", expected, "\n")
cat("rnorm(1) after set.seed(7) + unseeded fit:  ", after, "\n")
cat("global RNG stream disturbed:", !identical(expected, after), "\n")
