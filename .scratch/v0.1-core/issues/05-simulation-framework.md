# 05: Simulation framework

**What to build:** `simulate_direct_effect_network(n_drugs, n_comparisons, n_anchors, heterogeneity, seed)` generates true effects, a comparison graph, standard errors, observed comparisons, and anchor estimates — retaining `truth` — and `validate_recovery(fit, simulation)` reports bias, RMSE, interval coverage, and rank correlation. Demoable end to end: simulate a low-noise network, fit it with netmeta, see near-nominal recovery.

**Blocked by:** 03 (netmeta surface engine)

**Status:** done

- [x] Simulation returns network-ready comparisons and anchors tables plus retained truth; fully reproducible under a seed
- [x] Generated comparison graphs are connected by construction (or the generator reports components explicitly)
- [x] `validate_recovery()` works from the fit contract only — engine-agnostic by design
- [x] Recovery test in the suite: a low-noise simulation shows near-zero bias and near-nominal interval coverage of relative effects

## Comments

Done in R/simulate.R. Connectivity by construction via a random
spanning tree (drug k compares against a random earlier drug), then the
remaining comparisons between random pairs, repeats allowed. Generative
model matches the fitting model: y ~ N(delta, se² + tau²), anchors
a ~ N(theta, se²) vs placebo = 0; truth retains a placebo row at 0.
Seeded runs restore the caller's RNG state afterwards (tested). Beyond
the design signature, `effect_sd`, `se_range`, and `anchor_se_range`
are exposed so tests can dial noise. validate_recovery() re-centres
truth at the fit's reference, excludes the pinned reference row from
bias/RMSE/coverage, and returns a plain list (spec keeps the class
system minimal — no directeffect_validation until later versions).
Recovery test: 40 drugs, 250 low-noise comparisons → |bias| < 0.02,
RMSE < 0.05, coverage ≥ 0.85, rank correlation > 0.95.
