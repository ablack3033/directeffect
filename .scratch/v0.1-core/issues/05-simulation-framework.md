# 05: Simulation framework

**What to build:** `simulate_direct_effect_network(n_drugs, n_comparisons, n_anchors, heterogeneity, seed)` generates true effects, a comparison graph, standard errors, observed comparisons, and anchor estimates — retaining `truth` — and `validate_recovery(fit, simulation)` reports bias, RMSE, interval coverage, and rank correlation. Demoable end to end: simulate a low-noise network, fit it with netmeta, see near-nominal recovery.

**Blocked by:** 03 (netmeta surface engine)

**Status:** ready-for-agent

- [ ] Simulation returns network-ready comparisons and anchors tables plus retained truth; fully reproducible under a seed
- [ ] Generated comparison graphs are connected by construction (or the generator reports components explicitly)
- [ ] `validate_recovery()` works from the fit contract only — engine-agnostic by design
- [ ] Recovery test in the suite: a low-noise simulation shows near-zero bias and near-nominal interval coverage of relative effects
