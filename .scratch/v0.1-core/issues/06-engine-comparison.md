# 06: Engine comparison and equivalence CI test

**What to build:** `compare_engines(netmeta_fit, stan_fit)` returns per-drug estimates from both engines with differences and standardized differences; `plot_engine_comparison()` draws netmeta vs. Stan with the identity line. A permanent CI test simulates a network (Gaussian estimates, known standard errors, no heterogeneity, weak priors) and requires every |theta_stan − theta_netmeta| < ε. This makes the v0.1 goal — the same surface reconstructed independently by two engines — continuously enforced.

**Blocked by:** 03 (netmeta engine), 04 (Stan engine), 05 (Simulation framework)

**Status:** ready-for-agent

- [ ] Comparison table schema: `drug`, `netmeta`, `stan_mean`, `difference`, `standardized_difference`
- [ ] Plot includes the identity line y = x and works on any pair of fits from the two engines
- [ ] Reference and parameterization reconciled so the comparison is apples-to-apples
- [ ] Equivalence test runs in CI on every push and passes with a documented ε
