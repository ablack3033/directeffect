# 06: Engine comparison and equivalence CI test

**What to build:** `compare_engines(netmeta_fit, stan_fit)` returns per-drug estimates from both engines with differences and standardized differences; `plot_engine_comparison()` draws netmeta vs. Stan with the identity line. A permanent CI test simulates a network (Gaussian estimates, known standard errors, no heterogeneity, weak priors) and requires every |theta_stan − theta_netmeta| < ε. This makes the v0.1 goal — the same surface reconstructed independently by two engines — continuously enforced.

**Blocked by:** 03 (netmeta engine), 04 (Stan engine), 05 (Simulation framework)

**Status:** done

- [x] Comparison table schema: `drug`, `netmeta`, `stan_mean`, `difference`, `standardized_difference`
- [x] Plot includes the identity line y = x and works on any pair of fits from the two engines
- [x] Reference and parameterization reconciled so the comparison is apples-to-apples
- [x] Equivalence test runs in CI on every push and passes with a documented ε

## Comments

Done in R/compare-engines.R. compare_engines() accepts the two fits in
either order (detects engines, requires exactly one of each), requires
the same drug set and the same reference — differing references error
with instructions rather than being silently recentred, which keeps the
comparison apples-to-apples; the reference row's standardized
difference is NA (both engines report it as exact 0).
plot_engine_comparison() draws netmeta vs Stan posterior mean with the
y = x identity line and fixed aspect. The equivalence test
(test-engine-equivalence.R) simulates 8 drugs / 24 comparisons under
the model's own assumptions (no heterogeneity, known SEs), fits both
engines (Stan iter = 4000, fixed seed), and requires every
|theta_stan − theta_netmeta| < ε = 0.02 on the log scale — documented
in the test header: MCSE ~0.002 and prior shrinkage ~2e-4, so 0.02 has
an order of magnitude of slack yet fails on real model discrepancies.
Runs in CI on every push via the existing R-CMD-check workflow.
