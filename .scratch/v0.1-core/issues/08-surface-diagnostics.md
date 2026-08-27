# 08: Surface diagnostics — edge residuals and cycle consistency

**What to build:** `edge_residuals(fit)` reports observed, predicted, residual, and standardized residual for every comparison; `cycle_consistency(fit)` sums observed effects around a cycle basis of the network (never full cycle enumeration), with `plot_cycle_consistency()`. Together these answer "are my comparative estimates internally coherent?" on any fitted surface.

**Blocked by:** 03 (netmeta surface engine)

**Status:** ready-for-agent

- [ ] Edge residual output schema: `target`, `comparator`, `observed`, `predicted`, `residual`, `standardized_residual`
- [ ] Cycle sums computed over a cycle basis; a consistent network (hand-built or simulated) yields sums ≈ 0
- [ ] An injected inconsistent edge is flagged by both diagnostics in tests
- [ ] Both diagnostics work identically on netmeta and Stan fits, consuming only the fit contract
