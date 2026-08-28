# 08: Surface diagnostics — edge residuals and cycle consistency

**What to build:** `edge_residuals(fit)` reports observed, predicted, residual, and standardized residual for every comparison; `cycle_consistency(fit)` sums observed effects around a cycle basis of the network (never full cycle enumeration), with `plot_cycle_consistency()`. Together these answer "are my comparative estimates internally coherent?" on any fitted surface.

**Blocked by:** 03 (netmeta surface engine)

**Status:** done

- [x] Edge residual output schema: `target`, `comparator`, `observed`, `predicted`, `residual`, `standardized_residual`
- [x] Cycle sums computed over a cycle basis; a consistent network (hand-built or simulated) yields sums ≈ 0
- [x] An injected inconsistent edge is flagged by both diagnostics in tests
- [x] Both diagnostics work identically on netmeta and Stan fits, consuming only the fit contract

## Comments

Done in R/diagnostics.R. edge_residuals() predicts each comparison
from fit$effects and standardizes by the comparison's own std_error;
rows align with fit$comparisons. cycle_consistency() first pools
parallel comparisons of the same pair by precision (orienting to the
alphabetically first drug), then builds a fundamental cycle basis from
a spanning tree of the pooled simple graph — one basis cycle per
non-tree edge, never full enumeration; empty (zero-row, right-schema)
result for cycle-free networks. Cycle rows: cycle walk string, n_edges,
signed sum, std_error from pooled variances, z. Injected-inconsistency
test flips B−C from 0.4 to 0.9 and both diagnostics flag |z| > 1.96;
cycle_consistency agrees exactly across engines since it uses only
observed data. plot_cycle_consistency() (ggplot2, in Suggests) dots
z per cycle with ±1.96 guides; rlang enters Imports for the .data
pronoun (a top-level utils::globalVariables call trips a check NOTE).
