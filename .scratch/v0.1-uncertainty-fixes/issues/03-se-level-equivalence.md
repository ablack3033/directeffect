# 03: Engine equivalence extends to uncertainty

**What to build:** The package's central cross-validation claim —
two independent engines reconstructing the same surface — holds for
uncertainty, continuously: CI fails if the engines' standard errors
drift apart, on surface fits and anchored fits alike, the way it has
always failed on estimate drift.

**Blocked by:** 02 (anchoring covariance) — until the frequentist
anchored SEs are corrected, an SE-level equivalence test correctly
fails.

**Status:** ready-for-agent

- [ ] The engine-equivalence test additionally bounds the per-drug
      difference between frequentist standard errors and Stan posterior
      SDs on a simulated surface fit, with the tolerance justified
      against Monte Carlo error at the test's sampler settings (as the
      existing estimate EPSILON already is)
- [ ] The same SE-level bound is asserted for anchored fits
- [ ] The comparison-table documentation notes that
      `standardized_difference` is a yardstick, not a test statistic,
      because the two fits share the same data (no new statistic is
      introduced)
- [ ] The equivalence tests keep running in CI on every push
