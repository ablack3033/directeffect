# 03: Engine equivalence extends to uncertainty

**What to build:** The package's central cross-validation claim —
two independent engines reconstructing the same surface — holds for
uncertainty, continuously: CI fails if the engines' standard errors
drift apart, on surface fits and anchored fits alike, the way it has
always failed on estimate drift.

**Blocked by:** 02 (anchoring covariance) — until the frequentist
anchored SEs are corrected, an SE-level equivalence test correctly
fails.

**Status:** done

- [x] The engine-equivalence test additionally bounds the per-drug
      difference between frequentist standard errors and Stan posterior
      SDs on a simulated surface fit, with the tolerance justified
      against Monte Carlo error at the test's sampler settings (as the
      existing estimate EPSILON already is)
- [x] The same SE-level bound is asserted for anchored fits
- [x] The comparison-table documentation notes that
      `standardized_difference` is a yardstick, not a test statistic,
      because the two fits share the same data (no new statistic is
      introduced)
- [x] The equivalence tests keep running in CI on every push

## Comments

Done. The surface equivalence test now also bounds per-drug
|netmeta SE − Stan posterior SD| by SE_EPSILON = 0.01, justified in the
test header against the posterior-SD MCSE at iter = 4000
(sd/sqrt(2·ESS) ≈ 0.002) the same way EPSILON = 0.02 is for estimates;
a new anchored equivalence test asserts both bounds on anchored fits.
The anchored test uses a single anchor — the case where the frequentist
location shift and the anchored Stan model are mathematically the same
fit, so any drift is a real divergence; with several disagreeing
anchors the paths legitimately differ slightly (see 02), which the
anchoring oracle tests cover. `compare_engines()` docs and the formats
page now state that `standardized_difference` is a yardstick, not a
test statistic (the fits share data). Both tests run in the normal
suite on every push.
