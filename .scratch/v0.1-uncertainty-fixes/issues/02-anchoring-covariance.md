# 02: Correct covariance propagation in frequentist anchoring

**What to build:** A researcher who anchors a netmeta surface gets
absolute standard errors and intervals equal to the actual sampling
variability of the estimator — the algebra of their evidence reflected
honestly. With a single anchor, the anchored drug's absolute SE is
exactly its anchor's SE (the surface contribution cancels); in general
the two-stage answer matches the joint GLS of comparisons plus anchor
rows, which is what the anchored Stan model already computes. The
circular test that protected the old formula is replaced by oracles
that can actually fail.

**Blocked by:** 01 (fit-contract covariance).

**Status:** done

- [x] Absolute variances include the offset variance and the covariance
      between each drug's surface position and the offset, computed
      from the surface covariance; with multiple anchors, the precision
      weights come from the anchors' full proposal covariance, not just
      its diagonal
- [x] Mandatory example: anchored estimates and SEs match an
      independent joint-GLS oracle (anchor rows appended to the
      weighted comparison design) to tight numerical tolerance; the
      anchored drug's SE equals its anchor's SE exactly, as a
      closed-form assertion
- [x] Anchored frequentist SEs agree with the anchored Stan posterior
      SDs on the mandatory example within a Monte Carlo–justified
      tolerance (cross-engine agreement restored for uncertainty)
- [x] The circular "multiple anchors combine by precision, checked by
      hand" test is replaced by the joint-GLS oracle test, including a
      disagreeing-anchors case
- [x] A Monte Carlo calibration test (a few hundred replicates of the
      fast frequentist path on a small network) shows empirical
      coverage of nominal 95% anchored intervals within bounds derived
      from binomial Monte Carlo error
- [x] Point estimates are unchanged: the mandatory deterministic
      recovery (A = 0.7, B = 0.7, C = 0.3, null A–B comparison not
      implying null direct effects) still passes on both engines, and
      unanchored components are still refused
- [x] The anchoring documentation's propagation promise describes the
      corrected algebra, and the v0.1-core anchoring issue's
      "acceptable for v0.1" addendum is annotated as superseded by the
      review's measurements

## Comments

Done. `anchor_surface_offset()` now computes the offset by GLS over the
anchor proposals with covariance `diag(a_se²) + Σ[anchored, anchored]`
(full matrix, not its diagonal), and every absolute variance as
`Var(θ_d) + Var(offset) + 2·Cov(θ_d, offset)` from the fit-contract
covariance — engine internals untouched. On the mandatory example the
anchored fit now equals the joint GLS of comparisons + anchor rows to
1e-16 (estimates, SEs, and full covariance), C's absolute SE is exactly
its anchor's 0.04, and the anchored netmeta SEs agree with the anchored
Stan posterior SDs within 0.005. The circular hand test was replaced by
the `joint_gls()` oracle (helper), including the disagreeing-anchors
case, where a profiling identity makes the offset and its variance
equal the joint GLS at the reference drug exactly — asserted at 1e-10 —
while non-reference SEs are honestly bounded below by the joint's
(measured: the location-shift estimator, kept per the spec, is slightly
less efficient than a joint fit when anchors disagree; docs now say
so). Monte Carlo: 400 replicates of the two-stage path (hand-WLS
surface — proven equal to netmeta's elsewhere in the suite — feeding
the real `anchor_surface()`), single- and two-anchor configurations:
per-drug coverage within 4 binomial-MC SDs of 0.95 and reported SEs
within 4 MC SDs of the empirical sampling SD. The v0.1-core issue 07
addendum is annotated superseded with the review's measurements.
