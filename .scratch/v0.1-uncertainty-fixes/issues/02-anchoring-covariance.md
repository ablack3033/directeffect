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

**Status:** ready-for-agent

- [ ] Absolute variances include the offset variance and the covariance
      between each drug's surface position and the offset, computed
      from the surface covariance; with multiple anchors, the precision
      weights come from the anchors' full proposal covariance, not just
      its diagonal
- [ ] Mandatory example: anchored estimates and SEs match an
      independent joint-GLS oracle (anchor rows appended to the
      weighted comparison design) to tight numerical tolerance; the
      anchored drug's SE equals its anchor's SE exactly, as a
      closed-form assertion
- [ ] Anchored frequentist SEs agree with the anchored Stan posterior
      SDs on the mandatory example within a Monte Carlo–justified
      tolerance (cross-engine agreement restored for uncertainty)
- [ ] The circular "multiple anchors combine by precision, checked by
      hand" test is replaced by the joint-GLS oracle test, including a
      disagreeing-anchors case
- [ ] A Monte Carlo calibration test (a few hundred replicates of the
      fast frequentist path on a small network) shows empirical
      coverage of nominal 95% anchored intervals within bounds derived
      from binomial Monte Carlo error
- [ ] Point estimates are unchanged: the mandatory deterministic
      recovery (A = 0.7, B = 0.7, C = 0.3, null A–B comparison not
      implying null direct effects) still passes on both engines, and
      unanchored components are still refused
- [ ] The anchoring documentation's propagation promise describes the
      corrected algebra, and the v0.1-core anchoring issue's
      "acceptable for v0.1" addendum is annotated as superseded by the
      review's measurements
