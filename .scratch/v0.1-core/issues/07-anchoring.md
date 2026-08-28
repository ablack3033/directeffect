# 07: Anchor the surface (sea level), both engines

**What to build:** `anchor_surface(surface_fit, anchors)` positions a fitted relative surface against placebo = 0. The Bayesian path refits with the separate anchored Stan model — no arbitrary constraint; the anchors identify the absolute location while retaining their uncertainty. The frequentist path estimates the single precision-weighted location offset, propagating both surface and anchor uncertainty. Components without anchors are never silently anchored.

**Blocked by:** 03 (netmeta engine), 04 (Stan engine)

**Status:** done

- [x] The anchored Stan program is a separate model from the unanchored one — two simple deep modules, not one model full of conditionals
- [x] An anchored fit returns the same fit contract with `reference = "placebo"` and effects on the absolute scale
- [x] Mandatory deterministic test passes on BOTH engines: comparisons A−B = 0.0, A−C = 0.4, B−C = 0.4 plus a single anchor C = 0.3 recover A ≈ 0.7, B ≈ 0.7, C ≈ 0.3 — protecting the central property that a null active-comparator estimate does not imply a null direct effect
- [x] Anchoring a fit whose component has no anchor fails with an informative error — the package never silently picks a sea level
- [x] Anchor uncertainty propagates: inflating an anchor's standard error visibly widens the absolute intervals in a test

## Comments

Done. `anchor_surface(fit, anchors = NULL, ...)` in R/anchor-surface.R;
anchors default to the network's, an explicit table re-validates
against the network. Stan path (engine-stan.R) refits with
inst/stan/anchored_surface.stan — a separate program, vector[K] theta
with no constraint, anchors in the likelihood — so every drug is
sampled with real convergence diagnostics. Frequentist path
(anchor_surface_offset) implements the spec's precision-weighted
location offset: each anchor proposes offset a_m − surface_m with
variance a_se² + surface_se²; absolute se adds the offset variance to
the surface se. On the mandatory example netmeta recovers 0.7/0.7/0.3
to 1e-6 and Stan to 0.03; tests also assert the null A−B comparison
still yields A, B > 0.5. Bug found on the way: rstan drops the array
dimension of length-1 data (a single anchor), fixed with as.array() in
both adapters. Anchored-ness is signalled by reference = "placebo" —
the contract's component list is unchanged.

Code-review addendum (spec axis): the frequentist offset follows the
spec's wording literally — precision-weighted, adding offset variance
to surface variance — but ignores the covariance between the surface
estimate and the offset (both involve the anchored drug's surface
position). With a single anchor this double-counts that drug's surface
variance, making absolute intervals conservative rather than
miscalibrated-narrow. Acceptable for v0.1; revisit with the v0.2
robustness work if calibrated frequentist intervals matter.
