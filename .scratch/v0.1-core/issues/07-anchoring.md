# 07: Anchor the surface (sea level), both engines

**What to build:** `anchor_surface(surface_fit, anchors)` positions a fitted relative surface against placebo = 0. The Bayesian path refits with the separate anchored Stan model — no arbitrary constraint; the anchors identify the absolute location while retaining their uncertainty. The frequentist path estimates the single precision-weighted location offset, propagating both surface and anchor uncertainty. Components without anchors are never silently anchored.

**Blocked by:** 03 (netmeta engine), 04 (Stan engine)

**Status:** ready-for-agent

- [ ] The anchored Stan program is a separate model from the unanchored one — two simple deep modules, not one model full of conditionals
- [ ] An anchored fit returns the same fit contract with `reference = "placebo"` and effects on the absolute scale
- [ ] Mandatory deterministic test passes on BOTH engines: comparisons A−B = 0.0, A−C = 0.4, B−C = 0.4 plus a single anchor C = 0.3 recover A ≈ 0.7, B ≈ 0.7, C ≈ 0.3 — protecting the central property that a null active-comparator estimate does not imply a null direct effect
- [ ] Anchoring a fit whose component has no anchor fails with an informative error — the package never silently picks a sea level
- [ ] Anchor uncertainty propagates: inflating an anchor's standard error visibly widens the absolute intervals in a test
