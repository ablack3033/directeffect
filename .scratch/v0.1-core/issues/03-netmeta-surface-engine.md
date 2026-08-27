# 03: netmeta surface engine

**What to build:** `fit_surface(de, engine = "netmeta")` fits the relative direct-effect surface with netmeta behind an adapter and returns a `directeffect_fit` whose effects table has the stable tidy schema (`drug`, `estimate`, `std_error`, `lower`, `upper`, `scale`, `reference`, `engine`). Anchors are deliberately ignored; identification uses an arbitrary reference at 0. This is the tracer bullet: comparisons table → network → fit → tidy effects, verifiable by hand.

**Blocked by:** 02 (Direct-effect network object)

**Status:** ready-for-agent

- [ ] Column mapping into netmeta per the design (estimate→TE, std_error→seTE, target→treat1, comparator→treat2, study_id→studlab); common-effect estimates are the v0.1 result
- [ ] Returns a `directeffect_fit` with the full component set: effects, comparisons, anchors, heterogeneity, diagnostics, engine, engine_fit, network
- [ ] netmeta internals never leak — downstream code works from the fit contract only; `engine_fit` is the sole escape hatch
- [ ] On the spec's 3-drug example the relative surface matches the hand computation to tight numerical tolerance
- [ ] Fitting a multi-component network fails with an informative error pointing at `check_connectivity()` — never a silent wrong answer
