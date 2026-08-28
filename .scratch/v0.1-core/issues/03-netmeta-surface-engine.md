# 03: netmeta surface engine

**What to build:** `fit_surface(de, engine = "netmeta")` fits the relative direct-effect surface with netmeta behind an adapter and returns a `directeffect_fit` whose effects table has the stable tidy schema (`drug`, `estimate`, `std_error`, `lower`, `upper`, `scale`, `reference`, `engine`). Anchors are deliberately ignored; identification uses an arbitrary reference at 0. This is the tracer bullet: comparisons table → network → fit → tidy effects, verifiable by hand.

**Blocked by:** 02 (Direct-effect network object)

**Status:** done

- [x] Column mapping into netmeta per the design (estimate→TE, std_error→seTE, target→treat1, comparator→treat2, study_id→studlab); common-effect estimates are the v0.1 result
- [x] Returns a `directeffect_fit` with the full component set: effects, comparisons, anchors, heterogeneity, diagnostics, engine, engine_fit, network
- [x] netmeta internals never leak — downstream code works from the fit contract only; `engine_fit` is the sole escape hatch
- [x] On the spec's 3-drug example the relative surface matches the hand computation to tight numerical tolerance
- [x] Fitting a multi-component network fails with an informative error pointing at `check_connectivity()` — never a silent wrong answer

## Comments

Done. `fit_surface()` in R/fit-surface.R owns validation, engine
dispatch, the `reference` argument (default: first treatment), and the
`new_directeffect_fit()` contract constructor; the netmeta adapter in
R/engine-netmeta.R is the only file that touches netmeta objects
(common and random both requested; common used). The hand oracle in
tests solves the WLS normal equations directly and matches netmeta to
1e-6; reference invariance and anchors-are-ignored are also tested.
netmeta added to Suggests (engine deps stay optional). Local netmeta
3.7.0 was installed from the guido-s GitHub sources because the
sandbox's network policy blocks CRAN and apt has no netmeta for
R 4.3; CI installs it from CRAN normally.
