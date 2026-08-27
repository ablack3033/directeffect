# 10: Effect-surface and sea-level plots

**What to build:** `plot_effect_surface(fit)` renders the one-dimensional direct-effect surface (harmful → beneficial) with uncertainty intervals, selectable between `scale = "log"` and `scale = "natural"`; plus the signature surface + sea-level figure showing the relative structure positioned against the placebo = 0 line for anchored fits.

**Blocked by:** 07 (Anchoring)

**Status:** ready-for-agent

- [ ] Strictly one-dimensional rendering for a single outcome — no 2-D layout pretending to represent causal distance
- [ ] Both scales correct for multiplicative measures (log 0 ↔ HR = 1), with axis labels naming the effect measure
- [ ] Unanchored fits plot relative positions clearly labeled as relative to an arbitrary reference; anchored fits show the placebo = 0 sea-level line
- [ ] Engine-agnostic: identical calls work on netmeta and Stan fits
