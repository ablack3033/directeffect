# 10: Effect-surface and sea-level plots

**What to build:** `plot_effect_surface(fit)` renders the one-dimensional direct-effect surface (harmful → beneficial) with uncertainty intervals, selectable between `scale = "log"` and `scale = "natural"`; plus the signature surface + sea-level figure showing the relative structure positioned against the placebo = 0 line for anchored fits.

**Blocked by:** 07 (Anchoring)

**Status:** done

- [x] Strictly one-dimensional rendering for a single outcome — no 2-D layout pretending to represent causal distance
- [x] Both scales correct for multiplicative measures (log 0 ↔ HR = 1), with axis labels naming the effect measure
- [x] Unanchored fits plot relative positions clearly labeled as relative to an arbitrary reference; anchored fits show the placebo = 0 sea-level line
- [x] Engine-agnostic: identical calls work on netmeta and Stan fits

## Comments

Done in R/plot-surface.R, two functions on the shared
surface_plot_data()/surface_plot_labels() helpers. plot_effect_surface():
effect on x (harmful → beneficial by estimate order), drugs stacked on
a categorical axis — presentation only, the surface stays 1-D in the
effect; unanchored fits get a dashed reference line and a subtitle
"positions vs arbitrary reference X = 0 (not absolute effects)";
anchored fits get the solid placebo sea-level line. plot_sea_level()
is the signature figure — effect on y against the horizontal
placebo = 0 line with an inline label — and refuses relative fits ("a
relative surface has no sea level"). scale = "natural" back-transforms
via exp() on a log-spaced axis so log 0 ↔ HR = 1 (tested exactly);
axis labels name the effect measure. Both engines render both plots in
tests, via ggplot_build. Gotcha fixed: an annotate() with numeric x
placed before the point layer types the discrete axis as continuous —
annotation goes last.
