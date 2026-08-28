# 11: Vignettes and example dataset

**What to build:** The three v0.1 vignettes — the direct-effect model, netmeta vs. Stan, validating the surface — plus a small packaged example network, together walking the spec's full suggested workflow: construct → plot network → fit both engines → compare → diagnostics → anchor → plot surface.

**Blocked by:** 06 (Engine comparison), 07 (Anchoring), 08 (Surface diagnostics), 09 (Network plot), 10 (Surface plots)

**Status:** done

- [x] Example dataset ships with the package and powers the vignettes
- [x] All three vignettes build cleanly under `R CMD check`
- [x] The full suggested workflow from the spec is exercised end to end
- [x] Terminology consistent throughout: direct effects, surface, sea level

## Comments

Done. Example data: `example_comparisons` (16 head-to-head log-HR
estimates over 8 fictional drugs — abrelor, belvane, cordexa, delmiro,
envarin, folzane, gravix, harlent) and `example_anchors` (2 placebo
anchors), generated deterministically by data-raw/example-network.R
(seed 731, tau = 0.02, generating truth printed in the script);
documented in R/data.R, LazyData enabled, data-raw Rbuildignored.
Three vignettes (knitr/rmarkdown, engine-gated eval so they degrade
gracefully when Suggests are absent): direct-effect-model.Rmd opens
with the mandatory A/B/C story (null comparison, non-null direct
effects) then runs construct → check_connectivity → plot_network →
fit_surface → anchor_surface → plot_effect_surface → plot_sea_level;
netmeta-vs-stan.Rmd fits both engines, compare_engines +
plot_engine_comparison + anchored comparison, and notes the CI-enforced
ε = 0.02; validating-the-surface.Rmd covers edge_residuals,
cycle_consistency (+ injected inconsistency), plot_cycle_consistency,
and simulate/validate_recovery. Together they exercise the spec's full
suggested workflow. Surface / sea level / direct effects terminology
used throughout. Stan usage kept to netmeta-vs-stan.Rmd to bound
vignette build time.
