# 02: Direct-effect network object

**What to build:** `direct_effect_network()` turns a comparisons table (plus an optional anchors table and an effect measure) into a validated `directeffect_network` holding comparisons, anchors, treatments, graph, connected components, effect measure, and metadata — and no fitted model. `check_connectivity()` reports, per component, the drug count, comparison count, anchor count, and whether absolute effects are identifiable there. Malformed input fails loudly at construction; obviously incompatible estimand metadata warns. Demoable on the spec's 3-drug example.

**Blocked by:** 01 (Package skeleton with green CI)

**Status:** done

- [x] Constructor accepts the comparisons schema (`study_id`, `target`, `comparator`, `estimate`, `std_error`) and optional anchors schema (`study_id`, `drug`, `reference`, `estimate`, `std_error`), preserving any extra metadata columns
- [x] Validation errors on missing required columns, non-positive standard errors, self-comparisons, and missing/NA drug names
- [x] Warning when comparisons with differing estimand metadata (population, time_at_risk, outcome_definition, database, study_design, estimand) are combined
- [x] Graph and connected components are computed at construction; `check_connectivity()` prints the per-component identifiability report from the spec
- [x] A component with no anchor is reported as "relative effects identifiable only" — the package never implies it can be absolutely positioned
- [x] A print method summarises the network: drugs, comparisons, anchors, components, effect measure

## Comments

Done. `direct_effect_network()` in R/network.R (constructor +
validation + estimand warning + print method), `check_connectivity()`
in R/connectivity.R (per-component report data frame with a print
method rendering the design's exact format). Graph stored as an igraph
graph with study/estimate/std_error edge attributes so ticket 09 can
draw it without graph logic in plotting code; igraph added to Imports.
Anchors naming drugs absent from all comparisons error at construction
(spec user story 4). 52 tests pass; `R CMD check` Status: OK.
