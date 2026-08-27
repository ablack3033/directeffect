# directeffect v0.1 — surface and sea level

Status: ready-for-agent

## Problem Statement

Comparative-effectiveness studies mostly produce *relative* effects: drug A versus drug B, estimated from active-comparator designs. Each estimate `y_ij` constrains only the difference `theta_i - theta_j` between the latent direct effects of the two drugs. The researcher ultimately wants each drug's *direct* effect versus placebo or no treatment, but almost none of the evidence speaks to that directly — and a null active-comparator estimate does not mean either drug is safe, only that the two drugs sit at the same height.

Today there is no tool that reconstructs latent direct effects from a network of comparative estimates while keeping two questions separate:

1. Are my comparative-effect estimates internally coherent? (the *shape* of the surface)
2. Where does that coherent surface sit relative to placebo? (the *sea level*)

Standard network meta-analysis software answers relative questions against an arbitrary reference, offers no principled separate step for anchoring to absolute evidence, and provides no independent cross-implementation check that the reconstruction is correct.

## Solution

An R package, `directeffect`, in which the researcher:

1. Builds an explicit, validated network object from a table of comparative estimates, optionally alongside a table of absolute (placebo) anchors — with no model fitted at construction time.
2. Inspects the evidence: network plot, connected components, and per-component identifiability (relative-only vs. absolutely identifiable).
3. Fits the relative direct-effect **surface** with either of two independent engines — frequentist (netmeta) or Bayesian (Stan) — which are required to agree on the same statistical model. Anchors are deliberately ignored at this stage.
4. Applies absolute anchors as a distinct second step (**sea-level** estimation) that positions the whole surface relative to placebo = 0.
5. Receives one engine-agnostic tidy fit representation, and runs surface diagnostics (edge residuals, cycle consistency, engine comparison) and plots (network, effect surface, surface + sea level).
6. Can simulate networks with known truth and measure recovery, from the first release.

The v0.1 goal, verbatim from the design: *demonstrate that the same direct-effect surface is reconstructed independently by netmeta and Stan.*

## User Stories

1. As a comparative-effectiveness researcher, I want to construct a direct-effect network from a table of drug-versus-drug estimates, so that my evidence base is an explicit, validated object before any model is fitted.
2. As a researcher, I want to optionally attach a table of placebo-controlled (absolute) estimates to the network, so that absolute evidence is represented separately from comparative evidence.
3. As a researcher, I want the constructor to record the effect measure (e.g. HR, RR, OR), so that estimates are handled on the log scale internally and labeled correctly on output.
4. As a researcher, I want construction to fail loudly on malformed input (missing columns, non-positive standard errors, self-comparisons, drugs appearing only in anchors), so that data problems surface before fitting rather than as wrong answers.
5. As a researcher, I want optional study metadata (database, population, design, drug class, follow-up, outcome, estimand, …) preserved through the pipeline, so that downstream diagnostics and future bias models can use it.
6. As a researcher, I want a warning when comparisons with obviously incompatible estimand descriptions (population, time-at-risk, outcome definition, design) are combined, so that statistical connectivity is never mistaken for causal transportability.
7. As a researcher, I want a connectivity report listing each connected component with its drug count, comparison count, and anchor count, so that I know which components support absolute effects and which support relative effects only.
8. As a researcher, I want the package to refuse to silently pick an arbitrary sea level for an unanchored component, so that a relative-only quantity is never presented as an absolute effect.
9. As a researcher, I want to fit the relative surface with a frequentist engine (netmeta), so that I get standard network meta-analysis estimates without leaving the package's abstraction.
10. As a researcher, I want to fit the same surface with a Bayesian engine (Stan) under deliberately weak priors, so that an independent implementation of the same model is available.
11. As a researcher, I want surface fitting to ignore anchors entirely and use an arbitrary reference (theta_ref = 0), so that I can judge the internal coherence of comparative evidence on its own.
12. As a researcher, I want to anchor a fitted surface with one or more placebo anchors as a separate step, so that the decision about absolute position is explicit and revisable.
13. As a researcher, I want anchors to retain their uncertainty rather than pinning a drug's effect exactly, so that absolute conclusions honestly reflect the strength of the anchoring evidence.
14. As a researcher, I want the anchored Bayesian model to need no arbitrary identification constraint, so that the anchors, not a convention, determine the absolute location.
15. As a researcher, I want every fit returned in one common tidy representation (drug, estimate, std_error, lower, upper, scale, reference, engine), so that downstream code never needs to know which engine produced it.
16. As a researcher using Stan, I want posterior summaries (median, mean, sd, quantiles) and convergence diagnostics (rhat, ess_bulk, ess_tail) included in the fit, so that I can judge MCMC reliability without touching Stan internals.
17. As a methodologist, I want a first-class engine comparison (per-drug netmeta vs. Stan estimates, differences, standardized differences) and a comparison plot with the identity line, so that agreement between independent implementations is a routine check, not a research project.
18. As a methodologist, I want edge residuals (observed minus predicted for every comparison, raw and standardized), so that I can find comparisons that conflict with the fitted surface.
19. As a methodologist, I want cycle-consistency sums over a cycle basis of the network, so that inconsistency is quantified without enumerating every cycle in a large graph.
20. As a researcher, I want a network plot (drugs as nodes, comparisons as edges, with optional encodings for precision, anchor presence, and comparison counts), so that I can see the structure of my evidence before fitting.
21. As a researcher, I want a one-dimensional effect-surface plot with uncertainty intervals, orderable from harmful to beneficial, on either the log or the natural scale, so that results are readable without implying a fake 2-D geometry.
22. As a researcher, I want a surface + sea-level plot showing the relative structure positioned against the placebo = 0 line, so that the package's central conceptual decomposition is visible in one figure.
23. As a methodologist, I want to simulate networks with known true effects (configurable drugs, comparisons, anchors, heterogeneity, seed) where the simulation object retains truth, so that recovery can be measured rather than assumed.
24. As a methodologist, I want a recovery validation (bias, RMSE, interval coverage, rank correlation) comparing a fit against simulation truth, so that estimator quality is quantified continuously.
25. As a package developer, I want both engines hidden behind adapters returning the common fit class, so that engine internals (netmeta objects, Stan matrices) never leak into the rest of the package.
26. As a package developer, I want the deterministic worked example (A = 0.7, B = 0.7, C = 0.3 with a null A–B comparison and a single C–placebo anchor) encoded as a permanent test, so that the package's central claim — a null active-comparator estimate does not imply a null direct effect — is protected forever.
27. As a package developer, I want the engine-equivalence check to run in CI on every commit, so that the two implementations keep validating each other from day one.
28. As a paper author, I want the package's terminology (direct effects, surface, sea level) consistent across API, docs, and plots, so that methods sections read naturally.

## Implementation Decisions

- **Two-stage inference is the core abstraction.** Surface estimation (relative positions from comparisons) and sea-level estimation (absolute position from anchors) are separate operations with separate functions. The separation appears in the API: `direct_effect_network()` → `fit_surface()` → `anchor_surface()`. A combined convenience wrapper (`fit_direct_effects()`) may come later but must preserve the separation internally.
- **Statistical model.** Comparisons: `y_k ~ N(theta_target[k] - theta_comparator[k], se_k^2)` (common-effect; no heterogeneity term in v0.1). Anchors: `a_m ~ N(theta_drug[m], se_m^2)` with the convention `theta_placebo = 0`. Surface-only fits identify via `theta_ref = 0` for an arbitrary reference; anchored fits need no constraint because the anchors locate the surface.
- **Comparisons schema** (one row = one estimated drug-vs-drug effect): `study_id`, `target`, `comparator`, `estimate`, `std_error`. Estimates for multiplicative measures (HR/RR/OR) are supplied and stored on the log scale. Optional metadata columns are preserved: `database`, `population`, `design`, `drug_class`, `followup`, `outcome`, `analysis_id`, `negative_control`, `calibrated`, `n_target`, `n_comparator`.
- **Anchors schema**: `study_id`, `drug`, `reference` (normally `"placebo"`), `estimate`, `std_error`.
- **Network object.** `direct_effect_network()` returns a `directeffect_network` holding comparisons, anchors, treatments, the graph, connected components, effect measure, and metadata — and no fitted model. Its job is definition and validation of the evidence network.
- **Two engines, one model.** The frequentist engine wraps netmeta (mapping estimate→TE, std_error→seTE, target→treat1, comparator→treat2, study_id→studlab; common and random summaries requested, common used in v0.1). The Bayesian engine implements the identical likelihood in Stan with a weak `normal(0, 5)` prior on free effects so comparison against netmeta is meaningful. Each engine lives behind an adapter; netmeta and Stan internals are never exposed outside their adapters.
- **Two Stan models, not one.** The unanchored surface model and the anchored surface model are separate Stan programs rather than one program full of conditionals — two simple deep modules over one complicated universal model.
- **Frequentist anchoring** is a single location-offset estimate: the offset that best reconciles the fitted surface with the anchors, precision-weighted, propagating both surface and anchor uncertainty. (This decision was filled in during spec synthesis — the design showed anchoring only on the Stan fit, but the mandatory unit test requires both engines to recover absolute effects. See Further Notes.)
- **Common fit class.** Both engines return `directeffect_fit` with components: effects, comparisons, anchors, heterogeneity, diagnostics, engine, engine_fit, network. `fit$effects` has the stable schema `drug`, `estimate`, `std_error`, `lower`, `upper`, `scale`, `reference`, `engine`; Stan adds `median`, `mean`, `sd`, `q025`, `q975`, `rhat`, `ess_bulk`, `ess_tail`. Downstream functions (all plots and diagnostics) must work identically on either engine's fit.
- **Diagnostics in v0.1**: `check_connectivity()` (per-component identifiability report, run before fitting), `edge_residuals()` (observed, predicted, residual, standardized residual per comparison), `cycle_consistency()` over a cycle basis (never full cycle enumeration), and `compare_engines()` + `plot_engine_comparison()`.
- **Plots in v0.1**: `plot_network()` built on igraph/ggraph (graph logic outside plotting code), `plot_effect_surface()` (one-dimensional, uncertainty intervals, `scale = "log"` or `"natural"`), a surface + sea-level plot, and the engine-comparison plot.
- **Simulation module**: `simulate_direct_effect_network(n_drugs, n_comparisons, n_anchors, heterogeneity, seed)` generating truth → graph → contrasts → standard errors → observed comparisons → anchors, retaining `simulation$truth`; `validate_recovery(fit, simulation)` reports bias, RMSE, coverage, rank correlation.
- **Class system stays minimal**: `directeffect_network` and `directeffect_fit` only in v0.1 (`directeffect_validation` and `directeffect_error_model` arrive with later versions). No classes for intermediate operations.
- **No biological structure in the base model.** Drugs in the same class are not shrunk toward each other; class priors are a later, opt-in extension.
- **Estimand discipline.** The network object supports estimand-describing metadata and warns when obviously different specifications are combined; the package never implies that connectivity establishes transportability.
- **Naming.** Package `directeffect`; core abstraction `direct_effect_network`; the surface / sea-level metaphor is used in terminology, documentation, and plots but the public names stay scientifically sober.

## Testing Decisions

- **What makes a good test here:** exercise only the public seam — `direct_effect_network()` → `fit_surface()` / `anchor_surface()` → the `directeffect_fit` contract and the diagnostics that consume it. Never assert on netmeta or Stan internals; engines are implementation details behind adapters. Correctness comes from three oracles that need no mocking: hand-computable examples, cross-engine agreement, and simulation truth.
- **Mandatory deterministic test.** True effects placebo = 0, A = 0.7, B = 0.7, C = 0.3; comparisons A−B = 0.0, A−C = 0.4, B−C = 0.4; anchor C−placebo = 0.3. Both engines must recover A ≈ 0.7, B ≈ 0.7, C ≈ 0.3. This directly protects the central property: a null active-comparator estimate does not imply a null direct effect.
- **Engine equivalence test (CI, every commit).** Simulate under Gaussian estimates, known standard errors, no heterogeneity, weak priors; require `|theta_i_stan − theta_i_netmeta| < epsilon` for every treatment.
- **Identification test.** Surface fits are invariant to the choice of reference: shifting all true effects by a constant changes nothing observable; only differences are recovered from comparisons alone.
- **Anchor tests.** Anchoring positions the surface correctly and propagates anchor uncertainty; an unanchored component is reported relative-only and never silently anchored.
- **Cycle tests.** Consistent simulated networks give cycle sums ≈ 0; an injected inconsistent edge is flagged.
- **Connectivity tests.** Multi-component networks produce the correct per-component report and identifiability verdicts.
- **Simulation-recovery test.** `validate_recovery()` on simulated data shows near-zero bias and nominal interval coverage in the low-noise regime.
- **Prior art:** none — the repo is greenfield. Use testthat; Stan-dependent tests should tolerate MCMC noise via tolerances and fixed seeds, and CI must run the equivalence check continuously.

## Out of Scope

Everything below is designed (see `design.md` in this directory) but explicitly deferred:

- **v0.2 — robustness:** random-effects heterogeneity (global tau), leave-one-edge-out cross-validation, anchor roles (`fit` / `validate`) and `validate_anchors()`, leave-one-anchor-out (`loo_anchors()`), `comparator_robustness()`, `anchor_influence()` + its plot, `explain_effect()` evidence decomposition.
- **v0.3 — observational systematic error:** negative-control calibration, `fit_systematic_error()`, the `directeffect_error_model` class, database- and design-specific bias terms.
- **v0.4 — biological structure:** class-level hierarchical priors, mechanism-of-action priors, multiple outcomes, multivariate direct-effect surfaces.
- The `fit_direct_effects()` convenience wrapper (allowed later; not needed to prove the mathematics).
- Treatment-specific or comparison-specific heterogeneity (even v0.2 starts with one global tau).
- Any renaming — `directeffect` is the decided name; alternatives (causalsurface, effectmap, effectnet, drugposition, causalmap) were considered and rejected.

v0.3 and v0.4 must not be built until v0.1 and v0.2 validation demonstrates the basic model is useful.

## Further Notes

- **Design principle:** the package exposes the scientific abstraction (relative evidence → surface → + anchors → absolute direct effects) and hides netmeta objects, Stan matrices, igraph internals, and MCMC details unless explicitly requested via `engine_fit`.
- **Decision filled in during synthesis, flagged for review:** frequentist anchoring (the precision-weighted location offset) was not specified in the design document, which demonstrated `anchor_surface()` only on a Stan fit — but the mandatory unit test requires *both* engines to recover absolute effects, so some frequentist anchoring path must exist in v0.1. If you'd rather scope the netmeta engine to surface-only and run the deterministic absolute-recovery test against Stan alone, strike that decision and adjust the test.
- **Open question:** whether `direct_effect_network()` should also accept natural-scale estimates and log-transform them itself, or require callers to pass log-scale values (the design's examples pass `log(1.02)` etc. — the spec assumes caller-supplied log scale).
- **Vignettes** planned for v0.1 (direct-effect model, netmeta vs. Stan, validating the surface) plus a small example dataset are part of proving the mathematics publicly, but can trail the code within the version.
- The full design document — including v0.2–v0.4 details omitted here — is preserved as `design.md` alongside this spec.
