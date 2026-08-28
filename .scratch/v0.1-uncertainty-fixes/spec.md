# directeffect v0.1-uncertainty-fixes — honest evidence accounting

Status: ready-for-agent

Origin: statistical code review of 2026-08-28
(`.scratch/reviews/2026-08-28-stats-review.md`, commit `30a50af`).
Every defect this spec addresses was demonstrated numerically there;
the demonstration scripts live in `.scratch/reviews/demos/`.

## Problem Statement

The package's point estimates are right, but its uncertainty accounting
is not yet trustworthy — in a package whose stated purpose is honest
evidence accounting, and whose central claim is that two independent
engines continuously validate each other.

Concretely, as demonstrated in the review:

1. The frequentist `anchor_surface()` path reports absolute standard
   errors that are wrong for every drug except (coincidentally) the
   reference: on the package's own mandatory example, the anchored
   drug's SE is 76% too wide (0.0702 reported vs a true sampling SD of
   0.0400), and nominal 95% intervals cover at up to 99.9%. The
   anchored Stan engine gets the same quantities right, so the two
   engines silently disagree about anchored uncertainty.
2. A multi-arm study (several comparison rows sharing a `study_id`) is
   handled correctly by netmeta but treated as independent rows by the
   Stan likelihood and by the diagnostics — the "identical model /
   identical likelihood" claim in the README and docs is false for such
   inputs, precision is overstated by √(2/3) on a three-arm trial, and
   nothing detects the situation. The estimate-only engine-equivalence
   test stays green throughout.
3. `edge_residuals()` standardizes by the raw comparison SE, ignoring
   leverage: a bridge comparison (the only path between its endpoints)
   gets a standardized residual of exactly 0 whatever its estimate, so
   the least-corroborated comparisons are guaranteed to look perfectly
   consistent.
4. The test suite cannot catch any of the above: the only test of the
   frequentist anchoring formula re-derives that same formula (a
   circular oracle), no test compares engine standard errors, and no
   test measures anchored coverage.
5. Fitting with the Stan engine and no explicit seed consumes the
   caller's global RNG stream, and `validate_recovery()`'s rank
   correlation includes the pinned reference's exact (0, 0) pair while
   its other metrics exclude it.

## Solution

Make the uncertainty accounting as trustworthy as the point estimates,
without changing the package's two-stage abstraction or its public API:
propagate the surface's full covariance through frequentist anchoring
(the two-stage estimator then matches the joint GLS answer the anchored
Stan model already gives), detect multi-arm input loudly instead of
silently mishandling it, make edge residuals leverage-aware, upgrade the
test suite so every one of these properties is guarded by a
non-circular oracle, and clean up the RNG and metric hygiene. After
this version, the engines agree on uncertainties as well as estimates,
and every documented claim about uncertainty is either true or loudly
refused.

## User Stories

1. As a comparative-effectiveness researcher, I want the absolute
   standard errors from a frequentist anchored fit to equal the actual
   sampling variability of the estimator, so that my anchored intervals
   are calibrated rather than up to 76% too wide.
2. As a researcher, I want a drug carrying the network's single anchor
   to inherit exactly that anchor's uncertainty in the anchored fit, so
   that the algebra of my evidence (the surface contribution cancels)
   is reflected honestly in the output.
3. As a methodologist, I want the frequentist and Bayesian engines to
   agree on anchored standard errors, not just point estimates, so that
   the cross-engine validation the package is built on actually
   validates uncertainty.
4. As a researcher, I want the precision weights used to combine
   multiple anchors to account for the covariance between the anchored
   drugs' surface positions, so that disagreeing anchors are reconciled
   with the correct relative influence.
5. As a researcher, I want `plot_sea_level()` and printed anchored
   effects to inherit the corrected intervals automatically, so that
   every consumer of the fit contract benefits without new options.
6. As a researcher with multi-arm trials in my evidence base, I want
   network construction to tell me loudly that a study contributes more
   than one comparison, so that I know I am in territory the engines
   treat differently before I fit anything.
7. As a Bayesian user, I want the Stan engine to refuse a multi-arm
   network with an informative error naming the studies and pointing at
   the netmeta engine, so that I cannot unknowingly get √(2/3)-deflated
   standard errors from an independence likelihood.
8. As a methodologist, I want `cycle_consistency()` and
   `edge_residuals()` to carry the same multi-arm warning, so that
   diagnostics built on an independence assumption announce it when the
   assumption is violated.
9. As a paper author, I want the README, `fit_surface()` docs, and spec
   to state the model-identity claim precisely (identical likelihood
   for one-comparison-per-study networks in v0.1), so that the
   package's central claim is exactly true.
10. As a methodologist, I want standardized edge residuals divided by
    the residual's actual standard deviation (leverage-aware), so that
    the z-values I screen have unit variance under coherence.
11. As a methodologist, I want a bridge comparison's standardized
    residual reported as `NA` alongside its leverage rather than as 0,
    so that "uncheckable" is never displayed as "perfectly consistent".
12. As a package developer, I want the fit contract to carry the
    surface covariance from both engines, so that anchoring and
    leverage computations consume the contract rather than engine
    internals.
13. As a package developer, I want the engine-equivalence test to
    assert agreement of standard errors as well as estimates, so that a
    divergence in uncertainty between engines fails CI the way an
    estimate divergence always has.
14. As a package developer, I want the frequentist anchoring tested
    against an independent joint-GLS oracle (anchor rows appended to
    the design), so that the test can fail when the implementation is
    statistically wrong, unlike the current circular hand-derivation.
15. As a methodologist, I want a Monte Carlo coverage test of anchored
    intervals with bounds tied to the binomial Monte Carlo error, so
    that miscalibration of the package's headline output is caught
    continuously rather than hidden by a ±10-point tolerance.
16. As a package developer, I want a multi-arm fixture in the test
    suite, so that the input class the simulator cannot generate stops
    being structurally untestable.
17. As a researcher, I want `fit_surface()` and `anchor_surface()` to
    leave my global random-number state untouched when I do not pass a
    seed, so that fitting a model never silently changes my session's
    reproducibility — the standard `simulate_direct_effect_network()`
    already meets.
18. As a methodologist, I want `validate_recovery()`'s rank correlation
    computed over the same drugs as its bias, RMSE, and coverage
    (excluding the pinned reference), so that one phantom exact pair
    cannot move Spearman's ρ by ±0.3.
19. As a researcher, I want the mandatory deterministic example
    (A = 0.7, B = 0.7, C = 0.3) to keep passing on both engines with
    its corrected uncertainties, so that the package's founding
    property survives the fix.
20. As a package developer, I want the `directeffect_fit` effects
    schema to remain stable (corrections change values, not columns),
    so that downstream consumers of the contract need no changes.
21. As a researcher reading the docs, I want the anchoring
    documentation's promise — "their standard errors propagate into
    every absolute interval" — to describe the corrected propagation,
    so that prose and computation agree.
22. As a package developer, I want the v0.1-core issue record updated
    where it deferred the covariance problem as "acceptable for v0.1",
    so that the tracker reflects that the deferral was overturned by
    measurement.

## Implementation Decisions

- **The fix is variance algebra, not architecture.** The two-stage
  surface → sea-level abstraction, the public API, and the effects
  schema are unchanged. Only the numbers inside the anchored
  frequentist SEs, the standardized residuals, and the guarded input
  space change.
- **The fit contract gains a surface covariance.** Both engine adapters
  supply the covariance of the estimated effects versus the reference
  (frequentist: the common-effect covariance matrix the engine already
  computes; Bayesian: the posterior covariance of theta) as a new fit
  component alongside `effects`. It is an additive, non-breaking
  contract extension; the reference row is exact zeros. Anchoring and
  leverage-aware residuals consume it through the contract — engine
  internals stay behind the adapters, per the existing ADR-style
  decision that nothing outside an adapter touches an engine object.
- **Frequentist anchoring propagates the full covariance.** The offset
  remains a single location shift, but its variance and its covariance
  with every drug's surface position are computed from the surface
  covariance matrix, and each absolute variance is
  `Var(θ̂_d) + Var(offset) + 2·Cov(θ̂_d, offset)`. With one anchor this
  reduces to: anchored drug's SE = the anchor's SE, exactly; in general
  it agrees with the joint GLS of comparisons plus anchor rows, which
  is what the anchored Stan model computes — restoring cross-engine
  agreement. (The precision weights for multiple anchors likewise come
  from the anchors' proposal covariance, not just its diagonal.)
- **Multi-arm input is detected at the seam, not modeled.** Network
  construction warns when any `study_id` contributes more than one
  comparison, naming the studies. The Stan engine refuses such
  networks with an error that explains why (independence likelihood)
  and points to the netmeta engine, which handles them correctly.
  `cycle_consistency()` warns on multi-arm networks that its pooled
  edge variances assume independent comparisons. Modeling the
  within-study covariance in Stan is explicitly deferred (see Out of
  Scope), consistent with the package's fail-loudly philosophy: detect
  and refuse rather than silently approximate.
- **Docs state the model-identity claim precisely.** README,
  `fit_surface()` docs, and the spec language change from "identical
  likelihood" to identical likelihood for networks with one comparison
  per study, with multi-arm evidence supported by the frequentist
  engine only in this version.
- **Edge residuals become leverage-aware.** Standardization divides by
  the residual's actual standard deviation (comparison variance minus
  prediction variance, obtained from the surface covariance). The
  output gains a leverage column; where leverage is 1 within numerical
  tolerance (bridge comparisons), the standardized residual is `NA`
  and the documentation explains that such comparisons are
  uncorroborated rather than consistent.
- **RNG hygiene.** The Stan fitting functions' seed argument defaults
  to unset; when the caller supplies none, no draw is taken from R's
  global stream on the caller's behalf (the sampler's own seeding
  applies), and the fitting path is tested to leave `.Random.seed`
  untouched, matching the simulation module's existing standard.
- **Recovery metrics align.** Rank correlation in
  `validate_recovery()` is computed over the free drugs only, like
  bias, RMSE, and coverage.
- **Tracker correction.** The v0.1-core anchoring issue's addendum
  ("conservative rather than miscalibrated-narrow … acceptable for
  v0.1") is annotated as superseded: the review measured 76% SE
  inflation, up-to-99.9% coverage, and engine disagreement.

## Testing Decisions

- **What makes a good test here (unchanged philosophy):** exercise only
  the public seam — `direct_effect_network()` → `fit_surface()` /
  `anchor_surface()` → the fit contract and the diagnostics that
  consume it — and judge outputs against oracles that do not restate
  the implementation: hand-computable closed forms, an independent
  joint-GLS computation, cross-engine agreement, and simulation truth.
  The review's demonstration scripts are the seed material for the new
  oracles.
- **Replace the circular anchoring test** with the joint-GLS oracle:
  append anchor rows to the comparison design, solve the weighted
  normal equations independently, and require the anchored frequentist
  estimates and standard errors to match. Include the single-anchor
  special case (anchored drug's SE equals the anchor's SE exactly) as a
  closed-form assertion.
- **Extend engine equivalence to uncertainty:** the equivalence test
  additionally bounds the per-drug difference between frequentist SEs
  and Stan posterior SDs, on both surface and anchored fits, with the
  tolerance justified against Monte Carlo error at the test's sampler
  settings, as the existing EPSILON already is for estimates.
- **Anchored coverage test:** a Monte Carlo calibration test on a small
  network (a few hundred replicates of the fast frequentist path)
  requiring empirical coverage of nominal 95% intervals within a bound
  derived from binomial MC error — replacing reliance on the single
  loose `>= 0.85` unanchored assertion.
- **Multi-arm fixture tests:** a hand-built three-arm study asserting
  the construction warning, the Stan engine's refusal, the netmeta
  engine's correct SE (equal to the single trial's contrast SE), and
  the diagnostics' warning.
- **Residual tests:** the bridge fixture from the review (triangle plus
  bridge) asserting `NA` standardized residual and leverage 1 on the
  bridge and unit-variance standardization on cycle edges (checkable
  against the joint-GLS hat matrix).
- **RNG test:** mirror the existing "simulation does not disturb the
  global random-number state" test for an unseeded Stan fit.
- **Prior art:** the mandatory-example tests, the `wls_surface()`
  helper (extends naturally to anchor rows), `collect_warnings()`, the
  RNG-restoration test in the simulation tests, and the review's
  `demos/demo_fast.R` / `demos/demo_stan.R`.

## Out of Scope

- Modeling within-study covariance of multi-arm trials in the Stan
  likelihood (arm-based or contrast-based with covariance): a v0.2+
  feature; v0.1-uncertainty-fixes detects and refuses instead.
- Random-effects heterogeneity and everything else scheduled for
  v0.2–v0.4 in the design (leave-one-out machinery, anchor roles,
  systematic-error models, class priors).
- Reworking `compare_engines()`'s `standardized_difference` into a
  formal test statistic (the two fits share data; it remains a
  yardstick, and its documentation may say so, but no new statistic is
  introduced here).
- Any change to point-estimation, the effects schema, plotting
  geometry, or the package's two-stage abstraction.
- Back-transformation helpers or natural-scale SEs (bounds-only
  back-transformation stays, and stays correct).

## Further Notes

- The review confirmed the spec's own suspicion: v0.1-core's Further
  Notes flagged the frequentist anchoring as "filled in during
  synthesis, flagged for review". This spec is that review's verdict:
  keep the two-stage frequentist path, fix its covariance algebra.
- Overstated uncertainty is not "safe": it misranks drugs against
  decision thresholds, and because only one engine is wrong it
  contradicts the package's central cross-validation claim. Treat
  finding severity accordingly — this is the package's headline output.
- The reference-drug coincidence (its anchored SE is currently correct
  when anchors sit elsewhere) means spot checks at the reference will
  not reveal the defect; oracles must cover non-reference drugs.
- All numeric assertions proposed above have already been computed once
  in `.scratch/reviews/demos/`; porting them into testthat is
  transcription, not research.
