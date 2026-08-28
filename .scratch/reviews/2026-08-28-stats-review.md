# Statistical code review — directeffect v0.1

Reviewed with `.claude/skills/stats-review` on 2026-08-28, at commit
`30a50af`. Every finding below was demonstrated numerically with the
package's own code; the demonstration scripts are committed next to
this report under `demos/` (run from the repo root) and the key numbers
are reproduced inline.

## Verdict

The surface half of the package is statistically sound: both engines
reproduce the exact weighted-least-squares solution of the stated model,
orientation and scale conventions check out end to end, and the
identification story (arbitrary reference, pinned at exact zero) is
handled honestly. The sea-level half is not: the frequentist
`anchor_surface()` path reports standard errors that are wrong for every
drug except (coincidentally) the reference — up to 76% too wide on the
package's own mandatory example — because it drops every covariance in
the two-stage estimator. A cluster of secondary findings (multi-arm
studies silently violate the "identical likelihood" claim, bridge edges
get vacuous residuals, one circular test protects the broken formula)
follows the same theme: the point estimates are right, the *evidence
accounting* around them is not yet trustworthy.

## Findings (worst first)

### 1. [Wrong uncertainty] Frequentist anchoring misstates every absolute SE and interval

- Where: `R/anchor-surface.R:86-120` (`anchor_surface_offset`)
- What the code does: absolute variance = surface variance + offset
  variance, treating the offset as independent of each drug's surface
  estimate and the surface estimates as independent of each other. The
  offset is built *from* those same surface estimates
  (`offset = Σ w_m (a_m − θ̂_{drug(m)}) / Σ w_m`), so the dropped
  covariances are not small corrections — for the anchored drug the
  θ̂ terms cancel exactly and the code double-counts its surface
  variance twice over.
- What the mathematics gives: with a single anchor on drug C, the
  absolute estimate of C is algebraically `a` (the code's own point
  estimate confirms this), so its SE is the anchor's SE, exactly. In
  general the correct covariance is the GLS covariance of the joint
  Gaussian model (comparisons + anchor rows), which the package's own
  anchored Stan engine already implements.
- Demonstration (mandatory example: A−B = 0.0, A−C = 0.4, B−C = 0.4,
  all SE 0.05; anchor C = 0.3, SE 0.04; truth A = 0.7, B = 0.7,
  C = 0.3):

  | drug | reported SE | true sampling SD (20k-rep MC) | exact GLS SE | empirical coverage of the nominal 95% CI |
  |------|------------|------------------------------|--------------|------------------------------------------|
  | A    | 0.0572     | 0.0571                       | 0.0572       | 0.950 |
  | B    | 0.0702     | 0.0575                       | 0.0572       | 0.983 |
  | C    | 0.0702     | 0.0401                       | 0.0400       | 0.999 |

  The Monte Carlo redraws data from the model and reruns the package's
  own two-stage pipeline; the reported SE for C is 76% larger than the
  actual sampling SD of the package's own estimate. The package's
  Bayesian engine on the same data returns posterior SDs 0.0567 /
  0.0569 / 0.0396 (8000 iterations, matching the exact GLS column
  within Monte Carlo error) — the two engines disagree about the
  anchored uncertainties, and only the frequentist one is wrong. (A is right by coincidence: for the reference drug,
  surface variance 0 + offset variance happens to equal the correct
  formula when the anchors sit elsewhere.)
- Consequence: the anchored fit is the package's headline output.
  Everyone reading `absolute$effects` intervals, `plot_sea_level()`
  ribbons, or anchored `validate_recovery()` coverage sees overstated
  uncertainty — drugs appear less precisely located than the evidence
  warrants, in a package whose stated purpose is honest evidence
  accounting. The spec itself flagged this decision as "filled in
  during synthesis, flagged for review" (spec.md, Further Notes);
  this review confirms the flag was warranted.
- Fix: compute the offset and the absolute variances from the surface's
  full covariance matrix instead of its diagonal. `engine_fit` already
  carries it (`nm$Cov.common`, rows/columns vs the reference), and the
  two-stage estimator with correct propagation reduces to the same GLS
  answer the anchored Stan model gives. Keep the two-stage API;
  fix only the variance algebra.

### 2. [Silent scope] Multi-arm studies break the "identical likelihood" claim, undetected

- Where: `inst/stan/surface.stan:29`, `inst/stan/anchored_surface.stan:27`
  (rows independent); `R/engine-netmeta.R:12-22` (netmeta applies its
  multi-arm correction); `R/network.R` `validate_comparisons()` (no
  detection); claim in `README.md:30`, `R/fit-surface.R` docs
  ("identical likelihood"), spec.md ("Two engines, one model").
- What the code does: a study contributing several rows (a multi-arm
  trial, entered the standard way — one row per pairwise contrast,
  shared `study_id`) is treated as independent rows by the Stan
  likelihood, by `cycle_consistency()`'s pooling, and by
  `edge_residuals()`. netmeta, by design, detects the shared `studlab`
  and reweights for the within-study correlation. Nothing in the
  package detects the situation or warns.
- What the mathematics gives: a three-arm trial's three contrasts are
  linearly dependent (the third is determined by the other two);
  counting them as independent inflates information by 3/2.
- Demonstration: one 3-arm trial (B−A = 0.1, C−A = 0.3, C−B = 0.2, all
  contrast SEs 0.1). netmeta: SE(B−A) = 0.1000 (correct — all evidence
  is that one trial). Independent-rows WLS, which is exactly the Stan
  likelihood: SE = 0.0816 = 0.1·√(2/3), 18% overstated precision.
  Running the package's own engines on this network confirms it:
  Stan posterior SDs 0.0810 / 0.0819 vs netmeta 0.1000 (ratios 0.810
  and 0.819, i.e. √(2/3) within MC error), while the point-estimate
  differences are ≤ 2.5e-4 — so `compare_engines()`'s `difference`
  column shows nothing and the CI equivalence test (EPSILON = 0.02 on
  estimates only) passes while the engines fit different likelihoods.
- Consequence: any user with real multi-arm evidence — ubiquitous in
  drug networks — gets silently different answers from the two engines
  and overstated Bayesian precision, in the package whose central claim
  is that the engines validate each other. The diagnostics inherit the
  same independence assumption.
- Fix (v0.1-honest option): detect repeated `study_id` across rows in
  `validate_comparisons()` and warn or refuse — this matches the
  package's "fail loudly" philosophy and keeps v0.1's scope. The full
  fix (multi-arm covariance in the Stan likelihood) can wait for v0.2.

### 3. [Wrong uncertainty] Edge residuals ignore leverage; bridge edges are vacuously zero

- Where: `R/diagnostics.R:31-49` (`edge_residuals`)
- What the code does: `standardized_residual = (observed − predicted) /
  std_error`, where the prediction comes from a fit that used the
  observation itself.
- What the mathematics gives: Var(observed − predicted) =
  se² − Var(predicted) < se², so the z-values are conservative
  everywhere; in the limiting case of a bridge edge (the only path
  between its endpoints), the fit reproduces the observation exactly
  and the residual is identically zero whatever the data say.
- Demonstration: network A–B–C triangle plus bridge C–D. With the C–D
  estimate at 0.2 its standardized residual is 0 (3.8e-14); changing
  the C–D estimate to 5.0 — a hazard ratio of 148 — leaves it 0
  (2.0e-12).
- Consequence: the documented reading ("large standardized residuals
  mark comparisons that conflict with the surface implied by the rest
  of the network") invites exactly the wrong conclusion for
  uncorroborated comparisons: the least-checkable edges in the network
  are guaranteed to look perfectly consistent.
- Fix: standardize by √(se² − Var(predicted)) (the hat-matrix
  diagonal is available from the WLS/netmeta machinery), and report
  bridge edges as `NA` with an explanatory column rather than 0.

### 4. [Fragile evidence] The test suite cannot catch findings 1–2

- Where: `tests/testthat/test-anchor-surface.R:132-167`,
  `tests/testthat/test-engine-equivalence.R`,
  `tests/testthat/test-simulate.R:74-90`.
- What the tests do:
  - "multiple anchors combine by precision, checked by hand" re-derives
    the implementation's own offset formula line by line and asserts
    agreement — a circular oracle. It asserts the exact SEs finding 1
    proves wrong, so it *protects* the defect.
  - The engine-equivalence test compares point estimates only
    (`comparison$difference`), never `std_error`. It is blind to
    finding 1 (where the Stan engine is right and the netmeta path
    wrong) and to finding 2 (where the SEs differ by √(2/3)); its
    simulator also assigns a fresh `study_id` to every row, so
    multi-arm inputs are structurally unreachable by it.
  - Coverage is asserted once, unanchored, at `>= 0.85` for a nominal
    95% interval — loose enough to hide serious miscalibration; the
    anchored path's coverage (empirically ~0.99 under finding 1) is
    never tested.
- Fix: assert SE agreement in the equivalence test; add an
  anchored-fit Monte Carlo coverage test with bounds tied to the
  binomial MC error (e.g. 500 reps ⇒ 95% ± 2%); replace the circular
  anchor test with the GLS oracle from finding 1's demonstration.

### 5. [Housekeeping] RNG and metric hygiene

- `R/engine-stan.R:22,73`: the default
  `seed = sample.int(.Machine$integer.max, 1)` silently consumes the
  caller's RNG stream — after `set.seed(7)` followed by an unseeded
  `fit_surface(de, engine = "stan")`, the user's next `rnorm()` returns
  0.981 where the seeded sequence gives 2.287
  (demonstrated with a live fit). `simulate_direct_effect_network()`
  carefully restores `.Random.seed`; the fitting functions undo that
  care. Fix: leave the seed argument `NULL` and only pass one to rstan
  when the user supplied it, or save and restore `.Random.seed` around
  the draw.
- `R/simulate.R:204`: `rank_correlation` is computed over all drugs
  including the pinned reference's exact (0, 0) pair, while bias, RMSE
  and coverage exclude it. In noisy networks the phantom row moves
  Spearman's ρ materially (seed 4 of the demo: 0.657 reported vs 0.400
  among actually-estimated drugs; shifts of ±0.3 across seeds). Fix:
  `cor(effects$estimate[free], theta_true[free], ...)`.

## Claim inventory

| claim | where | verdict |
|---|---|---|
| Input validation fails loudly (columns, SEs > 0, self-comparisons, orphan anchors) | `R/network.R` | verified (code + tests) |
| Surface estimates = common-effect NMA of stated model, netmeta engine | `R/engine-netmeta.R` | verified — matches exact GLS to 1e-10 on mandatory example |
| Surface estimates, Stan engine, agree with netmeta (single-arm rows) | `R/engine-stan.R`, `inst/stan/surface.stan` | verified (equivalence test + review demos) |
| Sign/orientation conventions (target vs comparator → TE, theta indexing) | both engines | verified end to end (A = 0.7 recovery) |
| "Identical likelihood" for both engines | README, `fit_surface()` docs | **finding 2** — false for multi-arm studies, undetected |
| Reference choice arbitrary; differences invariant | `fit_surface()` docs | verified (tested; pinned row handled honestly) |
| Anchored Stan model: no constraint, anchor uncertainty propagates | `inst/stan/anchored_surface.stan` | verified — posterior SDs match exact GLS |
| Anchored frequentist point estimates | `R/anchor-surface.R` | verified (precision-weighted offset; MC unbiased) |
| Anchored frequentist SEs/intervals "propagate both uncertainties" | `R/anchor-surface.R` | **finding 1** — wrong for every non-reference drug |
| Refuses to position unanchored components | `R/anchor-surface.R`, `R/connectivity.R` | verified (code + tests) |
| Heterogeneity block (Q, tau, I²) passed through from netmeta | `R/engine-netmeta.R` | not checked — netmeta internals, out of scope |
| `edge_residuals` standardization | `R/diagnostics.R` | **finding 3** — leverage ignored; bridges vacuous |
| `cycle_consistency` sums, SEs, z over a cycle basis | `R/diagnostics.R` | verified against hand computation (0.05 / 0.0866 / 0.577); independence caveat inherits finding 2 |
| `compare_engines` difference / standardized difference | `R/compare-engines.R` | verified arithmetic; note: the two fits share data, so `standardized_difference` is a yardstick, not a z-test — and no SE comparison exists (finding 4) |
| `check_connectivity` per-component report | `R/connectivity.R` | verified (code + tests) |
| Simulation matches the fitting model; truth retained; RNG restored | `R/simulate.R` | verified (code + tests) |
| `validate_recovery` bias/RMSE/coverage | `R/simulate.R` | verified arithmetic; anchored-coverage miscalibration is finding 1's downstream effect |
| `validate_recovery` rank correlation | `R/simulate.R` | **finding 5b** — includes the pinned reference pair |
| Stan convergence surfaced (R-hat, ESS, warning) | `R/engine-stan.R` | verified (thresholds 1.01 / 400; reference row excluded correctly) |
| Natural-scale plots transform bounds, not SEs | `R/plot-surface.R` | verified (exp on estimate/lower/upper only) |
| Seeded paths reproducible; RNG hygiene | `R/simulate.R`, `R/engine-stan.R` | **finding 5a** for the fitting default seed |
| Example data regenerable from committed script | `data-raw/`, `test-example-data.R` | verified (test re-runs the generation) |

## What the tests do and do not establish

The suite genuinely establishes: the mandatory deterministic recovery
(both engines), surface-level engine agreement on point estimates for
single-arm-per-study networks, refusal behaviors, connectivity
reporting, simulation reproducibility and RNG restoration, and the fit
contract's schema stability. These are real, well-chosen oracles.

It does not establish: any property of the anchored frequentist
*uncertainties* (the only test of them is circular), any SE-level
agreement between engines, any behavior on multi-arm studies (the
simulator cannot generate them), or interval calibration to a
resolution better than ±10 percentage points. Those four gaps are
precisely where findings 1, 2 and 4 live: the package's declared
validation philosophy ("recovery measured, not assumed") is currently
implemented for point estimates but not for uncertainty.
