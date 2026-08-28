---
name: stats-review
description: Statistical code review by an expert statistical programmer — audits estimators, uncertainty, assumptions, randomness, and the tests that guard them, and demonstrates every finding numerically before reporting it. Use whenever the user asks to review, audit, check, or validate statistical code or an analysis: R packages, Stan/JAGS/PyMC models, meta-analysis, causal inference, simulation studies, ML evaluation, or any code that produces estimates, standard errors, intervals, p-values, or posteriors — even when they just say "review this" and the code is statistical.
---

# Statistical code review

Ordinary code review asks whether the code does what the author intended.
Statistical review asks a harder question: whether what the code does is
correct *inference* — and the author's tests usually cannot answer it,
because a wrong formula implemented faithfully passes every test written
from that same formula. Your leverage is independence: derive what the
computation should be from the model, then hold the code to it.

The unit of review is a **claim** — one estimator, interval, test
statistic, diagnostic, or promise in prose — never a file. Style,
naming, and architecture are out of scope unless they cause a wrong
number.

## Process

### 1. Inventory the claims

Read the prose first (README, spec/design docs, function docs,
vignettes), then the code. Build a table of every statistical claim the
software ships:

- every exported function that returns an estimate, standard error,
  interval, p-value, posterior summary, or diagnostic;
- every promise in documentation: "unbiased", "exact", "identical
  likelihood", "propagates uncertainty", "coverage", "reproducible";
- every generative claim in simulation code ("matches the fitting
  model").

Done when: every exported statistical function and every documented
statistical promise appears in the inventory. This table is the review's
backbone — each row must end the review with a verdict: **verified**,
**finding**, or **not checked** (with the reason).

### 2. Trace each claim from math to code

For each row, derive the correct computation independently — from the
stated model, on paper — before reading the implementation, so the code
is checked against the mathematics rather than against itself. Then read
the implementation line by line against your derivation. The checklist
below lists where statistical implementations characteristically bend;
walk all of it for each claim it touches.

### 3. Demonstrate every suspected finding

A finding you have not demonstrated is a hypothesis, and reviews full of
hypotheses get ignored. For each suspect, construct the smallest input
that separates the correct computation from the implemented one — a
closed-form case, a two- or three-point dataset you can compute by hand,
or a Monte Carlo check whose truth is known by construction. Run it with
the project's own runtime when one is available; otherwise present the
worked numbers by hand. Either way the finding must show two numbers
that disagree: what the code reports and what the mathematics gives.

Discard any suspect whose demonstration shows no difference. This step
kills the plausible-sounding false positives that make statistical
reviews unactionable.

### 4. Report

Use the report format at the end. Rank findings by the severity ladder,
worst first. For large codebases, steps 1–3 may be split across parallel
subagents by claim group (estimation, uncertainty, diagnostics,
simulation/tests) — but verification of each finding stays with whoever
reports it.

## Where statistical code bends

A flat checklist; apply every line that touches the claim under review.

**Estimators.**
- Sign and orientation conventions: target vs comparator, treat1 vs
  treat2, exposure vs referent — verify one concrete case end to end.
- Scale: log vs natural, and where the transform happens. Intervals
  back-transform by transforming the *bounds*; standard errors do not
  back-transform without the delta method.
- Weights: precision weights use variances, not standard errors; pooled
  variance is `1/sum(w)`, not `mean`.
- Identification: what pins the model (reference level, sum-to-zero,
  anchor)? Confirm reported quantities are invariant to arbitrary
  choices the docs call arbitrary.

**Uncertainty.** The most common serious defect class: point estimates
right, uncertainty wrong.
- Independence assumed where correlation exists. Any quantity built
  from several estimates that share data (same study, same network fit,
  estimate-plus-derived-offset) has covariance terms; adding variances
  silently drops them. Wrong in both directions: double-counting
  (conservative) and undercounting (anticonservative) — both are
  findings, because both misstate the evidence.
- Variance components dropped: heterogeneity, estimation error of a
  plug-in (a "known" SE, weight, or offset that was itself estimated),
  leverage in residual standardization (a residual of an observation
  that helped fit its own prediction has variance below its SE; an
  observation that *fully determines* its prediction — a bridge edge, a
  saturated cell — has residual identically 0, which reads as agreement
  but is vacuous).
- Interval construction: what the quantile multiplier assumes (normal
  vs t), one- vs two-sided, and whether the interval matches the
  reported std_error or comes from a different computation.

**Assumptions and data structure.**
- Rows treated as independent that share a unit: multi-arm trials
  contributing several correlated contrasts, clustered or repeated
  measurements, overlapping populations. If one code path corrects for
  this and another doesn't, any "same model" claim between them is
  false for exactly these inputs.
- Silent scope: code correct under an assumption the input can violate
  without triggering any refusal or warning. The fix to demand is
  detection, not necessarily a bigger model.
- Missing values: does every aggregation define its NA behaviour, or
  inherit whatever the language does?

**Randomness and reproducibility.**
- Seed discipline: seeded paths reproducible; unseeded paths honest
  about it. Functions must not silently consume or reset the caller's
  global RNG state (draw-a-default-seed idioms do).
- Monte Carlo error: every tolerance in a stochastic test needs a
  justification relating it to the MCSE at the settings used —
  otherwise it is either flaky or vacuous.
- Convergence handling: are R-hat/ESS/divergences checked and surfaced,
  and does anything downstream proceed on an unconverged fit?

**Numerics.**
- `solve(A) %*% b` vs `solve(A, b)`; explicit inverses of
  near-singular systems; catastrophic cancellation; sums of logs vs
  log-sum-exp; comparisons of floats with `==`.
- Weakly identified directions: what happens at the boundary
  (disconnected graph, zero variance, single observation)?

**Tests as evidence.** Review the tests as statistical instruments, not
as code.
- Circular oracles: a test that recomputes the implementation's own
  formula and compares — it verifies transcription, not correctness.
  Real oracles: hand-computed closed forms, an independent
  implementation, simulation truth, known invariances.
- Calibration claims (coverage, type-I error) need simulation tests
  with bounds tight enough to fail under the miscalibration the code
  could plausibly have — a 0.85 bound on nominal 95% coverage hides a
  lot.
- Agreement tests between implementations: check *what* is compared.
  Two engines agreeing on point estimates says nothing about their
  standard errors.

**Claims vs code.** Every quoted promise from step 1 is a claim to
verify literally: "identical likelihood" means identical for every
admissible input, "propagates uncertainty" means the reported variance
is the correct variance, not merely a larger one.

## Severity ladder

Rank findings by consequence for inference, not by code size:

1. **Wrong number** — a shipped estimate is incorrect for some
   admissible input.
2. **Wrong uncertainty** — estimates right; SEs, intervals, or
   diagnostics misstate the evidence (either direction).
3. **Silent scope** — correct under an assumption inputs can violate
   with no detection.
4. **Fragile evidence** — tests or docs claim more than they verify
   (circular oracles, untested promises, loose calibration bounds).
5. **Statistical housekeeping** — RNG hygiene, numerical robustness,
   inefficiency; anything that degrades trust or stability without yet
   producing a wrong number.

## Report format

Open with the verdict in two or three sentences: is the statistics
sound, and what is the worst finding. Then:

```
## Findings (worst first)
### N. [severity] One-line claim
- Where: file:line
- What the code does / what the mathematics gives
- Demonstration: the two disagreeing numbers and how they were produced
- Consequence: who is misled, and when
- Fix: the smallest correct change
## Claim inventory
| claim | where | verdict |
## What the tests do and do not establish
```

Verified claims earn a sentence too — a review that only lists defects
hides how much of the package is solid, and the inventory is what makes
the review's coverage auditable.
