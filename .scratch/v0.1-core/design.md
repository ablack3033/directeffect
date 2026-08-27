# directeffect — full design document

> Source: design discussion, 2026-08-27. This is the complete design covering v0.1–v0.4.
> The v0.1 spec cut from it lives at `spec.md` in this directory. Math is written in
> plain notation (`theta_i`, `y_ij ~ N(...)`) rather than LaTeX.

## Purpose

`directeffect` estimates latent direct drug effects from networks of observed comparative-effect estimates.

The package is designed for comparative-effectiveness studies where the primary observed quantities are relative effects between drugs:

    y_ij ≈ theta_i - theta_j

where:

- `y_ij` is an observed comparative effect for drug i versus drug j;
- `theta_i` is the latent direct effect of drug i;
- `theta_j` is the latent direct effect of drug j.

The package separates two orthogonal inference problems:

1. **Surface estimation** — determining the relative positions of drugs from drug-versus-drug comparisons.
2. **Sea-level estimation** — locating that surface on an absolute scale using placebo-controlled trials or other external absolute-effect anchors.

The package provides parallel frequentist and Bayesian implementations so that results can be independently reproduced and compared.

## 1. Conceptual model

### 1.1 Surface

For each drug i, define a latent direct effect `theta_i`.

For a comparative-effect estimate involving drugs i and j:

    y_ij ~ N(theta_i - theta_j, s_ij^2)

where `s_ij` is the standard error of the observed comparative estimate.

The set `{theta_1, ..., theta_K}` defines the **direct-effect surface**.

Only differences are identified from comparative-effect data: `theta_i - theta_j`. Therefore `theta_i' = theta_i + c` for every drug produces exactly the same comparative-effect likelihood.

The surface has a shape but no absolute vertical position.

## 2. Sea level

Absolute-effect evidence anchors the surface.

For an external estimate of drug i versus placebo/no treatment:

    a_i ~ N(theta_i - theta_P, s_ai^2)

The natural convention is `theta_P = 0`. Then:

    a_i ~ N(theta_i, s_ai^2)

These observations establish the **sea level** of the relative-effect surface.

Several anchors may be supplied. They should normally retain their uncertainty rather than fixing a drug effect exactly.

## 3. Fundamental separation

The package should maintain a strong conceptual distinction between:

    comparative evidence
            │
            ▼
     relative surface
     θA - θB
     θA - θC
     θB - θD
            │
            │ absolute anchors
            ▼
     positioned surface
     θA, θB, θC, θD
     relative to placebo = 0

This distinction should appear throughout the API.

A user should be able to estimate and inspect the surface without any absolute anchors. The user should then be able to apply one or more anchors separately.

This allows questions such as *"Are my comparative-effect estimates internally coherent?"* to remain separate from *"Where should this coherent surface sit relative to placebo?"*

## 4. Primary data structures

### Comparative effects

One row represents one estimated drug-versus-drug effect.

Minimum columns:

    comparisons
    #> study_id
    #> target
    #> comparator
    #> estimate
    #> std_error

For multiplicative effect measures such as HR, RR or OR, `estimate` should normally be represented internally on the log scale.

Example:

```r
comparisons <- tibble::tribble(
  ~study_id, ~target, ~comparator, ~estimate, ~std_error,
  "S1",      "A",     "B",         log(1.02), 0.07,
  "S2",      "A",     "C",         log(1.34), 0.09,
  "S3",      "B",     "C",         log(1.29), 0.08
)
```

Optional metadata should be preserved: `database`, `population`, `design`, `drug_class`, `followup`, `outcome`, `analysis_id`, `negative_control`, `calibrated`, `n_target`, `n_comparator`.

### Absolute anchors

A second table represents evidence about absolute/direct effects:

    anchors
    #> study_id
    #> drug
    #> reference
    #> estimate
    #> std_error

Usually `reference = "placebo"`.

Example:

```r
anchors <- tibble::tribble(
  ~study_id, ~drug, ~reference, ~estimate, ~std_error,
  "RCT1",    "C",   "placebo",  log(1.20), 0.04,
  "RCT2",    "D",   "placebo",  log(0.91), 0.06
)
```

## 5. Package object

Create an explicit problem specification before fitting:

```r
de <- direct_effect_network(
  comparisons,
  anchors = anchors,
  effect_measure = "HR"
)
```

Return:

    direct_effect_network
    ├── comparisons
    ├── anchors
    ├── treatments
    ├── graph
    ├── components
    ├── effect_measure
    └── metadata

This object should contain no fitted model. Its job is to define and validate the causal-effect network.

## 6. Primary API

The API should remain small.

- Construct network: `direct_effect_network()`
- Fit relative surface: `fit_surface()`

```r
surface_netmeta <- fit_surface(de, engine = "netmeta")
surface_stan    <- fit_surface(de, engine = "stan")
```

At this stage the model should deliberately ignore absolute anchors. Identification can use an arbitrary reference (`theta_reference = 0`) because only relative positions matter.

- Anchor the surface: `anchor_surface()`

```r
fit <- anchor_surface(surface_stan, anchors = anchors)
```

This makes the conceptual decomposition explicit:

```r
surface  <- fit_surface(de)
absolute <- anchor_surface(surface, anchors)
```

A convenience function may eventually combine them — `fit_direct_effects(de, engine = "stan")` — but internally it should preserve the separation.

## 7. Frequentist engine

Use **netmeta**.

Input mapping:

    directeffect       netmeta
    ---------------------------
    estimate      ->   TE
    std_error     ->   seTE
    target        ->   treat1
    comparator    ->   treat2
    study_id      ->   studlab

Conceptually:

```r
netmeta::netmeta(
  TE = estimate,
  seTE = std_error,
  treat1 = target,
  treat2 = comparator,
  studlab = study_id,
  data = comparisons,
  reference.group = reference,
  common = TRUE,
  random = TRUE
)
```

The wrapper should convert the resulting relative estimates back into the common `directeffect_fit` representation. Do not expose netmeta internals throughout the package. Use it behind an adapter: `fit_surface_netmeta()`. This is important information hiding.

## 8. Stan engine

Use the same statistical model explicitly.

Base model:

    y_k ~ N(theta_target[k] - theta_comparator[k], s_k^2)

For identifiability during surface fitting: `theta_1 = 0`. All other effects are relative to this arbitrary origin.

Stan pseudocode:

```stan
data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> target;
  array[N] int<lower=1, upper=K> comparator;
  vector[N] y;
  vector<lower=0>[N] se;
}
parameters {
  vector[K - 1] theta_free;
}
transformed parameters {
  vector[K] theta;
  theta[1] = 0;
  for (k in 2:K)
    theta[k] = theta_free[k - 1];
}
model {
  theta_free ~ normal(0, 5);
  for (n in 1:N)
    y[n] ~ normal(theta[target[n]] - theta[comparator[n]], se[n]);
}
```

The initial prior should be deliberately weak so that comparison against netmeta is meaningful.

## 9. Random-effects model

Add:

    y_ij,k ~ N(theta_i - theta_j, s_ij,k^2 + tau^2)

API:

```r
fit_surface(de, engine = "stan", heterogeneity = "common")
```

The initial implementation should have one global tau. Do not immediately add treatment-specific or comparison-specific heterogeneity.

## 10. Stan absolute model

For simultaneous surface + sea-level estimation:

    y_k ~ N(theta_target[k] - theta_comparator[k], s_k^2 + tau^2)

and

    a_m ~ N(theta_drug[m], s_am^2)

Once anchors exist, no arbitrary treatment constraint is needed — the anchors identify the absolute location.

This should be a **separate Stan model** from the unanchored model rather than filling one model with conditionals. Prefer two simple deep modules over one complicated universal model.

## 11. Common fit representation

Both engines should return the same user-facing class: `directeffect_fit`, with components:

    fit$effects
    fit$comparisons
    fit$anchors
    fit$heterogeneity
    fit$diagnostics
    fit$engine
    fit$engine_fit
    fit$network

`fit$effects` should have a stable tidy format: `drug`, `estimate`, `std_error`, `lower`, `upper`, `scale`, `reference`, `engine`.

Stan additionally supplies: `median`, `mean`, `sd`, `q025`, `q975`, `rhat`, `ess_bulk`, `ess_tail`.

The package should never require downstream plotting functions to know whether the originating fit came from Stan or netmeta.

## 12. Engine comparison

This deserves a first-class function:

```r
compare_engines(netmeta_fit, stan_fit)
```

Output: `drug`, `netmeta`, `stan_mean`, `difference`, `standardized_difference`.

Primary diagnostic plot: `plot_engine_comparison()` with x = theta_hat(netmeta), y = E[theta | Stan], including the identity line y = x.

For the base common-effect model, the two should agree to numerical tolerance once priors become effectively flat and parameterization is reconciled. This should become a package unit/integration test.

## 13. Surface diagnostics

These are more important than sophisticated estimation.

### Edge residuals

For every observed comparison:

    r_ij = y_ij - (theta_hat_i - theta_hat_j)

Function: `edge_residuals(fit)`.

Output: `target`, `comparator`, `observed`, `predicted`, `residual`, `standardized_residual`.

## 14. Cycle consistency

For every closed path i1 → i2 → ... → im → i1, calculate the signed sum C = Σ y_ij around the cycle. Under perfect consistency, C = 0.

Functions: `cycle_consistency(fit)`, `plot_cycle_consistency(fit)`.

Do not enumerate every possible cycle for large graphs. Use a **cycle basis**.

## 15. Direct versus indirect prediction

Given an observed edge A–B, temporarily exclude it. Estimate `theta_hat_A - theta_hat_B` from all remaining paths. Compare against the held-out edge.

API: `leave_one_edge_out(fit)`.

Output: `edge`, `observed`, `predicted`, `prediction_error`, `standardized_error`.

This is one of the central validation mechanisms.

## 16. Anchor validation

Absolute anchors should be usable either for fitting or validation.

Allow `anchor_role = c("fit", "validate")`:

```r
anchors <- tibble(
  drug = ...,
  estimate = ...,
  std_error = ...,
  role = c("fit", "fit", "validate")
)
```

Then `validate_anchors(fit)` performs held-out prediction.

The key question: **can relative observational evidence plus other placebo anchors recover a placebo-controlled effect that was hidden from the model?**

## 17. Leave-one-anchor-out validation

Provide `loo_anchors(de)`. For each anchor i:

1. remove it;
2. fit using all remaining anchors;
3. predict theta_i;
4. compare prediction with held-out RCT estimate.

Return: `drug`, `observed`, `observed_se`, `predicted`, `predicted_se`, `error`, `coverage`.

This should probably become the most important validation output of the package.

## 18. Connectivity diagnostics

Before fitting, calculate connected components: `check_connectivity(de)`.

A disconnected component without an anchor can have its relative surface estimated but cannot have an absolute sea level. Report this explicitly:

    Component 1:
      23 drugs
      104 comparisons
      3 absolute anchors
      absolute effects identifiable
    Component 2:
      7 drugs
      12 comparisons
      0 absolute anchors
      relative effects identifiable only

Never silently choose an arbitrary sea level and present it as an absolute effect.

## 19. Surface plots

For a single outcome the surface is one-dimensional: `plot_effect_surface(fit)`.

Conceptually:

    harmful <------------------------------------> beneficial
         A        C/D                F              B
    -----|---------||----------------|--------------|----

with uncertainty intervals. If using log HR, 0 corresponds to HR = 1.

The user should be able to select `scale = "log"` or `scale = "natural"`.

## 20. Network plot

`plot_network(de)`. Nodes: drugs. Edges: comparative-effect estimates.

Optional visual encodings: edge thickness = precision; node shape = presence of absolute anchor; node size = number of comparisons; drug-class grouping; edge label = observed effect.

Use igraph/ggraph rather than embedding graph logic inside plotting code.

## 21. Surface + sea-level plot

This should become a signature visualization.

Conceptually:

                            DIRECT EFFECT
                                 ↑
                                 Drug A
                                   ●
                                  / \
                                 /   \
                         Drug B ●     ● Drug C
                                 \   /
                                  \ /
                             Drug D ●
    --------------------- placebo = 0 ---------------------
                             SEA LEVEL

The surface is inferred predominantly from comparative evidence. Absolute anchors determine where the entire structure sits relative to the placebo line.

For a single outcome, render it one-dimensionally rather than pretending that a 2-D graph layout itself represents causal distance.

## 22. Anchor influence

Calculate how much each anchor determines absolute position: `anchor_influence(fit)`.

For each anchor a, refit without it and calculate shifts:

    Delta_i = theta_hat_i^(-a) - theta_hat_i

Plot: `plot_anchor_influence(fit)`.

This answers: *is our absolute conclusion dependent on one particular placebo trial?*

## 23. Comparator robustness

For a target drug: `comparator_robustness(fit, drug = "A")`.

Refit using: same-class comparators; outside-class comparators; random subsets; high-precision comparisons; low-precision comparisons; specified comparator groups.

This directly reflects the motivating scientific problem. The question is not *"do all comparisons agree?"* — it is *"do the different relative constraints lead to a stable inferred absolute position?"*

## 24. Evidence decomposition

For a particular drug: `explain_effect(fit, drug = "A")`.

Return the evidence paths that contribute most strongly to its inferred location. Possible output:

    Drug A estimated direct log-HR: 0.31
    Important constraints:
    A vs B        +0.03
    B vs placebo  +0.27
    A vs C        +0.17
    C vs placebo  +0.12
    A vs D        -0.08
    D vs B        +0.10
    B vs placebo  +0.27

This should emphasize interpretability. A latent estimate should remain traceable back to comparative evidence.

## 25. Negative controls

Negative-control calibration should be deliberately postponed until the core model is verified.

Initial extension:

```r
fit_surface(de, systematic_error = error_model)
```

with a separate object:

```r
error_model <- fit_systematic_error(negative_controls)
```

Possible model:

    y_ij ~ N(theta_i - theta_j + mu_b, s_ij^2 + sigma_b^2)

Later versions could allow `systematic_error = "database"` or `systematic_error = "design"`, but these should not enter v0.1. The causal surface model and bias model should remain separate modules.

## 26. Drug biology

The base model should not force drugs in the same class to have similar direct effects. Drug biology determines the inferred positions through the comparative evidence.

A later optional prior may introduce biological information:

    theta_i ~ N(mu_class(i), sigma_class^2)

API:

```r
fit_direct_effects(
  de,
  prior = hierarchical_class_prior(class = drug_class)
)
```

This must remain opt-in. The base estimator should make no assumption that drugs in the same pharmacological class have similar effects.

## 27. Important estimand restriction

A coherent scalar theta_i only exists if comparisons refer to sufficiently compatible causal estimands.

The package should therefore support metadata describing: `population`, `time_at_risk`, `outcome_definition`, `database`, `study_design`, `estimand` — and warn when obviously different specifications are combined.

The package should not imply that statistical connectivity establishes causal transportability.

## 28. Core classes

Keep the object system minimal:

- `directeffect_network`
- `directeffect_fit`
- `directeffect_validation`
- `directeffect_error_model`

Avoid creating separate classes for every intermediate operation.

## 29. Proposed package structure

    directeffect/
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── R/
    │   ├── network.R
    │   ├── fit.R
    │   ├── fit-netmeta.R
    │   ├── fit-stan.R
    │   ├── anchors.R
    │   ├── effects.R
    │   ├── diagnostics.R
    │   ├── validation.R
    │   ├── plotting.R
    │   └── utils.R
    ├── inst/
    │   └── stan/
    │       ├── surface.stan
    │       └── anchored_surface.stan
    ├── tests/
    │   └── testthat/
    │       ├── test-netmeta-stan-equivalence.R
    │       ├── test-identification.R
    │       ├── test-anchors.R
    │       ├── test-cycles.R
    │       ├── test-connectivity.R
    │       └── test-simulation-recovery.R
    ├── vignettes/
    │   ├── direct-effect-model.Rmd
    │   ├── netmeta-vs-stan.Rmd
    │   └── validating-the-surface.Rmd
    └── data/
        └── example_network.rda

## 30. Simulation framework

Simulation should exist from the first release.

```r
simulate_direct_effect_network(
  n_drugs = 20,
  n_comparisons = 100,
  n_anchors = 3,
  heterogeneity = 0.1,
  seed = 1
)
```

Procedure:

1. generate true theta_i;
2. generate a comparison graph;
3. calculate true contrasts delta_ij = theta_i - theta_j;
4. generate standard errors;
5. sample observed comparisons;
6. sample anchor estimates;
7. fit the model;
8. compare recovered versus known theta_i.

The simulation object should retain truth (`simulation$truth`) so `validate_recovery(fit, simulation)` can measure: bias; RMSE; interval coverage; rank correlation; edge prediction; anchor prediction.

## 31. Mandatory unit test

A deterministic low-noise example:

    true effects
    placebo = 0.0
    A       = 0.7
    B       = 0.7
    C       = 0.3

Comparisons:

    A - B = 0.0
    A - C = 0.4
    B - C = 0.4

Anchor:

    C - placebo = 0.3

Both engines must recover approximately A = 0.7, B = 0.7, C = 0.3.

This example directly protects the central conceptual property of the package: **a null active-comparator estimate does not imply a null direct drug effect.**

## 32. Engine equivalence test

Simulate a network under: Gaussian comparison estimates; known standard errors; no heterogeneity; weak Stan priors.

Require `|theta_i_stan - theta_i_netmeta| < epsilon` for every treatment.

This should run continuously in CI. It gives the project an independent implementation check from day one.

## 33. Suggested user workflow

```r
library(directeffect)

de <- direct_effect_network(
  comparisons = comparisons,
  anchors = placebo_trials,
  effect_measure = "HR"
)

plot_network(de)

freq  <- fit_surface(de, engine = "netmeta")
bayes <- fit_surface(de, engine = "stan")

compare_engines(freq, bayes)
cycle_consistency(freq)
edge_cv <- leave_one_edge_out(freq)

absolute <- anchor_surface(bayes, anchors = placebo_trials)

plot_effect_surface(absolute)
validate_anchors(absolute)
anchor_influence(absolute)
```

## 34. Version roadmap

### v0.1 — prove the mathematics

Implement only: network object; netmeta engine; Stan engine; common-effect model; absolute anchors; engine comparison; connectivity; edge residuals; cycle consistency; simulation; basic plots.

The goal: **demonstrate that the same direct-effect surface is reconstructed independently by netmeta and Stan.**

### v0.2 — prove robustness

Add: random effects; leave-one-edge-out; leave-one-anchor-out; comparator robustness; anchor influence; direct/indirect evidence diagnostics.

### v0.3 — model observational systematic error

Add: negative-control calibration; database-specific bias; design-specific bias.

### v0.4 — biological structure

Potentially add: class-level hierarchical priors; mechanism-of-action priors; multiple outcomes; multivariate direct-effect space.

Do not build v0.3 or v0.4 until validation of v0.1 and v0.2 demonstrates that the basic model is useful.

## 35. Package design principle

The package should expose the scientific abstraction:

    relative evidence
           ↓
    direct-effect surface
           +
    absolute anchors
           ↓
    absolute direct effects

It should not expose: netmeta objects, Stan matrices, igraph internals, MCMC implementation details — unless explicitly requested.

netmeta and Stan are implementation engines. The package abstraction is the direct-effect network.

## 36. Naming

`directeffect` is accurate and immediately understandable, but somewhat generic. Alternatives considered:

- **causalsurface** — captures the central model well (reconstruct the causal-effect surface from relative comparisons); leaves room for multiple outcomes later.
- **effectmap** — simple and memorable; maps relative treatment effects into a coherent set of direct effects. Downside: generic.
- **effectnet** — communicates the network structure clearly but sounds like generic network meta-analysis.
- **drugposition** — captures the triangulation idea but undersells the causal/statistical purpose.
- **causalmap** — good conceptually, although broader than the actual package.
- **directeffect** — probably the strongest scientific name. It says exactly what users want (`library(directeffect)`, `fit_direct_effects(...)`, `plot_direct_effects(...)`) and reads naturally in papers: *"Direct drug effects were estimated using the directeffect R package."*

Decision:

- **Package**: `directeffect`
- **Core abstraction**: `direct_effect_network`
- **Conceptual terminology**: *surface* and *sea level*

That combination keeps the public name scientifically sober while preserving the useful surveying/ocean metaphor internally and in explanatory material.
