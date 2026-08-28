# directeffect

<!-- badges: start -->
[![R-CMD-check](https://github.com/ablack3033/directeffect/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ablack3033/directeffect/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## What it is

Comparative-effectiveness studies mostly estimate *relative* effects:
drug A versus drug B. Each estimate constrains only the **difference**
between two drugs' latent direct effects — but the researcher wants
each drug's *direct* effect versus placebo, which almost none of the
evidence measures. A null A-vs-B comparison does not mean either drug
is safe; it only means the two drugs sit at the same height.

directeffect reconstructs direct effects from a network of comparative
estimates, in two deliberately separate steps:

1. **The surface** — `fit_surface()` estimates every drug's position
   *relative* to the others, from the comparisons alone. This answers:
   is my comparative evidence internally coherent, and what relative
   structure does it imply?
2. **The sea level** — `anchor_surface()` uses placebo-controlled
   estimates (anchors) to position the whole surface against
   placebo = 0, keeping the anchors' uncertainty. This answers: where
   does that structure sit in absolute terms?

Two independent estimation engines — frequentist
([netmeta](https://cran.r-project.org/package=netmeta)) and Bayesian
(Stan via rstan) — sit behind one interface and are required by the
test suite to reconstruct the same surface, estimates and standard
errors alike, so the implementations continuously validate each other.
For networks with one comparison per study the two engines fit the
identical likelihood; multi-arm trials (several comparisons sharing a
`study_id`) are supported by the netmeta engine only in this version —
netmeta models the within-study correlation, and the Stan engine
refuses such networks rather than mis-treat correlated rows as
independent. Diagnostics (leverage-aware edge residuals, cycle
consistency, connectivity/identifiability reports) and plots come
built in.

## Installation

``` r
# install.packages("pak")
pak::pak("ablack3033/directeffect")
```

For the frequentist engine also install `netmeta`; for the Bayesian
engine, `rstan`.

## How to use it

The package ships a realistic but entirely simulated statin evidence
base (real drug names, fake numbers — the generating truth ships as
`example_truth`): twelve head-to-head estimates among six statins plus
three placebo-controlled anchor trials.

``` r
library(directeffect)

# 1. Build a validated evidence network — no model fitted yet.
de <- direct_effect_network(
  example_comparisons,
  anchors = example_anchors,
  effect_measure = "HR"
)
de
#> <directeffect_network>
#>   Effect measure: HR (log scale)
#>   Drugs:          6
#>   Comparisons:    12
#>   Anchors:        3
#>   Components:     1

# 2. Check what is identifiable before fitting.
check_connectivity(de)

# 3. Fit the relative surface (anchors deliberately ignored),
#    then position it against placebo as a separate step.
absolute <- anchor_surface(fit_surface(de, engine = "netmeta"))
absolute
#> <directeffect_fit>
#>   Engine:         netmeta
#>   Effect measure: HR (log scale)
#>   Reference:      placebo = 0 (sea level; absolute direct effects)
#>
#>          drug estimate std_error  lower  upper
#>  atorvastatin   -0.324     0.046 -0.414 -0.233
#>   fluvastatin   -0.094     0.084 -0.259  0.071
#>    lovastatin   -0.129     0.088 -0.302  0.044
#>   pravastatin   -0.190     0.060 -0.308 -0.072
#>  rosuvastatin   -0.423     0.068 -0.556 -0.290
#>   simvastatin   -0.278     0.060 -0.396 -0.160

# 4. Plots and diagnostics.
plot_sea_level(absolute)
edge_residuals(fit_surface(de, engine = "netmeta"))
```

Only three of the six statins were ever measured against placebo; the
other three get absolute effects through the network. Because the
example data are simulated, this reconstruction is graded against the
known truth in the vignettes and the test suite: all six 95% intervals
cover the true values and the estimated ordering matches the true
ordering exactly.

## Data formats

All estimates enter and leave the package on the **log scale** for
multiplicative measures (HR, RR, OR): pass `log(1.02)`, not `1.02`. The
full reference, including validation rules, is at `?directeffect_formats`.

**Input — comparisons** (one row per drug-vs-drug estimate; required
by `direct_effect_network()`):

| column | type | meaning |
|---|---|---|
| `study_id` | character | study/analysis the estimate came from |
| `target` | character | drug the estimate is for |
| `comparator` | character | drug compared against (≠ `target`) |
| `estimate` | numeric | log-scale effect, target vs comparator |
| `std_error` | numeric > 0 | its standard error |

Extra columns are preserved; differing values in estimand-describing
columns (`population`, `database`, `design`, …) trigger a
transportability warning.

**Input — anchors** (optional; one row per placebo-controlled
estimate): `study_id`, `drug` (must appear in a comparison),
`reference` (normally `"placebo"`), `estimate`, `std_error`.

**Output — `fit$effects`** (both engines, anchored or not): one row per
drug with `drug`, `estimate`, `std_error`, `lower`, `upper`, `scale`
(always `"log"`), `reference` (a drug fixed at 0, or `"placebo"` after
anchoring), `engine`. Stan fits append `median`, `mean`, `sd`, `q025`,
`q975`, `rhat`, `ess_bulk`, `ess_tail`.

## Learn more

- `vignette("directeffect")` — getting started: the full workflow on
  the example data.
- `vignette("how-it-works")` — the method explained from scratch, and
  the package graded against the example data's known truth.
- `vignette("netmeta-vs-stan")` — the same surface reconstructed by
  two independent engines.
- `vignette("validating-the-surface")` — edge residuals, cycle
  consistency, and simulation-based recovery.

## Status

Early development — the v0.1 core (network object, dual estimation
engines, anchoring, and diagnostics) is complete. See
`.scratch/v0.1-core/` for the spec and ticket breakdown.
