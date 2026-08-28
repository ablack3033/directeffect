# directeffect

<!-- badges: start -->
[![R-CMD-check](https://github.com/ablack3033/directeffect/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ablack3033/directeffect/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

directeffect estimates latent direct drug effects from a network of
comparative-effect estimates. Comparative estimates identify only
*differences* between drugs — the effect surface. Absolute placebo anchors
position that surface against the placebo-defined sea level, turning
relative evidence into absolute direct effects. The package provides
frequentist and Bayesian estimation engines behind a single interface,
plus diagnostics for network coherence and identifiability.

## Installation

You can install the development version of directeffect from GitHub:

``` r
# install.packages("pak")
pak::pak("ablack3033/directeffect")
```

## Status

Early development — the v0.1 core (network object, dual estimation
engines, anchoring, and diagnostics) is being built. See
`.scratch/v0.1-core/` for the spec and ticket breakdown.
