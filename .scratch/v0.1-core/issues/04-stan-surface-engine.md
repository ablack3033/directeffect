# 04: Stan surface engine

**What to build:** `fit_surface(de, engine = "stan")` fits the identical common-effect model with the unanchored Stan program — arbitrary reference fixed to 0, deliberately weak normal(0, 5) priors — and returns the same `directeffect_fit` contract, extended with posterior summaries (`median`, `mean`, `sd`, `q025`, `q975`) and convergence diagnostics (`rhat`, `ess_bulk`, `ess_tail`).

**Blocked by:** 02 (Direct-effect network object)

**Status:** done

- [x] Unanchored Stan model implements y ~ N(theta_target − theta_comparator, se²) with the reference constraint and weak priors from the design
- [x] Same fit contract as the netmeta engine: downstream code cannot tell the engines apart except via the `engine` field and the extra Bayesian columns
- [x] Convergence diagnostics included in the effects table, with a warning when rhat/ess indicate trouble
- [x] On the spec's 3-drug example, posterior means match the hand-computed relative surface within MCMC tolerance under a fixed seed
- [x] Stan interface decided (cmdstanr vs. rstan) with CI feasibility and installation friction in mind; record the choice and rationale in this file's Comments

## Comments

Done. **Interface decision: rstan.** Rationale: rstan is a
self-contained CRAN package — `setup-r-dependencies` installs it as a
binary in CI with no extra workflow steps, whereas cmdstanr is off-CRAN
(extra repo) and needs a CmdStan toolchain download/build step at setup
time, which the sandbox's network policy cannot even reach. The engine
is fully behind the adapter in R/engine-stan.R, so swapping to cmdstanr
later touches one file.

Model at inst/stan/surface.stan, verbatim from the design except the
index arrays are `target_idx`/`comparator_idx` — `target` is a reserved
word in Stan and the design's pseudocode would not compile. The adapter
orders drugs so the chosen reference is index 1 (theta[1] = 0), then
reports effects in network order. Compiled models are cached per
session. `estimate`/`std_error` are posterior mean/sd; reference row is
pinned at exactly 0 with NA convergence diagnostics; rhat > 1.01 or
ess < 400 warns. Sandbox note: Debian's r-cran-bh strips Boost headers,
so system libboost-dev is symlinked into BH/include for local compiles;
CI uses the real CRAN BH and needs nothing.
