# 04: Stan surface engine

**What to build:** `fit_surface(de, engine = "stan")` fits the identical common-effect model with the unanchored Stan program — arbitrary reference fixed to 0, deliberately weak normal(0, 5) priors — and returns the same `directeffect_fit` contract, extended with posterior summaries (`median`, `mean`, `sd`, `q025`, `q975`) and convergence diagnostics (`rhat`, `ess_bulk`, `ess_tail`).

**Blocked by:** 02 (Direct-effect network object)

**Status:** ready-for-agent

- [ ] Unanchored Stan model implements y ~ N(theta_target − theta_comparator, se²) with the reference constraint and weak priors from the design
- [ ] Same fit contract as the netmeta engine: downstream code cannot tell the engines apart except via the `engine` field and the extra Bayesian columns
- [ ] Convergence diagnostics included in the effects table, with a warning when rhat/ess indicate trouble
- [ ] On the spec's 3-drug example, posterior means match the hand-computed relative surface within MCMC tolerance under a fixed seed
- [ ] Stan interface decided (cmdstanr vs. rstan) with CI feasibility and installation friction in mind; record the choice and rationale in this file's Comments
