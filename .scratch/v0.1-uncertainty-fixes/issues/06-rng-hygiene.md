# 06: Fitting never touches the caller's random-number state

**What to build:** A researcher who calls `set.seed()` and then fits
with the Stan engine — without passing a seed — gets exactly the same
subsequent random numbers as if the fit had never happened, matching
the standard the simulation module already meets. Explicitly seeded
fits stay reproducible.

**Blocked by:** None (can start immediately).

**Status:** done

- [x] An unseeded Stan surface or anchored fit leaves `.Random.seed`
      identical before and after (mirroring the existing "simulation
      does not disturb the global random-number state" test)
- [x] No draw is taken from R's global RNG on the caller's behalf when
      no seed is supplied — the sampler's own seeding applies
- [x] Two fits with the same explicit seed produce identical posterior
      summaries; the seed argument's documentation reflects the new
      default behavior
- [x] Existing seeded tests (mandatory Stan example, equivalence tests)
      pass unchanged

## Comments

Done. Both Stan fitters' `seed` defaults changed from
`sample.int(.Machine$integer.max, 1)` (a draw from the caller's
stream) to `NULL`; `stan_sampling_seed()` then derives the sampler
seed from the wall clock and process id, touching no R RNG state, so
successive unseeded fits still differ. Measurement during
implementation showed `rstan::sampling()` perturbs `.Random.seed` even
when handed an explicit seed, so `sampling_preserving_rng()` snapshots
and restores the caller's state around every sampling call —
mirroring the simulation module's standard. Tests: unseeded surface
and anchored fits leave `.Random.seed` identical (the mirror of the
simulate test); the same explicit seed reproduces identical effects
and covariance; existing seeded tests unchanged. `fit_surface()` and
`anchor_surface()` docs describe the new default.
