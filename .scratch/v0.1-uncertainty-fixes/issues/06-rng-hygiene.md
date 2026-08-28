# 06: Fitting never touches the caller's random-number state

**What to build:** A researcher who calls `set.seed()` and then fits
with the Stan engine — without passing a seed — gets exactly the same
subsequent random numbers as if the fit had never happened, matching
the standard the simulation module already meets. Explicitly seeded
fits stay reproducible.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] An unseeded Stan surface or anchored fit leaves `.Random.seed`
      identical before and after (mirroring the existing "simulation
      does not disturb the global random-number state" test)
- [ ] No draw is taken from R's global RNG on the caller's behalf when
      no seed is supplied — the sampler's own seeding applies
- [ ] Two fits with the same explicit seed produce identical posterior
      summaries; the seed argument's documentation reflects the new
      default behavior
- [ ] Existing seeded tests (mandatory Stan example, equivalence tests)
      pass unchanged
