# 01: Surface covariance in the fit contract, both engines

**What to build:** A fitted surface exposes not just each drug's
standard error but the full covariance of the estimated effects versus
the reference, from either engine, through the engine-agnostic fit
contract. A researcher (or downstream package code) can fit with
netmeta or Stan and read the same covariance component in the same
shape — the prefactoring that makes correct anchoring (02) and
leverage-aware residuals (05) easy.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] A netmeta surface fit carries the common-effect covariance of
      every drug versus the reference as a new fit component; the
      reference drug's row and column are exact zeros
- [ ] A Stan surface fit (and an anchored Stan fit) carries the
      posterior covariance of theta in the same component and shape
- [ ] The netmeta covariance matches an independent weighted
      normal-equations oracle (extend the existing hand WLS oracle to
      return the full covariance, not just its diagonal) to numerical
      tolerance; the Stan covariance matches within a Monte Carlo
      tolerance justified at the test's sampler settings
- [ ] The component's diagonal is consistent with the existing
      `std_error` column on both engines
- [ ] The extension is additive: the effects schema and existing fit
      components are unchanged, and the whole existing test suite
      passes untouched
- [ ] The formats reference documents the new component (shape, row
      order, reference-row convention, engine parity)
