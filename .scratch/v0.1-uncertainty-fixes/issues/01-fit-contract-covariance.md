# 01: Surface covariance in the fit contract, both engines

**What to build:** A fitted surface exposes not just each drug's
standard error but the full covariance of the estimated effects versus
the reference, from either engine, through the engine-agnostic fit
contract. A researcher (or downstream package code) can fit with
netmeta or Stan and read the same covariance component in the same
shape — the prefactoring that makes correct anchoring (02) and
leverage-aware residuals (05) easy.

**Blocked by:** None (can start immediately).

**Status:** done

- [x] A netmeta surface fit carries the common-effect covariance of
      every drug versus the reference as a new fit component; the
      reference drug's row and column are exact zeros
- [x] A Stan surface fit (and an anchored Stan fit) carries the
      posterior covariance of theta in the same component and shape
- [x] The netmeta covariance matches an independent weighted
      normal-equations oracle (extend the existing hand WLS oracle to
      return the full covariance, not just its diagonal) to numerical
      tolerance; the Stan covariance matches within a Monte Carlo
      tolerance justified at the test's sampler settings
- [x] The component's diagonal is consistent with the existing
      `std_error` column on both engines
- [x] The extension is additive: the effects schema and existing fit
      components are unchanged, and the whole existing test suite
      passes untouched
- [x] The formats reference documents the new component (shape, row
      order, reference-row convention, engine parity)

## Comments

Done. `new_directeffect_fit()` gains a `covariance` component (drug ×
drug, effects-table order, dimnames asserted) placed right after
`effects`. netmeta adapter derives it from `Lplus.matrix.common` (the
Moore–Penrose pseudoinverse of the weighted Laplacian):
`Cov(i-r, j-r) = L+_ij − L+_ir − L+_rj + L+_rr` — indexed by treatment,
so no sign bookkeeping, and it reflects netmeta's multi-arm variance
adjustments; matched the extended hand-WLS oracle to 1e-16 and the
three-arm fixture exactly. Stan adapter computes `cov()` over the
retained theta draws (`stan_theta_covariance()`), for surface and
anchored fits alike. Reference row/col forced to exact zeros on
unanchored fits. Anchored netmeta fits carry the propagated absolute
covariance too (needed by 02/05 and engine parity). The extension is
additive — effects schema untouched — but the three tests that assert
`names(fit)` verbatim were updated to include the new component, and
`wls_surface()` in the test helpers now returns the full covariance
(plus design and weights, reused by the hat-matrix oracle in 05).
Formats page documents shape, order, reference-row convention, and
engine parity.
