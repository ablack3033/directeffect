# 05: Leverage-aware edge residuals

**What to build:** A methodologist screening `edge_residuals()` sees
z-values with unit variance under coherence, a leverage column saying
how much each comparison checked itself, and `NA` — not a reassuring
zero — for a bridge comparison that nothing in the network can
corroborate. "Uncheckable" is never displayed as "perfectly
consistent".

**Blocked by:** 01 (fit-contract covariance).

**Status:** done

- [x] The residuals table gains a leverage column derived from the
      surface covariance (how much of each observation's variance is
      absorbed by its own prediction)
- [x] Standardized residuals divide by the residual's actual standard
      deviation (comparison variance minus prediction variance), not
      the raw comparison SE
- [x] The standardization is verified against an independent hat-matrix
      oracle from the weighted normal equations on a hand-computable
      network
- [x] Bridge fixture (triangle plus a pendant comparison): the bridge
      row reports leverage 1 and standardized residual `NA`, and
      changing the bridge's estimate can never make it display as
      consistent
- [x] Documentation explains the leverage column and that `NA` means
      uncorroborated, with the interpretive warning moved accordingly
- [x] Both engines' fits produce the same residual table shape, and the
      multi-arm caveat from ticket 04 (if landed) is preserved — no
      coordination needed beyond not clobbering each other's warnings

## Comments

Done. `edge_residuals()` consumes the fit-contract covariance:
`leverage = Var(pred)/se²` and
`standardized_residual = residual / sqrt(se² − Var(pred))`, verified
against an independent hat-matrix oracle (H = X(X'WX)⁻¹X'W on an
asymmetric 4-drug network: leverage = diag(H), residual SD =
diag((I−H)V(I−H)')) to 1e-8. Bridge rows are identified structurally —
`igraph::bridges()` on the comparisons multigraph — rather than by
thresholding leverage, so the NA is exact on either engine and immune
to Monte Carlo noise in a posterior covariance; parallel edges (which
corroborate each other) are correctly not bridges. The bridge fixture
asserts leverage 1 and NA at bridge estimates 0.2 and 5.0 alike. Table
gains `leverage` as the last column; docs (function + formats page +
validating-the-surface vignette) explain leverage and that NA means
uncorroborated. Multi-arm warning from 04 fires here too;
`assert_fit_covariance()` gives pre-covariance fit objects a clear
error.
