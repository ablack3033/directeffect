# 05: Leverage-aware edge residuals

**What to build:** A methodologist screening `edge_residuals()` sees
z-values with unit variance under coherence, a leverage column saying
how much each comparison checked itself, and `NA` — not a reassuring
zero — for a bridge comparison that nothing in the network can
corroborate. "Uncheckable" is never displayed as "perfectly
consistent".

**Blocked by:** 01 (fit-contract covariance).

**Status:** ready-for-agent

- [ ] The residuals table gains a leverage column derived from the
      surface covariance (how much of each observation's variance is
      absorbed by its own prediction)
- [ ] Standardized residuals divide by the residual's actual standard
      deviation (comparison variance minus prediction variance), not
      the raw comparison SE
- [ ] The standardization is verified against an independent hat-matrix
      oracle from the weighted normal equations on a hand-computable
      network
- [ ] Bridge fixture (triangle plus a pendant comparison): the bridge
      row reports leverage 1 and standardized residual `NA`, and
      changing the bridge's estimate can never make it display as
      consistent
- [ ] Documentation explains the leverage column and that `NA` means
      uncorroborated, with the interpretive warning moved accordingly
- [ ] Both engines' fits produce the same residual table shape, and the
      multi-arm caveat from ticket 04 (if landed) is preserved — no
      coordination needed beyond not clobbering each other's warnings
