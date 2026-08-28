# 07: Recovery metrics computed over the same drugs

**What to build:** A methodologist reading `validate_recovery()` gets
all four summaries — bias, RMSE, coverage, rank correlation — over the
same set of drugs: the free ones. The pinned reference's exact (0, 0)
pair can no longer move Spearman's ρ by ±0.3 while the other metrics
exclude it.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Rank correlation is computed over the free drugs only for surface
      fits; anchored fits (no pinned drug) are unchanged
- [ ] A test on a noisy simulated network asserts the free-drugs value
      and would fail under the old all-drugs computation (a case where
      the two demonstrably differ, as in the review)
- [ ] Documentation states that the reference row is excluded from all
      four metrics for surface fits
- [ ] The existing recovery tests (low-noise regime, anchored path)
      pass, with expectations updated only where the metric definition
      changed
